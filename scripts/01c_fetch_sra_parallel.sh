#!/usr/bin/env bash
# Purpose: Fetch the three validation .sra objects from the SRA S3 mirror with 8 parallel byte-range streams,
#   resumable; prefetch cannot resolve accessions through the network sandbox.
# Input: results/srr_order.txt (run accessions); SRA S3 mirror over HTTPS.
# Output: data/sra/<run>.final.sra plus .ready markers; logs/01c_fetch.log.
# 01c: fetch .sra via 8-way parallel range requests -> fasterq-dump -> kb count -> verify -> clean up
#  Reason for route change (self-judgment ③): ENA FASTQ measured at 0.35-0.47 MB/s/stream, totaling 2.0 MB/s with 3 parallel streams
#  (155.3 GB = 21.6 h). SRA S3 achieves 5.54 MB/s with 8-way parallel ranges, and .sra totals a smaller 132.9 GB
#  (= 6.7 h). Adopted S3 + fasterq-dump to reach the count stage within 24 hours.
set -uo pipefail
export KALLISTO_BIN="$(command -v kallisto)"
export BUSTOOLS_BIN="$(command -v bustools)"
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
CH=$W/data/chunks; SRA=$W/data/sra; FQ=$W/data/fastq; CNT=$W/data/counts; R=$W/ref
L=$W/logs/01c_fetch.log; mkdir -p $CH $SRA $FQ $CNT
export HOME=$W/.home; mkdir -p $HOME
S3=https://sra-pub-run-odp.s3.amazonaws.com/sra
NPAR=8
say(){ echo "[$(date -u +%FT%TZ)] $*" | tee -a $L; }

get_size(){ curl -sI --max-time 60 "$S3/$1/$1" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r' | tail -1; }

fetch_parallel(){ # $1=acc
  local acc=$1 total part i st en got tries
  [ -s "$SRA/$acc.sra.done" ] && { say "$acc: .sra already assembled"; return 0; }
  total=$(get_size $acc); [ -n "$total" ] || { say "$acc: サイズ取得失敗"; return 1; }
  part=$(( (total + NPAR - 1) / NPAR ))
  say "$acc: $total bytes, $NPAR chunks x $part"
  for tries in 1 2 3; do
    for i in $(seq 0 $((NPAR-1))); do
      st=$((i*part)); en=$((st+part-1)); [ $en -ge $total ] && en=$((total-1))
      # Note: curl rejects the combination of -C - and -r and exits immediately (cause of the initial failure).
      # Compute the range by summing the already-fetched partial byte count ourselves, and append with >>.
      ( got=$(stat -f%z "$CH/${acc}.p$i" 2>/dev/null || echo 0)
        want=$((en-st+1))
        [ "$got" -ge "$want" ] && exit 0
        curl -sS -r $((st+got))-${en} --retry 3 --retry-delay 10 --max-time 86400 \
             "$S3/$acc/$acc" >> "$CH/${acc}.p$i" ) &
    done
    wait
    got=$(cat $CH/${acc}.p* 2>/dev/null | wc -c | tr -d ' ')
    [ "$got" = "$total" ] && break
    say "$acc: 試行 $tries 未完 ($got/$total) 再開"
  done
  got=$(cat $CH/${acc}.p* | wc -c | tr -d ' ')
  if [ "$got" = "$total" ]; then
    cat $CH/${acc}.p0 $CH/${acc}.p1 $CH/${acc}.p2 $CH/${acc}.p3 $CH/${acc}.p4 $CH/${acc}.p5 $CH/${acc}.p6 $CH/${acc}.p7 > "$SRA/$acc.final.sra"
    [ "$(stat -f%z $SRA/$acc.final.sra)" = "$total" ] && { rm -f $CH/${acc}.p*; touch "$SRA/$acc.sra.done"; say "$acc: DOWNLOAD OK $total"; return 0; }
  fi
  say "$acc: DOWNLOAD FAILED ($got/$total)"; return 1
}

dump_and_count(){ # $1=acc $2=label
  local acc=$1 lab=$2 out=$CNT/$2 r1 r2
  if [ -s "$out/counts_unfiltered/cells_x_genes.mtx" ]; then say "$lab: counted already"; return 0; fi
  if [ ! -s "$FQ/${acc}_1.fastq" ] && [ ! -s "$FQ/${acc}_1.fastq.gz" ]; then
    say "$lab: fasterq-dump start"
    for a in 1 2; do
      fasterq-dump --split-files --threads 6 --temp "$FQ/t_$acc" -O "$FQ" "$SRA/$acc.final.sra" >> $L 2>&1 && break
      say "$lab: fasterq-dump 試行 $a 失敗"
    done
    rm -rf "$FQ/t_$acc"
  fi
  r1=$(ls $FQ/${acc}_1.fastq.gz $FQ/${acc}_1.fastq 2>/dev/null | head -1)
  r2=$(ls $FQ/${acc}_2.fastq.gz $FQ/${acc}_2.fastq 2>/dev/null | head -1)
  [ -s "$r1" ] && [ -s "$r2" ] || { say "$lab: FASTQ 生成失敗"; return 1; }
  say "$lab: R1 $(head -2 $r1 | tail -1 | wc -c | tr -d ' ')bp 相当 / kb count start ($(du -ch $r1 $r2 | tail -1 | cut -f1))"
  mkdir -p $out
  for a in 1 2; do
    kb count --kallisto "$KALLISTO_BIN" --bustools "$BUSTOOLS_BIN" -i $R/index.idx -g $R/t2g.txt -x 10xv3 -o $out -t 6 -m 24G --overwrite "$r1" "$r2" > $W/logs/kb_$lab.log 2>&1 && break
    say "$lab: kb count 試行 $a 失敗"
  done
  if [ -s "$out/counts_unfiltered/cells_x_genes.mtx" ]; then
    python - "$out" "$lab" <<'PY' 2>&1 | tee -a $L
import sys, json, os
out, lab = sys.argv[1], sys.argv[2]
p = os.path.join(out, "counts_unfiltered", "cells_x_genes.mtx")
with open(p) as f:
    f.readline(); f.readline(); nb, ng, nnz = map(int, f.readline().split())
    tot = sum(int(float(l.rsplit(" ",1)[1])) for l in f)
rec = dict(label=lab, barcodes=nb, genes=ng, nonzero=nnz, total_umi=tot)
ri = os.path.join(out, "run_info.json")
if os.path.exists(ri):
    d = json.load(open(ri)); rec.update({k: d.get(k) for k in ("n_processed","n_pseudoaligned","p_pseudoaligned","n_targets")})
json.dump(rec, open(os.path.join(out,"verify.json"),"w"), indent=1); print("VERIFY", json.dumps(rec))
PY
    say "$lab: COUNT OK"
    rm -f $FQ/${acc}_*.fastq "$SRA/$acc.final.sra"
    say "$lab: 掃除完了 空き $(df -h $W | tail -1 | awk '{print $4}')"
  else
    say "$lab: COUNT FAILED"; return 1
  fi
}

say "=== 01c start 空き $(df -h $W | tail -1 | awk '{print $4}') ==="
for pair in "SRR32119705 CAP1" "SRR32119703 CAP3" "SRR32119704 CAP2"; do
  set -- $pair
  fetch_parallel $1 && dump_and_count $1 $2 || say "$2: スキップ（前段失敗）"
done
say "=== 01c DONE counted=$(ls $CNT/*/verify.json 2>/dev/null | wc -l | tr -d ' ')/3 ==="
