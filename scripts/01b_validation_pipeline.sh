#!/usr/bin/env bash
# Purpose: Re-quantify validation cohort PRJNA1213680 (CAP1–3); fetch FASTQ from ENA or SRA S3, run kallisto |
#   bustools, validate, clean
# Input: none (fetches SRR32119705, SRR32119704, SRR32119703 from ENA/SRA S3)
# Output: per-sample matrix.ec, matrix.mtx, barcodes.tsv; validation logs; intermediate FASTQ deleted post-run
# 01b: Re-quantification of validation cohort PRJNA1213680
#   CAP1 = SRR32119705, CAP2 = SRR32119704, CAP3 = SRR32119703  (per requant_run_donor_map.csv)
#   CAP1/CAP3: fetch FASTQ directly from ENA (fasterq-dump not needed)
#   CAP2:      no FASTQ on ENA -> fetch .sra from SRA S3 mirror and fasterq-dump
#              (curl rather than prefetch. Reason: S3 reachability and range resume were verified in this session.
#                recorded in self-determination ② of docs/dev_log_ja.md)
#   Each run proceeds in the order: download -> kb count -> validation -> FASTQ deletion. Retries up to 2 times per stage.
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
FQ=$W/data/fastq; SRA=$W/data/sra; CNT=$W/data/counts; R=$W/ref
L=$W/logs/01b_validation.log
mkdir -p $FQ $SRA $CNT $W/logs
export HOME=$W/.home; mkdir -p $HOME
say(){ echo "[$(date -u +%FT%TZ)] $*" | tee -a $L; }

ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq
S3=https://sra-pub-run-odp.s3.amazonaws.com/sra

# --- Download function (range resume, size verification) ---
dl(){ # $1=url $2=dest
  local url=$1 dst=$2 want got
  want=$(curl -sI --max-time 60 "$url" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r' | tail -1)
  for a in 1 2 3; do
    got=$(stat -f%z "$dst" 2>/dev/null || echo 0)
    [ -n "$want" ] && [ "$got" = "$want" ] && { say "  $(basename $dst) already complete ($got)"; return 0; }
    curl -sS -C - --retry 3 --retry-delay 15 --max-time 86400 -o "$dst" "$url"
    got=$(stat -f%z "$dst" 2>/dev/null || echo 0)
    [ -n "$want" ] && [ "$got" = "$want" ] && { say "  $(basename $dst) ok $got bytes"; return 0; }
    say "  $(basename $dst) attempt $a incomplete ($got/${want:-?})"
  done
  say "  $(basename $dst) DOWNLOAD FAILED"; return 1
}

fetch_ena(){ # $1=SRR $2=subdir(3-digit)
  local acc=$1 sub=$2
  say "$acc: ENA fastq fetch"
  dl "$ENA/SRR321/$sub/$acc/${acc}_1.fastq.gz" "$FQ/${acc}_1.fastq.gz" || return 1
  dl "$ENA/SRR321/$sub/$acc/${acc}_2.fastq.gz" "$FQ/${acc}_2.fastq.gz" || return 1
}

fetch_sra(){ # $1=SRR  -> .sra acquisition + fasterq-dump
  local acc=$1
  say "$acc: SRA object fetch (ENA に fastq 無し)"
  dl "$S3/$acc/$acc" "$SRA/$acc.sra" || return 1
  say "$acc: fasterq-dump start"
  for a in 1 2; do
    fasterq-dump --split-files --include-technical --threads 6 --temp "$FQ/tmp_$acc" -O "$FQ" "$SRA/$acc.sra" >> $L 2>&1 \
      && break
    say "$acc: fasterq-dump attempt $a failed"
  done
  rm -rf "$FQ/tmp_$acc"
  # fasterq-dump outputs uncompressed -> kb accepts either gz/plain. Existence check only.
  ls -l $FQ/${acc}_*.fastq* 2>/dev/null | awk '{print "   ",$5,$9}' | tee -a $L
}

count_run(){ # $1=SRR $2=CAPlabel
  local acc=$1 lab=$2 r1 r2 out=$CNT/$2
  if [ -s "$out/counts_unfiltered/cells_x_genes.mtx" ]; then say "$lab: already counted — skip"; return 0; fi
  # Determination of R1/R2: for 10xv3, R1=28bp(CB+UMI), R2=cDNA. Use fasterq's _1/_2 as-is.
  r1=$(ls $FQ/${acc}_1.fastq.gz $FQ/${acc}_1.fastq 2>/dev/null | head -1)
  r2=$(ls $FQ/${acc}_2.fastq.gz $FQ/${acc}_2.fastq 2>/dev/null | head -1)
  [ -s "$r1" ] && [ -s "$r2" ] || { say "$lab: FASTQ 不足 (r1=$r1 r2=$r2) — count skip"; return 1; }
  say "$lab: kb count start ($(du -ch $r1 $r2 | tail -1 | cut -f1))"
  mkdir -p $out
  for a in 1 2; do
    kb count -i $R/index.idx -g $R/t2g.txt -x 10xv3 -o $out -t 6 -m 24G --overwrite "$r1" "$r2" \
      > $W/logs/kb_$lab.log 2>&1 && break
    say "$lab: kb count attempt $a failed (see logs/kb_$lab.log)"
  done
  if [ -s "$out/counts_unfiltered/cells_x_genes.mtx" ]; then
    say "$lab: COUNT OK"
    python - "$out" "$lab" <<'EOF' 2>&1 | tee -a $L
import sys, json, os
out, lab = sys.argv[1], sys.argv[2]
mtx = os.path.join(out, "counts_unfiltered", "cells_x_genes.mtx")
with open(mtx) as f:
    f.readline(); f.readline()
    nb, ng, nnz = map(int, f.readline().split())
    tot = 0
    for line in f: tot += int(float(line.rsplit(" ",1)[1]))
rec = dict(label=lab, barcodes=nb, genes=ng, nonzero=nnz, total_umi=tot)
ri = os.path.join(out, "run_info.json")
if os.path.exists(ri):
    d = json.load(open(ri))
    rec.update(n_processed=d.get("n_processed"), n_pseudoaligned=d.get("n_pseudoaligned"),
               p_pseudoaligned=d.get("p_pseudoaligned"))
json.dump(rec, open(os.path.join(out, "verify.json"), "w"), indent=1)
print("VERIFY", json.dumps(rec))
EOF
    # Delete FASTQ and .sra after passing validation
    rm -f $FQ/${acc}_1.fastq.gz $FQ/${acc}_2.fastq.gz $FQ/${acc}_1.fastq $FQ/${acc}_2.fastq $FQ/${acc}_3.fastq $SRA/$acc.sra
    say "$lab: FASTQ/.sra 削除、空き $(df -h $W | tail -1 | awk '{print $4}')"
  else
    say "$lab: COUNT FAILED"; return 1
  fi
}

# ---- Main body: download 3 runs in parallel (sharing bandwidth), count in order of completion ----
say "=== 01b start | 空き $(df -h $W | tail -1 | awk '{print $4}') ==="
( fetch_ena SRR32119705 005 && touch $FQ/.rdy_CAP1 ) &
( fetch_ena SRR32119703 003 && touch $FQ/.rdy_CAP3 ) &
( fetch_sra SRR32119704      && touch $FQ/.rdy_CAP2 ) &
wait
say "ダウンロード段 終了: $(ls $FQ/.rdy_* 2>/dev/null | wc -l | tr -d ' ')/3 完了"

# Wait for index (01a may be running concurrently)
w=0
while [ ! -s $R/index.idx ]; do sleep 60; w=$((w+60)); [ $w -ge 14400 ] && { say "index 待機タイムアウト — count 段を中止"; exit 1; }; done
say "index 利用可"

[ -f $FQ/.rdy_CAP1 ] && count_run SRR32119705 CAP1
[ -f $FQ/.rdy_CAP3 ] && count_run SRR32119703 CAP3
[ -f $FQ/.rdy_CAP2 ] && count_run SRR32119704 CAP2
say "=== 01b DONE | counted: $(ls -d $CNT/CAP*/ 2>/dev/null | wc -l | tr -d ' ')/3 | 空き $(df -h $W | tail -1 | awk '{print $4}') ==="
