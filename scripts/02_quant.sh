#!/usr/bin/env bash
# Purpose: Quantify RNA-seq data (capacity-prioritized, run-by-run) with staged downloads and SRA conversion
# Input: none (fetches SRR IDs from ENA; resumes partial .fastq.gz, converts .sra where needed)
# Output: data/fastq/{CAP1,CAP2,CAP3}/ directories with .fastq.gz files for kb count
# 02: Quantification for the restarted version (capacity-prioritized, one run at a time). Replaces 01c.
#
# Addressing the two causes that made 01c fail:
#   (a) Because the old 01b was writing to the same data/fastq/, it misjudged the state as "both present" when only R1 was present, and
#       Skipped fasterq-dump and it failed immediately. -> Per-run subdirectory + md5 verification of both mates.
#   (b) fasterq-dump failed with 'disk-limit exeeded!'. -> Avoid the uncompressed extraction itself.
#
# Route A (CAP1=SRR32119705, CAP3=SRR32119703): fastq.gz exists in ENA.
#   Existing files are verified against ENA's md5; reuse if matching, resume fetching missing parts with curl -C -.
#   kb count reads .gz directly, so no extraction is performed.
# Route B (CAP2=SRR32119704): fastq does not exist in ENA (the corresponding field in filereport is empty).
#   Output gz directly from .sra using fastq-dump --split-files --gzip (no temporary space used).
#   Use fasterq-dump only when free space exceeds .sra size ×8 (explicitly specifying --temp,
#   Set --disk-limit / --disk-limit-tmp to 80% of free space).
#
# Mate order is _1 -> _2. The previous run record (frozen kb_info.json), in the same order,
# obtained -x 10xv3 / p_pseudoaligned 69.0%, so this order is correct.
# For kallisto/bustools, specify the full path to the conda binaries (the bundled kb binaries do not work in this environment).
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
R=$W/ref; FQ=$W/data/fq2; SRA=$W/data/sra; CNT=$W/data/counts; L=$W/logs/02_quant.log
ENA=$R/ena_filereport.tsv
export HOME=$W/.home; mkdir -p "$HOME" "$CNT" "$FQ"
export KALLISTO_BIN="$(command -v kallisto)"
export BUSTOOLS_BIN="$(command -v bustools)"
PY=${APSC_PYTHON:-python}
say(){ echo "[$(date -u +%FT%TZ)] $*" >> $L; }
free_bytes(){ df -k "$W" | tail -1 | awk '{print $4*1024}'; }
progress(){ date -u +%FT%TZ > $W/logs/quant_progress.txt; }   # For monitoring driver progress

[ -s "$ENA" ] || curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1213680&result=read_run&fields=run_accession,fastq_md5,fastq_bytes,fastq_ftp&format=tsv" -o "$ENA"

ena_field(){ awk -F'\t' -v a="$1" -v c="$2" '$1==a{print $c}' "$ENA"; }
ena_part(){  echo "$(ena_field "$1" "$2")" | cut -d';' -f"$3"; }

# ---- Fetch a single mate from ENA (resumable). Up to 3 attempts until md5 matches ----
fetch_mate(){ # $1=acc $2=mate(1|2)
  local acc=$1 m=$2 d=$FQ/$acc f url want wmd5 got
  d=$FQ/$acc; mkdir -p "$d"; f=$d/${acc}_${m}.fastq.gz
  url=https://$(ena_part "$acc" 4 "$m"); want=$(ena_part "$acc" 3 "$m"); wmd5=$(ena_part "$acc" 2 "$m")
  [ -z "$want" ] && { say "$acc mate$m: ENA に記録なし"; return 1; }
  for a in 1 2 3; do
    got=$(stat -f%z "$f" 2>/dev/null || echo 0)
    if [ "$got" = "$want" ]; then
      if [ "$(md5 -q "$f")" = "$wmd5" ]; then say "$acc mate$m: md5 一致 ($want B)"; return 0; fi
      say "$acc mate$m: サイズ一致だが md5 不一致 — 破棄して再取得"; rm -f "$f"
    fi
    say "$acc mate$m: 取得 ($got/$want B) 試行 $a"
    curl -sS -C - --retry 4 --retry-delay 15 --max-time 86400 -o "$f" "$url"
    progress
  done
  say "$acc mate$m: 取得失敗"; return 1
}

# ---- Route B: directly to gz from .sra ----
dump_gz(){ # $1=acc
  local acc=$1 d=$FQ/$acc s=$SRA/$1.final.sra sz fr
  mkdir -p "$d"
  [ -s "$d/${acc}_1.fastq.gz" ] && [ -s "$d/${acc}_2.fastq.gz" ] && { say "$acc: gz 既存"; return 0; }
  sz=$(stat -f%z "$s"); fr=$(free_bytes)
  if [ "$fr" -gt $((sz*8)) ]; then
    local lim=$((fr*8/10/1000000000))
    say "$acc: 空き $((fr/1000000000))G > .sra×8 ($((sz*8/1000000000))G) -> fasterq-dump (--disk-limit ${lim}G)"
    mkdir -p "$d/tmp"
    fasterq-dump --split-files --threads 6 --temp "$d/tmp" --disk-limit ${lim}G --disk-limit-tmp ${lim}G \
                 -O "$d" "$s" >> $L 2>&1 && { gzip -f "$d"/${acc}_*.fastq; rm -rf "$d/tmp"; return 0; }
    say "$acc: fasterq-dump 失敗 -> fastq-dump にフォールバック"; rm -rf "$d/tmp" "$d"/${acc}_*.fastq
  else
    say "$acc: 空き $((fr/1000000000))G <= .sra×8 ($((sz*8/1000000000))G) -> fastq-dump --gzip（一時領域なし）"
  fi
  fastq-dump --split-files --gzip -O "$d" "$s" >> $L 2>&1
  progress
  [ -s "$d/${acc}_1.fastq.gz" ] && [ -s "$d/${acc}_2.fastq.gz" ]
}

# ---- kb count + verification + delete FASTQ for that run ----
count_run(){ # $1=acc $2=label
  local acc=$1 lab=$2 d=$FQ/$acc out=$CNT/$lab r1=$FQ/$1/${1}_1.fastq.gz r2=$FQ/$1/${1}_2.fastq.gz
  [ -s "$out/verify.json" ] && { say "$lab: 定量済み — skip"; return 0; }
  [ -s "$r1" ] && [ -s "$r2" ] || { say "$lab: FASTQ が揃っていない"; return 1; }
  say "$lab: kb count start（R1 $(stat -f%z $r1) B / R2 $(stat -f%z $r2) B、空き $(($(free_bytes)/1000000000))G）"
  mkdir -p "$out"
  kb count --kallisto "$KALLISTO_BIN" --bustools "$BUSTOOLS_BIN" \
     -i $R/index.idx -g $R/t2g.txt -x 10xv3 -o "$out" -t 6 -m 24G --overwrite "$r1" "$r2" \
     > $W/logs/kb_$lab.log 2>&1
  progress
  [ -s "$out/counts_unfiltered/cells_x_genes.mtx" ] || { say "$lab: kb count 失敗（logs/kb_$lab.log）"; return 1; }
  $PY - "$out" "$lab" <<'PY' >> $L 2>&1
import sys, json, os
out, lab = sys.argv[1], sys.argv[2]
p = os.path.join(out, "counts_unfiltered", "cells_x_genes.mtx")
with open(p) as f:
    f.readline(); f.readline(); nb, ng, nnz = map(int, f.readline().split())
    tot = sum(int(float(l.rsplit(" ", 1)[1])) for l in f)
rec = dict(label=lab, barcodes=nb, genes=ng, nonzero=nnz, total_umi=tot)
ri = os.path.join(out, "run_info.json")
if os.path.exists(ri):
    d = json.load(open(ri)); rec.update({k: d.get(k) for k in ("n_processed","n_pseudoaligned","p_pseudoaligned","n_targets")})
json.dump(rec, open(os.path.join(out, "verify.json"), "w"), indent=1)
print("VERIFY", json.dumps(rec))
PY
  say "$lab: COUNT OK -> $(cat $out/verify.json | tr -d '\n ' | cut -c1-200)"
  rm -rf "$d"                       # Delete only this run's FASTQ (keep .sra until all runs complete)
  say "$lab: FASTQ 削除、空き $(($(free_bytes)/1000000000))G"
  progress
}

say "=== 02 start 空き $(($(free_bytes)/1000000000))G ==="
progress
# Do not touch the old data/fastq/ since an old curl that cannot be stopped keeps writing to it.
# Isolate output to data/fq2/<acc>/ (the corruption in 01c was caused by concurrent writes to the same path).

# Order: CAP1 -> CAP3 -> CAP2 (ENA route first, CAP2 which consumes capacity last)
for pair in "SRR32119705 CAP1" "SRR32119703 CAP3" "SRR32119704 CAP2"; do
  set -- $pair; acc=$1; lab=$2
  [ -s "$CNT/$lab/verify.json" ] && { say "$lab: 完了済み — skip"; continue; }
  # Route selection (self-judgment 9, changed at 2026-09-19 09:45Z):
  #   Initially instructed to treat CAP1/CAP3 as route A (reuse ENA's fastq.gz and resume-fetch the missing portion),
  #   but that premise (that the local _1 is complete) broke down. _1 matched ENA byte-for-byte yet had an md5 mismatch and was corrupted
  #   (caused by concurrent writes from the old 01b that could not be stopped), requiring a full re-fetch of all data.
  #   Measured ENA throughput is 0.32 MB/s; _1 alone for CAP1 takes 29.3 hours. Meanwhile, all 3 runs of .sra
  #   are already on hand, and vdb-validate returns md5 ok (internal checksum) for CAP1.
  #   Therefore, as long as .sra is available, use route B (generate locally). ENA is demoted to a fallback for when .sra is unavailable.
  if [ -s "$SRA/$acc.final.sra" ]; then
    say "$lab ($acc): 経路B = ローカル .sra から gz 生成（ENA は 0.32 MB/s のため不使用）"
    dump_gz "$acc" || { say "$lab: gz 生成失敗 — 次へ"; continue; }
  elif [ -n "$(ena_part "$acc" 4 1)" ]; then
    say "$lab ($acc): 経路A = ENA fastq.gz（.sra が無いため）"
    fetch_mate "$acc" 1 && fetch_mate "$acc" 2 || { say "$lab: 取得失敗 — 次へ"; continue; }
  else
    say "$lab: .sra も ENA fastq も無い — 次へ"; continue
  fi
  count_run "$acc" "$lab" || say "$lab: 定量失敗 — 次へ"
done
say "=== 02 DONE counted=$(ls $CNT/*/verify.json 2>/dev/null | wc -l | tr -d ' ')/3 空き $(($(free_bytes)/1000000000))G ==="
progress
