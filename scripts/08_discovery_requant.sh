#!/usr/bin/env bash
# Purpose: Re-quantify 5 discovery lesions across 20 runs (5 samples × 4 lanes); aggregate lane FASTQs per
#   sample for unified kb count
# Input: ref/discovery_runs.tsv, ref/chemistry_decision.json, SRA archives
# Output: kb_info.json, abundance matrices per sample; FASTQ and SRA cleaned up; stops when free disk < 60 GB
# 08: Re-quantification of the 5 discovery lesions. 20 runs in ref/discovery_runs.tsv (5 samples × 4 runs).
#
# Self-judgment ⑫: acquisition and dump are performed **per run**, kb count is performed **per sample** (FASTQ for 4 runs is
#   passed to a single kb count call). Since the 20 runs are lane splits of the same library, counting per run would
#   cause the same cell barcode and same UMI to be counted as separate molecules, breaking UMI deduplication.
#   The previous validation-side requant also passed the 44 splits to a single kb count call (the call in kb_info.json).
# chemistry uses the value from ref/chemistry_decision.json (determined by the measurement in 07).
# After each sample completes, rm the FASTQ; after each run completes, rm the .sra. Stop once df falls below 60 GB.
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
R=$W/ref; SRA=$W/data/dsra; FQ=$W/data/dfq; CNT=$W/data/counts_discovery
L=$W/logs/08_requant.log; LOCK=$W/logs/08.lock
export HOME=$W/.home; mkdir -p "$HOME" "$SRA" "$FQ" "$CNT"
KA=$(command -v kallisto); BU=$(command -v bustools)
PY=${APSC_PYTHON:-python}
S3=https://sra-pub-run-odp.s3.amazonaws.com/sra
MINFREE_G=60
say(){ echo "[$(date -u +%FT%TZ)] $*" >> $L; }
beat(){ date -u +%FT%TZ > $W/logs/requant_heartbeat.txt; }
free_g(){ df -k "$W" | tail -1 | awk '{print int($4/1048576)}'; }
sz(){ stat -f%z "$1" 2>/dev/null || echo 0; }

mkdir "$LOCK" 2>/dev/null || { say "別インスタンスが稼働中（$LOCK）— 起動しない"; exit 0; }
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

CHEM=$(python3 -c "import json;print(json.load(open('$R/chemistry_decision.json'))['chemistry'])" 2>/dev/null)
[ -z "$CHEM" ] && { say "chemistry 未決定（ref/chemistry_decision.json 無し）— 中止"; exit 1; }
say "=== 08 start chemistry=$CHEM 空き $(free_g)G ==="; beat

fetch_run(){ # $1=acc  8-way parallel range retrieval
  local acc=$1 f=$SRA/$acc.sra want got
  [ -f "$SRA/$acc.ok" ] && return 0
  want=$(curl -sI "$S3/$acc/$acc" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r' | tail -1)
  [ -z "$want" ] && { say "$acc: サイズ取得失敗"; return 1; }
  local part=$((want/8+1))
  for try in 1 2 3; do
    for i in $(seq 0 7); do
      local st=$((i*part)) en=$(( (i+1)*part - 1 )); [ $en -ge $want ] && en=$((want-1))
      ( g=$(sz "$SRA/$acc.p$i"); w=$((en-st+1)); [ "$g" -ge "$w" ] && exit 0
        curl -sS -r $((st+g))-${en} --retry 3 --retry-delay 10 --max-time 86400 "$S3/$acc/$acc" >> "$SRA/$acc.p$i" ) &
    done
    wait; beat
    local tot=0; for i in $(seq 0 7); do tot=$((tot + $(sz "$SRA/$acc.p$i"))); done
    [ "$tot" = "$want" ] && break
    say "$acc: 試行 $try 未完 ($tot/$want)"
  done
  cat $SRA/$acc.p0 $SRA/$acc.p1 $SRA/$acc.p2 $SRA/$acc.p3 $SRA/$acc.p4 $SRA/$acc.p5 $SRA/$acc.p6 $SRA/$acc.p7 > "$f"
  rm -f $SRA/$acc.p*
  got=$(sz "$f")
  [ "$got" = "$want" ] || { say "$acc: 結合後サイズ不一致 ($got/$want)"; return 1; }
  touch "$SRA/$acc.ok"; say "$acc: 取得 $got B"; beat
}

dump_run(){ # $1=acc -> gz into $FQ/<sample>/, .sra is deleted
  local acc=$1 d=$2
  [ -s "$d/${acc}_1.fastq.gz" ] && [ -s "$d/${acc}_2.fastq.gz" ] && return 0
  say "$acc: fastq-dump --gzip"
  fastq-dump --split-files --gzip -O "$d" "$SRA/$acc.sra" >> $L 2>&1
  for m in 1 2; do
    [ -s "$d/${acc}_$m.fastq.gz" ] || { say "$acc: mate$m 無し"; return 1; }
  done
  rm -f "$SRA/$acc.sra"; beat
}

for SAMP in A B C D E; do
  out=$CNT/Sample$SAMP
  [ -s "$out/verify.json" ] && { say "Sample$SAMP: 定量済み — skip"; continue; }
  ACCS=$(awk -F'\t' -v s="Sample $SAMP" '$2==s{print $1}' $R/discovery_runs.tsv)
  [ -z "$ACCS" ] && { say "Sample$SAMP: ラン不明 — skip"; continue; }
  d=$FQ/Sample$SAMP; mkdir -p "$d"
  ok=1
  for acc in $ACCS; do
    if [ "$(free_g)" -lt $MINFREE_G ]; then say "空き $(free_g)G < ${MINFREE_G}G — 停止"; exit 2; fi
    fetch_run "$acc" && dump_run "$acc" "$d" || { say "Sample$SAMP/$acc: 失敗 — このサンプルを飛ばす"; ok=0; break; }
  done
  [ $ok -eq 1 ] || continue
  FILES=""
  for acc in $ACCS; do FILES="$FILES $d/${acc}_1.fastq.gz $d/${acc}_2.fastq.gz"; done
  say "Sample$SAMP: kb count -x $CHEM（4ラン分を1回に渡す）空き $(free_g)G"
  kb count --kallisto "$KA" --bustools "$BU" -i $R/index.idx -g $R/t2g.txt -x $CHEM \
     -o "$out" -t 6 -m 24G --overwrite $FILES > $W/logs/kb_Sample$SAMP.log 2>&1
  beat
  if [ -s "$out/counts_unfiltered/cells_x_genes.mtx" ]; then
    $PY - "$out" "Sample$SAMP" "$CHEM" <<'PY' >> $L 2>&1
import sys, json, os
out, lab, chem = sys.argv[1], sys.argv[2], sys.argv[3]
p = os.path.join(out, "counts_unfiltered", "cells_x_genes.mtx")
with open(p) as f:
    f.readline(); f.readline(); nb, ng, nnz = map(int, f.readline().split())
    tot = sum(int(float(l.rsplit(" ", 1)[1])) for l in f)
rec = dict(label=lab, chemistry=chem, barcodes=nb, genes=ng, nonzero=nnz, total_umi=tot)
ri = os.path.join(out, "run_info.json")
if os.path.exists(ri):
    d = json.load(open(ri)); rec.update({k: d.get(k) for k in ("n_processed","n_pseudoaligned","p_pseudoaligned","n_targets")})
json.dump(rec, open(os.path.join(out, "verify.json"), "w"), indent=1)
print("VERIFY", json.dumps(rec))
PY
    say "Sample$SAMP: COUNT OK $(tr -d '\n ' < $out/verify.json | cut -c1-200)"
    rm -rf "$d"
    say "Sample$SAMP: FASTQ 解放、空き $(free_g)G"
  else
    say "Sample$SAMP: kb count 失敗（logs/kb_Sample$SAMP.log）"
  fi
  beat
done
say "=== 08 DONE counted=$(ls $CNT/*/verify.json 2>/dev/null | wc -l | tr -d ' ')/5 空き $(free_g)G ==="
beat
