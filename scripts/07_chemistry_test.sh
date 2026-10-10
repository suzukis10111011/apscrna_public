#!/usr/bin/env bash
# Purpose: Decide the 10x chemistry for the discovery samples by running kb count under both v2 and v3 on the
#   first 20M reads of one run.
# Input: one discovery SRA run (GSE181688 = PRJNA753049, GSE197680 = PRJNA811571); ref/index.idx, ref/t2g.txt.
# Output: ref/chemistry_decision.json (verdict with pseudoalignment rate and whitelist-matching barcode
#   counts); logs/07_chemistry.log.
# 07: 10x chemistry determination for the discovery side. Using the first 20 million reads of one run, try both 10xv2 and 10xv3,
# and decide based on the pseudoalignment rate and the number of barcodes matching the whitelist.
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
R=$W/ref; T=$W/data/chem; L=$W/logs/07_chemistry.log
export HOME=$W/.home; mkdir -p "$HOME" "$T"
KA=$(command -v kallisto); BU=$(command -v bustools)
say(){ echo "[$(date -u +%FT%TZ)] $*" >> $L; }
say "=== 07 start 空き $(df -k $W | tail -1 | awk '{print int($4/1048576)}')G ==="

# List of runs in the discovery cohort (GSE181688=PRJNA753049, GSE197680=PRJNA811571)
for P in PRJNA753049 PRJNA811571; do
  curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$P&result=read_run&fields=run_accession,sample_title,read_count,fastq_bytes&format=tsv" \
    | tail -n +2 >> $R/discovery_runs.tsv
done
sort -u $R/discovery_runs.tsv -o $R/discovery_runs.tsv
say "探索ラン数: $(wc -l < $R/discovery_runs.tsv | tr -d ' ')"
ACC=$(head -1 $R/discovery_runs.tsv | cut -f1)
say "chemistry 判定に使うラン: $ACC"

if [ ! -s $T/$ACC.sra ]; then
  say "$ACC: .sra 取得開始"
  curl -sS --retry 4 --retry-delay 15 --max-time 86400 -o $T/$ACC.sra \
    "https://sra-pub-run-odp.s3.amazonaws.com/sra/$ACC/$ACC" || { say "$ACC: 取得失敗"; exit 1; }
  say "$ACC: 取得完了 $(stat -f%z $T/$ACC.sra) B"
fi
if [ ! -s $T/${ACC}_1.fastq.gz ]; then
  say "$ACC: 先頭2,000万スポットを dump"
  fastq-dump --split-files --gzip -X 20000000 -O $T $T/$ACC.sra >> $L 2>&1
fi
R1=$T/${ACC}_1.fastq.gz; R2=$T/${ACC}_2.fastq.gz
[ -s "$R1" ] && [ -s "$R2" ] || { say "dump 失敗"; exit 1; }
say "R1 長 $(gzcat $R1 | head -2 | tail -1 | wc -c | tr -d ' ') / R2 長 $(gzcat $R2 | head -2 | tail -1 | wc -c | tr -d ' ')"

for X in 10xv2 10xv3; do
  o=$T/out_$X
  [ -s "$o/counts_unfiltered/cells_x_genes.mtx" ] && { say "$X: 既存 — skip"; continue; }
  say "$X: kb count"
  kb count --kallisto "$KA" --bustools "$BU" -i $R/index.idx -g $R/t2g.txt -x $X \
     -o "$o" -t 6 -m 24G --overwrite "$R1" "$R2" > $W/logs/kb_chem_$X.log 2>&1
  if [ -s "$o/run_info.json" ]; then
    say "$X: $(python3 -c "
import json;d=json.load(open('$o/run_info.json'))
bc=sum(1 for _ in open('$o/counts_unfiltered/cells_x_genes.barcodes.txt'))
print(f\"p_pseudoaligned={d.get('p_pseudoaligned')} n_processed={d.get('n_processed')} barcodes={bc}\")")"
  else
    say "$X: 失敗（logs/kb_chem_$X.log）"
  fi
done
say "=== 07 DONE ==="
