#!/usr/bin/env bash
# Purpose: Build kallisto index from Ensembl release-112 human cDNA and map transcript IDs to gene symbols
# Input: none (fetches Ensembl GRCh38 release-112 cdna.all.fa.gz)
# Output: kallisto index and transcript-to-gene mapping table (tid, gid, symbol)
# 01a: Rebuild kallisto index. Same conditions as frozen/scripts/mkidx.sh.
#   Reference = Ensembl release-112 Homo_sapiens.GRCh38.cdna.all.fa.gz, k=31 (kallisto default)
#   t2g = transcript -> gene -> symbol from cDNA headers
set -uo pipefail
W="${APSC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
R=$W/ref; L=$W/logs/01a_mkidx.log
mkdir -p $R
say(){ echo "[$(date -u +%FT%TZ)] $*" | tee -a $L; }

if [ -s $R/index.idx ] && [ -s $R/t2g.txt ]; then say "index already present: $(du -h $R/index.idx | cut -f1)"; exit 0; fi

CD=$R/cdna.fa.gz
if [ ! -s $CD ]; then
  say "cDNA fetch start (Ensembl release-112)"
  for a in 1 2 3; do
    curl -sS --retry 4 --retry-delay 10 --max-time 3600 -o $CD \
      "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" && break
    say "cDNA fetch attempt $a failed"; sleep 15
  done
  [ -s $CD ] || { say "cDNA FETCH FAILED"; exit 1; }
fi
say "cDNA ok $(du -h $CD | cut -f1)"

say "build t2g"
python - <<'EOF' 2>&1 | tee -a $L
import gzip, re, os
R = f"{W}/ref"
n = 0
with gzip.open(f"{R}/cdna.fa.gz", "rt") as f, open(f"{R}/t2g.txt", "w") as o:
    for line in f:
        if not line.startswith(">"): continue
        h = line[1:].strip()
        tid = h.split()[0]
        g = re.search(r"gene:(\S+)", h)
        s = re.search(r"gene_symbol:(\S+)", h)
        gid = g.group(1) if g else tid
        sym = s.group(1) if s else gid
        o.write(f"{tid}\t{gid}\t{sym}\n"); n += 1
print("t2g transcripts:", n)
EOF
[ -s $R/t2g.txt ] || { say "t2g FAILED"; exit 1; }

say "kallisto index start (t=6)"
kallisto index -i $R/index.new.idx -t 6 $CD >> $L 2>&1
if [ -s $R/index.new.idx ]; then
  mv $R/index.new.idx $R/index.idx
  say "index built $(python3 -c "import os;print('%.2f GB'%(os.path.getsize('$R/index.idx')/1e9))")"
else
  say "INDEX BUILD FAILED"; exit 1
fi
say "01a DONE"
