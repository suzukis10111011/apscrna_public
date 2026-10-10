#!/usr/bin/env python
"""Purpose: Recover the base composition of positions 1-26 of read 1 of the discovery
  run SRR15373534, which the saved record does not cover, and check it against the
  published positions 27-32.
Input: a FASTQ of read 1 of SRR15373534 (the first 200,000 reads; see the log for how it
  was fetched), results/tables/23_r1_base_composition.csv
Output: 08_external_review/tables/E16_r1_base_composition_recomputed.csv,
  08_external_review/logs/e16_r1_base_composition.log, and
  results/tables/23_r1_base_composition.csv extended to positions 1-32

The published columns A, C, G, T are left exactly as they were and stay empty for the
positions the first run did not cover; the recomputed values are carried in their own
columns, so no published number changes. Positions 27-32 are measured twice, on
different subsamples of the same run, and the two are compared in the log.

Fetch (recorded here so the subsample can be rebuilt):
  curl -s -r 0-69999999 -o part1.gz \\
    https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/034/SRR15373534/SRR15373534_1.fastq.gz
  gzip -dc part1.gz | head -800000 > r1_200k.fastq      # 200,000 reads
The NCBI SRA toolkit could not be used (fastq-dump reports "Failed to call external
services" in this sandbox), so the run was taken from its ENA mirror.
"""
import datetime
import os
import sys

import numpy as np
import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
FQ = sys.argv[1] if len(sys.argv) > 1 else "/tmp/sracomp/r1_200k.fastq"
OUT = f"{W}/08_external_review/tables"
LOGD = f"{W}/08_external_review/logs"
PUB = f"{W}/results/tables/23_r1_base_composition.csv"
NPOS = 32
os.makedirs(LOGD, exist_ok=True)

seqs = []
with open(FQ) as fh:
    for i, line in enumerate(fh):
        if i % 4 == 1:
            seqs.append(line.rstrip("\n")[:NPOS])
n_reads = len(seqs)
arr = np.frombuffer("".join(s.ljust(NPOS, "\0") for s in seqs).encode(),
                    dtype="S1").reshape(n_reads, NPOS)
rows = []
for p in range(NPOS):
    col = arr[:, p]
    n = int((col != b"\0").sum())
    r = dict(pos=p + 1, n_reads=n)
    for b in ("A", "C", "G", "T", "N"):
        r[f"{b}_recomputed"] = round(100 * float((col == b.encode()).sum()) / n, 2)
    rows.append(r)
NEW = pd.DataFrame(rows)
NEW["subsample"] = f"first {n_reads:,} reads of read 1, from the ENA mirror"
NEW["run"] = "SRR15373534"
NEW.to_csv(f"{OUT}/E16_r1_base_composition_recomputed.csv", index=False)

OLD = pd.read_csv(PUB)
pub_cols = [c for c in OLD.columns if c != "pos"]
EXT = NEW.merge(OLD, on="pos", how="left")
EXT = EXT[["pos"] + pub_cols + [f"{b}_recomputed" for b in "ACGTN"]
          + ["n_reads", "subsample", "run"]]
EXT.to_csv(PUB, index=False)

both = EXT.dropna(subset=["T"])
dmax = float((both["T"] - both["T_recomputed"]).abs().max())
rand = EXT[(EXT.pos >= 20) & (EXT.pos <= 26)]
lines = [
    f"[{datetime.datetime.utcnow():%Y-%m-%dT%H:%M:%SZ}] "
    f"read 1 base composition recomputed for SRR15373534",
    f"fastq: {FQ}  reads used: {n_reads:,}  positions: 1-{NPOS}",
    "fetch: curl -s -r 0-69999999 -o part1.gz "
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/034/SRR15373534/SRR15373534_1.fastq.gz"
    " ; gzip -dc part1.gz | head -800000 > r1_200k.fastq",
    "the NCBI SRA toolkit is unreachable from this sandbox "
    '(fastq-dump: "Failed to call external services"), so the ENA mirror was used',
    "",
    "positions 27-32, published (first 2,000,000 reads) against recomputed "
    f"(first {n_reads:,} reads): largest difference in T = {dmax:.2f} points",
    both[["pos", "T", "T_recomputed"]].to_string(index=False),
    "",
    "positions 20-26, the range the manuscript calls random:",
    rand[["pos", "A_recomputed", "C_recomputed", "G_recomputed",
          "T_recomputed"]].to_string(index=False),
    f"T over positions 20-26: {rand.T_recomputed.min():.2f}% to "
    f"{rand.T_recomputed.max():.2f}%  (the manuscript states 33-36%)",
]
open(f"{LOGD}/e16_r1_base_composition.log", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
