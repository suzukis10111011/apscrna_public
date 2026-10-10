#!/usr/bin/env python
"""Chemistry check for the three validation runs (SRR32119703-05).

The source publication (Zhang et al., Front Immunol 2025, doi:10.3389/fimmu.2025.1559220)
states 10x Chromium Single Cell 3' v2. Every kb count call in this repository used
-x 10xv3 (scripts/01b_validation_pipeline.sh:71), and the kb logs confirm the
10x_version3 on-list was applied to CAP1-CAP3. ref/chemistry_decision.json records a
10xv2 verdict, but scripts/07_chemistry_test.sh ran that test on a DISCOVERY run
(the first row of ref/discovery_runs.tsv), never on a validation run.

This script takes the first 200,000 reads of each validation run and reports:
  read 1 length
  base composition at positions 27 and 28 (1-based)
  the fraction of the top 5,000 barcodes (first 16 bases of read 1, ranked by read
  count) present in the 10x v2 (737K-august-2016) and v3 (3M-february-2018) on-lists,
  and the same fraction over all reads

Discriminator: with v2 chemistry the cell barcode is 16 bases and the UMI 10, so the
poly-T of the capture oligo begins at position 27; with v3 the UMI is 12 and positions
27-28 are still UMI. The published discovery-run test used exactly this rule.

Read source. SRR32119703 and SRR32119705 have ENA fastq, which is streamed and cut at
200,000 reads. SRR32119704 has no ENA fastq (read_count 0 in the ENA file report), so
its .sra object is fetched from the NCBI open-data mirror and dumped with
fastq-dump -X 200000, the command named in the request. The route is recorded per run.

Output: 06_submission_prep/tables/C6_validation_chemistry_check.csv
        06_submission_prep/tables/C6_base_composition_by_position.csv
"""
import gzip
import os
import subprocess
import sys
import time
from collections import Counter

import pandas as pd

W = "/path/to/Documents/apscrna_v2"
OUT = f"{W}/06_submission_prep/tables"
WORK = f"{W}/data/chem2"
N_READS = 200_000
TOPN = 5_000
WL = {"v2_737K": f"{W}/data/counts_discovery/SampleA/10x_version2_whitelist.txt",
      "v3_3M": f"{W}/data/counts/CAP1/10x_version3_whitelist.txt"}
ENA = {"SRR32119703": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR321/003/"
                      "SRR32119703/SRR32119703_1.fastq.gz",
       "SRR32119705": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR321/005/"
                      "SRR32119705/SRR32119705_1.fastq.gz"}
DONOR = {"SRR32119703": "CAP1", "SRR32119704": "CAP2", "SRR32119705": "CAP3"}
NOTE = ("pre-submission check; no frozen file, gene set, candidate factor or threshold "
        "is involved")


def log(m):
    s = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {m}"
    print(s, flush=True)
    open(f"{W}/06_submission_prep/logs/run_submission_prep.log", "a").write(s + "\n")


def reads_from_ena(url, n):
    """Stream the gzip fastq and stop after n reads."""
    p = subprocess.Popen(["curl", "-sS", "--retry", "4", url], stdout=subprocess.PIPE)
    out = []
    try:
        with gzip.open(p.stdout, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    out.append(line.rstrip("\n"))
                    if len(out) >= n:
                        break
    finally:
        p.stdout.close()
        p.terminate()
    return out


def reads_from_sra(acc, n):
    """fastq-dump -X n on the locally fetched .sra, exactly as the request names."""
    sra = f"{WORK}/{acc}.sra"
    fq = f"{WORK}/{acc}_1.fastq"
    if not os.path.exists(fq):
        env = dict(os.environ, HOME=f"{W}/.home")
        r = subprocess.run(["fastq-dump", "-X", str(n), "--split-files", "-O", WORK,
                            sra], capture_output=True, text=True, env=env)
        log(f"{acc}: fastq-dump rc={r.returncode} {r.stderr.strip()[:200]}")
    out = []
    with open(fq) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out.append(line.rstrip("\n"))
                if len(out) >= n:
                    break
    return out


wl = {}
for k, p in WL.items():
    with open(p) as fh:
        wl[k] = {ln.strip() for ln in fh if ln.strip()}
    log(f"on-list {k}: {len(wl[k])} barcodes")

rows, comp_rows = [], []
for acc in ["SRR32119703", "SRR32119704", "SRR32119705"]:
    if acc in ENA:
        seqs = reads_from_ena(ENA[acc], N_READS)
        route = ("ENA fastq, first 200,000 reads "
                 "(ftp.sra.ebi.ac.uk/vol1/fastq/...)")
    else:
        EXPECTED = {"SRR32119704": 40563234654}
        sra = f"{WORK}/{acc}.sra"
        if not os.path.exists(sra) or os.path.getsize(sra) < EXPECTED.get(acc, 1):
            have = os.path.getsize(sra) if os.path.exists(sra) else 0
            log(f"{acc}: .sra incomplete ({have} of {EXPECTED.get(acc)} bytes), skipped")
            continue
        seqs = reads_from_sra(acc, N_READS)
        route = "NCBI open-data .sra, fastq-dump -X 200000 --split-files"
    lens = Counter(len(s) for s in seqs)
    L = lens.most_common(1)[0][0]
    log(f"{acc}: {len(seqs)} reads, read-1 length {dict(lens)}")

    for pos in range(25, 33):
        if L < pos:
            continue
        c = Counter(s[pos - 1] for s in seqs if len(s) >= pos)
        tot = sum(c.values())
        comp_rows.append(dict(run_accession=acc, donor_label=DONOR[acc], position=pos,
                              **{b: round(100 * c.get(b, 0) / tot, 2)
                                 for b in "ACGTN"},
                              n_reads=tot, note=NOTE))

    bc = Counter(s[:16] for s in seqs if len(s) >= 16)
    top = [b for b, _ in bc.most_common(TOPN)]
    r = dict(run_accession=acc, donor_label=DONOR[acc], read_source=route,
             n_reads_examined=len(seqs), read1_length=L,
             read1_length_counts="; ".join(f"{k}bp x{v}" for k, v in lens.most_common(3)),
             n_distinct_barcodes=len(bc), n_top_barcodes=len(top), note=NOTE)
    for k in WL:
        r[f"top{TOPN}_in_{k}_pct"] = round(100 * sum(b in wl[k] for b in top) / len(top), 2)
        r[f"all_reads_in_{k}_pct"] = round(
            100 * sum(n for b, n in bc.items() if b in wl[k]) / sum(bc.values()), 2)
    for pos in (27, 28):
        if L >= pos:
            c = Counter(s[pos - 1] for s in seqs if len(s) >= pos)
            tot = sum(c.values())
            r[f"pos{pos}_T_pct"] = round(100 * c.get("T", 0) / tot, 2)
    r["expected_if_v2"] = ("positions 27 and 28 are the start of the capture-oligo "
                           "poly-T, so both are ~98% T; barcodes match the 737K "
                           "on-list")
    r["expected_if_v3"] = ("positions 27 and 28 are still inside the 12-base UMI, so "
                           "no T enrichment; barcodes match the 3M on-list")
    rows.append(r)

T = pd.DataFrame(rows)
T.to_csv(f"{OUT}/C6_validation_chemistry_check.csv", index=False)
pd.DataFrame(comp_rows).to_csv(f"{OUT}/C6_base_composition_by_position.csv", index=False)
cols = ["run_accession", "donor_label", "read1_length", "pos27_T_pct", "pos28_T_pct",
        "top5000_in_v2_737K_pct", "top5000_in_v3_3M_pct",
        "all_reads_in_v2_737K_pct", "all_reads_in_v3_3M_pct"]
print(T[[c for c in cols if c in T.columns]].to_string(index=False))
log("r5 done")
