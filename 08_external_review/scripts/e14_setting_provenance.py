#!/usr/bin/env python
"""Purpose: Check the 25 worksheet rows the student marked as absent from the tables -
  software versions, analysis settings and read-structure figures - against the code,
  the run logs and the files under results/, and record where each one lives.
Input: 07_submission_bundle/tables/manuscript_number_agreement.csv, and the scripts,
  logs, ref/ records and results/tables files of the repository
Output: 08_external_review/tables/E14_setting_provenance.csv

Every value is read out of the named file at run time - a log line, a cell of a result
table, a field of a JSON record, or the version of the installed distribution - so the
recorded value is what the file says, not what the manuscript says. A row whose value is
in no saved file is recorded as such instead of being asserted.
"""
import importlib.metadata as md
import json
import os
import re
import subprocess

import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
BUN = f"{W}/07_submission_bundle/tables"
OUT = f"{W}/08_external_review/tables"

A = pd.read_csv(f"{BUN}/manuscript_number_agreement.csv")
A = A[A.version == sorted(v for v in A.version.unique()
                          if v.startswith("clean_"))[-1]].reset_index(drop=True)


def grep(path, pattern):
    """(file:line, line text) of the first line of `path` matching `pattern`"""
    full = f"{W}/{path}"
    if not os.path.exists(full):
        return "", f"{path} not found"
    with open(full, errors="replace") as fh:
        for i, line in enumerate(fh, start=1):
            if re.search(pattern, line):
                return f"{path}:{i}", line.strip()[:260]
    return "", f"no line of {path} matches {pattern}"


def cell(path, query, col):
    d = pd.read_csv(f"{W}/{path}")
    v = d.query(query)[col]
    return path, ("" if v.empty else v.iloc[0])


def jsonf(path, *keys):
    d = json.load(open(f"{W}/{path}"))
    for k in keys:
        d = d[k]
    return path, d


def dist(name):
    try:
        return f"the environment apsc", md.version(name)
    except md.PackageNotFoundError:
        return "the environment apsc", ""


CHEM = "results/tables/23_chemistry_umi_comparison.csv"
BASE = "results/tables/23_r1_base_composition.csv"
DEC = "ref/chemistry_decision.json"
PRE = "scripts/03_preprocess_discovery.py"
VAL = "scripts/05_validation_and_projection.py"

# No -> (what the manuscript states, kind, spec, category)
ITEMS = {
    16: ("the kit is the 10x Single Cell 3' Reagent Kit v2", "json",
         (DEC, "chemistry"), "read structure / kit"),
    20: ("Cell Ranger 3.1.0 for lesions B-E", "external", None, "software version"),
    39: ("the comparison run used 10x v3 settings", "grep",
         (CHEM, r"^10xv3"), "read structure / kit"),
    43: ("kallisto index k = 31", "grep",
         ("logs/01a_mkidx.log", r"k-mer length"), "pipeline argument"),
    44: ("Ensembl release 112 GRCh38 cDNA", "grep",
         ("logs/01a_mkidx.log", r"release-112"), "reference build"),
    48: ("33% of the counts of MYH11 fall on the secondary identifier", "csv",
         (CHEM.replace("23_chemistry_umi_comparison",
                       "20_duplicate_symbol_umi_split"),
          "symbol == 'MYH11'", "minor_frac"), "measured value"),
    49: ("the symbol table was supplemented for 122 genes", "grep",
         ("docs/private/dev_log_ja.md", r"GEO features で 122 件補完"), "measured value"),
    55: ("435,676,976 reads, 69.6% pseudoaligned for SRR32119705", "grep",
         ("logs/02b_quant.log", r'"p_pseudoaligned": 69\.6'), "log value"),
    57: ("the discovery runs were quantified with kb count -x 10xv2", "grep",
         ("logs/08_requant.log", r"kb count -x 10xv2"), "pipeline argument"),
    58: ("the base composition of read 1", "grep",
         (DEC, r"位置27-28"), "read structure / kit"),
    60: ("position 28 of read 1", "csv", (BASE, "pos == 28", "T"),
         "read structure / kit"),
    61: ("position 27 is thymine in 97.9% of reads", "csv", (BASE, "pos == 27", "T"),
         "measured value"),
    62: ("position 28 is thymine in 97.8% of reads", "csv", (BASE, "pos == 28", "T"),
         "measured value"),
    64: ("the 16 + 10 nt barcode-UMI structure of v2", "json",
         (DEC, "decided_by"), "read structure / kit"),
    66: ("positions 20-26 were random", "absent", None, "measured value"),
    67: ("T 33-36% over positions 20-26", "absent", None, "measured value"),
    68: ("T 33-36% over positions 20-26", "absent", None, "measured value"),
    69: ("1.15 million deduplicated UMIs under v3", "csv",
         (CHEM, "chemistry == '10xv3'", "total_umi_dedup"), "measured value"),
    70: ("median 16 UMIs per barcode among the top 5,000 under v3", "csv",
         (CHEM, "chemistry == '10xv3'", "median_umi_top5k"), "measured value"),
    71: ("9.06 million deduplicated UMIs under v2", "csv",
         (CHEM, "chemistry == '10xv2'", "total_umi_dedup"), "measured value"),
    87: ("counts normalized to 10,000 UMIs per cell", "grep",
         (PRE, r"normalize_total"), "pipeline argument"),
    97: ("Leiden 0.2-1.5 in steps of 0.1", "grep",
         (VAL, r"np\.arange\(0\.2, 1\.51, 0\.1\)"), "pipeline argument"),
    116: ("score_genes with a control set of 50 genes", "grep",
          (VAL, r"ctrl_size=50"), "pipeline argument"),
    120: ("omnipath 1.0.12", "pkg", "omnipath", "software version"),
    205: ("decoupler 2.2.0", "pkg", "decoupler", "software version"),
}

rows = []
for no, (claim, kind, spec, cat) in sorted(ITEMS.items()):
    printed = str(A.iloc[no - 1].printed)
    if kind == "grep":
        where, value = grep(*spec)
    elif kind == "csv":
        path, v = cell(*spec)
        where, value = f"{path} ({spec[1]})", v
    elif kind == "json":
        path, v = jsonf(spec[0], spec[1])
        where, value = f"{path} ({spec[1]})", str(v)[:160]
    elif kind == "pkg":
        where, value = dist(spec)
    elif kind == "external":
        where, value = ("GEO series metadata of GSE181688 and GSE197680 "
                        "(the depositors' processing, not run here)", "")
    else:
        where, value = ("", "")
    ok = ""
    if kind == "absent":
        ok = "not in a saved file"
        where = (f"{BASE} and {DEC} record positions 27-32 only; the subsampled fastq "
                 "of SRR15373534 is no longer on disk")
    elif kind == "external":
        ok = "not verifiable in this repository"
    elif not where:
        ok = "not found in the file named"
    elif cat in ("read structure / kit",):
        ok = "recorded in the named file"
    else:
        t = printed.replace(",", "").replace("%", "").strip()
        try:
            pv = float(t)
            nums = [abs(float(x)) for x in
                    re.findall(r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", str(value))]
            dp = len(t.split(".")[1]) if "." in t else 0
            pv = abs(pv)
            ok = ("confirmed" if any(round(c, dp) == round(pv, dp)
                                     for b in nums
                                     for c in (b, b * 100, b / 100, b / 1e6, b * 1e6))
                  else "differs from the file")
        except ValueError:
            ok = "confirmed" if t in str(value) else "differs from the file"
    rows.append(dict(No=no, printed=printed, category=cat, claim=claim,
                     value_in_the_file=value, where_it_lives=where, status=ok))
E14 = pd.DataFrame(rows)
E14.to_csv(f"{OUT}/E14_setting_provenance.csv", index=False)
pd.set_option("display.max_colwidth", 58)
print(E14[["No", "printed", "category", "value_in_the_file", "where_it_lives",
           "status"]].to_string(index=False))
print("\n", E14.status.value_counts().to_dict())
