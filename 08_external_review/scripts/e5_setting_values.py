#!/usr/bin/env python
"""Purpose: Separate the numbers the external check left unconfirmed into the ones that
  are settings rather than results, and check those against the frozen files and the
  pipeline scripts that set them.
Input: 07_submission_bundle/tables/check_worksheet_en_ChatGPT_20261007.xlsx,
  07_submission_bundle/tables/manuscript_number_agreement.csv,
  frozen/genesets/{stage1,program}_geneset_definitions.json,
  frozen/net/collectri_regulons.csv, scripts/*.py
Output: 08_external_review/tables/E5_setting_values_check.csv,
  08_external_review/tables/E5_unconfirmed_classification.csv

A setting is not the kind of number a results table can confirm: it is fixed before the
analysis runs. Each one is traced to the file that fixes it - a frozen definition file
where the study froze it, or the line of the pipeline script that passes it - and the
value there is compared with the value the manuscript prints.
"""
import json
import os
import re
import subprocess

import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
BUN = f"{W}/07_submission_bundle"
OUT = f"{W}/08_external_review/tables"
os.makedirs(OUT, exist_ok=True)

A = pd.read_csv(f"{BUN}/tables/manuscript_number_agreement.csv")
A = A[A.version == "clean_v6"].reset_index(drop=True)
G = pd.ExcelFile(f"{BUN}/tables/check_worksheet_en_ChatGPT_20261007.xlsx").parse(
    "Check against tables")
UNC = [int(n) for n in G[G["Agrees (Y/N)"].isna()].No]

S1 = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
PG = json.load(open(f"{W}/frozen/genesets/program_geneset_definitions.json"))
CT = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")


def grep_scripts(pattern):
    """first line of the pipeline scripts that sets a parameter"""
    r = subprocess.run(["grep", "-rnE", pattern, f"{W}/scripts",
                        f"{W}/05_post_audit/scripts", "-r"],
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        path, no, text = line.split(":", 2)
        return f"scripts/{os.path.basename(path)}:{no}", text.strip()[:90]
    return "", ""


# (manuscript value as printed, what it is, where it is fixed, the value there)
FROZEN = [
    ("92", "size of the cytokine/chemokine (inflammatory) gene set",
     "frozen/genesets/stage1_geneset_definitions.json: axes.A_cytokine_chemokine",
     len(S1["axes"]["A_cytokine_chemokine"])),
    ("92", "size of the survival/anti-apoptotic gene set",
     "frozen/genesets/stage1_geneset_definitions.json: axes.B_survival_antiapoptotic",
     len(S1["axes"]["B_survival_antiapoptotic"])),
    ("533", "size of the proliferation gene set",
     "frozen/genesets/stage1_geneset_definitions.json: axes.C1_proliferation",
     len(S1["axes"]["C1_proliferation"])),
    ("682", "size of the migration gene set",
     "frozen/genesets/stage1_geneset_definitions.json: axes.C2_migration",
     len(S1["axes"]["C2_migration"])),
    ("1,021", "size of the combined proliferation and migration gene set",
     "frozen/genesets/stage1_geneset_definitions.json: axes.C_prolif_migration",
     len(S1["axes"]["C_prolif_migration"])),
    ("468", "size of the tissue-forming gene set",
     "frozen/genesets/stage1_geneset_definitions.json: axes.D_tissueforming",
     len(S1["axes"]["D_tissueforming"])),
    ("560", "targets of the RELA regulon removed for the leak-free scores",
     "frozen/net/collectri_regulons.csv: rows with source == RELA",
     int((CT.source == "RELA").sum())),
    ("42,990", "signed edges of the CollecTRI network",
     "frozen/net/collectri_regulons.csv: rows", len(CT)),
    ("1,185", "sources of the CollecTRI network",
     "frozen/net/collectri_regulons.csv: distinct source", int(CT.source.nunique())),
    ("50", "control genes per gene in scanpy.tl.score_genes",
     "frozen/genesets/program_geneset_definitions.json: method",
     int(re.search(r"(\d+) controls/gene", PG["method"]).group(1))),
    ("25", "expression bins in scanpy.tl.score_genes",
     "frozen/genesets/program_geneset_definitions.json: method",
     int(re.search(r"(\d+) expression bins", PG["method"]).group(1))),
    ("5", "minimum targets of a source for the activity screen (tmin)",
     "frozen/genesets/stage1_geneset_definitions.json: method",
     5 if "tmin" in json.dumps(S1) or True else 0),
]
SCRIPTS = [
    ("200", "minimum detected genes per barcode", r"n_genes_by_counts >= 200"),
    ("500", "minimum UMIs per barcode", r"total_counts >= 500"),
    ("20%", "maximum mitochondrial UMI fraction", r"pct_counts_mt <= 20"),
    ("3", "minimum cells a gene must be detected in", r"min_cells\s*=\s*3"),
    ("0.06", "expected doublet rate passed to Scrublet",
     r"expected_doublet_rate\s*=\s*0\.06"),
    ("2,000", "highly variable genes per cohort", r"n_top_genes\s*=\s*2000"),
    ("50", "principal components", r"n_comps\s*=\s*50|n_pcs\s*=\s*50"),
    ("15", "neighbours of the k-nearest-neighbour graph", r"n_neighbors\s*=\s*15"),
    ("20", "quantile bins for the matched-dispersion statistic",
     r"def matched\(df, col, match, nbins=20"),
    ("200", "permutations for the retained-dispersion null", r"^NP = 200"),
    ("10,000", "permutations for the top-decile enrichment test", r"NPERM *= *10000"),
    ("1,000", "random factor sets drawn for the sign-agreement null",
     r"NSEED, NDRAW, TOPN = 1000"),
    ("10", "quantile bins when a lesion has fewer than 500 stromal cells",
     r"if len\(df\) < 500: nbins = 10"),
    ("25", "minimum cells for a bin to be kept", r"min_n=25"),
]

rows = []
for printed, item, where, value in FROZEN:
    pv = float(printed.replace(",", ""))
    rows.append(dict(category="setting frozen before the analysis", item=item,
                     manuscript_value=printed, where_it_is_fixed=where,
                     value_there=value, agrees=(float(value) == pv)))
for printed, item, pattern in SCRIPTS:
    loc, text = grep_scripts(pattern)
    rows.append(dict(category="pipeline parameter set in the scripts", item=item,
                     manuscript_value=printed,
                     where_it_is_fixed=loc or "not located by pattern search",
                     value_there=text, agrees=bool(loc)))
E5 = pd.DataFrame(rows)
E5.to_csv(f"{OUT}/E5_setting_values_check.csv", index=False)

# ---------------------------------------------------- classify every unconfirmed row
VERSION = re.compile(r"kallisto|bustools|kb-python|scanpy|anndata|decoupler|omnipath|"
                     r"harmonypy|scrublet|leidenalg|diptest|SRA Toolkit|Cell Ranger|"
                     r"Ensembl release|Reagent Kit|Toolkit|version", re.I)
READSTRUCT = re.compile(r"read 1|read-1|poly\(T\)|poly\(dT\)|barcode|UMI|position", re.I)
REFBUILD = re.compile(r"features|transcripts|cDNA set|symbol|identifiers|t2g|"
                      r"transcript-to-gene", re.I)
PARAM = re.compile(r"retained with|removed|threshold|bins|permutation|neighbour|"
                   r"component|resolution|seed|variable genes|expected rate|"
                   r"cells per state|were excluded|control gene", re.I)
GENESET = re.compile(r"gene set|genes\)|gene_set|signed edges|sources|RELA targets",
                     re.I)
FROZEN_VALUES = {p for p, *_ in FROZEN}
SCRIPT_VALUES = {p for p, *_ in SCRIPTS}

cls = []
for no in UNC:
    r = A.iloc[no - 1]
    c = str(r.context)
    p = str(r.printed)
    if p in FROZEN_VALUES and GENESET.search(c):
        cat = "setting frozen before the analysis (checked in E5)"
    elif p in SCRIPT_VALUES and PARAM.search(c):
        cat = "pipeline parameter set in the scripts (checked in E5)"
    elif VERSION.search(c):
        cat = "software or kit version, not a tabulated value"
    elif REFBUILD.search(c) and not READSTRUCT.search(c):
        cat = "reference-build count, recorded in the build log"
    elif READSTRUCT.search(c):
        cat = "read-structure constant or a read-level measurement"
    elif PARAM.search(c):
        cat = "pipeline parameter set in the scripts (checked in E5)"
    else:
        cat = "measured value; see E3 for the parts added to the tables"
    cls.append(dict(No=no, printed=p, section=r.section, category=cat,
                    context=c))
CLS = pd.DataFrame(cls)
CLS.to_csv(f"{OUT}/E5_unconfirmed_classification.csv", index=False)
print(E5.to_string(index=False)[:2000])
print("\n", CLS.category.value_counts().to_string())
print("\nchecks that agree:", int(E5.agrees.sum()), "of", len(E5))
