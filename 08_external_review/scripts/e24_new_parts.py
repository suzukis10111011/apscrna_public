#!/usr/bin/env python
"""Purpose: Build the three parts the v15 legends promise but the workbook did not hold:
  the members of the six programme gene sets, the validation values of every evaluable
  factor, and the factor-by-factor correlation matrix of Figure 5D.
Input: frozen/genesets/stage1_geneset_definitions.json,
  data/processed/all_cells_{discovery,validation}.h5ad (gene symbols only),
  results/tables/33_tf_validation_all730.csv,
  results/tables/32_tf_activity_spearman_matrix.csv
Output: 08_external_review/tables/E24_geneset_members.csv,
  08_external_review/tables/E25_validation_all_factors.csv,
  08_external_review/tables/E26_tf_activity_matrix.csv

The gene spaces are read straight out of the two matrices (the var symbols), so the
presence flags are what the analysis actually saw rather than a recount of the coverage
table; the recount is checked against 20b_geneset_coverage_before_after.csv.
"""
import json
import os

import h5py
import numpy as np
import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
TAB = f"{W}/results/tables"
OUT = f"{W}/08_external_review/tables"


def symbols(path):
    """the gene symbols of an .h5ad, without reading the matrix"""
    with h5py.File(path, "r") as h:
        v = h["var"]
        key = "symbol" if "symbol" in v else "_index"
        d = v[key]
        if isinstance(d, h5py.Group):          # categorical
            cats = np.array(d["categories"]).astype(str)
            return set(cats[np.array(d["codes"])])
        return set(np.array(d).astype(str))


DISC = symbols(f"{W}/data/processed/all_cells_discovery.h5ad")
VAL = symbols(f"{W}/data/processed/all_cells_validation.h5ad")
AX = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))["axes"]

rows = []
for s, genes in AX.items():
    for g in sorted(genes):
        rows.append(dict(programme_gene_set=s, n_genes_in_set=len(genes), gene_symbol=g,
                         in_discovery_gene_space=g in DISC,
                         in_validation_gene_space=g in VAL,
                         source_file_of_part="frozen/genesets/"
                                             "stage1_geneset_definitions.json"))
MEM = pd.DataFrame(rows)
MEM.to_csv(f"{OUT}/E24_geneset_members.csv", index=False)

COV = pd.read_csv(f"{TAB}/20b_geneset_coverage_before_after.csv")
chk = (MEM.groupby("programme_gene_set")
          .agg(n_def=("gene_symbol", "size"),
               found_validation=("in_validation_gene_space", "sum"),
               found_discovery=("in_discovery_gene_space", "sum")).reset_index())
chk["gene_set_in_coverage_table"] = "axis " + chk.programme_gene_set
chk = chk.merge(COV[["gene_set", "n_def", "found_after"]],
                left_on="gene_set_in_coverage_table", right_on="gene_set",
                how="left", suffixes=("", "_cov"))

V = pd.read_csv(f"{TAB}/33_tf_validation_all730.csv")
V = V.assign(source_file_of_part="33_tf_validation_all730.csv",
             figure_panel="Figure 3A, B")
V.to_csv(f"{OUT}/E25_validation_all_factors.csv", index=False)

M = pd.read_csv(f"{TAB}/32_tf_activity_spearman_matrix.csv")
M = M.assign(source_file_of_part="32_tf_activity_spearman_matrix.csv",
             figure_panel="Figure 5D")
M.to_csv(f"{OUT}/E26_tf_activity_matrix.csv", index=False)

print(chk[["programme_gene_set", "n_def", "found_discovery", "found_validation",
           "n_def_cov", "found_after"]].to_string(index=False))
print(f"\nmembers {len(MEM)} rows | validation factors {len(V)} rows, "
      f"{int(V.in_validation.sum())} evaluable, {V.axis.nunique()} programmes | "
      f"matrix {M.shape[0]} rows x {M.shape[1]} columns")
