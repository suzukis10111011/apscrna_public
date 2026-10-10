#!/usr/bin/env python
"""Purpose: Trace the "117/120 and 101/103 agreements" of the Methods back to the script
  output they come from, and report what the same rule gives under each estimator.
Input: results/tables/18_tf_sign_by_cellset.csv (written by scripts/18_stromal_markerdefined.py)
Output: 08_external_review/tables/E13_candidate_sign_agreement.csv,
  08_external_review/tables/E13_sign_mismatches.csv,
  08_external_review/tables/E15_sign_agreement_part.csv (the same counts as a part of
  Table S5, so that each count can be read from one cell)

The count is the one scripts/18_stromal_markerdefined.py prints at the end: for each
cohort and each of the six axes, the candidate TFs whose pooled rho has the same sign in
the marker-defined and in the cluster-based stromal set, over the TFs evaluable in both
(the script's cmp_.dropna()). The file carries two estimators of that rho, so the rule is
applied to each.
"""
import os

import numpy as np
import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
OUT = f"{W}/08_external_review/tables"
SRC = "results/tables/18_tf_sign_by_cellset.csv"

S = pd.read_csv(f"{W}/{SRC}")
EST = {"pooled_rho_published_plain_spearman":
       "plain Spearman, the estimator the published file was written with",
       "pooled_rho":
       "depth-controlled rho on [1, log10(total_counts), n_genes_by_counts], the "
       "estimator the column was recomputed with after the audit"}

rows, mism = [], []
for col, what in EST.items():
    p = S.pivot_table(index=["cohort", "axis", "TF"], columns="cellset",
                      values=col).dropna()
    same = np.sign(p.marker_defined) == np.sign(p.primary_cluster_based)
    for (coh, ax), g in same.groupby(level=[0, 1]):
        rows.append(dict(estimator=col, estimator_note=what, cohort=coh, axis=ax,
                         n_agree=int(g.sum()), n_evaluable=int(g.size),
                         agreement=f"{int(g.sum())}/{int(g.size)}"))
    for coh, g in same.groupby(level=0):
        rows.append(dict(estimator=col, estimator_note=what, cohort=coh,
                         axis="ALL SIX AXES", n_agree=int(g.sum()),
                         n_evaluable=int(g.size),
                         agreement=f"{int(g.sum())}/{int(g.size)}"))
    # the two programmes of the paper, which is what "the candidates" means elsewhere
    two = same[same.index.get_level_values("axis").isin(
        ["lf_A_cytokine_chemokine", "lf_D_tissueforming"])]
    for coh, g in two.groupby(level=0):
        rows.append(dict(estimator=col, estimator_note=what, cohort=coh,
                         axis="THE TWO PROGRAMMES (A and D)", n_agree=int(g.sum()),
                         n_evaluable=int(g.size),
                         agreement=f"{int(g.sum())}/{int(g.size)}"))
    bad = p[~same]
    for (coh, ax, tf), r in bad.iterrows():
        mism.append(dict(estimator=col, cohort=coh, axis=ax, TF=tf,
                         rho_marker_defined=round(float(r.marker_defined), 4),
                         rho_primary_cluster_based=round(
                             float(r.primary_cluster_based), 4)))
E13 = pd.DataFrame(rows)
E13["source_file"] = SRC
E13["where_in_the_tables"] = (
    "Table S5, part 'candidate sign by cell set' (rows 24-487): cohort in column L, "
    "cellset in M, axis in AD, TF in AE, pooled_rho in AF, "
    "pooled_rho_published_plain_spearman in AI")
E13.to_csv(f"{OUT}/E13_candidate_sign_agreement.csv", index=False)
PART = E13.rename(columns={"source_file": "source_file_of_part",
                           "where_in_the_tables": "where_the_counts_come_from",
                           # S5 already has a column called estimator, belonging to the
                           # per-TF part; this part gets its own name so that no
                           # existing column changes meaning
                           "estimator": "estimator_of_the_pooled_rho",
                           "estimator_note": "estimator_of_the_pooled_rho_note"})
PART.to_csv(f"{OUT}/E15_sign_agreement_part.csv", index=False)
pd.DataFrame(mism).to_csv(f"{OUT}/E13_sign_mismatches.csv", index=False)

tot = E13[E13.axis == "ALL SIX AXES"]
print(tot[["estimator", "cohort", "agreement"]].to_string(index=False))
print("\nmanuscript prints 117/120 (discovery) and 101/103 (validation)")
for _, r in tot.iterrows():
    print(f"  {r.cohort:14s} {r.estimator:38s} {r.agreement}")
print("\nmismatching factors, plain Spearman:")
print(pd.DataFrame(mism).query("estimator == 'pooled_rho_published_plain_spearman'")
      .to_string(index=False))
