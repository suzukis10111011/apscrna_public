#!/usr/bin/env python
"""Purpose: Assess robustness of gene-set scoring methods across parameters
Input: {TAB}/14_robustness_genesets_scoring.csv
Output: figS3_genesets, {TAB}/_paneln_figS3.csv

Fig. S3: Robustness of gene set x scoring method.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
B = pd.read_csv(f"{TAB}/14_robustness_genesets_scoring.csv")
P = B[B.lesion == "POOLED"].copy()
P["combo"] = P.inflammatory_set.str[:10] + " / " + P.tissueforming_set.str[:10] + "\n" + P.scoring
P = P.sort_values("depth_controlled_rho").reset_index(drop=True)
f, gs = gridfig(180, 86, nrows=1, ncols=2, width_ratios=[1.25, 0.85], wspace=0.44)
a1 = f.add_subplot(gs[0, 0])
col = [DISC if v < 0 else VAL for v in P.depth_controlled_rho]
a1.barh(range(len(P)), P.depth_controlled_rho, color=col, height=0.66, lw=0)
a1.axvline(0, color="0.5", lw=0.7)
a1.set_yticks(range(len(P))); a1.set_yticklabels(P.combo, fontsize=4.9)
a1.set_xlabel("Depth-controlled Spearman $\\rho$ (pooled discovery)", fontsize=6.6)
neg = int((P.depth_controlled_rho < 0).sum())
a1.set_title(f"Negative in {neg}/{len(P)} gene-set x scoring combinations", loc="left", fontsize=7.4)
a2 = f.add_subplot(gs[0, 1])
a2.scatter(P.defA_frac_retained, range(len(P)), s=14, c="0.35", lw=0)
a2.axvline(1.0, color="0.7", lw=0.7, ls="--")
a2.set_yticks([]); a2.set_xlabel("retained dispersion\n(definition A)", fontsize=6.6)
a2.set_title(f"{P.defA_frac_retained.min():.3f}-{P.defA_frac_retained.max():.3f}", loc="left", fontsize=7.4)
panel(a1, "A", dx=-0.52); panel(a2, "B", dx=-0.16)
save(f, "figS3_genesets")
pd.DataFrame([dict(figure="figS3_genesets", panel="A+B", n_combinations=len(P),
                   n_negative=neg, cohort="discovery pooled")]).to_csv(f"{TAB}/_paneln_figS3.csv", index=False)
