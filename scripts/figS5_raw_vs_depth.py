#!/usr/bin/env python
"""Purpose: Compare Spearman rho between raw and depth-controlled inflammation-tissue-forming correlations
  across discovery and validation cohorts
Input: 04_rho_inflam_vs_tissueforming_discovery.csv, 05_rho_inflam_vs_tissueforming_validation.csv
Output: figS5_raw_vs_depth, _paneln_figS5.csv

Fig. S5: Comparison of raw vs depth-controlled Spearman rho (both cohorts, both definitions).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
D = pd.read_csv(f"{TAB}/04_rho_inflam_vs_tissueforming_discovery.csv")
V = pd.read_csv(f"{TAB}/05_rho_inflam_vs_tissueforming_validation.csv")
D["cohort_lab"] = "discovery"; V["cohort_lab"] = "validation"
A = pd.concat([D, V], ignore_index=True)
A = A[A.cellset == "stromal"]
f, gs = gridfig(180, 76, nrows=1, ncols=2, width_ratios=[1, 1], wspace=0.34)
MK = {"discovery": "o", "validation": "s"}
NP = []
for i, defn in enumerate(["axis", "composite"]):
    a = f.add_subplot(gs[0, i])
    d = A[A.definition == defn]
    for coh, g in d.groupby("cohort_lab"):
        pooled = g.lesion.astype(str).str.startswith("POOLED")
        a.scatter(g.raw_rho[~pooled], g.depth_controlled_rho[~pooled], s=26, marker=MK[coh],
                  c=DISC if coh == "discovery" else VAL, edgecolor="white", linewidth=0.5,
                  label=f"{coh} (per lesion)")
        a.scatter(g.raw_rho[pooled], g.depth_controlled_rho[pooled], s=60, marker="*",
                  c=DISC if coh == "discovery" else VAL, edgecolor="0.2", linewidth=0.4,
                  label=f"{coh} (pooled)")
    lo = min(d.raw_rho.min(), d.depth_controlled_rho.min()) - 0.05
    hi = max(d.raw_rho.max(), d.depth_controlled_rho.max()) + 0.05
    a.plot([lo, hi], [lo, hi], color="0.6", lw=0.7, ls="--")
    a.axhline(0, color="0.8", lw=0.6); a.axvline(0, color="0.8", lw=0.6)
    a.set_xlim(lo, hi); a.set_ylim(lo, hi)
    a.set_xlabel("raw Spearman $\\rho$", fontsize=6.6)
    a.set_ylabel("depth-controlled Spearman $\\rho$", fontsize=6.6)
    a.set_title(f"{defn} definition", loc="left", fontsize=7.4)
    if i == 0: a.legend(fontsize=5.2, frameon=False, handlelength=0.9, loc="upper left")
    panel(a, "AB"[i], dx=-0.22)
    NP.append(dict(figure="figS5_raw_vs_depth", panel="AB"[i], definition=defn, n_points=len(d),
                   n_lesions=int((~d.lesion.astype(str).str.startswith("POOLED")).sum())))
save(f, "figS5_raw_vs_depth")
pd.DataFrame(NP).to_csv(f"{TAB}/_paneln_figS5.csv", index=False)
