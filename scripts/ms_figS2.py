#!/usr/bin/env python
"""Purpose: Compare integration methods (Harmony vs uncorrected PCA) on TF discovery and robustness
Input: {TAB}/04_tf_screen_discovery_candidates.csv, {TAB}/14_robustness_integration.csv,
  {TAB}/14_tf_rank_by_integration.csv
Output: FigS2_integration

Supplementary Figure S2: integration method (Harmony vs uncorrected PCA).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np
I = pd.read_csv(f"{TAB}/14_robustness_integration.csv").iloc[0]
R = pd.read_csv(f"{TAB}/14_tf_rank_by_integration.csv")
S04 = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
n_screen = int(S04.TF.nunique())
f, gs = gridfig(180, 80, nrows=1, ncols=2, width_ratios=[0.8, 1.1], wspace=0.36)
a1 = f.add_subplot(gs[0, 0])
vals = [int(I.stromal_n_harmony), int(I.stromal_n_uncorrected), int(I.stromal_intersection)]
a1.bar(range(3), vals, color=[DISC, "#9970ab", "0.62"], lw=0, width=0.62)
# Value labels are removed in the plain figure and moved to a pptx text box (numbox)
for i, v in enumerate(vals): a1.text(i, v * 1.02, f"{v:,}", ha="center", fontsize=7, numbox=True)
a1.set_xticks(range(3))
a1.set_xticklabels(["Harmony", "uncorrected\nPCA", "intersection"], fontsize=6.2)
a1.set_ylabel(comp("stromal") + " identified"); a1.set_ylim(0, max(vals) * 1.16)

panel(a1, "A", dx=-0.24)
a2 = f.add_subplot(gs[0, 1])
c = [x for x in R.columns if "harmony" in x.lower()][0]
u = [x for x in R.columns if "uncorrect" in x.lower()][0]
a2.scatter(R[c], R[u], s=5, c="0.45", lw=0, rasterized=True)
m = float(np.nanmax(np.abs(np.r_[R[c].values, R[u].values]))) * 1.06
a2.plot([-m, m], [-m, m], color=VAL, lw=0.8, ls="--")
a2.axhline(0, color="0.85", lw=0.6); a2.axvline(0, color="0.85", lw=0.6)
a2.set_xlim(-m, m); a2.set_ylim(-m, m)
a2.set_xlabel("pooled Spearman $\\rho$, Harmony", fontsize=6.6)
a2.set_ylabel("pooled Spearman $\\rho$, uncorrected PCA", fontsize=6.6)


panel(a2, "B", dx=-0.20)
save(f, "FigS2_integration")
paneln([dict(figure="FigS2", panel="A", content="stromal set by integration method",
             source_table="14_robustness_integration.csv", n_cells=int(I.stromal_n_harmony),
             n_groups=3, group_unit="methods", cohort="discovery_n5"),
        dict(figure="FigS2", panel="B", content="TF-level rho, Harmony vs uncorrected",
             source_table="14_tf_rank_by_integration.csv", n_cells=int(I.stromal_n_harmony),
             n_groups=int(I.n_tf), group_unit="TFs", cohort="discovery_n5")], "FigS2")
