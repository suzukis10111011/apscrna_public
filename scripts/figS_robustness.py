#!/usr/bin/env python
"""Purpose: Assess robustness of gene set, scoring method, and integration approaches
Input: {TAB}/14_robustness_genesets_scoring.csv
Output: figS_robustness

Supplementary figure: Robustness of gene sets, scoring methods, and integration approaches.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
R = pd.read_csv(f"{TAB}/14_robustness_genesets_scoring.csv")
P = R[R.lesion == "POOLED"]
f, gs = gridfig(180, 78, nrows=1, ncols=2, width_ratios=[1.45, 1], wspace=0.40)
a1 = f.add_subplot(gs[0, 0])
meths = list(P.scoring.unique())
combos = P[["inflammatory_set", "tissueforming_set"]].drop_duplicates().apply(
    lambda r: f"{r.inflammatory_set[:26]}\n× {r.tissueforming_set[:26]}", axis=1).tolist()
Mx = np.full((len(combos), len(meths)), np.nan)
keys = P[["inflammatory_set", "tissueforming_set"]].drop_duplicates().values.tolist()
for i, (a, b) in enumerate(keys):
    for j, m in enumerate(meths):
        v = P[(P.inflammatory_set == a) & (P.tissueforming_set == b) & (P.scoring == m)].depth_controlled_rho
        if len(v): Mx[i, j] = float(v.iloc[0])
im = a1.imshow(Mx, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto")
a1.set_xticks(range(len(meths))); a1.set_xticklabels(meths, fontsize=6.4)
a1.set_yticks(range(len(combos))); a1.set_yticklabels(combos, fontsize=5.6)
for i in range(Mx.shape[0]):
    for j in range(Mx.shape[1]):
        if not np.isnan(Mx[i, j]):
            a1.text(j, i, f"{Mx[i,j]:+.3f}", ha="center", va="center", fontsize=5.8,
                    color="white" if abs(Mx[i, j]) > 0.18 else "0.15")
cb = f.colorbar(im, ax=a1, fraction=0.045, pad=0.03); cb.set_label("Depth-controlled $\\rho$ (pooled)", fontsize=6.2); cb.ax.tick_params(labelsize=6)
a2 = f.add_subplot(gs[0, 1])
d = P[["inflammatory_set", "tissueforming_set", "scoring", "defA_frac_retained"]].dropna()
a2.scatter(range(len(d)), d.defA_frac_retained, s=16, color=DISC, linewidths=0)
a2.set_ylim(0.97, 1.005); a2.set_xticks([])
a2.set_xlabel("Gene-set $\\times$ scoring-method combinations")
a2.set_ylabel("retained dispersion (pooled)")
a2.axhline(1.0, color="0.55", lw=0.6, ls="--")
I = f"{TAB}/14_robustness_integration.csv"
if os.path.exists(I):
    r = pd.read_csv(I).iloc[0]
    a2.text(0.02, 0.06, f"Integration: Jaccard {r.stromal_jaccard:.3f} / TF-rank rho {r.tf_rank_spearman:.3f}",
            transform=a2.transAxes, fontsize=6, color="0.3")
panel(a1, "A", dx=-0.30); panel(a2, "B", dx=-0.22)
save(f, "figS_robustness")
