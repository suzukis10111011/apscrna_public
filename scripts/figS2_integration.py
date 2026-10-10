#!/usr/bin/env python
"""Purpose: Compare integration methods (Harmony vs uncorrected PCA) on batch-corrected RNA-seq data
Input: 14_robustness_integration.csv, 14_tf_rank_by_integration.csv
Output: figS2_integration (plot), _paneln_figS2.csv (panel metrics)

Figure S2: Comparison of integration methods (Harmony vs uncorrected PCA).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
I = pd.read_csv(f"{TAB}/14_robustness_integration.csv").iloc[0]
R = pd.read_csv(f"{TAB}/14_tf_rank_by_integration.csv")
f, gs = gridfig(180, 74, nrows=1, ncols=2, width_ratios=[0.8, 1.1], wspace=0.34)
a1 = f.add_subplot(gs[0, 0])
vals = [int(I.stromal_n_harmony), int(I.stromal_n_uncorrected), int(I.stromal_intersection)]
a1.bar(range(3), vals, color=[DISC, "#9970ab", "0.6"], lw=0, width=0.62)
for i, v in enumerate(vals): a1.text(i, v * 1.02, f"{v:,}", ha="center", fontsize=6)
a1.set_xticks(range(3)); a1.set_xticklabels(["Harmony", "uncorrected", "intersection"], fontsize=6.3)
a1.set_ylabel("stromal cells"); a1.set_ylim(0, max(vals) * 1.18)
a1.set_title(f"Stromal set: Jaccard {I.stromal_jaccard:.3f}", loc="left", fontsize=7.4)
a2 = f.add_subplot(gs[0, 1])
# Note: the 'harmony'/'uncorrected' columns in 14_tf_rank_by_integration.csv are pooled rho, not ranks
# (Range -0.42 to +0.51). If the axis starts at 0, about half of the negative side disappears from the figure, so use a symmetric axis.
c = [x for x in R.columns if "harmony" in x.lower()][0]
u = [x for x in R.columns if "uncorrect" in x.lower()][0]
a2.scatter(R[c], R[u], s=5, c="0.45", lw=0, rasterized=True)
m = float(np.nanmax(np.abs(np.r_[R[c].values, R[u].values]))) * 1.06
a2.plot([-m, m], [-m, m], color=VAL, lw=0.8, ls="--")
a2.axhline(0, color="0.85", lw=0.6); a2.axvline(0, color="0.85", lw=0.6)
a2.set_xlim(-m, m); a2.set_ylim(-m, m)
a2.set_xlabel("pooled Spearman $\\rho$, Harmony", fontsize=6.6)
a2.set_ylabel("pooled Spearman $\\rho$, uncorrected PCA", fontsize=6.6)
a2.set_title(f"TF-level agreement: Spearman {I.tf_rank_spearman:.3f} over {int(I.n_tf)} TFs;\ntop-20 overlap {int(I.top20_overlap)}/20",
             loc="left", fontsize=7.2)
panel(a1, "A", dx=-0.22); panel(a2, "B", dx=-0.20)
save(f, "figS2_integration")
pd.DataFrame([dict(figure="figS2_integration", panel="A", n_stromal_harmony=int(I.stromal_n_harmony),
                   n_stromal_uncorrected=int(I.stromal_n_uncorrected)),
              dict(figure="figS2_integration", panel="B", n_TFs=len(R), axis_quantity="pooled Spearman rho vs lf_D_tissueforming")]).to_csv(f"{TAB}/_paneln_figS2.csv", index=False)
