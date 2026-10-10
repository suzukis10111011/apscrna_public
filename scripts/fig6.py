#!/usr/bin/env python
"""Purpose: Plotting inter-axis coupling across compartments
Input: {TAB}/12_coupling_by_compartment.csv
Output: fig6_compartments

Figure 6: Inter-axon coupling in other compartments.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
C = pd.read_csv(f"{TAB}/12_coupling_by_compartment.csv")
C = C[C.excluded.fillna("") == ""]
f, gs = gridfig(180, 76, nrows=1, ncols=2, width_ratios=[1.3, 1], wspace=0.36)
a1 = f.add_subplot(gs[0, 0])
piv = C[C.lesion != "POOLED"].pivot_table(index="compartment", columns="lesion", values="depth_controlled_rho")
order = C[C.lesion == "POOLED"].set_index("compartment").depth_controlled_rho.sort_values().index
piv = piv.reindex([o for o in order if o in piv.index])
for j, comp in enumerate(piv.index):
    v = piv.loc[comp].dropna()
    a1.scatter(v.values, [j] * len(v), s=14, color=DISC, alpha=0.8, linewidths=0)
    p = C[(C.compartment == comp) & (C.lesion == "POOLED")].depth_controlled_rho
    if len(p): a1.plot([float(p.iloc[0])] * 2, [j - 0.26, j + 0.26], color="0.15", lw=1.3)
a1.axvline(0, color="0.5", lw=0.6)
a1.set_yticks(range(len(piv))); a1.set_yticklabels(piv.index, fontsize=6.4)
a1.set_xlabel("Depth-controlled Spearman $\\rho$ (inflammatory vs tissue-forming axis)")
a1.set_title("dots = lesions, bar = pooled", loc="left", fontsize=7)
a2 = f.add_subplot(gs[0, 1])
d = C[(C.compartment == "immune_myeloid") & C.myeloid_M1_vs_M2_depth_controlled_rho.notna()]
dl = d[d.lesion != "POOLED"]
a2.bar(range(len(dl)), dl.myeloid_M1_vs_M2_depth_controlled_rho, color=OKABE[2], width=0.6, lw=0)
p = d[d.lesion == "POOLED"]
if len(p): a2.axhline(float(p.myeloid_M1_vs_M2_depth_controlled_rho.iloc[0]), color="0.25", lw=0.9, ls="--")
a2.axhline(0, color="0.5", lw=0.6)
a2.set_xticks(range(len(dl))); a2.set_xticklabels(dl.lesion, fontsize=6.5)
a2.set_ylabel("Depth-controlled $\\rho$\n(myeloid M1 vs M2/repair)")
a2.set_title("Myeloid-specific axes", loc="left", fontsize=7)
panel(a1, "A", dx=-0.24); panel(a2, "B", dx=-0.22)
save(f, "fig6_compartments")
