#!/usr/bin/env python
"""Purpose: Generate cohort composition and interstitial enrichment visualization
Input: 03c_compartment_by_lesion_discovery.csv, 03c_enrichment_by_lesion_discovery.csv
Output: fig1_cohort_composition

Figure 1: Cohort composition and interstitial enrichment.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
C = pd.read_csv(f"{TAB}/03c_compartment_by_lesion_discovery.csv", index_col=0)
E = pd.read_csv(f"{TAB}/03c_enrichment_by_lesion_discovery.csv")
E = E[E.score == "axis_D_tissueforming"]
f, gs = gridfig(180, 72, nrows=1, ncols=2, width_ratios=[1.25, 1], wspace=0.34)
a1 = f.add_subplot(gs[0, 0])
pct = C[[c for c in C.columns if c.startswith("pct_")]].rename(columns=lambda s: s[4:])
pct = pct[pct.sum().sort_values(ascending=False).index]
b = np.zeros(len(pct))
for i, c in enumerate(pct.columns):
    a1.bar(range(len(pct)), pct[c], bottom=b, color=OKABE[i % len(OKABE)], label=c, width=0.72, lw=0)
    b += pct[c].values
a1.set_xticks(range(len(pct))); a1.set_xticklabels(pct.index, fontsize=6.5)
a1.set_ylabel("Fraction of cells (%)"); a1.set_ylim(0, 100)
a1.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.30), handlelength=1.1, columnspacing=0.9)
a2 = f.add_subplot(gs[0, 1])
d = E[E.lesion != "POOLED_discovery"]
a2.bar(range(len(d)), d.enrichment_fold, color=DISC, width=0.62, lw=0)
p = E[E.lesion == "POOLED_discovery"]
if len(p): a2.axhline(float(p.enrichment_fold.iloc[0]), color="0.25", lw=0.8, ls="--")
a2.axhline(6.42, color=VAL, lw=0.8, ls=":")
a2.text(len(d) - 0.4, 6.42, " grant 6.42", color=VAL, fontsize=6, va="bottom", ha="right")
a2.axhline(1, color="0.6", lw=0.5)
a2.set_xticks(range(len(d))); a2.set_xticklabels(d.lesion, fontsize=6.5)
a2.set_ylabel("Stromal enrichment in\ntissue-forming top decile (fold)")
panel(a1, "A", dx=-0.10); panel(a2, "B", dx=-0.19)
save(f, "fig1_cohort_composition")
