#!/usr/bin/env python
"""Purpose: Draft figure comparing this cohort with the independent periodontitis dataset GSE164241.
Input: none (reads results/tables/13_gse164241_summary.csv and 13_gse164241_tf_sign_reproduction.csv)
Output: fig7_gse164241

Figure 7: Comparison with the periodontitis dataset GSE164241. If the output of 13 is absent, nothing is created.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
p = f"{TAB}/13_gse164241_summary.csv"
if not os.path.exists(p):
    print("13 の出力が無い — 図7 はスキップ"); sys.exit(0)
R = pd.read_csv(p); Rd = R[(R.group == "ALL") & (R.donor != "POOLED")]
f, gs = gridfig(180, 72, nrows=1, ncols=2, wspace=0.38)
a1 = f.add_subplot(gs[0, 0])
a1.bar(range(len(Rd)), Rd.depth_controlled_rho, color=OKABE[3], width=0.62, lw=0)
a1.axhline(0, color="0.5", lw=0.6)
a1.axhline(-0.0337, color=DISC, lw=0.9, ls=":")
a1.text(len(Rd) - 0.4, -0.0337, " discovery pooled", color=DISC, fontsize=6, va="top", ha="right")
a1.set_xticks(range(len(Rd))); a1.set_xticklabels(Rd.donor, rotation=90, fontsize=5.8)
a1.set_ylabel("Depth-controlled $\\rho$ (inflammatory vs tissue-forming axis)")
a1.set_title("GSE164241 gingival stromal cells, per donor", loc="left", fontsize=7, pad=8)
a2 = f.add_subplot(gs[0, 1])
q = f"{TAB}/13_gse164241_tf_sign_reproduction.csv"
if os.path.exists(q):
    T = pd.read_csv(q); T = T[(T.group == "ALL") & (T.present == True)]
    cnt = T.groupby("axis").agg(n=("TF", "size"), m=("sign_matches_discovery", "sum")).sort_values("m")
    a2.barh(range(len(cnt)), cnt.m / cnt.n * 100, color=OKABE[3], height=0.6, lw=0)
    a2.set_yticks(range(len(cnt))); a2.set_yticklabels([i.replace("lf_", "") for i in cnt.index], fontsize=6.2)
    a2.set_xlabel("Candidates matching discovery sign (%)"); a2.set_xlim(0, 105)
    for j, (_, r) in enumerate(cnt.iterrows()):
        a2.text(r.m / r.n * 100 + 1.5, j, f"{int(r.m)}/{int(r.n)}", fontsize=6, va="center")
panel(a1, "A", dx=-0.17, dy=1.14); panel(a2, "B", dx=-0.46, dy=1.14)
save(f, "fig7_gse164241")
