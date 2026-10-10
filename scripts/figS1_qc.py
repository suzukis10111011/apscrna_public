#!/usr/bin/env python
"""Purpose: Plot cell counts per QC stage across both cohorts and samples (Figure S1)
Input: {TAB}/28_qc_all_stages.csv
Output: figS1_qc, {TAB}/_paneln_figS1.csv

Figure S1: Cell counts by QC stage (both cohorts, by sample).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
Q = pd.read_csv(f"{TAB}/28_qc_all_stages.csv")
f, gs = gridfig(180, 72, nrows=1, ncols=2, width_ratios=[1, 1], wspace=0.30)
NP = []
for i, (coh, lab) in enumerate([("discovery_n5", "Discovery (n=5)"), ("validation_n3", "Validation (n=3)")]):
    a = f.add_subplot(gs[0, i])
    d = Q[Q.cohort == coh]
    stages = list(dict.fromkeys(d.stage)); les = sorted(d.lesion.unique())
    x = np.arange(len(les)); w = 0.8 / max(len(stages), 1)
    for j, s in enumerate(stages):
        v = [float(d[(d.lesion == l) & (d.stage == s)].n_cells.sum()) for l in les]
        a.bar(x + j * w - 0.4 + w / 2, v, width=w * 0.92, lw=0,
              color=plt.cm.Blues(0.35 + 0.28 * j) if i == 0 else plt.cm.Oranges(0.35 + 0.28 * j),
              label=s.replace("_", " "))
    a.set_xticks(x); a.set_xticklabels(les, fontsize=6.5)
    a.set_ylabel("cells"); a.set_title(lab, loc="left", fontsize=7.5)
    a.legend(fontsize=5.4, handlelength=0.9, frameon=False, ncol=1, loc="upper right")
    panel(a, "AB"[i], dx=-0.17)
    NP.append(dict(figure="figS1_qc", panel="AB"[i], cohort=coh, n_samples=len(les),
                   n_cells_final=int(d[d.stage == stages[-1]].n_cells.sum()), stages=";".join(stages)))
save(f, "figS1_qc")
pd.DataFrame(NP).to_csv(f"{TAB}/_paneln_figS1.csv", index=False)
