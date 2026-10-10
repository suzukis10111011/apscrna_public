#!/usr/bin/env python
"""Purpose: Plot cell counts across QC filtering stages
Input: {TAB}/28_qc_all_stages.csv
Output: FigS1_qc

Supplementary Figure S1: cell counts at each QC stage.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np
Q = pd.read_csv(f"{TAB}/28_qc_all_stages.csv")
f, gs = gridfig(180, 76, nrows=1, ncols=2, wspace=0.30)
STAGE = {"0_raw": "raw", "0_raw_ge100umi": "raw (\u2265100 UMIs)",
         "1_qc_pass": "QC pass", "2_singlet": "singlets"}
NP = []
for i, (coh, lab, cm) in enumerate([("discovery_n5", "Discovery (5 lesions)", plt.cm.Blues),
                                    ("validation_n3", "Validation (3 lesions)", plt.cm.Oranges)]):
    a = f.add_subplot(gs[0, i]); d = Q[Q.cohort == coh]
    stages = list(dict.fromkeys(d.stage)); les = [l for l in LESION_ORDER if l in set(d.lesion)]
    x = np.arange(len(les)); w = 0.8 / max(len(stages), 1)
    for j, s in enumerate(stages):
        v = [float(d[(d.lesion == l) & (d.stage == s)].n_cells.sum()) for l in les]
        a.bar(x + j * w - 0.4 + w / 2, v, width=w * 0.92, lw=0,
              color=cm(0.34 + 0.26 * j), label=STAGE.get(s, s.replace("_", " ")))
    a.set_xticks(x); a.set_xticklabels(les, fontsize=6.4)
    a.set_ylabel("cells"); a.set_title(lab, loc="left", fontsize=7.4)
    a.legend(fontsize=5.3, handlelength=0.9, loc="upper right")
    panel(a, "AB"[i], dx=-0.17)
    NP.append(dict(figure="FigS1", panel="AB"[i], content="cells per QC stage",
                   source_table="28_qc_all_stages.csv",
                   n_cells=int(d[d.stage == stages[-1]].n_cells.sum()), n_groups=len(les),
                   group_unit="lesions", cohort=coh))
save(f, "FigS1_qc")
paneln(NP, "FigS1")
