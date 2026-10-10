#!/usr/bin/env python
"""Purpose: Rank ten named factors by leave-one-lesion-out cross-validation performance
Input: {TAB}/34_lolo_rulebased.csv
Output: FigS6_leave_one_out

Supplementary Figure S6: leave-one-lesion-out ranks of the ten named factors.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np
L = pd.read_csv(f"{TAB}/34_lolo_rulebased.csv")
L["is_candidate"] = L.is_candidate.astype(str).str.lower().isin(["true", "1", "1.0"])
ORD = ["RELB", "NFKB2", "RELA", "NFKB1", "REL", "SP7", "RUNX2", "DLX5", "MYRF", "SMAD3"]
held = sorted(L.held_out_lesion.unique())
# Height 118 mm: setting in-cell values to 7 pt causes 10 rows x 5 columns to collapse at 80 mm
f, gs = gridfig(180, 118, nrows=1, ncols=2, width_ratios=[1.15, 1.0], wspace=0.40)
a1 = f.add_subplot(gs[0, 0])
M = L.pivot_table(index="TF", columns="held_out_lesion", values="rank_within_axis",
                  observed=True).reindex(ORD)
CM1, V0, V1 = "Blues", 1, 25          # Pale monochrome palette (instruction [3])
im = a1.imshow(M.values, cmap=CM1, vmin=V0, vmax=V1, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M.values[i, j]
        a1.text(j, i, "\u2014" if not np.isfinite(v) else f"{int(v)}", ha="center", va="center",
                fontsize=7, color=cell_text_colour(CM1, v, V0, V1), numbox=True)
a1.set_xticks(range(len(held))); a1.set_xticklabels([f"\u2212{h}" for h in held], fontsize=6.2)
a1.set_yticks(range(len(ORD))); a1.set_yticklabels(ORD, fontsize=6.2)
a1.set_xlabel("lesion left out", fontsize=6.6)
a1.axhline(4.5, color="0.35", lw=1.0)
cb = f.colorbar(im, ax=a1, fraction=0.045, pad=0.03)
cb.set_label("rank within programme", fontsize=6); cb.ax.tick_params(labelsize=5.4)

a2 = f.add_subplot(gs[0, 1])
cnt = L.groupby("TF").is_candidate.sum().reindex(ORD)
n_all = int((cnt == len(held)).sum())
a2.barh(range(len(cnt)), cnt.values, color="#1b7837", height=0.62, lw=0)
a2.invert_yaxis()
a2.set_yticks(range(len(cnt))); a2.set_yticklabels(ORD, fontsize=6.2)
a2.set_xlim(0, len(held) + 0.4); a2.set_xticks(range(len(held) + 1))
a2.axvline(len(held), color="0.6", lw=0.7, ls="--")
a2.set_xlabel("runs in which the factor remains a candidate", fontsize=6.4)

panel(a1, "A", dx=-0.16); panel(a2, "B", dx=-0.22)
save(f, "FigS6_leave_one_out")
paneln([dict(figure="FigS6", panel="A+B", content="leave-one-lesion-out ranks and candidacy",
             source_table="34_lolo_rulebased.csv", n_cells=6215, n_groups=len(ORD),
             group_unit="factors", cohort="discovery_n5")], "FigS6")
