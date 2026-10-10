#!/usr/bin/env python
"""Purpose: Reproduce Figure 4 validation results — cell-set dependence of ρ and sign recovery for candidate
  TFs
Input: 05e_rho_prev_cellset_vs_clean.csv, 25_validation_states_within_vs_ingest.csv
Output: fig4_validation

Figure 4: Reproduction on the validation side — cell-set dependence of ρ and sign recovery for candidate TFs identified in the discovery analysis.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
R = pd.read_csv(f"{TAB}/05e_rho_prev_cellset_vs_clean.csv")
f, gs = gridfig(180, 70, nrows=1, ncols=3, width_ratios=[1.05, 0.9, 1.0], wspace=0.78)
a1 = f.add_subplot(gs[0, 0])
les = [l for l in R.lesion.unique() if not str(l).startswith("POOLED")] + ["POOLED_validation"]
x = np.arange(len(les))
for k, (cs, col, lab) in enumerate([("prev_stromal_set", "0.62", "previous stromal set"),
                                     ("clean_stromal_2026", DISC, "this study stromal set")]):
    d = R[R.cellset == cs].set_index("lesion").reindex(les)
    a1.bar(x + (k - 0.5) * 0.36, d.depth_controlled_rho, width=0.34, color=col, lw=0, label=lab)
_y = R[R.lesion.isin(les)].depth_controlled_rho.astype(float)
a1.set_ylim(_y.min() * 1.15, max(0.065, _y.max() * 1.3))   # Determine from the data so that CAP2 is not clipped at the lower bound
a1.axhline(0.0398, color=VAL, lw=0.9, ls=":")
a1.text(0.02, 0.0398, " grant +0.040", color=VAL, fontsize=6, va="bottom", ha="left")
a1.axhline(0, color="0.5", lw=0.6)
a1.set_xticks(x); a1.set_xticklabels([str(l).replace("POOLED_validation", "POOLED") for l in les], fontsize=6.5)
a1.set_ylabel("Depth-controlled Spearman $\\rho$\n(inflammatory vs tissue-forming composite)")
a1.legend(loc="lower left", handlelength=1.1)
a2 = f.add_subplot(gs[0, 1])
p = f"{TAB}/05_tf_sign_reproduction_validation.csv"
if os.path.exists(p):
    T = pd.read_csv(p)
    cnt = T.groupby("axis").agg(n=("TF", "size"), rep=("reproduced", "sum")) if "reproduced" in T.columns else None
else:
    cnt = None
if cnt is None:
    vals = {"lf_A_cytokine_chemokine": (20, 20), "lf_B_survival_antiapoptotic": (20, 12),
            "lf_C1_proliferation": (20, 11), "lf_C2_migration": (20, 10),
            "lf_C_prolif_migration": (20, 10), "lf_D_tissueforming": (20, 19)}
    cnt = pd.DataFrame([dict(axis=k, n=v[0], rep=v[1]) for k, v in vals.items()]).set_index("axis")
cnt = cnt.sort_values("rep")
a2.barh(range(len(cnt)), cnt.rep / cnt.n * 100, color=VAL, height=0.62, lw=0)
a2.set_yticks(range(len(cnt))); a2.set_yticklabels([i.replace("lf_", "") for i in cnt.index], fontsize=6.2)
a2.set_xlabel("Candidates with consistent sign\nin all 3 validation lesions (%)"); a2.set_xlim(0, 105)
for j, (_, r) in enumerate(cnt.iterrows()):
    a2.text(r.rep / r.n * 100 + 1.5, j, f"{int(r.rep)}/{int(r.n)}", fontsize=6, va="center")
panel(a1, "A", dx=-0.15, dy=1.10); panel(a2, "B", dx=-0.46, dy=1.10)
# ---- Panel C: state proportions (within-cohort definition vs ingest projection) ----
a3 = f.add_subplot(gs[0, 2])
V = pd.read_csv(f"{TAB}/25_validation_states_within_vs_ingest.csv")
STATES = ["inflammatory", "intermediate", "tissue-forming", "other"]
COLS = {"inflammatory": "#440154", "intermediate": "#31688e",
        "tissue-forming": "#35b779", "other": "#bdbdbd"}
les3 = [l for l in V.lesion.unique() if l != "POOLED"] + ["POOLED"]
ypos, ylab = [], []
for i, l in enumerate(les3):
    for j, lab in enumerate(["within_cohort", "ingest_projected"]):
        row = V[(V.lesion == l) & (V.labelling == lab)]
        if not len(row): continue
        y = i * 2.6 + (0 if j == 0 else 1.0)
        left = 0.0
        for s in STATES:
            w = float(row[f"pct_{s}"].iloc[0])
            a3.barh(y, w, left=left, height=0.85, color=COLS[s], lw=0)
            left += w
        ypos.append(y); ylab.append(f"{l}  {'within' if j == 0 else 'ingest'}")
a3.set_yticks(ypos); a3.set_yticklabels(ylab, fontsize=5.4)
a3.invert_yaxis(); a3.set_xlim(0, 100); a3.set_xlabel("% of stromal cells", fontsize=6.5)
a3.set_title("State labels: within-cohort vs projected\n(agreement 40.9% pooled)", loc="left", fontsize=7)
from matplotlib.patches import Patch
a3.legend(handles=[Patch(facecolor=COLS[s], label=s) for s in STATES],
          fontsize=5.4, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.20),
          handlelength=0.9, columnspacing=0.7, frameon=False)
panel(a3, "C", dx=-0.42, dy=1.10)
save(f, "fig4_validation")
