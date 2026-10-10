#!/usr/bin/env python
"""Purpose: Reproduce supplementary Figure S8 showing candidate validation for four remaining programmes across
  two datasets
Input: 05_tf_screen_validation_reproduction.csv, 13_gse164241_tf_sign_reproduction.csv
Output: FigS8_other_programs

Supplementary Figure S8: candidate reproduction for the four remaining programmes.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np
V = pd.read_csv(f"{TAB}/05_tf_screen_validation_reproduction.csv")
G = pd.read_csv(f"{TAB}/13_gse164241_tf_sign_reproduction.csv")
for d in (V, G):
    for c in d.columns:
        if d[c].dtype == object and set(d[c].astype(str).str.lower().unique()) <= {"true", "false", "nan"}:
            d[c] = d[c].astype(str).str.lower().eq("true")
MAIN = ["lf_A_cytokine_chemokine", "lf_D_tissueforming"]
others = [a for a in V.axis.unique() if a not in MAIN]
f, gs = gridfig(180, 84, nrows=1, ncols=2, width_ratios=[1, 1], wspace=0.52)
rows = []
a1 = f.add_subplot(gs[0, 0])
Rv = []
for a in others:
    d = V[(V.axis == a) & V.in_validation]
    Rv.append(dict(axis=a, n=len(d), rep=int(d.reproduced.sum())))
Rv = pd.DataFrame(Rv)
y = np.arange(len(Rv))
a1.barh(y, 100 * Rv.rep / Rv.n, color=VAL, height=0.6, lw=0)
for j, r in Rv.iterrows():
    a1.text(100 * r.rep / r.n + 1.5, j, f"{r.rep}/{r.n}", fontsize=5.8, va="center", bare_keep=True)
a1.set_yticks(y); a1.set_yticklabels([prog(a) for a in Rv.axis], fontsize=5.8)
move_ticklabels(a1, "y")        # Move the 4 tick labels to the text box
a1.set_xlim(0, 122); a1.set_xlabel("candidates with a consistent sign (%)", fontsize=6.4)
a1.set_title("Validation cohort (3 lesions)", loc="left", fontsize=7.2)
panel(a1, "A", dx=-0.62)
a2 = f.add_subplot(gs[0, 1])
Rg = []
for a in others:
    d = G[(G.group == "ALL") & (G.axis == a) & G.present]
    if not len(d): continue
    Rg.append(dict(axis=a, n=len(d), rep=int(d.sign_matches_discovery.sum())))
Rg = pd.DataFrame(Rg)
if len(Rg):
    y = np.arange(len(Rg))
    a2.barh(y, 100 * Rg.rep / Rg.n, color="#8c6bb1", height=0.6, lw=0)
    for j, r in Rg.iterrows():
        a2.text(100 * r.rep / r.n + 1.5, j, f"{r.rep}/{r.n}", fontsize=5.8, va="center", bare_keep=True)
    a2.set_yticks(y); a2.set_yticklabels([prog(a) for a in Rg.axis], fontsize=5.8)
    move_ticklabels(a2, "y")    # Also 4 lines
    a2.set_xlim(0, 122); a2.set_xlabel("candidates with the discovery sign (%)", fontsize=6.4)
a2.set_title("Independent cohort GSE164241", loc="left", fontsize=7.2)
panel(a2, "B", dx=-0.62)
save(f, "FigS8_other_programs")
paneln([dict(figure="FigS8", panel="A", content="candidate reproduction, four remaining programmes",
             source_table="05_tf_screen_validation_reproduction.csv", n_cells=1779,
             n_groups=len(Rv), group_unit="programmes", cohort="validation_n3"),
        dict(figure="FigS8", panel="B", content="same, independent cohort",
             source_table="13_gse164241_tf_sign_reproduction.csv", n_cells=19652,
             n_groups=len(Rg), group_unit="programmes", cohort="GSE164241")], "FigS8")
