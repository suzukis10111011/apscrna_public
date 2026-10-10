#!/usr/bin/env python
"""Purpose: Reproduce Figure 7 coupling analysis in independent periodontitis cohort (GSE164241)
Input: {TAB}/13_gse164241_summary.csv, {TAB}/13_gse164241_tf_sign_reproduction.csv
Output: Fig7_independent_cohort

Figure 7: the same coupling in an independent periodontitis cohort (GSE164241).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np

G = pd.read_csv(f"{TAB}/13_gse164241_summary.csv")
S = pd.read_csv(f"{TAB}/13_gse164241_tf_sign_reproduction.csv")
for c in ("present", "all_donors_same_sign", "sign_matches_discovery"):
    if c in S: S[c] = S[c].astype(str).str.lower().isin(["true", "1", "1.0"])
AX2 = ["lf_A_cytokine_chemokine", "lf_D_tissueforming"]
# Keep the groups as 3, as in the original table. If Periodontitis-LAD1 is merged into Periodontitis,
# the plotted points (9 donors) and the pooled line (calculated from 8 donors) become inconsistent.
GRPC = {"Healthy": "#1b7837", "Periodontitis": "#b2182b", "Periodontitis-LAD1": "#8c6bb1"}
GRPN = {"Healthy": "Healthy", "Periodontitis": "Periodontitis",
        "Periodontitis-LAD1": "Periodontitis\nwith LAD1"}

def grp_label(g):
    return str(g)

f = plt.figure(figsize=(180 * MM, 78 * MM))
gs = f.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.40)

# ---- A: per-donor rho, coloured by disease group ----
a1 = f.add_subplot(gs[0, 0])
per = G[(G.donor != "POOLED") & G.depth_controlled_rho.notna() & (G.group != "ALL")].copy()
per["grp"] = per.group.map(grp_label)
groups = [g for g in ["Healthy", "Periodontitis", "Periodontitis-LAD1"] if g in set(per.grp)]
xoff, handles = {}, []
for k, g in enumerate(groups):
    d = per[per.grp == g].sort_values("depth_controlled_rho")
    prev = sum(len(per[per.grp == g2]) for g2 in groups[:k])
    xs = np.arange(len(d)) + prev + 1.4 * k
    a1.scatter(xs, d.depth_controlled_rho, s=16, color=GRPC[g], lw=0, zorder=3)
    xoff[g] = (xs.min(), xs.max())
    pool = G[(G.donor == "POOLED") & (G.group.map(grp_label) == g)]
    if len(pool):
        y = float(pool.depth_controlled_rho.iloc[0])
        a1.hlines(y, xs.min() - 0.45, xs.max() + 0.45, color=GRPC[g], lw=1.4, zorder=4)
        _v = f"{y:+.2f}".replace("-", "\u2212")     # 2 decimal places; minus sign is U+2212
        a1.text(xs.mean(), y + 0.012, ("pooled " if len(d) > 1 else "") + _v,
                ha="center", fontsize=5.4, color=GRPC[g])
    a1.text(xs.mean(), a1.get_ylim()[0], "", ha="center")
    handles.append((g, len(d)))
a1.axhline(0, color="0.5", lw=0.7)
a1.set_xticks([np.mean(xoff[g]) for g in groups])
a1.set_xticklabels([f"{GRPN.get(g, g)}\n({n} donor{'s' if n > 1 else ''})" for g, n in handles],
                   fontsize=5.9)
a1.set_ylabel("Depth-controlled Spearman $\\rho$\n(inflammatory vs tissue-forming)", fontsize=6.3)
a1.set_title("GSE164241 gingiva", loc="left", fontsize=7.0)
panel(a1, "A", dx=-0.20)

# ---- B: candidate sign reproduction, two programmes only ----
a2 = f.add_subplot(gs[0, 1])
d = S[(S.group == "ALL") & S.present & S.axis.isin(AX2)]
rows = [dict(axis=a, n=len(g), rep=int(g.sign_matches_discovery.sum()),
             cons=int(g.all_donors_same_sign.sum())) for a, g in d.groupby("axis")]
R = pd.DataFrame(rows)
y = np.arange(len(R))
a2.barh(y - 0.19, 100 * R.rep / R.n, height=0.36, color="#8c6bb1", lw=0,
        label="same sign as discovery")
a2.barh(y + 0.19, 100 * R.cons / R.n, height=0.36, color="0.62", lw=0,
        label="consistent across all donors")
for j, r in R.iterrows():
    a2.text(100 * r.rep / r.n + 1.5, j - 0.19, f"{r.rep}/{r.n}", fontsize=5.6, va="center", bare_keep=True)
    a2.text(100 * r.cons / r.n + 1.5, j + 0.19, f"{r.cons}/{r.n}", fontsize=5.6, va="center", bare_keep=True)
a2.set_yticks(y); a2.set_yticklabels([prog(a).replace(" programme", "") for a in R.axis], fontsize=6.3)
a2.set_xlim(0, 122); a2.set_xlabel("Discovery candidates reproducing (%)", fontsize=6.4)
a2.legend(fontsize=5.4, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=1, handlelength=0.9)
a2.set_title("Candidate reproduction", loc="left", fontsize=7.0)
panel(a2, "B", dx=-0.34)

save(f, "Fig7_independent_cohort")
paneln([dict(figure="Fig7", panel="A", content="per-donor rho by disease group",
             source_table="13_gse164241_summary.csv", n_cells=int(per.n_cells.sum()),
             n_groups=int(per.donor.nunique()), group_unit="donors", cohort="GSE164241"),
        dict(figure="Fig7", panel="B", content="candidate sign reproduction, two programmes",
             source_table="13_gse164241_tf_sign_reproduction.csv", n_cells=int(per.n_cells.sum()),
             n_groups=len(R), group_unit="programmes", cohort="GSE164241")], "Fig7")
