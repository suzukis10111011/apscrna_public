#!/usr/bin/env python
# Pre-submission redraw (request 5). Copy of scripts/ms_fig3.py with one panel changed.
# Reads published tables from results/tables and the post-audit baselines from
# 05_post_audit/tables; writes figures to 06_submission_prep/figures and the panel
# record to 06_submission_prep/tables. The published figure files are not touched.
import os as _o
_o.environ.setdefault("APSC_MSFIGDIR", "../06_submission_prep/figures")
_o.environ.setdefault("APSC_TABDIR", "/path/to/Documents/apscrna_v2/06_submission_prep/tables")
import sys as _s; _s.path.insert(0, "/path/to/Documents/apscrna_v2/scripts")
AUD = "/path/to/Documents/apscrna_v2/05_post_audit/tables"
TABIN = "/path/to/Documents/apscrna_v2/results/tables"
from _msstyle import *
import pandas as pd, numpy as np

A730 = pd.read_csv(f"{TABIN}/33_tf_validation_all730.csv")
# Via CSV, NaN gets mixed into the boolean column, making it object type (~ returns -1/-2 and gets interpreted as a column name).
for _c in ("in_validation", "discovery_is_candidate", "same_sign_pooled"):
    A730[_c] = A730[_c].fillna(False).astype(str).str.lower().isin(["true", "1", "1.0"])
L = pd.read_csv(f"{TABIN}/34_lolo_rulebased.csv")   # Rule-based 10 factors (28 is the old set and does not include RELA/NFKB1)
L["is_candidate"] = L.is_candidate.fillna(False).astype(str).str.lower().isin(["true", "1", "1.0"])
S04 = pd.read_csv(f"{TABIN}/04_tf_screen_discovery_candidates.csv")
S04["is_candidate"] = S04.is_candidate.fillna(False).astype(str).str.lower().isin(["true", "1", "1.0"])
INFL5 = ["RELB", "NFKB2", "RELA", "NFKB1", "REL"]
TISS5 = ["SP7", "RUNX2", "DLX5", "MYRF", "SMAD3"]
KEY10 = INFL5 + TISS5
LABEL = {"lf_D_tissueforming": TISS5, "lf_A_cytokine_chemokine": INFL5}
AX2 = ["lf_D_tissueforming", "lf_A_cytokine_chemokine"]

# Height 224 mm: place c and d side by side (85 mm width each), changing 4 rows -> 3 rows,
# and tightened the margins at the top and between rows. The scatter in a/b and the bars in e/f are unchanged in both dimensions and content.
# Since saving uses bbox_inches='tight' (pad 2.54 mm x2), saving at canvas 224 mm results in approximately 229 mm
f = plt.figure(figsize=(180 * MM, 224 * MM))
gs = f.add_gridspec(3, 2, height_ratios=[1.00, 1.45, 0.72], hspace=0.40, wspace=0.42,
                    top=0.985, bottom=0.070, left=0.125, right=0.945)

# ---- A, B: discovery vs validation pooled rho ----
for k, axis in enumerate(AX2):
    ax = f.add_subplot(gs[0, k])
    d = A730[(A730.axis == axis) & A730.in_validation]
    bg = d[~d.discovery_is_candidate]; cd = d[d.discovery_is_candidate]
    ax.scatter(bg.discovery_pooled_rho, bg.validation_pooled_rho, s=3, c="0.78", lw=0, rasterized=True)
    ax.scatter(cd.discovery_pooled_rho, cd.validation_pooled_rho, s=13, c=DISC,
               edgecolor="white", linewidth=0.3, zorder=3)
    lab = d[d.TF.isin(LABEL[axis])].sort_values("discovery_pooled_rho")
    # Gene name overlaps are resolved with adjustText (fixed offsets would overlap depending on the number of points).
    # Leader lines are attached only to the labeled plate (in the plain plate, they are hidden along with the text).
    from adjustText import adjust_text
    txts = [ax.text(r.discovery_pooled_rho, r.validation_pooled_rho, r.TF, fontsize=7,
                    color="#b2182b", zorder=4, ha="left", va="center")
            for _, r in lab.iterrows()]
    if txts:
        adjust_text(txts, ax=ax, expand=(1.6, 2.0), force_text=(0.7, 1.0),
                    arrowprops=dict(arrowstyle="-", color="#b2182b", lw=0.35, shrinkA=0, shrinkB=1.5))
        if BARE:
            # Since adjust_text adds leader lines itself as a FancyArrowPatch,
            # it does not go through the annotate patch of _msstyle. They are removed here for the plain plate.
            from matplotlib.patches import FancyArrowPatch
            for _ch in ax.get_children():
                if isinstance(_ch, FancyArrowPatch):
                    _ch.set_visible(False)
    m = 1.06 * max(abs(d.discovery_pooled_rho).max(), abs(d.validation_pooled_rho).max())
    ax.plot([-m, m], [-m, m], color="0.6", lw=0.6, ls="--")
    ax.axhline(0, color="0.85", lw=0.6); ax.axvline(0, color="0.85", lw=0.6)
    ax.set_xlim(-m, m); ax.set_ylim(-m, m)
    ax.set_xlabel("Discovery pooled $\\rho$", fontsize=6.6)
    ax.set_ylabel("Validation pooled $\\rho$", fontsize=6.6)
    same = int(cd.same_sign_pooled.sum())
    ax.set_title(prog(axis), loc="left", fontsize=7.0)
    panel(ax, "AB"[k], dx=-0.22)

# ---- C, D: candidate TF activity rho per lesion, both cohorts, one panel per programme ----
DL = [f"rho_{c}" for c in ["A", "B", "C", "D", "E"]]
VL = [f"val_rho_{c}" for c in ["CAP1", "CAP2", "CAP3"]]


def cand_heatmap(ax, axis, letter):
    cand = S04[(S04.axis == axis) & S04.is_candidate].nsmallest(20, "rank_within_axis").TF.tolist()
    Md = S04[(S04.axis == axis) & S04.TF.isin(cand)].set_index("TF").reindex(cand)[DL]
    Mv = A730[(A730.axis == axis) & A730.TF.isin(cand)].set_index("TF").reindex(cand)[VL]
    M = pd.concat([Md, Mv], axis=1)
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels([c.replace("rho_", "").replace("val_", "") for c in M.columns], fontsize=6.2)
    ax.set_yticks(range(len(cand))); ax.set_yticklabels(cand, fontsize=5.6)
    ax.axvline(4.5, color="0.15", lw=1.1)
    ax.text(2.0, -1.1, "Discovery", ha="center", fontsize=6.4, color=DISC)
    ax.text(6.0, -1.1, "Validation", ha="center", fontsize=6.4, color=VAL)
    cb = f.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("depth-controlled $\\rho$", fontsize=6); cb.ax.tick_params(labelsize=5.6)
    ax.set_title(f"Top-20 candidates, {prog(axis)}", loc="left", fontsize=7.2, pad=22)
    panel(ax, letter, dx=-0.075, dy=1.13)
    return M, cand


a3 = f.add_subplot(gs[1, 0])
M, cand = cand_heatmap(a3, AX2[0], "C")
a3b = f.add_subplot(gs[1, 1])
Mi, cand_i = cand_heatmap(a3b, AX2[1], "D")

# ---- D: fraction of candidates reproducing sign in validation ----
a4 = f.add_subplot(gs[2, 0])
rows = []
for axis in AX2:
    d = A730[(A730.axis == axis) & A730.in_validation & A730.discovery_is_candidate]
    rows.append(dict(axis=axis, n=len(d), all3=int((d.n_val_lesions_same_sign == 3).sum()),
                     pooled=int(d.same_sign_pooled.sum())))
R = pd.DataFrame(rows)
y = np.arange(len(R))
a4.barh(y - 0.19, 100 * R.pooled / R.n, height=0.36, color=VAL, lw=0, label="pooled $\\rho$ same sign")
a4.barh(y + 0.19, 100 * R.all3 / R.n, height=0.36, color="#8c6bb1", lw=0, label="same sign in all 3 lesions")
for j, r in R.iterrows():
    a4.text(100 * r.pooled / r.n + 1.5, j - 0.19, f"{r.pooled}/{r.n}", fontsize=5.6, va="center", bare_keep=True)
    a4.text(100 * r.all3 / r.n + 1.5, j + 0.19, f"{r.all3}/{r.n}", fontsize=5.6, va="center", bare_keep=True)
# POST-AUDIT baseline, added on top of the published panel: the mean sign agreement over
# 1,000 draws of 20 non-candidate factors of the same axis that are themselves
# sign-consistent across the five discovery lesions (05_post_audit/tables/
# Q1b_null_summary.csv, as_screened activities, all-3-lesion rule). One dashed line per
# axis; the values are in the caption and in Table S11, not on the figure.
NQ = pd.read_csv(f"{AUD}/Q1b_null_summary.csv")
NQ = NQ[(NQ.version == "as_screened")
        & (NQ.pool == "sign_consistent_non_candidates")].set_index("axis")
for j, r in R.iterrows():
    xn = 100 * float(NQ.loc[r.axis, "null_all3_mean"]) / r.n
    a4.plot([xn, xn], [j + 0.012, j + 0.368], color="0.25", lw=0.9,
            ls=(0, (2.2, 1.3)), zorder=5, solid_capstyle="butt",
            label="random 20 factors, mean" if j == 0 else None)
a4.set_yticks(y); a4.set_yticklabels([prog(a).replace(" programme", "") for a in R.axis], fontsize=6.2)
a4.set_xlim(0, 118); a4.set_xlabel("Candidates reproducing in validation (%)", fontsize=6.4)
a4.legend(fontsize=5.0, loc="upper center", bbox_to_anchor=(0.5, -0.26), ncol=3,
          handlelength=1.2, columnspacing=0.7)
panel(a4, "E", dx=-0.34)

# ---- E: leave-one-lesion-out ----
a5 = f.add_subplot(gs[2, 1])
cnt = L.groupby("TF").is_candidate.sum().reindex(KEY10)
assert cnt.notna().all(), f"LOLO 表に無い因子: {list(cnt[cnt.isna()].index)}"
nr = int(L.held_out_lesion.nunique())
col = ["#1b7837" if v == nr else "#b2182b" for v in cnt]
a5.barh(range(len(cnt)), cnt.values, color=col, height=0.62, lw=0)
a5.invert_yaxis()
a5.set_yticks(range(len(cnt))); a5.set_yticklabels(KEY10, fontsize=6.0)
a5.set_xlim(0, nr + 0.4); a5.set_xticks(range(nr + 1))
a5.axvline(nr, color="0.6", lw=0.7, ls="--")
a5.set_xlabel("Leave-one-lesion-out runs retaining the factor", fontsize=6.4)
n_all = int((cnt == nr).sum())

panel(a5, "F", dx=-0.30)

save(f, "Fig3_tf_screen_v3")
paneln([dict(figure="Fig3", panel="AB", content="discovery vs validation pooled rho, all screened factors",
             source_table="33_tf_validation_all730.csv", n_cells=6215 + 1779,
             n_groups=int(A730[A730.in_validation == True].TF.nunique()), group_unit="TFs", cohort="both"),
        dict(figure="Fig3", panel="C", content=f"per-lesion rho of top-20 candidates, {prog(AX2[0])}",
             source_table="04_tf_screen_discovery_candidates.csv; 33_tf_validation_all730.csv",
             n_cells=6215 + 1779, n_groups=8, group_unit="lesions", cohort="both"),
        dict(figure="Fig3", panel="D", content=f"per-lesion rho of top-20 candidates, {prog(AX2[1])}",
             source_table="04_tf_screen_discovery_candidates.csv; 33_tf_validation_all730.csv",
             n_cells=7994, n_groups=8, group_unit="lesions", cohort="both"),
        dict(figure="Fig3", panel="E", content="candidate sign reproduction, with the "
                                              "random-20-factor null mean marked",
             source_table="33_tf_validation_all730.csv; "
                          "05_post_audit/tables/Q1b_null_summary.csv",
             n_cells=1779, n_groups=2, group_unit="programmes", cohort="validation_n3"),
        dict(figure="Fig3", panel="F", content="leave-one-lesion-out stability",
             source_table="34_lolo_rulebased.csv", n_cells=6215, n_groups=len(KEY10),
             group_unit="factors", cohort="discovery_n5")], "Fig3")
