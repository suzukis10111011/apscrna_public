#!/usr/bin/env python
# Pre-submission redraw (request 5). Copy of scripts/ms_fig4.py with one panel changed.
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

RD = pd.read_csv(f"{TABIN}/04_rho_inflam_vs_tissueforming_discovery.csv")
RV = pd.read_csv(f"{TABIN}/05_rho_inflam_vs_tissueforming_validation.csv")
PN = pd.read_csv(f"{TABIN}/33_defA_permutation_null.csv")
SP = pd.read_csv(f"{TABIN}/33_state_pct_8lesions.csv")
SP = SP[SP.lesion.astype(str).isin(LESION_ORDER)]   # Exclude POOLED rows
SC = pd.read_csv(f"{TABIN}/33_stromal_program_scores.csv")
DP = pd.read_csv(f"{TABIN}/33_dip_by_lesion.csv")

f = plt.figure(figsize=(180 * MM, 196 * MM))
gs = f.add_gridspec(4, 2, height_ratios=[1.0, 1.0, 1.35, 0.8], hspace=0.66, wspace=0.34)

# ---- A: per-lesion depth-controlled rho (composite score) ----
a1 = f.add_subplot(gs[0, 0])
def comp_rows(R, coh):
    d = R[(R.cellset == "stromal") & (R.definition == "composite")]
    return [(str(r.lesion), float(r.depth_controlled_rho), int(r.n)) for _, r in d.iterrows()]
cd = comp_rows(RD, "discovery"); cv = comp_rows(RV, "validation")
per_d = [r for r in cd if not str(r[0]).upper().startswith("POOL")]
per_v = [r for r in cv if not str(r[0]).upper().startswith("POOL")]
pool_d = [r for r in cd if str(r[0]).upper().startswith("POOL")]
pool_v = [r for r in cv if str(r[0]).upper().startswith("POOL")]
for k, (per, pool, col, lab) in enumerate([(per_d, pool_d, DISC, "Discovery"), (per_v, pool_v, VAL, "Validation")]):
    xs = np.full(len(per), k) + np.linspace(-0.13, 0.13, len(per))
    a1.scatter(xs, [p[1] for p in per], s=17, color=col, lw=0, label=f"{lab} lesions")
    for j, (x, p) in enumerate(zip(xs, per)):
        a1.annotate(p[0], (x, p[1]), fontsize=4.8, xytext=(0, 5 if j % 2 == 0 else -9),
                    textcoords="offset points", ha="center", color="0.3")
    if pool:
        a1.scatter([k], [pool[0][1]], s=46, marker="D", facecolor="none", edgecolor=col,
                   linewidth=1.1, label=f"{lab} pooled")
a1.axhline(0, color="0.5", lw=0.7)
a1.set_xticks([0, 1]); a1.set_xticklabels(["Discovery (n=5)", "Validation (n=3)"], fontsize=6.4)
a1.set_ylabel("Depth-controlled Spearman $\\rho$\n(inflammatory vs tissue-forming)", fontsize=6.2, labelpad=1.5)
a1.set_xlim(-0.45, 1.45)
a1.legend(fontsize=5.2, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.16),
          handlelength=0.9, columnspacing=0.7)
panel(a1, "A", dx=-0.24)

# ---- B: retained dispersion with permutation null band ----
# The published panel: bars of the retained dispersion against the within-lesion RELA
# permutation band. One post-audit baseline is added on top of it
# (05_post_audit/tables/Q2_definitionA_baselines.csv) - the retained dispersion when the
# matching uses the mean activity of 10 factors uncorrelated with the tissue-forming
# axis, marked per lesion. The lowest of those (CAP1, 0.883) sits below the published
# y range, so the axis starts lower; no axis break and no floor line. The numbers are in
# the caption and in Tables S3 and S11, not on the bars.
a2 = f.add_subplot(gs[0, 1])
PN["lesion"] = pd.Categorical(PN.lesion, [l for l in LESION_ORDER if l in set(PN.lesion)], ordered=True)
PN = PN.sort_values("lesion")
x = np.arange(len(PN))
BL = pd.read_csv(f"{AUD}/Q2_definitionA_baselines.csv")
nullv = BL[BL.conditioning_group == "null_tf_mean"].set_index("lesion").frac_retained
YLO = 0.875
a2.fill_between([-0.6, len(PN) - 0.4], PN.null_lo.min(), PN.null_hi.max(),
                color="0.86", lw=0, zorder=0, label="permutation null (95%)")
a2.axhline(float(PN.null_mean.mean()), color="0.45", lw=0.7, ls=":", zorder=1, label="null mean")
cols = [DISC if c == "discovery_n5" else VAL for c in PN.cohort]
a2.bar(x, PN.frac_retained, width=0.62, color=cols, lw=0, zorder=2)
a2.scatter(x, [float(nullv[l]) for l in PN.lesion], marker="_", s=90, lw=1.1,
           color="#762a83", zorder=5, label="unrelated factors (mean of 10)")
for j, r in enumerate(PN.itertuples()):
    if r.p_permutation > 0.05:
        a2.text(j, r.null_hi + 0.004, "ns", ha="center", fontsize=5.4, color="#b2182b", bare_keep=True)
a2.set_xticks(x); a2.set_xticklabels(PN.lesion.astype(str), fontsize=6.2)
a2.set_ylim(YLO, 1.022); a2.set_xlim(-0.6, len(PN) - 0.4)
a2.set_ylabel("Retained dispersion (SD ratio)", fontsize=6.2, labelpad=1.5)
a2v = a2.twinx()
a2v.set_ylim(YLO ** 2, 1.022 ** 2)
a2v.set_ylabel("variance ratio", fontsize=6.4); a2v.tick_params(labelsize=5.8)
a2v.spines["right"].set_visible(True)
a2.legend(fontsize=5.0, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3,
          handlelength=0.9, columnspacing=0.7)

panel(a2, "B", dx=-0.22)

# ---- C: state proportions, cohort-internal definition ----
a3 = f.add_subplot(gs[1, :])
SP["lesion"] = pd.Categorical(SP.lesion, [l for l in LESION_ORDER if l in set(SP.lesion)], ordered=True)
SP = SP.sort_values("lesion")
x = np.arange(len(SP)); bottom = np.zeros(len(SP))
for s in ["inflammatory", "intermediate", "tissue-forming", "other"]:
    v = SP[s].values
    a3.bar(x, v, bottom=bottom, width=0.66, color=STATE[s], lw=0, label=s)
    bottom += v
a3.axvline(4.5, color="0.35", lw=0.8, ls="--")
a3.text(2.0, 104, "Discovery", ha="center", fontsize=6.4, color=DISC)
a3.text(6.0, 104, "Validation", ha="center", fontsize=6.4, color=VAL)
a3.set_xticks(x)
a3.set_xticklabels([f"{r.lesion}\n({int(r.n_cells):,})" for r in SP.itertuples()], fontsize=6.0)
a3.set_ylim(0, 110); a3.set_yticks([0, 25, 50, 75, 100])
a3.set_ylabel("Stromal cells (%)")
a3.legend(fontsize=5.6, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.52), handlelength=0.9)
a3.set_title("Cell states", loc="left", fontsize=7.2)
panel(a3, "C", dx=-0.062, dy=1.10)

# ---- D: per-lesion 2D density ----
gsd = gs[2, :].subgridspec(2, 4, hspace=0.50, wspace=0.26)
les = [l for l in LESION_ORDER if l in set(SC.lesion)]
xlim = (np.quantile(SC.inflammatory, 0.002), np.quantile(SC.inflammatory, 0.998))
ylim = (np.quantile(SC.tissueforming, 0.002), np.quantile(SC.tissueforming, 0.998))
for i, l in enumerate(les):
    ax = f.add_subplot(gsd[i // 4, i % 4])
    g = SC[SC.lesion == l]
    ax.hexbin(g.inflammatory, g.tissueforming, gridsize=17, cmap="Greys", mincnt=1, linewidths=0)
    ax.axhline(0, color="0.75", lw=0.5); ax.axvline(0, color="0.75", lw=0.5)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.tick_params(labelsize=5.2)
    col = DISC if g.cohort.iloc[0] == "discovery_n5" else VAL
    ax.set_title(f"{l}  n={len(g):,}", loc="left", fontsize=6.2, color=col)
    if i % 4: ax.set_yticklabels([])
    if i < 4: ax.set_xticklabels([])
    if i == 4: ax.set_xlabel(prog("lf_A_cytokine_chemokine"), fontsize=5.8)
    if i == 0: ax.set_ylabel(prog("lf_D_tissueforming"), fontsize=5.8)
    if i == 0: panel(ax, "D", dx=-0.40, dy=1.30)

# ---- E: dip test ----
a5 = f.add_subplot(gs[3, :])
D2 = DP[DP.lesion != "POOLED"].copy()
D2["lesion"] = pd.Categorical(D2.lesion, les, ordered=True)
D2 = D2.sort_values(["program", "lesion"])
# The value in the table's program column ("Inflammatory program") is the schema, so only align the spelling for display
def _prg(s): return str(s).replace("program", "programme")


for k, prg in enumerate(sorted(D2.program.unique())):
    d = D2[D2.program == prg]
    a5.scatter(np.arange(len(d)) + (k - 0.5) * 0.18, d.dip_p, s=20,
               marker="o" if k == 0 else "^",
               color=["#440154" if k == 0 else "#35b779"], lw=0, label=_prg(prg))
a5.axhline(0.05, color="#b2182b", lw=0.7, ls="--")
a5.text(len(les) - 0.5, 0.075, "p = 0.05", fontsize=5.4, color="#b2182b", ha="right", bare_keep=True)

a5.set_xticks(range(len(les))); a5.set_xticklabels(les, fontsize=6.2)
a5.set_ylim(0, 1.08); a5.set_ylabel("Hartigan dip test, p", fontsize=6.4)
a5.axvline(4.5, color="0.35", lw=0.8, ls="--")
a5.legend(fontsize=5.4, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.40), handlelength=0.9)

panel(a5, "E", dx=-0.062, dy=1.16)

save(f, "Fig4_coupling_coexistence_v3")
paneln([dict(figure="Fig4", panel="A", content="per-lesion depth-controlled rho (composite scores)",
             source_table="04_rho_inflam_vs_tissueforming_discovery.csv; 05_rho_inflam_vs_tissueforming_validation.csv",
             n_cells=int(SC.groupby('cohort').size().sum()), n_groups=len(per_d) + len(per_v),
             group_unit="lesions", cohort="both"),
        dict(figure="Fig4", panel="B", content="retained dispersion vs permutation null, "
                                              "with the unrelated-factor mean marked "
                                              "per lesion",
             source_table="33_defA_permutation_null.csv; "
                          "05_post_audit/tables/Q2_definitionA_baselines.csv",
             n_cells=int(PN.n.sum()), n_groups=len(PN), group_unit="lesions", cohort="both"),
        dict(figure="Fig4", panel="C", content="state proportions, cohort-internal definition",
             source_table="33_state_pct_8lesions.csv", n_cells=int(SP.n_cells.sum()),
             n_groups=len(SP), group_unit="lesions", cohort="both"),
        dict(figure="Fig4", panel="D", content="per-lesion 2D density of the two programme scores",
             source_table="33_stromal_program_scores.csv", n_cells=len(SC),
             n_groups=len(les), group_unit="lesions", cohort="both"),
        dict(figure="Fig4", panel="E", content="Hartigan dip test per lesion and programme",
             source_table="33_dip_by_lesion.csv", n_cells=len(SC), n_groups=len(D2),
             group_unit="lesion x programme", cohort="both")], "Fig4")
