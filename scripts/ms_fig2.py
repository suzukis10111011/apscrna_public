#!/usr/bin/env python
"""Purpose: Generate Figure 2 visualizations showing tissue-forming programme enrichment in stromal compartment
Input: stromal/epithelial state definitions, enrichment permutation results, UMAP coordinates
  (discovery/validation cohorts)
Output: Fig2_enrichment

Figure 2: the tissue-forming programme concentrates in the stromal compartment.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np

Ud = pd.read_csv(f"{TAB}/33_umap_discovery.csv")
Uv = pd.read_csv(f"{TAB}/33_umap_validation.csv")
PM = pd.read_csv(f"{TAB}/28_enrichment_permutation.csv")
IN = pd.read_csv(f"{TAB}/33_ingest_umap_validation.csv")
Sd = pd.read_csv(f"{TAB}/33_umap_discovery.csv")          # stromal subset taken below
St = pd.read_csv(f"{TAB}/04b_state_definition_discovery.csv") if os.path.exists(f"{TAB}/04b_state_definition_discovery.csv") else None

f = plt.figure(figsize=(180 * MM, 132 * MM))
gs = f.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.25, 1.0],
                    hspace=0.55, wspace=0.34)

# ---- A: tissue-forming score by compartment, both cohorts ----
a1 = f.add_subplot(gs[0, :])
order = ["stromal", "pericyte_smc", "endothelial", "epithelial", "immune"]
pos, labs = [], []
for i, c in enumerate(order):
    for k, (U, coh, col) in enumerate([(Ud, "discovery", DISC), (Uv, "validation", VAL)]):
        v = U[U.compartment == c].tissueforming.values
        if len(v) < 20: continue
        x = i + (-0.19 if k == 0 else 0.19)
        pv = a1.violinplot([v], positions=[x], widths=0.33, showextrema=False, showmedians=True)
        for b in pv["bodies"]:
            b.set_facecolor(col); b.set_alpha(0.62); b.set_edgecolor("none")
        pv["cmedians"].set_color("0.15"); pv["cmedians"].set_linewidth(0.8)
    pos.append(i); labs.append(comp(c))
a1.set_xticks(pos); a1.set_xticklabels(labs, fontsize=6.4)
a1.set_ylabel(prog("lf_D_tissueforming") + "\nscore", fontsize=6.6)
a1.axhline(0, color="0.8", lw=0.6)
from matplotlib.patches import Patch
a1.legend(handles=[Patch(facecolor=DISC, alpha=0.62, label="Discovery (5 lesions)"),
                   Patch(facecolor=VAL, alpha=0.62, label="Validation (3 lesions)")],
          fontsize=5.8, loc="upper right", handlelength=0.9)
a1.set_title("Programme score by compartment", loc="left", fontsize=7.4)
panel(a1, "A", dx=-0.065)

# ---- B: composition of all cells vs top decile, per lesion ----
a2 = f.add_subplot(gs[1, 0])
les = [l for l in LESION_ORDER if l in set(Ud.lesion)] + ["POOLED_discovery"]
x = np.arange(len(les))
st_all = [100 * (Ud[Ud.lesion == l].compartment == "stromal").mean() if l != "POOLED_discovery"
          else 100 * (Ud.compartment == "stromal").mean() for l in les]
st_top, folds = [], []
for l in les:
    r = PM[(PM.lesion == l) & (PM.compartment == "stromal")]
    st_top.append(float(r.pct_top_decile.iloc[0]) if len(r) else np.nan)
    folds.append(float(r.enrichment_fold.iloc[0]) if len(r) else np.nan)
a2.bar(x - 0.19, st_all, width=0.36, color="0.72", lw=0, label="all cells")
a2.bar(x + 0.19, st_top, width=0.36, color=DISC, lw=0, label="top decile of\ntissue-forming score")
a2.set_xticks(x); a2.set_xticklabels([l.replace("POOLED_discovery", "POOL") for l in les], fontsize=6.3)
a2.set_ylabel(comp("stromal") + " (%)", fontsize=6.6); a2.set_ylim(0, 112)
# Move legend outside the axis (below): inside the axis it overlaps the bars and the enrichment line
a2.legend(fontsize=5.6, loc="upper center", bbox_to_anchor=(0.5, -0.30),
          ncol=3, handlelength=0.9, columnspacing=0.8)
a2b = a2.twinx()
a2b.plot(x, folds, "o-", color="#b2182b", lw=0.9, ms=3.4, label="enrichment fold")
a2b.set_ylabel("enrichment (fold)", color="#b2182b", fontsize=6.6)
a2b.tick_params(axis="y", colors="#b2182b", labelsize=6)
a2b.set_ylim(0, max(folds) * 1.25); a2b.spines["right"].set_visible(True)
a2b.spines["right"].set_color("#b2182b")

panel(a2, "B", dx=-0.17)

# ---- C: cohort-internal stromal UMAPs side by side, shared colour scale ----
# The old 2c (ingest-projecting the validation onto the discovery embedding) was moved to FigS7c. The two clusters were cohort boundaries
# themselves, and were an artifact of the projection (diag_fig2c.png / 34_fig2c_lobes.csv).
# Here, the stromal UMAP (Harmony-corrected) computed within each cohort is placed side by side,
# and colored by the tissue-forming score using the same p_D_tissueforming as panel a (to keep the same quantity within the same figure).
SU = pd.read_csv(f"{TAB}/33_stromal_umap_by_cohort.csv")
SCORE = "p_D_tissueforming"
vmin, vmax = float(SU[SCORE].min()), float(SU[SCORE].max())
sub = gs[1, 1].subgridspec(1, 3, width_ratios=[1.0, 1.0, 0.07], wspace=0.12)
axc = []
for k, (coh, lab) in enumerate([("discovery_n5", "Discovery"), ("validation_n3", "Validation")]):
    ax = f.add_subplot(sub[0, k])
    d = SU[SU.cohort == coh]
    sm = ax.scatter(d.umap1, d.umap2, c=d[SCORE], cmap="viridis", vmin=vmin, vmax=vmax,
                    s=1.6, lw=0, alpha=0.7, rasterized=True)
    # Cohort names are placed as "tick labels" (specified to remain even in the plain figure version)
    ax.set_xticks([float(d.umap1.mean())]); ax.set_xticklabels([lab], fontsize=6.6)
    ax.tick_params(axis="x", length=0)
    ax.set_yticks([])
    ax.set_aspect("equal")
    ax.set_xlabel("UMAP 1", fontsize=6.4)
    if k == 0:
        ax.set_ylabel("UMAP 2", fontsize=6.4)
    axc.append((ax, sm, len(d)))
cax = f.add_subplot(sub[0, 2])
cb = f.colorbar(axc[0][1], cax=cax)
# 7 pt は _msstyle の下限。下限未満で置くと save() の引き上げで回転ラベルが伸び、
# tight 保存領域（引き上げ前に決まる）から右にはみ出して文字が切れる。
cb.ax.tick_params(labelsize=7.0)
# ラベルはカラーバーの左（図の内側）に置く。右に置くと tight 保存で端が切れる。
cb.ax.yaxis.set_label_position("left")
cb.set_label("Tissue-forming programme score", fontsize=7.0, labelpad=3)
a3 = axc[0][0]
panel(a3, "C", dx=-0.26)

save(f, "Fig2_enrichment")
paneln([dict(figure="Fig2", panel="A", content="programme score by compartment",
             source_table="33_umap_discovery.csv; 33_umap_validation.csv",
             n_cells=len(Ud) + len(Uv), n_groups=len(order), group_unit="compartments", cohort="both"),
        dict(figure="Fig2", panel="B", content="stromal share of all cells vs top decile + fold",
             source_table="28_enrichment_permutation.csv", n_cells=len(Ud),
             n_groups=len(les), group_unit="lesions (+pooled)", cohort="discovery_n5"),
        dict(figure="Fig2", panel="C", content="cohort-internal stromal UMAPs coloured by tissue-forming score",
             source_table="33_stromal_umap_by_cohort.csv", n_cells=int(len(SU)),
             n_groups=2, group_unit="cohorts", cohort="both")], "Fig2")
