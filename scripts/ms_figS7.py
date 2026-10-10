#!/usr/bin/env python
"""Purpose: Generate Supplementary Figure S7 comparing within-cohort state labels to cross-cohort UMAP
  projections
Input: {TAB}/25_state_crosstab_within_vs_ingest.csv, {TAB}/25_validation_states_within_vs_ingest.csv,
  {TAB}/33_ingest_umap_validation.csv, {TAB}/33_umap_discovery.csv
Output: FigS7_state_labelling

Supplementary Figure S7: cohort-internal state labels vs cross-cohort projection.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np
V = pd.read_csv(f"{TAB}/25_validation_states_within_vs_ingest.csv")
CT = pd.read_csv(f"{TAB}/25_state_crosstab_within_vs_ingest.csv", index_col=0)
STATES = ["inflammatory", "intermediate", "tissue-forming", "other"]
les = [l for l in LESION_ORDER if l in set(V.lesion)]
# Height 82 -> 150 mm: because c (former Fig2c) was added to the bottom row. Width remains 180 mm
f, gs = gridfig(180, 168, nrows=2, ncols=2, width_ratios=[1.3, 1.0],
                height_ratios=[1.0, 1.45], wspace=0.42, hspace=0.40)
a1 = f.add_subplot(gs[0, 0])
w = 0.38
for k, (lab, key) in enumerate([("within cohort", "within_cohort"), ("projected", "ingest_projected")]):
    d = V[V.labelling == key].set_index("lesion").reindex(les)
    bottom = np.zeros(len(les))
    for s in STATES:
        v = d[f"pct_{s}"].values.astype(float)
        a1.bar(np.arange(len(les)) + (k - 0.5) * w, v, bottom=bottom, width=w * 0.9,
               color=STATE[s], lw=0, label=s if k == 0 else None)
        bottom += v
for j, l in enumerate(les):
    a1.text(j - 0.5 * w, 102, "within", ha="center", fontsize=4.6, rotation=90, color="0.35")
    a1.text(j + 0.5 * w, 102, "projected", ha="center", fontsize=4.6, rotation=90, color="0.35")
    ag = V[(V.labelling == "within_cohort") & (V.lesion == l)].agreement_pct
    if len(ag) and pd.notna(ag.iloc[0]):
        a1.text(j, -7, f"{float(ag.iloc[0]):.1f}% agree", ha="center", fontsize=5.0, color="#b2182b")
a1.set_xticks(range(len(les))); a1.set_xticklabels(les, fontsize=6.4)
a1.set_ylim(-10, 122); a1.set_yticks([0, 25, 50, 75, 100])
a1.set_ylabel(comp("stromal") + " (%)")
a1.legend(fontsize=5.4, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.30), handlelength=0.9)

panel(a1, "A", dx=-0.16)
a2 = f.add_subplot(gs[0, 1])
M = CT.copy()
CM2 = "Blues"
im = a2.imshow(M.values, cmap=CM2, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        a2.text(j, i, f"{int(M.values[i, j]):,}", ha="center", va="center", fontsize=7,
                color=cell_text_colour(CM2, M.values[i, j], M.values.min(), M.values.max()),
                numbox=True)
a2.set_xticks(range(M.shape[1])); a2.set_xticklabels(M.columns, rotation=40, ha="right", fontsize=5.8)
a2.set_yticks(range(M.shape[0])); a2.set_yticklabels(M.index, fontsize=5.8)
a2.set_xlabel("projected label", fontsize=6.4); a2.set_ylabel("cohort-internal label", fontsize=6.4)
a2.set_title("Cell-level confusion", loc="left", fontsize=7.2)
panel(a2, "B", dx=-0.34)
# ---- C: former Fig2c — validation projected (ingested) into the discovery embedding (colored by cohort) ----
# The two clusters are exactly the cohort boundary itself; kept to show that this is an artifact of the projection
# (Diagnostics are in results/figures_ms/diag_fig2c.png and results/tables/34_fig2c_lobes.csv).
Ud = pd.read_csv(f"{TAB}/33_umap_discovery.csv")
IN = pd.read_csv(f"{TAB}/33_ingest_umap_validation.csv")
sd = Ud[Ud.compartment == "stromal"]
a3 = f.add_subplot(gs[1, :])
a3.scatter(sd.umap1, sd.umap2, s=1.6, c="#DDDDDD", lw=0, alpha=0.9, rasterized=True,
           zorder=1, label=f"Discovery ({len(sd):,})")
a3.scatter(IN.umap1, IN.umap2, s=2.4, c="#E69F00", lw=0, alpha=0.85, rasterized=True,
           zorder=2, label=f"Validation ({len(IN):,})")
a3.set_aspect("equal")
a3.set_xticks([]); a3.set_yticks([])
a3.set_xlabel("UMAP 1", fontsize=6.4); a3.set_ylabel("UMAP 2", fontsize=6.4)
a3.legend(fontsize=6.0, markerscale=4, loc="lower left", handletextpad=0.3)
a3.set_title("Validation projected onto the discovery embedding", loc="left", fontsize=7.2)
panel(a3, "C", dx=-0.09)

save(f, "FigS7_state_labelling")
paneln([dict(figure="FigS7", panel="A", content="state proportions by labelling route",
             source_table="25_validation_states_within_vs_ingest.csv", n_cells=1779,
             n_groups=len(les), group_unit="lesions", cohort="validation_n3"),
        dict(figure="FigS7", panel="B", content="cell-level confusion between routes",
             source_table="25_state_crosstab_within_vs_ingest.csv", n_cells=1779,
             n_groups=int(M.size), group_unit="label pairs", cohort="validation_n3"),
        dict(figure="FigS7", panel="C", content="validation cells projected onto the discovery embedding",
             source_table="33_umap_discovery.csv; 33_ingest_umap_validation.csv",
             n_cells=len(sd) + len(IN), n_groups=2, group_unit="cohorts", cohort="both")], "FigS7")
