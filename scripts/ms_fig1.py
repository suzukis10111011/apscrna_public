#!/usr/bin/env python
"""Purpose: Generate Figure 1 showing study design and cohort composition across discovery and validation
  cohorts
Input: 03c_compartment_by_lesion_discovery.csv, 05_compartment_by_lesion_validation.csv,
  33_compartment_marker_dotplot.csv, 33_umap_discovery.csv, 33_umap_validation.csv
Output: Fig1_design_composition

Figure 1: study design and cohort composition.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np
import matplotlib.patches as mp

Ud = pd.read_csv(f"{TAB}/33_umap_discovery.csv")
Uv = pd.read_csv(f"{TAB}/33_umap_validation.csv")
DP = pd.read_csv(f"{TAB}/33_compartment_marker_dotplot.csv")
Cd = pd.read_csv(f"{TAB}/03c_compartment_by_lesion_discovery.csv", index_col=0)
Cv = pd.read_csv(f"{TAB}/05_compartment_by_lesion_validation.csv", index_col=0)

# Height 240 mm: if the legend is arranged in a single column, c becomes 6 rows (including Schwann / neural).
# Dimension that ensures an absolute margin so that these 6 rows do not encroach on e (row ratios and order unchanged)
# Height 226 mm: tightened the margins at the top edge and between rows (row ratios and panel layout unchanged)
f = plt.figure(figsize=(180 * MM, 226 * MM))
gs = f.add_gridspec(3, 3, height_ratios=[0.72, 1.0, 0.9], width_ratios=[1, 1, 1],
                    hspace=0.46, wspace=0.42, top=0.965, bottom=0.075,
                    left=0.085, right=0.965)
a1 = f.add_subplot(gs[0, :])
a1.set_xlim(0, 10); a1.set_ylim(-0.35, 2.3); a1.axis("off")
def box(x, y, w, h, txt, fc, ec):
    a1.add_patch(mp.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06",
                                   linewidth=0.7, facecolor=fc, edgecolor=ec))
    a1.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=6.4)
box(0.1, 1.28, 2.0, 0.72, "Discovery\n5 lesions", "#dbe6f2", DISC)
box(2.6, 1.28, 2.3, 0.72, "Programme scores\nTF activity screen", "#dbe6f2", DISC)
box(5.4, 1.28, 2.2, 0.72, "Candidate factors\nand thresholds", "#dbe6f2", DISC)
box(0.1, 0.16, 2.0, 0.72, "Validation\n3 lesions", "#fbe4cd", VAL)
box(2.6, 0.16, 2.3, 0.72, "Same pipeline,\ncohort-internal cuts", "#fbe4cd", VAL)
box(5.4, 0.16, 2.2, 0.72, "Sign confirmation\nonly", "#fbe4cd", VAL)
box(8.1, 0.72, 1.8, 0.72, "Independent\ncohort\n(GSE164241)", "#f0e6f5", "#8c6bb1")
for x0, y0 in [(2.1, 1.64), (4.9, 1.64), (2.1, 0.52), (4.9, 0.52)]:
    a1.annotate("", xy=(x0 + 0.5, y0), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color="0.35"))
a1.annotate("", xy=(3.7, 1.24), xytext=(3.7, 0.92),
            arrowprops=dict(arrowstyle="-|>", lw=0.9, color="0.2"))

a1.annotate("", xy=(8.05, 1.08), xytext=(7.65, 1.08),
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color="0.35"))
panel(a1, "A", dx=-0.03, dy=1.16)

# ---- B, C: all-cell UMAPs ----
for k, (U, lab, ax) in enumerate([(Ud, "Discovery, 5 lesions", f.add_subplot(gs[1, 0])),
                                  (Uv, "Validation, 3 lesions", f.add_subplot(gs[1, 1]))]):
    for c, g in U.groupby("compartment"):
        ax.scatter(g.umap1, g.umap2, s=0.8, c=COMPC.get(c, "0.6"), lw=0,
                   rasterized=True, label=f"{comp(c)} ({len(g):,})")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1", fontsize=6.4); ax.set_ylabel("UMAP 2", fontsize=6.4)
    ax.set_title(f"{lab}  n={len(U):,}", loc="left", fontsize=7.2)
    ax.legend(fontsize=4.9, markerscale=6, handletextpad=0.2, borderpad=0.2,
              loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1,
              columnspacing=0.7, handlelength=0.8)
    panel(ax, "BC"[k], dx=-0.13)

# ---- D: marker dot plot ----
a4 = f.add_subplot(gs[1, 2])
d = DP[DP.cohort == "discovery"].copy()
genes = list(dict.fromkeys(d.gene))
comps = [c for c in ["stromal", "immune", "endothelial", "pericyte_smc", "epithelial"] if c in set(d.compartment)]
for i, c in enumerate(comps):
    for j, g in enumerate(genes):
        r = d[(d.compartment == c) & (d.gene == g)]
        if not len(r): continue
        a4.scatter(i, j, s=2 + 0.42 * float(r.detection_pct.iloc[0]),
                   c=[float(r.mean_lognorm.iloc[0])], cmap="viridis",
                   vmin=0, vmax=float(d.mean_lognorm.max()), edgecolor="none")
a4.set_xticks(range(len(comps))); a4.set_xticklabels([comp(c) for c in comps], rotation=40,
                                                     ha="right", fontsize=5.6)
a4.set_yticks(range(len(genes))); a4.set_yticklabels(genes, fontsize=5.6, style="italic")
a4.set_xlim(-0.6, len(comps) - 0.4); a4.set_ylim(-0.8, len(genes) - 0.2)
sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, float(d.mean_lognorm.max())))
cb = f.colorbar(sm, ax=a4, fraction=0.045, pad=0.04); cb.set_label("mean log-normalized", fontsize=7.0)
cb.ax.tick_params(labelsize=5.4)
a4.set_title("Compartment markers", loc="left", fontsize=6.8)
panel(a4, "D", dx=-0.42)

# ---- E: compartment fractions per lesion ----
a5 = f.add_subplot(gs[2, :])
pcd = [c for c in Cd.columns if c.startswith("pct_")]
pcv = [c for c in Cv.columns if c.startswith("pct_")]
order = [c[4:] for c in pcd] + [c[4:] for c in pcv if c[4:] not in [x[4:] for x in pcd]]
rows = []
for les in Cd.index: rows.append(("discovery", str(les), {c[4:]: float(Cd.loc[les, c]) for c in pcd}))
for les in Cv.index: rows.append(("validation", str(les), {c[4:]: float(Cv.loc[les, c]) for c in pcv}))
x = np.arange(len(rows))
bottom = np.zeros(len(rows))
for c in order:
    v = np.array([r[2].get(c, 0.0) for r in rows])
    a5.bar(x, v, bottom=bottom, width=0.68, color=COMPC.get(c, "0.6"), lw=0, label=comp(c))
    bottom += v
a5.axvline(4.5, color="0.35", lw=0.8, ls="--")
a5.text(2.0, 104, "Discovery", ha="center", fontsize=6.6, color=DISC)
a5.text(6.0, 104, "Validation", ha="center", fontsize=6.6, color=VAL)
a5.set_xticks(x); a5.set_xticklabels([r[1] for r in rows], fontsize=6.4)
a5.set_ylim(0, 110); a5.set_yticks([0, 25, 50, 75, 100])
a5.set_ylabel("Fraction of cells (%)")
a5.legend(fontsize=5.4, ncol=6, loc="lower center", bbox_to_anchor=(0.5, -0.30),
          handlelength=0.9, columnspacing=0.8)
panel(a5, "E", dx=-0.055, dy=1.12)

save(f, "Fig1_design_composition")
paneln([dict(figure="Fig1", panel="A", content="study design schematic", source_table="-",
             n_cells=None, n_groups=8, group_unit="lesions", cohort="both"),
        dict(figure="Fig1", panel="B", content="all-cell UMAP", source_table="33_umap_discovery.csv",
             n_cells=len(Ud), n_groups=Ud.lesion.nunique(), group_unit="lesions", cohort="discovery_n5"),
        dict(figure="Fig1", panel="C", content="all-cell UMAP", source_table="33_umap_validation.csv",
             n_cells=len(Uv), n_groups=Uv.lesion.nunique(), group_unit="lesions", cohort="validation_n3"),
        dict(figure="Fig1", panel="D", content="compartment marker dot plot",
             source_table="33_compartment_marker_dotplot.csv", n_cells=len(Ud),
             n_groups=len(comps), group_unit="compartments", cohort="discovery_n5"),
        dict(figure="Fig1", panel="E", content="compartment fractions per lesion",
             source_table="03c_compartment_by_lesion_discovery.csv; 05_compartment_by_lesion_validation.csv",
             n_cells=len(Ud) + len(Uv), n_groups=len(rows), group_unit="lesions", cohort="both")],
       "Fig1")
