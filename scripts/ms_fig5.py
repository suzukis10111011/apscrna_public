#!/usr/bin/env python
"""Purpose: Generate Figure 5 panels comparing receptor states and TF coactivity with SMAD3 sensitivity
  analysis
Input: {TAB}/10_receptors_by_state.csv, {TAB}/10_receptors_dotplot_input.csv,
  {TAB}/32_coactivity_ruleset_sensitivity.csv, {TAB}/32_tf_activity_spearman_matrix.csv
Output: Fig5_receptors_coactivity

Fig. 5: A/B comparison of receptor states, C/D bilateral TF coactivity (SMAD3 sensitivity).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np
R = pd.read_csv(f"{TAB}/10_receptors_by_state.csv")
Dp = pd.read_csv(f"{TAB}/10_receptors_dotplot_input.csv")
C = pd.read_csv(f"{TAB}/32_coactivity_ruleset_sensitivity.csv")   # 32: Sets defined by the rules
M = pd.read_csv(f"{TAB}/32_tf_activity_spearman_matrix.csv")
f, gs = gridfig(180, 150, nrows=2, ncols=2, hspace=0.62, wspace=0.42,
                height_ratios=[1, 1.08])
a1 = f.add_subplot(gs[0, 0])
d = R[R.cohort == "discovery_n5"].sort_values("delta_infl_minus_tiss")
col = [DISC if p < 0.1 else "0.72" for p in d.pvalue.fillna(1)]
a1.barh(range(len(d)), d.delta_infl_minus_tiss, color=col, height=0.66, lw=0)
a1.set_yticks(range(len(d))); a1.set_yticklabels(d.gene, fontsize=6.2, style="italic")
a1.axvline(0, color="0.5", lw=0.6)
a1.set_xlabel("$\\Delta$ log2CPM (inflammatory $-$ tissue-forming)")
a1.set_title("Discovery, pseudobulk (5 lesions)", loc="left", fontsize=7.0)
a2 = f.add_subplot(gs[0, 1])
dd = Dp[Dp.cohort == "discovery_n5"].groupby(["gene", "state"], as_index=False).agg(
    pct=("pct_expressing", "mean"), ex=("mean_expr_log2cpm", "mean"))
genes = list(d.gene)[::-1]; states = ["inflammatory", "tissue-forming"]
for i, st in enumerate(states):
    s = dd[dd.state == st].set_index("gene").reindex(genes)
    sizes = np.clip(s.pct.values, 0, 100) * 0.9 + 2
    sc_ = a2.scatter([i] * len(genes), range(len(genes)), s=sizes, c=s.ex.values,
                     cmap="viridis", vmin=0, vmax=float(np.nanmax(dd.ex)), linewidths=0.3, edgecolors="white")
a2.set_xticks(range(len(states))); a2.set_xticklabels(["inflammatory", "tissue-forming"], fontsize=6.2)
a2.set_yticks(range(len(genes))); a2.set_yticklabels(genes, fontsize=6.2, style="italic")
a2.set_xlim(-0.6, 1.6)
cb = f.colorbar(sc_, ax=a2, fraction=0.05, pad=0.04); cb.set_label("mean log2CPM", fontsize=6.2); cb.ax.tick_params(labelsize=6)

a3 = f.add_subplot(gs[1, 0])
# 32: On the inflammation side, inflammation-axis candidates x NF-kB family; on the tissue-formation side, tissue-formation-axis candidates x osteoblast lineage.
# Also list the old set (including ATF1/MTA1) alongside, to show that the conclusion does not depend on the choice of set.
cc = C[(C.group == "lesion") & (C.label != "POOLED")]
ISETS = [("infl_new_NFkB", "NF-\u03baB family set"), ("infl_old", "alternative set (ATF1, MTA1)")]
for k, (coh, col2) in enumerate([("discovery", DISC), ("validation", VAL)]):
    for j, (iset, ilab) in enumerate(ISETS):
        for v, mk, lab in [("tiss_new5_incl_SMAD3", "o", "with SMAD3"),
                           ("tiss_new4_excl_SMAD3", "^", "without SMAD3")]:
            s = cc[(cc.cohort == coh) & (cc.inflammatory_set == iset) & (cc.tissueforming_set == v)]
            a3.scatter([k + (-0.24 + 0.16 * j) + (0.16 if "excl" in v else 0.0)] * len(s), s.spearman_rho,
                       marker=mk, s=15, color=col2, alpha=0.55 if j else 0.9, linewidths=0,
                       label=f"{ilab} / {lab}" if k == 0 else None)
a3.axhline(0, color="0.5", lw=0.6)
a3.set_xticks([0, 1]); a3.set_xticklabels(["Discovery (n=5)", "Validation (n=3)"], fontsize=6.5)
a3.set_ylabel("Per-lesion Spearman $\\rho$\n(inflammatory vs tissue-forming TF means)", fontsize=6.4)
a3.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.36), handlelength=1.0,
          columnspacing=0.8, fontsize=5.2, frameon=False)
a4 = f.add_subplot(gs[1, 1])
md = M[M.cohort == "discovery"].set_index("TF").drop(columns=["cohort"])
# Sort in set order (alphabetical order would split the NF-kB group and osteoblast group blocks)
ORDER = ["RELB", "NFKB2", "RELA", "NFKB1", "REL", "ATF1", "MTA1",
         "SP7", "RUNX2", "DLX5", "MYRF", "SMAD3", "MYC"]
ordr = [g for g in ORDER if g in md.index] + [g for g in md.index if g not in ORDER]
md = md.loc[ordr, ordr]
im = a4.imshow(md.values, cmap="RdBu_r", vmin=-0.8, vmax=0.8)
a4.set_xticks(range(len(md))); a4.set_xticklabels(md.index, rotation=90, fontsize=5.8)
a4.axvline(6.5, color="0.25", lw=0.9); a4.axhline(6.5, color="0.25", lw=0.9)
a4.set_yticks(range(len(md))); a4.set_yticklabels(md.index, fontsize=5.8)
cb2 = f.colorbar(im, ax=a4, fraction=0.046, pad=0.04); cb2.set_label("Spearman ρ", fontsize=6.2); cb2.ax.tick_params(labelsize=6)
a4.set_title("TF activity correlation (discovery)", loc="left", fontsize=7.0)
for a, L, dx in [(a1, "A", -0.30), (a2, "B", -0.34), (a3, "C", -0.20), (a4, "D", -0.26)]:
    panel(a, L, dx=dx)
save(f, "Fig5_receptors_coactivity")

paneln([dict(figure="Fig5", panel="A", content="receptor difference between states (pseudobulk)",
             source_table="10_receptors_by_state.csv", n_cells=6215, n_groups=int(R[R.cohort=="discovery_n5"].gene.nunique()),
             group_unit="receptor genes", cohort="discovery_n5"),
        dict(figure="Fig5", panel="B", content="receptor dot plot by state",
             source_table="10_receptors_dotplot_input.csv", n_cells=6215,
             n_groups=int(Dp.gene.nunique()) if "gene" in Dp.columns else None,
             group_unit="receptor genes", cohort="discovery_n5"),
        dict(figure="Fig5", panel="C", content="per-lesion co-activity, rule-based and alternative sets",
             source_table="32_coactivity_ruleset_sensitivity.csv", n_cells=7994, n_groups=8,
             group_unit="lesions", cohort="both"),
        dict(figure="Fig5", panel="D", content="TF activity correlation matrix",
             source_table="32_tf_activity_spearman_matrix.csv", n_cells=6215,
             n_groups=int(md.shape[0]), group_unit="TFs", cohort="discovery_n5")], "Fig5")
