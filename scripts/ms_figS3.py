#!/usr/bin/env python
"""Purpose: Supplementary Figure S3 - gene-set and scoring-method robustness as two heatmaps (6 gene-set pairings x 3 scoring methods).
Input: {TAB}/14_robustness_genesets_scoring.csv
Output: FigS3_genesets_scoring, {TAB}/_figS3_setcodes.csv

Supplementary Figure S3: gene-set and scoring-method robustness, as two heatmaps.

We abandoned 18 vertical bars (if set names are not abbreviated, they become 5-line labels that collapse at 7 pt),
and instead made two heatmap panels with rows = 3 inflammation-side sets x 2 tissue-formation-side sets (6 rows), columns = score_genes / ULM / AUCell (3 columns).
Row labels are abbreviation codes, and the correspondence table is also listed in bare/colour_keys.csv.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _msstyle import *
import pandas as pd, numpy as np

B = pd.read_csv(f"{TAB}/14_robustness_genesets_scoring.csv")
P = B[B.lesion == "POOLED"].copy()

INF = [("current_lf_A_cytokine_chemokine", "I"), ("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "H"),
       ("GOBP_INFLAMMATORY_RESPONSE", "G")]
TIS = [("current_lf_D_tissueforming", "T"), ("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "R")]
SCO = [("score_genes", "score_genes"), ("ulm", "ULM"), ("aucell", "AUCell")]
ROWS = [(i, t) for i, _ in INF for t, _ in TIS]
RLAB = [f"{ic} \u00d7 {tc}" for _, ic in INF for _, tc in TIS]
assert len(P) == len(ROWS) * len(SCO), (len(P), len(ROWS), len(SCO))


def grid(col):
    A = np.full((len(ROWS), len(SCO)), np.nan)
    for r, (i, t) in enumerate(ROWS):
        for c, (s, _) in enumerate(SCO):
            m = P[(P.inflammatory_set == i) & (P.tissueforming_set == t) & (P.scoring == s)]
            assert len(m) == 1, (i, t, s, len(m))
            A[r, c] = float(m[col].iloc[0])
    return A


RHO, DSP = grid("depth_controlled_rho"), grid("defA_frac_retained")
# wspace 0.72: spacing so that the colorbar label of a and the y-axis label of b do not overlap
f, gs = gridfig(180, 86, nrows=1, ncols=2, wspace=0.72)


def heat(ax, A, cmap, vmin, vmax, fmt, cblab):
    im = ax.imshow(A, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    for r in range(A.shape[0]):
        for c in range(A.shape[1]):
            ax.text(c, r, format(A[r, c], fmt), ha="center", va="center", fontsize=7,
                    color=cell_text_colour(cmap, A[r, c], vmin, vmax), bare_keep=True)
    ax.set_xticks(range(len(SCO))); ax.set_xticklabels([s for _, s in SCO], fontsize=7)
    ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(RLAB, fontsize=7)
    ax.set_xlabel("scoring method", fontsize=7)
    ax.set_ylabel("inflammatory \u00d7 tissue-forming set", fontsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    cb = f.colorbar(im, ax=ax, fraction=0.048, pad=0.04)
    cb.set_label(cblab, fontsize=7); cb.ax.tick_params(labelsize=7)
    return im


a1 = f.add_subplot(gs[0, 0])
m = float(np.nanmax(np.abs(RHO)))
heat(a1, RHO, "RdBu_r", -m, m, "+.3f", "depth-controlled Spearman $\\rho$")
a2 = f.add_subplot(gs[0, 1])
heat(a2, DSP, "Blues", float(np.nanmin(DSP)), 1.0, ".4f", "retained dispersion (SD ratio)")

neg = int((RHO < 0).sum())
panel(a1, "A", dx=-0.30); panel(a2, "B", dx=-0.30)
save(f, "FigS3_genesets_scoring")
paneln([dict(figure="FigS3", panel="a", content="depth-controlled rho, gene-set x scoring",
             source_table="14_robustness_genesets_scoring.csv", n_cells=6215,
             n_groups=RHO.size, group_unit="combinations", cohort="discovery_n5"),
        dict(figure="FigS3", panel="b", content="retained dispersion, same combinations",
             source_table="14_robustness_genesets_scoring.csv", n_cells=6215,
             n_groups=DSP.size, group_unit="combinations", cohort="discovery_n5")], "FigS3")
print(f"rho: {neg}/{RHO.size} が負 | 範囲 {RHO.min():+.4f} .. {RHO.max():+.4f}")
print(f"retained dispersion: {DSP.min():.4f} .. {DSP.max():.4f}")
# Correspondence table of abbreviation codes (material also listed in colour_keys.csv)
pd.DataFrame([dict(code=c, gene_set=k, side="inflammatory") for k, c in INF] +
             [dict(code=c, gene_set=k, side="tissue-forming") for k, c in TIS]
             ).to_csv(f"{TAB}/_figS3_setcodes.csv", index=False)
