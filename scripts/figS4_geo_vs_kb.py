#!/usr/bin/env python
"""Purpose: Compare GEO matrix versus kb requantification robustness across stromal Jaccard and TF agreement
  metrics
Input: {TAB}/09_geo_vs_kb.csv, {TAB}/09_stromal_jaccard.csv, {TAB}/09_tf_agreement_geo_vs_kb.csv
Output: figS4_geo_vs_kb, {TAB}/_paneln_figS4.csv

Figure S4: GEO matrix version vs kb requantification version (supplementary robustness check).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figstyle import *
import pandas as pd, numpy as np
P = pd.read_csv(f"{TAB}/09_geo_vs_kb.csv")
J = pd.read_csv(f"{TAB}/09_stromal_jaccard.csv")
T = pd.read_csv(f"{TAB}/09_tf_agreement_geo_vs_kb.csv")
KB = [c for c in ("kb_allbarcodes", "kb_geobarcodes") if c in P.columns]
CK = {"kb_allbarcodes": "#9970ab", "kb_geobarcodes": "#1b7837"}
f, gs = gridfig(180, 78, nrows=1, ncols=3, width_ratios=[1, 1, 1], wspace=0.52)
a1 = f.add_subplot(gs[0, 0])
E = P[P.quantity == "enrichment_fold_tissueforming"].dropna(subset=["geo"])
x = np.arange(len(E))
a1.bar(x - 0.22, E.geo, width=0.2, color=DISC, lw=0, label="GEO matrix")
for k, c in enumerate(KB):
    a1.bar(x + 0.02 + k * 0.22, pd.to_numeric(E[c], errors="coerce"), width=0.2, color=CK[c], lw=0,
           label=c.replace("kb_", "kb, ").replace("allbarcodes", "all barcodes").replace("geobarcodes", "GEO barcodes"))
a1.set_xticks(x); a1.set_xticklabels([str(l).replace("POOLED_discovery", "POOL") for l in E.lesion], fontsize=5.8)
a1.set_ylabel("enrichment fold\n(tissue-forming top decile)", fontsize=6.4)
a1.legend(fontsize=5.1, frameon=False, handlelength=0.9, loc="upper left")
a1.set_title("Enrichment", loc="left", fontsize=7.4)
a2 = f.add_subplot(gs[0, 1])
from matplotlib.patches import Patch
cmps = list(J.comparison.unique())
les_b = list(dict.fromkeys(J.lesion))
for k, cmp_ in enumerate(cmps):
    d = J[J.comparison == cmp_].set_index("lesion").reindex(les_b)
    col = CK["kb_geobarcodes"] if "geobarcodes" in cmp_ else CK["kb_allbarcodes"]
    a2.bar(np.arange(len(les_b)) + (k - 0.5) * 0.38, d.jaccard.values, width=0.36, color=col, lw=0)
a2.set_xticks(range(len(les_b)))
a2.set_xticklabels([str(l).replace("POOLED", "POOL") for l in les_b], fontsize=6.0)
a2.set_ylim(0, 1.12); a2.set_ylabel("Jaccard of stromal barcodes", fontsize=6.4)
a2.legend(handles=[Patch(facecolor=CK["kb_allbarcodes"], label="vs kb, all barcodes"),
                   Patch(facecolor=CK["kb_geobarcodes"], label="vs kb, GEO barcodes")],
          fontsize=5.2, frameon=False, handlelength=0.9, loc="upper center",
          bbox_to_anchor=(0.5, -0.16), ncol=1)
a2.set_title("Stromal-set agreement", loc="left", fontsize=7.4)
a3 = f.add_subplot(gs[0, 2])
T2 = T[T.axis.isin(["lf_A_cytokine_chemokine", "lf_D_tissueforming"])]
y = np.arange(len(T2))
a3.barh(y, T2.same_sign_pct, color=[CK.get(s, "0.6") for s in T2.source], height=0.6, lw=0)
a3.set_yticks(y); a3.set_yticklabels([f"{r.axis.replace('lf_','')}\n{r.source.replace('kb_allbarcodes','kb all BC').replace('kb_geobarcodes','kb GEO BC')}" for r in T2.itertuples()], fontsize=5.2)
a3.set_xlim(0, 105); a3.set_xlabel("TFs with same sign of pooled $\\rho$ (%)", fontsize=6.4)
for j, r in enumerate(T2.itertuples()):
    a3.text(r.same_sign_pct + 1.2, j, f"{int(r.overlap)}/{int(r.n_candidates_geo)} cand.", fontsize=5, va="center")
a3.set_title("TF screen agreement", loc="left", fontsize=7.4)
panel(a1, "A", dx=-0.30); panel(a2, "B", dx=-0.26); panel(a3, "C", dx=-0.52)
save(f, "figS4_geo_vs_kb")
pd.DataFrame([dict(figure="figS4_geo_vs_kb", panel="A", n_lesions=len(E)),
              dict(figure="figS4_geo_vs_kb", panel="B", n_comparisons=J.comparison.nunique(),
                   n_geo_stromal=int(J[J.lesion == "POOLED"].n_geo.iloc[0]) if (J.lesion == "POOLED").any() else None),
              dict(figure="figS4_geo_vs_kb", panel="C", n_axes=T2.axis.nunique(), n_TFs_common=int(T2.n_TF_common.max()))
              ]).to_csv(f"{TAB}/_paneln_figS4.csv", index=False)
