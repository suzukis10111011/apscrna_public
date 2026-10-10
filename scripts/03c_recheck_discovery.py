#!/usr/bin/env python
"""Purpose: Recheck discovery-side compartment annotations using per-cell score_genes and cluster-wise majority
  vote; verify stromal ratios and enrichment fold-changes against 03b.
Input: {OUT}/all_cells_discovery.h5ad
Output: {OUT}/all_cells_discovery.h5ad, {TAB}/03c_enrichment_by_lesion_discovery.csv,
  {TAB}/03c_vs_03b_confusion_discovery.csv, {TAB}/03c_vs_03b_counts_discovery.csv

03c: Recheck the discovery-side compartment annotations using the common function (per-cell score_genes + cluster-wise majority vote).
Since 03b uses cross-cluster z-scores, when a single cell type is scattered across multiple clusters the effects cancel out,
resulting in 33.5% unassigned in discovery lesion B, with two plasma cell clusters (6,962 cells) left unclassified.
Confirm that the stromal ratio and enrichment fold-change do not change substantially from 03b.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".numba_cache"))
import scanpy as sc, pandas as pd, numpy as np
from _utils import annotate_compartment

W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables"); OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
D = sc.read_h5ad(f"{OUT}/all_cells_discovery.h5ad")
old = D.obs["compartment"].astype(str).copy() if "compartment" in D.obs else None
D, votes, pct = annotate_compartment(D, "leiden_0.5", "discovery", TAB)
print("03c 区画比率 (%):"); print(pct.to_string())

cmp_rows = []
if old is not None:
    ct = pd.crosstab(old, D.obs["compartment"].astype(str))
    ct.to_csv(f"{TAB}/03c_vs_03b_confusion_discovery.csv")
    print("\n03b -> 03c の移動:"); print(ct.to_string())
    for lab in sorted(set(old) | set(D.obs["compartment"].astype(str))):
        cmp_rows.append(dict(compartment=lab, n_03b=int((old == lab).sum()),
                             n_03c=int((D.obs["compartment"].astype(str) == lab).sum())))
    pd.DataFrame(cmp_rows).to_csv(f"{TAB}/03c_vs_03b_counts_discovery.csv", index=False)

pooled = 100 * float((D.obs.compartment == "stromal").mean())
print(f"\n探索 pooled 間質比率: 03c {pooled:.2f}%  (03b 13.99%, 申請書 14.8%)")

# Enrichment (stromal share among the top decile of tissue-formation scores): recompute if the 04 score column remains
sc_path = f"{OUT}/all_cells_discovery_obs_scored.parquet"
rows = []
if os.path.exists(sc_path):
    O = pd.read_parquet(sc_path)
    O = O.loc[D.obs_names.intersection(O.index)]
    comp = D.obs.loc[O.index, "compartment"].astype(str)
    for col, lab in [("lf_D_tissueforming", "axis_D_tissueforming"),
                     ("tissueforming_composite", "composite_tissueforming")]:
        if col not in O.columns: continue
        for les in list(sorted(O.lesion.unique())) + ["POOLED_discovery"]:
            m = np.ones(len(O), bool) if les.startswith("POOLED") else (O.lesion.values == les)
            g, c = O[m], comp[m]
            thr = np.quantile(g[col].values, 0.90); top = g[col].values >= thr
            base = float((c == "stromal").mean()); topf = float((c[top] == "stromal").mean())
            rows.append(dict(lesion=les, score=lab, n_cells=int(m.sum()),
                             stromal_pct_overall=round(100*base, 2), stromal_pct_top_decile=round(100*topf, 2),
                             enrichment_fold=round(topf/base, 3) if base > 0 else np.nan))
    E = pd.DataFrame(rows); E.to_csv(f"{TAB}/03c_enrichment_by_lesion_discovery.csv", index=False)
    print("\n濃縮（03c 区画）:")
    print(E[E.score == "axis_D_tissueforming"][["lesion","stromal_pct_overall","stromal_pct_top_decile","enrichment_fold"]].to_string(index=False))

D.write(f"{OUT}/all_cells_discovery.h5ad")
json.dump(dict(stromal_pct_pooled_03c=round(pooled, 2), stromal_pct_pooled_03b=13.99, grant=14.8),
          open(f"{TAB}/03c_stromal_pct_discovery.json", "w"), indent=1)
print("\n03c DONE")
