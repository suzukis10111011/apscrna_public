#!/usr/bin/env python
"""Purpose: Define stromal cell states by inflammation and tissue-formation scores using Leiden 0.5 clusters
  from discovery cohort
Input: data/processed/stromal_discovery.h5ad
Output: data/processed/stromal_discovery_labelled.h5ad, results/tables/04b_state_definition_discovery.csv,
  results/tables/04b_state_pct_discovery.csv

04b (first half of step 5): Define state labels for the stromal cells of the 5 discovery lesions.
For each Leiden 0.5 cluster, take the median inflammation score and tissue-formation score,
and name the 4 states by high/low relative to the median across the whole discovery set:
    high inflammation and low tissue-formation -> inflammatory
    low inflammation and high tissue-formation -> tissue-forming
    both high               -> intermediate
    both low               -> other
Input : data/processed/stromal_discovery.h5ad
Output: data/processed/stromal_discovery_labelled.h5ad
       results/tables/04b_state_definition_discovery.csv
       results/tables/04b_state_pct_discovery.csv
Validation-side values are not used (this script only reads discovery data).
"""
import os, sys, numpy as np, pandas as pd, scanpy as sc, warnings
warnings.filterwarnings("ignore")
import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
INFL, TISS = "p_A_cytokine_chemokine", "p_D_tissueforming"

St = sc.read_h5ad(f"{OUT}/stromal_discovery.h5ad")
print("stromal_discovery:", St.shape)
CL = "leiden_0.5"
assert CL in St.obs, f"{CL} が無い"
m_i = float(np.median(St.obs[INFL])); m_t = float(np.median(St.obs[TISS]))
print(f"探索全体の中央値: 炎症 {m_i:.4f} / 組織形成 {m_t:.4f}")

g = St.obs.groupby(CL, observed=True)
tab = pd.DataFrame({"n_cells": g.size(), "median_inflammatory": g[INFL].median(), "median_tissueforming": g[TISS].median()})
tab["inflammatory_high"] = tab.median_inflammatory > m_i
tab["tissueforming_high"] = tab.median_tissueforming > m_t
def name(r):
    if r.inflammatory_high and not r.tissueforming_high: return "inflammatory"
    if r.tissueforming_high and not r.inflammatory_high: return "tissue-forming"
    if r.inflammatory_high and r.tissueforming_high:     return "intermediate"
    return "other"
tab["state"] = tab.apply(name, axis=1)
tab["reference_median_inflammatory"] = m_i
tab["reference_median_tissueforming"] = m_t
tab.to_csv(f"{TAB}/04b_state_definition_discovery.csv")
print("\nクラスタ -> 状態:")
print(tab[["n_cells","median_inflammatory","median_tissueforming","state"]].round(4).to_string())

St.obs["state"] = St.obs[CL].map(tab.state).astype("category")
ct = pd.crosstab(St.obs.lesion, St.obs.state)
pct = (100 * (ct.T / ct.sum(1)).T).round(2)
o = ct.copy(); o.columns = [f"n_{c}" for c in ct.columns]
for c in pct.columns: o[f"pct_{c}"] = pct[c]
o["n_total"] = ct.sum(1); o["cohort"] = "discovery_n5"
o.to_csv(f"{TAB}/04b_state_pct_discovery.csv")
print("\n状態×病変 (%):"); print(pct.to_string())
St.write(f"{OUT}/stromal_discovery_labelled.h5ad")
print("\n04b DONE")
