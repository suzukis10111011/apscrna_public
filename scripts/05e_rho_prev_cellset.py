#!/usr/bin/env python
"""Purpose: Validate whether ρ sign discrepancies explained by prior stromal contamination; extract previous-
  session barcodes from current QC-pass cells, compute composite ρ with current normalization/scoring,
  compare to clean stromal set
Input: validation_covariates_stromal.csv (prior stromal barcodes),
  {W}/data/processed/all_cells_validation.h5ad (current QC-pass cells)
Output: {TAB}/05e_rho_prev_cellset_vs_clean.csv

05e: Validate whether the sign discrepancy in the validation ρ can be explained by "contamination of the previous stromal set."
The previous stromal 5,459 barcodes (validation_covariates_stromal.csv) are extracted from this round's
QC-pass cells, the composite ρ is computed with this round's normalization/scoring/estimator, and compared side by side with this round's clean stromal set.
The estimator is identical to 05 (raw = Spearman; depth-controlled = Spearman after ranking and
OLS residualization on [1, logdepth, n_genes]).
Note: this file is not under frozen/percell/; it is in the artifact store (proj_5d34a2e101d1).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".numba_cache"))
import scanpy as sc, pandas as pd, numpy as np
from scipy.stats import spearmanr, rankdata

W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
A, B = "inflammatory_composite", "tissueforming_composite"

def dctrl(x, y, dep, ng):
    M = np.column_stack([np.ones(len(x)), dep, ng])
    def res(v):
        b, *_ = np.linalg.lstsq(M, rankdata(v), rcond=None); return rankdata(v) - M @ b
    return spearmanr(res(x), res(y)).statistic

V = sc.read_h5ad(f"{W}/data/processed/all_cells_validation.h5ad")
P = pd.read_csv(f"{W}/frozen/percell/validation_covariates_stromal.csv")
prev_bare = pd.DataFrame({"donor": P.iloc[:, 0].astype(str).str.split("_").str[0],
                          "bc": P.iloc[:, 0].astype(str).str.split("_").str[1]})
# The previous donor names are correct (what was wrong was my run→donor mapping). Match using raw barcode + donor.
# The obs_names format differs by step: "CAP1_<bc>" (old) and "<bc>_CAP1" (current, concat's index_unique)
# Both cases are possible, so make it possible to extract the raw barcode either way.
import re as _re
def _bare(s):
    s = str(s)
    m = _re.match(r"^(CAP\d+)[_:](.+)$", s)
    if m: return m.group(2)
    m = _re.match(r"^(.+?)[_-](CAP\d+)(-\d+)?$", s)
    if m: return m.group(1)
    return s
cur = pd.DataFrame({"lesion": V.obs.lesion.astype(str).values,
                    "bc": [_bare(s) for s in V.obs_names]}, index=V.obs_names)
key_prev = set(prev_bare.donor + "|" + prev_bare.bc)
in_prev = (cur.lesion + "|" + cur.bc).isin(key_prev).values
print(f"前回間質 {len(P)} のうち今回の QC 通過細胞に存在: {int(in_prev.sum())}")
assert in_prev.sum() > 1000, f"前回集合の突合が壊れている: {int(in_prev.sum())} 件"
print("病変別:", cur.lesion[in_prev].value_counts().to_dict())

O = V.obs.copy(); O["prev_set"] = in_prev
rows = []
for tag, mask in [("prev_stromal_set", in_prev), ("clean_stromal_2026", (O.compartment == "stromal").values)]:
    sub_all = O[mask]
    for les in list(sorted(sub_all.lesion.astype(str).unique())) + ["POOLED_validation"]:
        s = sub_all if les.startswith("POOLED") else sub_all[sub_all.lesion.astype(str) == les]
        if len(s) < 20: 
            rows.append(dict(cellset=tag, lesion=les, n=len(s), raw_rho=np.nan,
                             depth_controlled_rho=np.nan, note="not computed (n<20)")); continue
        rows.append(dict(cellset=tag, lesion=les, n=len(s),
                         raw_rho=round(spearmanr(s[A], s[B]).statistic, 4),
                         depth_controlled_rho=round(dctrl(s[A].values, s[B].values,
                                                          s.logdepth.values, s.n_genes.values), 4),
                         note=""))
R = pd.DataFrame(rows)
R["definition"] = "composite (inflammatory vs tissueforming)"
# The pooled values in the application form are not public (treated the same way as docs/private/grant_factors.txt).
# Only when a reference value is needed, check the full private-side version 05e_rho_prev_cellset_vs_clean_with_grant.csv.
R["source"] = np.where(R.cellset == "prev_stromal_set",
                       "pre-correction stromal barcode set, recomputed with the current matrix, normalisation and scoring",
                       "stromal aggregation by compartment judgment in this study")
R.to_csv(f"{TAB}/05e_rho_prev_cellset_vs_clean.csv", index=False)
print("\n=== composite ρ: 前回集合 vs 今回集合 ===")
print(R[["cellset","lesion","n","raw_rho","depth_controlled_rho","note"]].to_string(index=False))
