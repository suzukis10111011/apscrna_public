#!/usr/bin/env python
"""Purpose: Does the co-activity result depend on which tissue-forming factors are in the
  fixed list? The NF-kB five are held fixed and the tissue-forming side is replaced by two
  alternatives, per lesion and pooled, in both cohorts.
Input: data/processed/stromal_discovery_labelled.h5ad, stromal_validation.h5ad and the
  matching TF-activity parquets, results/tables/04_tf_screen_discovery_candidates.csv
Output: 08_external_review/tables/E40_coactivity_set_variants.csv

The estimator is the one of scripts/11_coactivity.py and 32_coactivity_ruleset.py: the
mean of the raw ULM activities of the set members per cell, Spearman rho between the two
set means, and Hartigan's dip test on each marginal. The published values of the current
sets are recomputed here from the stored activities and must reproduce
32_coactivity_ruleset_sensitivity.csv; the script asserts that before going on.

Variants of the tissue-forming side:
  current_5_incl_SMAD3   SP7, RUNX2, DLX5, MYRF, SMAD3  (the fixed list, as published)
  current_4_excl_SMAD3   the same without SMAD3
  top5_by_rank           SP7, RUNX2, DLX5, HAND2, SP3   (the five best-ranked candidates)
  all_positive_candidates every candidate of the tissue-forming axis with pooled rho > 0
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/nbc")
import numpy as np
import pandas as pd
import anndata as ad
import diptest
from scipy.stats import spearmanr

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
P, T = f"{W}/data/processed", f"{W}/results/tables"
OUT = f"{W}/08_external_review/tables"
NOTE = ("post hoc, in answer to a question about the fixed tissue-forming list; frozen "
        "gene sets, candidate factors and thresholds are unchanged")

INFL = ["RELB", "NFKB2", "RELA", "NFKB1", "REL"]
S04 = pd.read_csv(f"{T}/04_tf_screen_discovery_candidates.csv")
cand_D = S04[(S04.axis == "lf_D_tissueforming") & S04.is_candidate]
POSITIVE = sorted(cand_D[cand_D.pooled_rho > 0].TF)
VARIANTS = {
    "current_5_incl_SMAD3": ["SP7", "RUNX2", "DLX5", "MYRF", "SMAD3"],
    "current_4_excl_SMAD3": ["SP7", "RUNX2", "DLX5", "MYRF"],
    "top5_by_rank": ["SP7", "RUNX2", "DLX5", "HAND2", "SP3"],
    "all_positive_candidates": POSITIVE,
}

rows = []
for coh, h5, pq in (("discovery_n5", "stromal_discovery_labelled.h5ad",
                     "stromal_discovery_tf_activity.parquet"),
                    ("validation_n3", "stromal_validation.h5ad",
                     "stromal_validation_tf_activity.parquet")):
    obs = ad.read_h5ad(f"{P}/{h5}", backed="r").obs.copy()
    ACT = pd.read_parquet(f"{P}/{pq}")
    D = pd.DataFrame({"lesion": obs.lesion.astype(str).values})
    for name, members in VARIANTS.items():
        cols = [t for t in members if t in ACT.columns]
        assert cols, (coh, name)
        D[name] = ACT[cols].mean(1).values
        D.attrs[name] = cols
    D["inflammatory"] = ACT[[t for t in INFL if t in ACT.columns]].mean(1).values
    for name, members in VARIANTS.items():
        used = [t for t in members if t in ACT.columns]
        for les, g in list(D.groupby("lesion", observed=True)) + [("POOLED", D)]:
            x, y = g.inflammatory.values, g[name].values
            r = spearmanr(x, y)
            rows.append(dict(
                cohort=coh, lesion=str(les), tissueforming_set=name,
                n_factors_defined=len(members), n_factors_used=len(used),
                factors_used="; ".join(used), n_cells=len(g),
                spearman_rho=round(float(r.statistic), 4), spearman_p=float(r.pvalue),
                dip_p_inflammatory=round(float(diptest.diptest(x)[1]), 4),
                dip_p_tissueforming=round(float(diptest.diptest(y)[1]), 4)))
E40 = pd.DataFrame(rows)
E40["inflammatory_set"] = "RELB; NFKB2; RELA; NFKB1; REL (the NF-kB five, held fixed)"
E40["estimator"] = ("mean of the raw ULM activities of the set members per cell; Spearman "
                    "rho between the two set means; Hartigan dip test on each marginal")
E40["note"] = NOTE

# the two published variants must come out as in the published sensitivity table
PUB = pd.read_csv(f"{T}/32_coactivity_ruleset_sensitivity.csv")
PUB = PUB[PUB.inflammatory_set == "infl_new_NFkB"]
pairs = {("current_5_incl_SMAD3", "tiss_new5_incl_SMAD3"),
         ("current_4_excl_SMAD3", "tiss_new4_excl_SMAD3")}
for mine, theirs in pairs:
    for coh, short in (("discovery_n5", "discovery"), ("validation_n3", "validation")):
        a = E40[(E40.cohort == coh) & (E40.tissueforming_set == mine)
                & (E40.lesion == "POOLED")].spearman_rho.iloc[0]
        b = PUB[(PUB.cohort == short) & (PUB.tissueforming_set == theirs)
                & (PUB.group == "lesion") & (PUB.label == "POOLED")].spearman_rho.iloc[0]
        assert abs(a - b) < 1e-4, (mine, coh, a, b)

E40.to_csv(f"{OUT}/E40_coactivity_set_variants.csv", index=False)
print("positive candidates of the tissue-forming axis:", len(POSITIVE), POSITIVE)
print(E40[E40.lesion == "POOLED"][["cohort", "tissueforming_set", "n_factors_used",
                                   "n_cells", "spearman_rho", "spearman_p"]].to_string(index=False))
print("E40 DONE", E40.shape)
