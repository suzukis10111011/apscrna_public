#!/usr/bin/env python
"""Purpose: Compute the correlation between MYC activity and the mean activity of each
  of the two TF sets, the quantity the manuscript reports as +0.45 and -0.06, which no
  saved table carried.
Input: data/processed/stromal_discovery_labelled.h5ad,
  data/processed/stromal_validation.h5ad, frozen/net/collectri_regulons.csv,
  results/tables/32_set_definition.csv
Output: 08_external_review/tables/E6_myc_vs_set_means.csv

The activity estimator, the network, the minimum regulon size and the set membership
are taken from the co-activity analysis (scripts/32_coactivity_ruleset.py) without
change: decoupleR ULM over the CollecTRI network with tmin=5 on the stromal cells, the
set activity being the mean of the member activities, and Spearman as the estimator.
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
sys.path.insert(0, f"{W}/scripts")
from _utils import ulm_on_symbols  # noqa: E402  the estimator used throughout

TAB = f"{W}/results/tables"
OUT = f"{W}/08_external_review/tables"
os.makedirs(OUT, exist_ok=True)

NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")[["source", "target",
                                                             "weight"]]
SD = pd.read_csv(f"{TAB}/32_set_definition.csv")
SETS = {name: SD[SD.set_name == name].TF.tolist() for name in SD.set_name.unique()}
PATHS = [("discovery_n5", f"{W}/data/processed/stromal_discovery_labelled.h5ad"),
         ("validation_n3", f"{W}/data/processed/stromal_validation.h5ad")]

rows = []
for cohort, path in PATHS:
    A = sc.read_h5ad(path)
    ACT = ulm_on_symbols(A, NET, tmin=5)
    assert "MYC" in ACT.columns, "MYC has no activity estimate in this cohort"
    les = A.obs["lesion"].astype(str).values
    for set_name, members in SETS.items():
        have = [t for t in members if t in ACT.columns]
        mean_act = ACT[have].mean(axis=1).values
        myc = ACT["MYC"].values
        # per-factor correlations as well, so the set-mean value can be compared with
        # the mean and the median of the individual correlations
        per = {t: float(spearmanr(myc, ACT[t].values).statistic) for t in have}
        rows.append(dict(cohort=cohort, lesion="POOLED", set_name=set_name,
                         factors_in_set=";".join(have), n_factors=len(have),
                         n_cells=int(len(myc)), summary="mean of per-factor rho",
                         spearman_rho_MYC_vs_set_mean=round(np.mean(list(per.values())), 4),
                         spearman_p=np.nan))
        rows.append(dict(cohort=cohort, lesion="POOLED", set_name=set_name,
                         factors_in_set=";".join(have), n_factors=len(have),
                         n_cells=int(len(myc)), summary="median of per-factor rho",
                         spearman_rho_MYC_vs_set_mean=round(np.median(list(per.values())), 4),
                         spearman_p=np.nan))
        for t, v in per.items():
            rows.append(dict(cohort=cohort, lesion="POOLED", set_name=set_name,
                             factors_in_set=t, n_factors=1, n_cells=int(len(myc)),
                             summary="single factor vs MYC",
                             spearman_rho_MYC_vs_set_mean=round(v, 4),
                             spearman_p=float(spearmanr(myc, ACT[t].values).pvalue)))
        for group in ["POOLED"] + sorted(set(les)):
            m = np.ones(len(myc), bool) if group == "POOLED" else (les == group)
            r = spearmanr(myc[m], mean_act[m])
            rows.append(dict(cohort=cohort, lesion=group, set_name=set_name,
                             factors_in_set=";".join(have), n_factors=len(have),
                             n_cells=int(m.sum()),
                             summary="rho with the set mean activity",
                             spearman_rho_MYC_vs_set_mean=round(float(r.statistic), 4),
                             spearman_p=float(r.pvalue)))
    del A, ACT

E6 = pd.DataFrame(rows)
E6.to_csv(f"{OUT}/E6_myc_vs_set_means.csv", index=False)
pool = E6[(E6.lesion == "POOLED") & (E6.summary != "single factor vs MYC")]
print(pool[["cohort", "set_name", "summary", "n_cells",
            "spearman_rho_MYC_vs_set_mean"]].to_string(index=False))
print("\nper factor, discovery:")
print(E6[(E6.cohort == "discovery_n5") & (E6.summary == "single factor vs MYC")]
      [["set_name", "factors_in_set", "spearman_rho_MYC_vs_set_mean"]]
      .to_string(index=False))
