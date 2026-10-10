#!/usr/bin/env python
"""POST HOC (audit item 2).

Definition A (results/tables/04_definitionA_discovery_v2.csv,
05_definitionA_validation_v2.csv) reports, per lesion, the within-bin SD of the
tissue-forming score after matching on RELA activity, divided by its total SD. The
ratio was never put next to a baseline, so a reader cannot tell a conditioning that
removes variance from a conditioning that removes none.

All eight lesions are recomputed here with the identical binning rule (20 quantile
bins, 10 when the lesion has fewer than 500 cells, bins with fewer than 25 cells
dropped, cell-count-weighted within-bin SD) under four conditioning variables:

  RELA_activity        what the manuscript conditions on (04, 05 and 12)
  inflammation_score   lf_A_cytokine_chemokine, the second conditioning named in the
                       audit request; not used anywhere in the manuscript
  null_tf_mean         mean over 10 factors whose depth-controlled rho with the
                       tissue-forming axis is closest to zero in that cohort
                       (selected post hoc; names in Q2_null_tf_selection.csv)
  floor_self           the tissue-forming score itself, i.e. the lowest ratio the
                       procedure can return

Both the SD ratio (frac_retained, as reported) and the variance ratio (its square)
are given.

Output: 05_post_audit/tables/Q2_definitionA_baselines.csv
        05_post_audit/tables/Q2_null_tf_selection.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q0_common import AUD, OUT, POSTHOC_NOTE, PROC, log  # noqa: E402

TIS, INF = "lf_D_tissueforming", "lf_A_cytokine_chemokine"
N_NULL = 10


def matched_residual(df, col, by, min_n=25):
    """Identical to scripts/05 and scripts/12 matched_residual."""
    nbins = 10 if len(df) < 500 else 20
    d = df.dropna(subset=[col, by]).copy()
    d["bin"] = pd.qcut(d[by], nbins, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin", observed=True)
         if len(g) >= min_n]
    if not B:
        return dict(n_cells=len(d), n_bins=0, bins_used=nbins, median_bin_n=np.nan,
                    total_sd=tot, matched_sd=np.nan, frac_retained=np.nan)
    n = np.array([b[0] for b in B]); s = np.array([b[1] for b in B])
    within = float(np.sqrt(np.average(s ** 2, weights=n)))
    return dict(n_cells=len(d), n_bins=len(B), bins_used=nbins,
                median_bin_n=int(np.median(n)), total_sd=tot, matched_sd=within,
                frac_retained=within / tot)


rows, sel_rows = [], []
for cohort, coh_lab, h5, act in [
        ("discovery", "discovery_n5", f"{PROC}/stromal_discovery.h5ad",
         f"{PROC}/stromal_discovery_tf_activity.parquet"),
        ("validation", "validation_n3", f"{PROC}/stromal_validation.h5ad",
         f"{PROC}/stromal_validation_tf_activity.parquet")]:
    A = sc.read_h5ad(h5)
    O = A.obs.copy()
    O["lesion"] = O.lesion.astype(str)
    ACT = pd.read_parquet(act)
    O["RELA_activity"] = ACT["RELA"].values
    O["inflammation_score"] = O[INF].values
    O["floor_self"] = O[TIS].values

    R = pd.read_parquet(f"{OUT}/rho_{cohort}_{TIS}_as_screened.parquet")
    null_tfs = list(R.POOLED.abs().sort_values().head(N_NULL).index)
    for tf in null_tfs:
        O[f"nullTF_{tf}"] = ACT[tf].values
        sel_rows.append(dict(cohort=coh_lab, TF=tf,
                             pooled_depth_controlled_rho_with_tissueforming=round(
                                 float(R.POOLED[tf]), 6),
                             selection="10 smallest |pooled rho| with "
                                       "lf_D_tissueforming in this cohort",
                             note=POSTHOC_NOTE))
    log(f"{cohort}: null factors {null_tfs}")

    for les in sorted(O.lesion.unique()):
        g = O[O.lesion == les]
        per_null = []
        for tf in null_tfs:
            r = matched_residual(g, TIS, f"nullTF_{tf}")
            per_null.append(r["frac_retained"])
            rows.append(dict(cohort=coh_lab, lesion=les, conditioning=f"null_tf_{tf}",
                             conditioning_group="null_tf_individual", **r))
        for name in ["RELA_activity", "inflammation_score", "floor_self"]:
            r = matched_residual(g, TIS, name)
            rows.append(dict(cohort=coh_lab, lesion=les, conditioning=name,
                             conditioning_group=name, **r))
        pn = np.array(per_null, float)
        rows.append(dict(cohort=coh_lab, lesion=les, conditioning="null_tf_mean",
                         conditioning_group="null_tf_mean", n_cells=len(g),
                         n_bins=np.nan, bins_used=10 if len(g) < 500 else 20,
                         median_bin_n=np.nan, total_sd=g[TIS].std(ddof=1),
                         matched_sd=np.nan, frac_retained=float(np.nanmean(pn)),
                         null_tf_frac_retained_sd=float(np.nanstd(pn, ddof=1)),
                         null_tf_frac_retained_min=float(np.nanmin(pn)),
                         null_tf_frac_retained_max=float(np.nanmax(pn))))

T = pd.DataFrame(rows)
T["variance_retained"] = T.frac_retained ** 2
T["sd_reduction_pct"] = 100 * (1 - T.frac_retained)
T["variance_reduction_pct"] = 100 * (1 - T.variance_retained)
T["score_conditioned"] = TIS
T["note"] = POSTHOC_NOTE + (" Only RELA_activity is the conditioning used in the "
                            "manuscript; the other three are baselines added here.")
for c in ["total_sd", "matched_sd", "frac_retained", "variance_retained",
          "sd_reduction_pct", "variance_reduction_pct", "null_tf_frac_retained_sd",
          "null_tf_frac_retained_min", "null_tf_frac_retained_max"]:
    if c in T:
        T[c] = T[c].astype(float).round(6)
T.to_csv(f"{AUD}/tables/Q2_definitionA_baselines.csv", index=False)
pd.DataFrame(sel_rows).to_csv(f"{AUD}/tables/Q2_null_tf_selection.csv", index=False)

M = T[T.conditioning_group.isin(["RELA_activity", "inflammation_score",
                                 "null_tf_mean", "floor_self"])]
print(M.pivot_table(index=["cohort", "lesion"], columns="conditioning",
                    values="frac_retained").round(4).to_string())
log("Q2 done")
