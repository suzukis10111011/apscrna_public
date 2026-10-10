#!/usr/bin/env python
"""POST HOC (audit item 1, step 1 of 2).

Builds the full TF x lesion matrices of depth-controlled rho that the null test
needs, for two activity versions:

  as_screened       TF activity exactly as used in the manuscript
                    (ULM on frozen/net/collectri_regulons.csv, tmin=5)
  overlap_removed   the same ULM after deleting, from the network, every target
                    that is a member of the axis gene set being correlated
                    (frozen/genesets/stage1_geneset_definitions.json -> axes[k]).
                    Factors left with fewer than 5 targets are dropped by tmin=5
                    and are recorded by name.

The as_screened matrices are checked against the values already in
results/tables (04_tf_screen_discovery_candidates.csv, 33_tf_validation_all730.csv),
which were produced by the manuscript scripts; the check is a correctness test of
the vectorised dctrl used here, not a re-derivation of the reported numbers.

Output: data/post_audit/rho_{cohort}_{axis}_{version}.parquet
        data/post_audit/act_{cohort}_{axis}_overlap_removed.parquet
        05_post_audit/tables/Q1a_activity_version_check.csv
        05_post_audit/tables/Q1a_tmin_dropped_factors.csv
        05_post_audit/tables/Q1a_sign_disagreement_detail.csv
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q0_common import (AUD, AXES, AXIS_SOURCE, OUT, PROC, TAB, W, log,  # noqa: E402
                       POSTHOC_NOTE, rho_by_lesion)

sys.path.insert(0, f"{W}/scripts")
from _utils import ulm_on_symbols  # noqa: E402

os.makedirs(OUT, exist_ok=True)
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
SCREEN = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
A730 = pd.read_csv(f"{TAB}/33_tf_validation_all730.csv")

checks, dropped, sgn_detail = [], [], []
for cohort, h5, act_saved in [("discovery", f"{PROC}/stromal_discovery.h5ad",
                               f"{PROC}/stromal_discovery_tf_activity.parquet"),
                              ("validation", f"{PROC}/stromal_validation.h5ad",
                               f"{PROC}/stromal_validation_tf_activity.parquet")]:
    log(f"{cohort}: reading {os.path.basename(h5)}")
    A = sc.read_h5ad(h5)
    obs = A.obs.copy()
    obs["lesion"] = obs["lesion"].astype(str)
    ACT0 = pd.read_parquet(act_saved)
    assert len(ACT0) == A.n_obs, (len(ACT0), A.n_obs)
    log(f"{cohort}: cells {A.n_obs}, saved activity {ACT0.shape}")

    for axis in AXES:
        # ---- as screened ----
        R = rho_by_lesion(ACT0, obs, axis)
        R.to_parquet(f"{OUT}/rho_{cohort}_{axis}_as_screened.parquet")
        if cohort == "discovery":
            S = SCREEN[SCREEN.axis == axis].set_index("TF")
            ref = S[[f"rho_{l}" for l in ["A", "B", "C", "D", "E"]] + ["pooled_rho"]]
            ref.columns = ["A", "B", "C", "D", "E", "POOLED"]
        else:
            S = A730[A730.axis == axis].set_index("TF")
            ref = S[[f"val_rho_{l}" for l in ["CAP1", "CAP2", "CAP3"]]
                    + ["validation_pooled_rho"]]
            ref.columns = ["CAP1", "CAP2", "CAP3", "POOLED"]
        shared = [t for t in ref.index if t in R.index]
        d = (R.loc[shared, ref.columns] - ref.loc[shared]).abs()
        dis = np.sign(R.loc[shared, ref.columns].values) != np.sign(ref.loc[shared].values)
        for i, j in zip(*np.where(dis)):
            sgn_detail.append(dict(cohort=cohort, axis=axis, TF=shared[i],
                                   column=ref.columns[j],
                                   saved_rounded_value=float(ref.values[i, j]),
                                   recomputed_value=float(
                                       R.loc[shared, ref.columns].values[i, j]),
                                   saved_rounds_to_zero=bool(
                                       abs(float(ref.values[i, j])) == 0.0),
                                   note=POSTHOC_NOTE))
        checks.append(dict(cohort=cohort, axis=axis, version="as_screened",
                           n_tf_compared=len(shared), n_tf_recomputed=len(R),
                           max_abs_diff_vs_saved=float(d.values.max()),
                           n_sign_disagreements=int((np.sign(R.loc[shared, ref.columns])
                                                     != np.sign(ref.loc[shared])).values.sum()),
                           reference_table=("04_tf_screen_discovery_candidates.csv"
                                            if cohort == "discovery"
                                            else "33_tf_validation_all730.csv"),
                           note="saved values are rounded to 4 dp; " + POSTHOC_NOTE))

        # ---- overlap removed ----
        drop = set(AXJ["axes"][AXIS_SOURCE[axis]])
        sub = NET[~NET.target.astype(str).isin(drop)][["source", "target", "weight"]]
        log(f"{cohort}/{axis}: network {len(NET)} -> {len(sub)} edges "
            f"({len(drop)} axis genes removed as targets)")
        ACT1 = ulm_on_symbols(A, sub, tmin=5)
        ACT1.to_parquet(f"{OUT}/act_{cohort}_{axis}_overlap_removed.parquet")
        lost = sorted(set(ACT0.columns) - set(ACT1.columns))
        cand = set(SCREEN[(SCREEN.axis == axis) & SCREEN.is_candidate].TF)
        for tf in lost:
            n_all = int((NET.source == tf).sum())
            n_kept = int((sub.source == tf).sum())
            dropped.append(dict(cohort=cohort, axis=axis, TF=tf,
                                is_frozen_candidate=tf in cand,
                                regulon_size_full=n_all, regulon_size_after_removal=n_kept,
                                reason=f"tmin=5: {n_kept} targets left of {n_all} after "
                                       f"removing the {AXIS_SOURCE[axis]} gene set",
                                note=POSTHOC_NOTE))
        log(f"{cohort}/{axis}: activity {ACT0.shape[1]} -> {ACT1.shape[1]} factors "
            f"({len(lost)} dropped, {len(cand & set(lost))} of them frozen candidates)")
        R1 = rho_by_lesion(ACT1, obs, axis)
        R1.to_parquet(f"{OUT}/rho_{cohort}_{axis}_overlap_removed.parquet")
        checks.append(dict(cohort=cohort, axis=axis, version="overlap_removed",
                           n_tf_compared=0, n_tf_recomputed=len(R1),
                           max_abs_diff_vs_saved=np.nan, n_sign_disagreements=0,
                           reference_table="(none: recomputed for this audit)",
                           note=f"network edges {len(NET)} -> {len(sub)}; " + POSTHOC_NOTE))

C = pd.DataFrame(checks)
C.to_csv(f"{AUD}/tables/Q1a_activity_version_check.csv", index=False)
D = pd.DataFrame(dropped)
D.to_csv(f"{AUD}/tables/Q1a_tmin_dropped_factors.csv", index=False)
SD = pd.DataFrame(sgn_detail)
SD.to_csv(f"{AUD}/tables/Q1a_sign_disagreement_detail.csv", index=False)
log(f"Q1a done. max |diff| vs saved = "
    f"{C.max_abs_diff_vs_saved.max():.2e}; dropped rows {len(D)}")
print(C.to_string(index=False))
