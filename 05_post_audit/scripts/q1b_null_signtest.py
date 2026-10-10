#!/usr/bin/env python
"""POST HOC (audit item 1, step 2 of 2).

The pre-registered test counts how many of the 20 frozen candidates of an axis keep,
in all three validation lesions, the sign they had in discovery
(results/tables/05_tf_screen_validation_reproduction.csv: 20/20 on
lf_A_cytokine_chemokine, 19/20 on lf_D_tissueforming). That count was never compared
with anything. Here the same count is computed for 1,000 random draws of 20
non-candidate factors of the same axis, carrying the discovery sign into validation
through the identical residualisation (depth-controlled rho per lesion), the identical
decision rule (all three validation lesions) and the identical cell set
(data/processed/stromal_validation.h5ad, CAP1/CAP2/CAP3, n = 860/566/353).

Two draw pools are reported:
  all_non_candidates             every factor of the axis that is not a candidate
  sign_consistent_non_candidates the same, restricted to factors that are already
                                 sign-consistent across the five discovery lesions
                                 (the screen's own first filter), which is the pool a
                                 candidate was actually drawn from

Two specificity statistics, same 1,000 draws:
  (a) effect size  median |validation pooled rho| over the 20 factors
  (b) ranking      Jaccard(discovery top 20 by |pooled rho|, validation top 20 by
                   |pooled rho|), against Jaccard(random 20, validation top 20)

Both activity versions are covered; for overlap_removed the discovery sign and the
validation rho both come from the overlap-removed activities, and candidates lost to
tmin=5 are named and removed from the denominator.

Output: 05_post_audit/tables/Q1b_null_summary.csv
        05_post_audit/tables/Q1b_candidate_detail.csv
        05_post_audit/tables/Q1b_null_draws.csv.gz
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q0_common import (AUD, AXES, LESIONS, OUT, POSTHOC_NOTE, TAB, jaccard,  # noqa: E402
                       log)

NSEED, NDRAW, TOPN = 1000, 20, 20
SCREEN = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
REPORTED = pd.read_csv(f"{TAB}/05_tf_screen_validation_reproduction.csv")
DROPPED = pd.read_csv(f"{AUD}/tables/Q1a_tmin_dropped_factors.csv")

rows, draws, detail = [], [], []
for axis in AXES:
    frozen_cand = list(SCREEN[(SCREEN.axis == axis) & SCREEN.is_candidate].TF)
    rep = REPORTED[REPORTED.axis == axis]
    reported_k = int(rep.reproduced.sum()) if len(rep) else np.nan
    for version in ["as_screened", "overlap_removed"]:
        Rd = pd.read_parquet(f"{OUT}/rho_discovery_{axis}_{version}.parquet")
        Rv = pd.read_parquet(f"{OUT}/rho_validation_{axis}_{version}.parquet")
        shared = [t for t in Rd.index if t in Rv.index]
        Rd, Rv = Rd.loc[shared], Rv.loc[shared]
        dsign = np.sign(Rd["POOLED"])
        vles = LESIONS["validation"]
        agree = pd.DataFrame({l: np.sign(Rv[l]) == dsign for l in vles}, index=shared)
        n_agree = agree.sum(axis=1)
        strict = n_agree == len(vles)
        major = n_agree >= (len(vles) // 2 + 1)
        disc_consistent = (np.sign(Rd[LESIONS["discovery"]]).nunique(axis=1) == 1)

        cand = [t for t in frozen_cand if t in shared]
        lost = [t for t in frozen_cand if t not in shared]
        for tf in frozen_cand:
            d = DROPPED[(DROPPED.axis == axis) & (DROPPED.TF == tf)]
            detail.append(dict(
                axis=axis, version=version, TF=tf, in_denominator=tf in cand,
                discovery_pooled_rho=float(Rd.POOLED[tf]) if tf in shared else np.nan,
                validation_pooled_rho=float(Rv.POOLED[tf]) if tf in shared else np.nan,
                n_validation_lesions_same_sign=int(n_agree[tf]) if tf in shared else np.nan,
                reproduced_reported_rule=bool(strict[tf]) if tf in shared else np.nan,
                dropped_reason=("" if tf in cand else
                                (d.reason.iloc[0] if len(d) else
                                 "absent from the validation activity matrix")),
                note=POSTHOC_NOTE))

        obs_k = int(strict[cand].sum())
        obs_kmaj = int(major[cand].sum())
        obs_eff = float(Rv.POOLED[cand].abs().median())
        disc_top = list(Rd.POOLED[disc_consistent].abs()
                        .sort_values(ascending=False).head(TOPN).index)
        val_top = list(Rv.POOLED.abs().sort_values(ascending=False).head(TOPN).index)
        obs_j = jaccard(disc_top, val_top)

        for pool_name, pool in [("all_non_candidates",
                                 [t for t in shared if t not in frozen_cand]),
                                ("sign_consistent_non_candidates",
                                 [t for t in shared if t not in frozen_cand
                                  and bool(disc_consistent[t])])]:
            K, KM, EF, JA = [], [], [], []
            for seed in range(NSEED):
                d = list(np.random.default_rng(seed).choice(
                    pool, size=min(NDRAW, len(pool)), replace=False))
                k, km = int(strict[d].sum()), int(major[d].sum())
                e, j = float(Rv.POOLED[d].abs().median()), jaccard(d, val_top)
                K.append(k); KM.append(km); EF.append(e); JA.append(j)
                draws.append(dict(axis=axis, version=version, pool=pool_name, seed=seed,
                                  n_drawn=len(d), n_sign_agree_all3=k,
                                  n_sign_agree_majority=km,
                                  median_abs_validation_rho=round(e, 6),
                                  jaccard_with_validation_top20=round(j, 4)))
            K, KM, EF, JA = map(np.array, (K, KM, EF, JA))
            p = lambda arr, o: round((int((arr >= o).sum()) + 1) / (NSEED + 1), 4)
            rows.append(dict(
                axis=axis, version=version, pool=pool_name,
                n_candidates_frozen=len(frozen_cand), n_candidates_in_denominator=len(cand),
                candidates_dropped=",".join(lost), pool_size=len(pool), n_seeds=NSEED,
                reported_reproduced_count=reported_k if version == "as_screened" else np.nan,
                candidates_sign_agree_all3=obs_k,
                candidates_sign_agree_majority=obs_kmaj,
                null_all3_mean=round(float(K.mean()), 3), null_all3_sd=round(float(K.std(ddof=1)), 3),
                null_all3_q025=float(np.quantile(K, 0.025)),
                null_all3_q975=float(np.quantile(K, 0.975)),
                null_all3_max=int(K.max()),
                n_draws_ge_candidates_all3=int((K >= obs_k).sum()),
                empirical_p_sign_all3=p(K, obs_k),
                null_majority_mean=round(float(KM.mean()), 3),
                n_draws_ge_candidates_majority=int((KM >= obs_kmaj).sum()),
                empirical_p_sign_majority=p(KM, obs_kmaj),
                candidates_median_abs_validation_rho=round(obs_eff, 6),
                null_median_abs_rho_mean=round(float(EF.mean()), 6),
                null_median_abs_rho_q975=round(float(np.quantile(EF, 0.975)), 6),
                empirical_p_effect=p(EF, obs_eff),
                candidate_percentile_in_null_effect=round(float((EF < obs_eff).mean() * 100), 1),
                observed_jaccard_disc_top20_vs_val_top20=round(obs_j, 4),
                null_jaccard_mean=round(float(JA.mean()), 4),
                null_jaccard_q975=round(float(np.quantile(JA, 0.975)), 4),
                empirical_p_jaccard=p(JA, obs_j),
                discovery_top20=",".join(sorted(disc_top)),
                validation_top20=",".join(sorted(val_top)),
                note=POSTHOC_NOTE + " The frozen candidate list is not altered; "
                     "empirical p = (draws at least as extreme + 1) / (n_seeds + 1)."))
            log(f"{axis}/{version}/{pool_name}: candidates {obs_k}/{len(cand)} "
                f"vs null {K.mean():.2f} (max {K.max()}), p={p(K, obs_k)}")

pd.DataFrame(rows).to_csv(f"{AUD}/tables/Q1b_null_summary.csv", index=False)
pd.DataFrame(detail).to_csv(f"{AUD}/tables/Q1b_candidate_detail.csv", index=False)
pd.DataFrame(draws).to_csv(f"{AUD}/tables/Q1b_null_draws.csv.gz", index=False,
                           compression="gzip")
log("Q1b done")
