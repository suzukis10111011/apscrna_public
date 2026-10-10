#!/usr/bin/env python
"""POST HOC (audit item 4, the 05 vs 18 part).

On the same clean stromal cells the repository holds two values for the same two
quantities:

  depth-controlled rho   04/05 (-0.0653 pooled discovery, -0.0562/-0.2318/-0.1617
                         CAP1-3) vs 18 (0.0349/-0.1510/-0.1822 CAP1-3)
  definition A ratio     04/05 (0.9759/0.9880/1.0042) vs 18 (0.9513/0.9859/1.0083)

Four candidate causes are separated here by recomputing every combination on the
identical cells:

  score columns     p_A x p_D (whole axis, 04/05) vs lf_A x lf_D (RELA targets
                    removed, 12/13/14/18)
  depth covariate   log10(total_counts) (04/05) vs log10(total_counts + 1) (12/13/
                    14/18)
  RELA input        ULM run on the stromal subset (04/05) vs on the whole all-cell
                    object (12/18)
  binning           20 fixed quantiles (04) vs 10 when the lesion has < 500 cells
                    (05/12/18)

Cell-set identity is checked first: if the two scripts disagree on which cells are
stromal, none of the rest is interpretable.

Output: 05_post_audit/tables/Q4a_rho_decomposition.csv
        05_post_audit/tables/Q4a_defA_decomposition.csv
        05_post_audit/tables/Q4a_cellset_identity.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q0_common import AUD, POSTHOC_NOTE, PROC, TAB, W, dctrl_matrix, log  # noqa: E402

sys.path.insert(0, f"{W}/scripts")
from _utils import ulm_on_symbols  # noqa: E402

NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
T18 = pd.read_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv")


def matched(df, col, by, nbins, min_n=25):
    d = df.dropna(subset=[col, by]).copy()
    d["bin"] = pd.qcut(d[by], nbins, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin", observed=True)
         if len(g) >= min_n]
    if not B:
        return np.nan, 0
    n = np.array([b[0] for b in B]); s = np.array([b[1] for b in B])
    return float(np.sqrt(np.average(s ** 2, weights=n)) / tot), len(B)


rho_rows, def_rows, ident = [], [], []
for cohort, coh_lab, st_h5, all_h5, st_act, ref_rho_file, ref_def_file in [
        ("discovery", "discovery_n5", f"{PROC}/stromal_discovery.h5ad",
         f"{PROC}/all_cells_discovery.h5ad",
         f"{PROC}/stromal_discovery_tf_activity.parquet",
         "04_rho_inflam_vs_tissueforming_discovery.csv",
         "04_definitionA_discovery_v2.csv"),
        ("validation", "validation_n3", f"{PROC}/stromal_validation.h5ad",
         f"{PROC}/all_cells_validation.h5ad",
         f"{PROC}/stromal_validation_tf_activity.parquet",
         "05_rho_inflam_vs_tissueforming_validation.csv",
         "05_definitionA_validation_v2.csv")]:
    St = sc.read_h5ad(st_h5)
    O = St.obs.copy()
    O["lesion"] = O.lesion.astype(str)
    ACT_st = pd.read_parquet(st_act)
    log(f"{cohort}: stromal {St.shape}")

    A = sc.read_h5ad(all_h5)
    prim = (A.obs.compartment.astype(str) == "stromal").values
    log(f"{cohort}: all cells {A.shape}, compartment==stromal {int(prim.sum())}")
    ACT_all = ulm_on_symbols(A, NET[["source", "target", "weight"]], tmin=5)
    rela_all = ACT_all["RELA"].values[prim]

    # ---- cell-set identity between the 04/05 stromal object and 18's primary set
    t18 = T18[(T18.cohort == coh_lab) & (T18.cellset == "primary_cluster_based")]
    for les in sorted(O.lesion.unique()) + ["POOLED"]:
        n05 = len(O) if les == "POOLED" else int((O.lesion == les).sum())
        r = t18[t18.lesion == les]
        n18 = int(r.n_cells.iloc[0]) if len(r) else np.nan
        ident.append(dict(cohort=coh_lab, lesion=les,
                          n_cells_stromal_object_04_05=n05,
                          n_cells_primary_cluster_based_18=n18,
                          identical=bool(n05 == n18),
                          note="18 selects compartment=='stromal' from the all-cell "
                               "object; 04/05 write that same selection to "
                               f"{os.path.basename(st_h5)}. " + POSTHOC_NOTE))
    assert int(prim.sum()) == len(O), (int(prim.sum()), len(O))

    # ---- rho: 2 score versions x 2 depth forms
    for scores, slab in [(("p_A_cytokine_chemokine", "p_D_tissueforming"),
                          "p_ (whole axis)"),
                         (("lf_A_cytokine_chemokine", "lf_D_tissueforming"),
                          "lf_ (RELA targets removed)")]:
        for plus1, dlab in [(False, "log10(total_counts)"),
                            (True, "log10(total_counts + 1)")]:
            dep = np.log10(O.total_counts.values.astype(float) + (1 if plus1 else 0))
            for les in sorted(O.lesion.unique()) + ["POOLED"]:
                m = np.ones(len(O), bool) if les == "POOLED" else (O.lesion == les).values
                x = O[scores[0]].values[m]
                y = O[scores[1]].values[m]
                rho_rows.append(dict(
                    cohort=coh_lab, lesion=les, n_cells=int(m.sum()),
                    score_columns=slab, depth_covariate=dlab,
                    raw_rho=round(float(spearmanr(x, y).statistic), 4),
                    depth_controlled_rho=round(float(dctrl_matrix(
                        x.reshape(-1, 1), y, dep[m], O.n_genes.values[m])[0]), 4),
                    matches_script=("04/05" if slab.startswith("p_")
                                    and not plus1 else
                                    ("18 (and 12/13/14)" if slab.startswith("lf_")
                                     and plus1 else "neither (intermediate)")),
                    note=POSTHOC_NOTE))

    # ---- definition A: 2 RELA inputs x 2 binning rules
    for rela, rlab in [(ACT_st["RELA"].values, "ULM on the stromal subset (04/05)"),
                       (rela_all, "ULM on the all-cell object (12/18)")]:
        O2 = O.copy()
        O2["RELA"] = rela
        for rule, blab in [("fixed20", "20 quantiles always (04)"),
                           ("adaptive", "10 quantiles when n < 500 (05/12/18)")]:
            for les in sorted(O2.lesion.unique()):
                g = O2[O2.lesion == les]
                nb = 20 if rule == "fixed20" else (10 if len(g) < 500 else 20)
                fr, used = matched(g, "lf_D_tissueforming", "RELA", nb)
                def_rows.append(dict(
                    cohort=coh_lab, lesion=les, n_cells=len(g), RELA_input=rlab,
                    binning_rule=blab, bins_requested=nb, bins_used=used,
                    frac_retained=round(fr, 4) if fr == fr else np.nan,
                    variance_retained=round(fr ** 2, 4) if fr == fr else np.nan,
                    note=POSTHOC_NOTE))
    del A, ACT_all

R = pd.DataFrame(rho_rows)
R.to_csv(f"{AUD}/tables/Q4a_rho_decomposition.csv", index=False)
D = pd.DataFrame(def_rows)
D.to_csv(f"{AUD}/tables/Q4a_defA_decomposition.csv", index=False)
pd.DataFrame(ident).to_csv(f"{AUD}/tables/Q4a_cellset_identity.csv", index=False)
print(R[R.cohort == "validation_n3"].pivot_table(
    index="lesion", columns=["score_columns", "depth_covariate"],
    values="depth_controlled_rho").to_string())
print(D[D.binning_rule.str.startswith("10")].pivot_table(
    index=["cohort", "lesion"], columns="RELA_input", values="frac_retained").to_string())
log("Q4a done")
