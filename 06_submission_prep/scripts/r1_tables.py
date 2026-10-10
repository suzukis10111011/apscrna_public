#!/usr/bin/env python
"""Pre-submission tidy-up, table side (requests 1, 3, 4).

Nothing in results/tables, data/processed or frozen/ is modified. New versions are
written to 06_submission_prep/tables; r4 records which new file supersedes which
published one.

(1) S3 rebuilt so that every rho carries its score convention in the column name:
    composite, plain (p_, whole axis) and leak-free (lf_, RELA targets removed) side
    by side, for both cell sets and both cohorts. The depth covariates are the ones
    04/05 use, so the composite and plain columns reproduce the published values
    exactly; the leak-free columns are new.

(3) Definition A recomputed with the RELA activity taken from a ULM run on the cell
    set being reported, which is what 04/05 do. The values currently in S5 (18) and
    S8 (12) come from a ULM on the whole all-cell object; both are given with the
    difference.

(4) S2 regenerated at full precision (no 4-dp rounding), all six axes.
    S5's TF sign table regenerated with the depth-controlled estimator that S2 uses,
    instead of the plain Spearman of 18:111-112.

Output (06_submission_prep/tables):
    S3_rho_by_score_convention.csv
    S3_reproduction_check_vs_published.csv
    S2_tf_screen_discovery_full_precision.csv
    S2_tf_screen_validation_reproduction_full_precision.csv
    S5_stromal_cellsets.csv
    S5_tf_sign_by_cellset_depth_controlled.csv
    definitionA_RELA_input_unification.csv
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

W = "/path/to/Documents/apscrna_v2"
sys.path.insert(0, f"{W}/05_post_audit/scripts")
sys.path.insert(0, f"{W}/scripts")
from q0_common import dctrl_matrix  # noqa: E402
from _utils import ulm_on_symbols  # noqa: E402

TAB = f"{W}/results/tables"
OUT = f"{W}/06_submission_prep/tables"
PROC = f"{W}/data/processed"
SEED, Q = 0, 80
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
AXES6 = [f"lf_{k}" for k in AXJ["axes"]]
POS = ["COL1A1", "DCN", "LUM"]
NEG = {"immune": ["PTPRC"], "endothelial": ["PECAM1"], "pericyte": ["RGS5"],
       "epithelial": ["KRT14"]}
SCORE_SETS = {
    "composite": ("inflammatory_composite", "tissueforming_composite",
                  "equal-weighted composite of the component programme scores"),
    "plain_p": ("p_A_cytokine_chemokine", "p_D_tissueforming",
                "whole-axis score (p_), RELA target genes retained"),
    "leakfree_lf": ("lf_A_cytokine_chemokine", "lf_D_tissueforming",
                    "leak-free score (lf_), RELA target genes removed from the axis"),
}
NOTE = ("pre-submission regeneration; frozen gene sets, candidate factors and "
        "thresholds unchanged")


def log(m):
    s = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {m}"
    print(s, flush=True)
    open(f"{W}/06_submission_prep/logs/run_submission_prep.log", "a").write(s + "\n")


def dctrl1(x, y, dep, ng):
    return float(dctrl_matrix(np.asarray(x, float).reshape(-1, 1), y, dep, ng)[0])


def matched(df, col, by, min_n=25):
    nbins = 10 if len(df) < 500 else 20
    d = df.dropna(subset=[col, by]).copy()
    if len(d) < 60:
        return np.nan, 0, nbins
    d["bin"] = pd.qcut(d[by], nbins, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin", observed=True)
         if len(g) >= min_n]
    if not B:
        return np.nan, 0, nbins
    n = np.array([b[0] for b in B]); s = np.array([b[1] for b in B])
    return float(np.sqrt(np.average(s ** 2, weights=n)) / tot), len(B), nbins


OBS_PARQUET = {"discovery_n5": f"{PROC}/all_cells_discovery_obs_scored.parquet",
               "validation_n3": None}

COH = [("discovery_n5", "discovery", f"{PROC}/all_cells_discovery.h5ad",
        f"{PROC}/stromal_discovery.h5ad",
        f"{PROC}/stromal_discovery_tf_activity.parquet",
        "04_discovery_scores_and_screen.py"),
       ("validation_n3", "validation", f"{PROC}/all_cells_validation.h5ad",
        f"{PROC}/stromal_validation.h5ad",
        f"{PROC}/stromal_validation_tf_activity.parquet",
        "05_validation_and_projection.py")]

s3_rows, s5_rows, sign_rows, screen_rows, rep_rows = [], [], [], [], []
CAND = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
disc_cand_rho = {}

for coh, tag, all_h5, st_h5, st_act, script in COH:
    A = sc.read_h5ad(all_h5)
    if OBS_PARQUET[coh]:
        # 04 writes the scored obs to a parquet instead of rewriting the h5ad;
        # it is row-aligned with the object (checked on the barcode index).
        P = pd.read_parquet(OBS_PARQUET[coh])
        assert (P.index == A.obs.index).all(), "scored obs parquet is not aligned"
        for c in P.columns:
            if c not in A.obs.columns:
                A.obs[c] = P[c].values
    O = A.obs.copy()
    O["lesion"] = O.lesion.astype(str)
    O.index = pd.RangeIndex(len(O))
    prim = (O.compartment.astype(str) == "stromal").values
    log(f"{coh}: all cells {A.shape}, stromal {int(prim.sum())}")

    for cellset, mask in [("all_cells", np.ones(len(O), bool)), ("stromal", prim)]:
        sub = O[mask]
        for les in sorted(sub.lesion.unique()) + [f"POOLED_{tag}"]:
            g = sub if les.startswith("POOLED") else sub[sub.lesion == les]
            row = dict(cohort=coh, lesion=les, cell_set=cellset, n_cells=len(g),
                       depth_covariate="log10(total_counts), n_genes_by_counts "
                                       "(as in 04/05)",
                       source_script=script,
                       source_file="06_submission_prep/tables/"
                                   "S3_rho_by_score_convention.csv",
                       note=NOTE)
            for lab, (a, b, desc) in SCORE_SETS.items():
                row[f"raw_rho_{lab}"] = float(spearmanr(g[a], g[b]).statistic)
                row[f"depth_controlled_rho_{lab}"] = dctrl1(
                    g[a].values, g[b].values, g.logdepth.values, g.n_genes.values)
                row[f"score_columns_{lab}"] = f"{a} x {b} - {desc}"
            s3_rows.append(row)

    s2i = {}
    for gid, s in zip(A.var_names, A.var["symbol"].astype(str)):
        s2i.setdefault(s, gid)

    def sco(name, gs):
        ids = [s2i[x] for x in gs if x in s2i]
        if not ids:
            return 0
        sc.tl.score_genes(A, ids, ctrl_size=50, n_bins=25, random_state=SEED,
                          score_name=name)
        return len(ids)

    sco("m_stromal", POS)
    neg_ok = [k for k, gs in NEG.items() if sco(f"m_{k}", gs) > 0]
    M = A.obs.copy()
    M["lesion"] = M.lesion.astype(str)
    M.index = pd.RangeIndex(len(M))
    sel = np.zeros(len(M), bool)
    for les, g in M.groupby("lesion"):
        thr = np.percentile(g.m_stromal, Q)
        m = g.m_stromal >= thr
        for k in neg_ok:
            m &= g[f"m_{k}"] < g[f"m_{k}"].median()
        sel[M.index.get_indexer(g.index[m])] = True
    inter = int((sel & prim).sum()); uni = int((sel | prim).sum())
    log(f"{coh}: marker-defined {int(sel.sum())} / primary {int(prim.sum())} "
        f"Jaccard {inter / uni:.4f}")

    ACT_all = ulm_on_symbols(A, NET[["source", "target", "weight"]], tmin=5)
    O["logdepth_plus1"] = np.log10(O.total_counts.astype(float) + 1)

    for cellset, mask in [("marker_defined", sel), ("primary_cluster_based", prim)]:
        if cellset == "primary_cluster_based":
            ACT_sub = pd.read_parquet(st_act)
        else:
            Asub = A[mask].copy()
            ACT_sub = ulm_on_symbols(Asub, NET[["source", "target", "weight"]], tmin=5)
            del Asub
        assert len(ACT_sub) == int(mask.sum()), (len(ACT_sub), int(mask.sum()))
        sub = O[mask].copy()
        sub["RELA_cellset_ulm"] = ACT_sub["RELA"].values
        sub["RELA_allcell_ulm"] = ACT_all["RELA"].values[mask]

        for les in sorted(sub.lesion.unique()) + ["POOLED"]:
            g = sub if les == "POOLED" else sub[sub.lesion == les]
            if len(g) < 60:
                continue
            fr_new, nb_new, bins_req = matched(g, "lf_D_tissueforming",
                                               "RELA_cellset_ulm")
            fr_old, nb_old, _ = matched(g, "lf_D_tissueforming", "RELA_allcell_ulm")
            s5_rows.append(dict(
                cohort=coh, cell_set=cellset, lesion=les, n_cells=len(g),
                jaccard_marker_vs_primary=(round(inter / uni, 4)
                                           if cellset == "marker_defined" else np.nan),
                score_columns="lf_A_cytokine_chemokine x lf_D_tissueforming "
                              "(leak-free: RELA target genes removed from the axis)",
                raw_rho=float(spearmanr(g.lf_A_cytokine_chemokine,
                                        g.lf_D_tissueforming).statistic),
                depth_controlled_rho_as_in_04_05=dctrl1(
                    g.lf_A_cytokine_chemokine.values, g.lf_D_tissueforming.values,
                    g.logdepth.values, g.n_genes.values),
                depth_controlled_rho_as_in_18=dctrl1(
                    g.lf_A_cytokine_chemokine.values, g.lf_D_tissueforming.values,
                    g.logdepth_plus1.values, g.n_genes.values),
                defA_frac_retained_RELA_from_cellset_ulm=fr_new,
                defA_frac_retained_RELA_from_allcell_ulm=fr_old,
                defA_difference=(fr_new - fr_old
                                 if (fr_new == fr_new and fr_old == fr_old) else np.nan),
                defA_n_bins_used=nb_new, defA_bins_requested=bins_req,
                defA_conditioning="RELA activity; the first column uses a ULM run on "
                                  "this cell set (the 04/05 convention), the second a "
                                  "ULM run on the whole all-cell object (18 and 12)",
                source_script="06_submission_prep/scripts/r1_tables.py",
                note=NOTE))

        for axis in sorted(CAND.axis.unique()):
            for tf in CAND[(CAND.axis == axis) & CAND.is_candidate].TF:
                if tf not in ACT_sub.columns:
                    continue
                x = ACT_sub[tf].values
                per_p, per_d = [], []
                for les in sorted(sub.lesion.unique()):
                    m2 = (sub.lesion == les).values
                    if m2.sum() < 30:
                        continue
                    per_p.append(float(spearmanr(x[m2], sub[axis].values[m2]).statistic))
                    per_d.append(dctrl1(x[m2], sub[axis].values[m2],
                                        sub.logdepth.values[m2], sub.n_genes.values[m2]))
                sign_rows.append(dict(
                    cohort=coh, cell_set=cellset, axis=axis, TF=tf,
                    n_lesions=len(per_d),
                    pooled_rho_plain_spearman=float(
                        spearmanr(x, sub[axis].values).statistic),
                    pooled_rho_depth_controlled=dctrl1(
                        x, sub[axis].values, sub.logdepth.values, sub.n_genes.values),
                    all_same_sign_plain_spearman=bool(
                        per_p and (all(v > 0 for v in per_p) or all(v < 0 for v in per_p))),
                    all_same_sign_depth_controlled=bool(
                        per_d and (all(v > 0 for v in per_d) or all(v < 0 for v in per_d))),
                    estimator="depth-controlled rho on [1, log10(total_counts), "
                              "n_genes_by_counts], the estimator used in S2; the plain "
                              "Spearman column is the one 18:111-112 published",
                    score_columns=f"TF activity x {axis} (leak-free axis score)",
                    note=NOTE))
        del ACT_sub

    St = sc.read_h5ad(st_h5)
    SO = St.obs.copy()
    SO["lesion"] = SO.lesion.astype(str)
    ACT_st = pd.read_parquet(st_act)
    les_list = sorted(SO.lesion.unique())
    for axis in AXES6:
        y = SO[axis].values
        Rm = {}
        for les in les_list:
            m = (SO.lesion == les).values
            Rm[les] = dctrl_matrix(ACT_st.values[m], y[m], SO.logdepth.values[m],
                                   SO.n_genes.values[m])
        Rm["POOLED"] = dctrl_matrix(ACT_st.values, y, SO.logdepth.values,
                                    SO.n_genes.values)
        RM = pd.DataFrame(Rm, index=ACT_st.columns)
        if coh == "discovery_n5":
            for tf, r in RM.iterrows():
                per = r[les_list].values
                q = CAND[(CAND.axis == axis) & (CAND.TF == tf)]
                screen_rows.append(dict(
                    axis=axis, TF=tf, cohort=coh, cell_set="stromal",
                    pooled_rho=float(r["POOLED"]), abs_pooled=abs(float(r["POOLED"])),
                    n_pos=int((per > 0).sum()),
                    all5_same_sign=bool((per > 0).all() or (per < 0).all()),
                    **{f"rho_{l}": float(v) for l, v in zip(les_list, per)},
                    is_candidate=bool(q.is_candidate.iloc[0]) if len(q) else False,
                    score_columns=f"TF activity x {axis} (leak-free axis score)",
                    depth_covariate="log10(total_counts), n_genes_by_counts",
                    note=NOTE + "; full precision, no rounding"))
                disc_cand_rho[(axis, tf)] = float(r["POOLED"])
        else:
            for tf in CAND[(CAND.axis == axis) & CAND.is_candidate].TF:
                if tf not in RM.index:
                    rep_rows.append(dict(axis=axis, TF=tf, in_validation=False,
                                         note=NOTE))
                    continue
                r = RM.loc[tf]
                d = disc_cand_rho.get((axis, tf), np.nan)
                per = r[les_list].values
                rep_rows.append(dict(
                    axis=axis, TF=tf, in_validation=True,
                    discovery_pooled_rho=d, validation_pooled_rho=float(r["POOLED"]),
                    n_val_same_sign_as_discovery=int(
                        sum(np.sign(v) == np.sign(d) for v in per)),
                    reproduced=bool(all(np.sign(v) == np.sign(d) for v in per)),
                    **{f"rho_{l}": float(v) for l, v in zip(les_list, per)},
                    score_columns=f"TF activity x {axis} (leak-free axis score)",
                    depth_covariate="log10(total_counts), n_genes_by_counts",
                    cell_set="stromal",
                    note=NOTE + "; full precision, no rounding"))
    del A, ACT_all, St

S3 = pd.DataFrame(s3_rows)
S3.to_csv(f"{OUT}/S3_rho_by_score_convention.csv", index=False)
pd.DataFrame(s5_rows).to_csv(f"{OUT}/S5_stromal_cellsets.csv", index=False)
pd.DataFrame(sign_rows).to_csv(f"{OUT}/S5_tf_sign_by_cellset_depth_controlled.csv",
                               index=False)
pd.DataFrame(screen_rows).to_csv(f"{OUT}/S2_tf_screen_discovery_full_precision.csv",
                                 index=False)
pd.DataFrame(rep_rows).to_csv(
    f"{OUT}/S2_tf_screen_validation_reproduction_full_precision.csv", index=False)

chk = []
for T, coh in [(pd.read_csv(f"{TAB}/04_rho_inflam_vs_tissueforming_discovery.csv"),
                "discovery_n5"),
               (pd.read_csv(f"{TAB}/05_rho_inflam_vs_tissueforming_validation.csv"),
                "validation_n3")]:
    for _, r in T.iterrows():
        lab = "plain_p" if r.definition == "axis" else "composite"
        n = S3[(S3.cohort == coh) & (S3.lesion == r.lesion)
               & (S3.cell_set == r.cellset)]
        if not len(n):
            continue
        new = round(float(n.iloc[0][f"depth_controlled_rho_{lab}"]), 4)
        chk.append(dict(cohort=coh, lesion=r.lesion, cell_set=r.cellset,
                        published_definition_label=r.definition,
                        score_convention=lab, published=r.depth_controlled_rho,
                        regenerated_rounded_to_4dp=new,
                        abs_difference=abs(r.depth_controlled_rho - new), note=NOTE))
C = pd.DataFrame(chk)
C.to_csv(f"{OUT}/S3_reproduction_check_vs_published.csv", index=False)

S5 = pd.read_csv(f"{OUT}/S5_stromal_cellsets.csv")
D4 = pd.read_csv(f"{TAB}/04_definitionA_discovery_v2.csv")
D5 = pd.read_csv(f"{TAB}/05_definitionA_validation_v2.csv")
T18 = pd.read_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv")
C12 = pd.read_csv(f"{TAB}/12_coupling_by_compartment.csv")
rows = []
for _, r in S5.iterrows():
    pub, pub_src = np.nan, ""
    if r.cell_set == "primary_cluster_based" and r.lesion != "POOLED":
        ref = D4 if r.cohort == "discovery_n5" else D5
        q = ref[ref.lesion == r.lesion]
        if len(q):
            pub = float(q.frac_retained.iloc[0])
            pub_src = ("04_definitionA_discovery_v2.csv" if r.cohort == "discovery_n5"
                       else "05_definitionA_validation_v2.csv")
    q18 = T18[(T18.cohort == r.cohort) & (T18.cellset == r.cell_set)
              & (T18.lesion == r.lesion)]
    q12 = C12[(C12.compartment == "stromal") & (C12.lesion == r.lesion)
              & (C12.cohort == r.cohort)]
    d18 = (r.defA_frac_retained_RELA_from_cellset_ulm
           - float(q18.defA_frac_retained.iloc[0])) if len(q18) else np.nan
    rows.append(dict(
        cohort=r.cohort, cell_set=r.cell_set, lesion=r.lesion, n_cells=r.n_cells,
        unified_value_RELA_from_cellset_ulm=r.defA_frac_retained_RELA_from_cellset_ulm,
        published_04_05_stromal_ulm=pub, published_04_05_source=pub_src,
        published_18_allcell_ulm=(float(q18.defA_frac_retained.iloc[0])
                                  if len(q18) else np.nan),
        published_12_allcell_ulm=(float(q12.defA_frac_retained.iloc[0])
                                  if len(q12) and r.cell_set == "primary_cluster_based"
                                  else np.nan),
        difference_vs_18=d18,
        replacement_needed=bool(d18 == d18 and abs(d18) >= 5e-5),
        note=NOTE))
pd.DataFrame(rows).to_csv(f"{OUT}/definitionA_RELA_input_unification.csv", index=False)
log(f"r1 done: S3 {len(S3)} rows (max |diff| vs published "
    f"{C.abs_difference.max():.1e}), S5 {len(S5)}, sign {len(sign_rows)}, "
    f"S2 discovery {len(screen_rows)}, S2 validation {len(rep_rows)}")
