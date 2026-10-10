import os
#!/usr/bin/env python
"""Purpose: Redefine TF co-activity rulesets from screen candidates (lf_A inflammation, lf_D tissue formation
  axes) intersected with NF-κB and osteoblast/mesenchymal families; compare new ruleset against prior set on
  sign consistency, Spearman correlation, PC1 variance, and Hartigan dip.
Input: 04_tf_screen_discovery_candidates.csv, 11_coactivity_smad3_sensitivity.csv, collectri_regulons.csv
Output: 32_coactivity_ruleset_sensitivity.csv, 32_crossside_pairs_sign_consistency.csv,
  32_pc1_is_function_of_r.csv,
  32_old_vs_new_comparison.csv, 32_regulon_overlap.csv, 32_set_definition.csv,
  32_tf_activity_spearman_matrix.csv

32: 共活性解析の TF セットを本論文のスクリーン結果だけから規則で定義し直して再実行。

規則（ハードコードせず 04 の候補表から交差を取る）:
  炎症側     = 炎症軸 lf_A の候補（上位20・5病変同符号）∩ NF-kB ファミリー
               {RELB, NFKB2, RELA, NFKB1, REL}
  組織形成側 = 組織形成軸 lf_D の候補 ∩ 骨芽細胞系・間葉系 TF
               {SP7, RUNX2, DLX5, MYRF, SMAD3}、および SMAD3 を除いた4因子
旧セット（NFKB2/RELB/REL/ATF1/MTA1）と並べて、結論が変わるかを判定する。
推定量は 11_coactivity.py と同一（Spearman、PC1 分散比、Hartigan dip、群間ペアの符号一致）。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("32")
_VAL_STATE = os.environ.get("APSC_VAL_STATE", "state_within")
import scanpy as sc, pandas as pd, numpy as np
from scipy.stats import spearmanr
from _utils import ulm_on_symbols
import diptest

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
DISC = f"{W}/data/processed/stromal_discovery_labelled.h5ad"
VAL = f"{W}/data/processed/stromal_validation.h5ad"
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")[["source", "target", "weight"]]

# ---------- Set definition by rule ----------
NFKB_FAMILY = ["RELB", "NFKB2", "RELA", "NFKB1", "REL"]
OSTEO_MESEN = ["SP7", "RUNX2", "DLX5", "MYRF", "SMAD3"]
S04 = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")


def candidates(axis):
    return set(S04[(S04.axis == axis) & (S04.is_candidate)].TF)


cand_A, cand_D = candidates("lf_A_cytokine_chemokine"), candidates("lf_D_tissueforming")
INFL_NEW = [t for t in NFKB_FAMILY if t in cand_A]
TISS_NEW5 = [t for t in OSTEO_MESEN if t in cand_D]
TISS_NEW4 = [t for t in TISS_NEW5 if t != "SMAD3"]
INFL_OLD = ["NFKB2", "RELB", "REL", "ATF1", "MTA1"]
assert INFL_NEW and TISS_NEW5, (INFL_NEW, TISS_NEW5)

rows = []
for lab, axis, fam, sel in [("inflammatory_NFkB_family", "lf_A_cytokine_chemokine", NFKB_FAMILY, INFL_NEW),
                            ("tissueforming_osteo_mesenchymal", "lf_D_tissueforming", OSTEO_MESEN, TISS_NEW5)]:
    for tf in fam:
        r = S04[(S04.axis == axis) & (S04.TF == tf)]
        rows.append(dict(set_name=lab, axis=axis, TF=tf,
                         pooled_rho=round(float(r.pooled_rho.iloc[0]), 4) if len(r) else np.nan,
                         all5_same_sign=bool(r.all5_same_sign.iloc[0]) if len(r) else False,
                         rank_within_axis=float(r.rank_within_axis.iloc[0]) if len(r) else np.nan,
                         is_candidate=bool(r.is_candidate.iloc[0]) if len(r) else False,
                         included_by_rule=tf in sel))
for tf in INFL_OLD:
    r = S04[(S04.axis == "lf_A_cytokine_chemokine") & (S04.TF == tf)]
    rows.append(dict(set_name="inflammatory_OLD_set", axis="lf_A_cytokine_chemokine", TF=tf,
                     pooled_rho=round(float(r.pooled_rho.iloc[0]), 4), all5_same_sign=bool(r.all5_same_sign.iloc[0]),
                     rank_within_axis=float(r.rank_within_axis.iloc[0]), is_candidate=bool(r.is_candidate.iloc[0]),
                     included_by_rule=tf in INFL_NEW))
SD = pd.DataFrame(rows)
SD.to_csv(f"{TAB}/32_set_definition.csv", index=False)
print("=== 規則によるセット ===")
print(f"  炎症側（新）    : {INFL_NEW}")
print(f"  炎症側（旧）    : {INFL_OLD}  → 新で外れる: {[t for t in INFL_OLD if t not in INFL_NEW]} / 新で入る: {[t for t in INFL_NEW if t not in INFL_OLD]}")
print(f"  組織形成側(5)   : {TISS_NEW5}")
print(f"  組織形成側(4)   : {TISS_NEW4}")

# ---------- Regulon overlap (interpretive assumption) ----------
reg = {t: set(NET.loc[NET.source == t, "target"]) for t in set(INFL_NEW + TISS_NEW5 + INFL_OLD) if (NET.source == t).any()}
ov = []
for i, a in enumerate(sorted(reg)):
    for b in sorted(reg):
        if a >= b: continue
        inter = len(reg[a] & reg[b]); uni = len(reg[a] | reg[b])
        ov.append(dict(TF_a=a, TF_b=b, n_a=len(reg[a]), n_b=len(reg[b]), n_shared=inter,
                       jaccard=round(inter / uni, 4) if uni else np.nan,
                       frac_of_smaller=round(inter / min(len(reg[a]), len(reg[b])), 4) if min(len(reg[a]), len(reg[b])) else np.nan))
OV = pd.DataFrame(ov).sort_values("jaccard", ascending=False)
OV.to_csv(f"{TAB}/32_regulon_overlap.csv", index=False)
print("\n=== regulon 重なり上位6（活性は独立でない） ===")
print(OV.head(6)[["TF_a", "TF_b", "n_a", "n_b", "n_shared", "jaccard", "frac_of_smaller"]].to_string(index=False))


def shape2d(x, y):                                    # Verbatim from 11 (does not change the estimator)
    Z = np.column_stack([(x - x.mean()) / x.std(ddof=1), (y - y.mean()) / y.std(ddof=1)])
    ev, evec = np.linalg.eigh(np.cov(Z.T)); o = np.argsort(ev)[::-1]; ev, evec = ev[o], evec[:, o]
    r = spearmanr(x, y)
    return dict(n=len(x), spearman_rho=round(r.statistic, 4), spearman_p=r.pvalue,
                pc1_var_frac=round(float(ev[0] / ev.sum()), 4),
                pc1_angle_deg=round(float(np.degrees(np.arctan2(evec[1, 0], evec[0, 0])) % 180), 1),
                dip_p_infl=round(float(diptest.diptest(x)[1]), 4),
                dip_p_tiss=round(float(diptest.diptest(y)[1]), 4))


SETS = {"infl_new_NFkB": INFL_NEW, "infl_old": INFL_OLD,
        "tiss_new5_incl_SMAD3": TISS_NEW5, "tiss_new4_excl_SMAD3": TISS_NEW4}
sens, mats, pairs = [], [], []
for path, coh, statecol in [(DISC, "discovery", "state"), (VAL, "validation", _VAL_STATE)]:
    if not os.path.exists(path):
        print(f"{coh}: {path} 無し — スキップ"); continue
    A = sc.read_h5ad(path)
    ACT = ulm_on_symbols(A, NET, tmin=5)
    need = sorted(set(sum(SETS.values(), [])))
    have = [t for t in need if t in ACT.columns]
    missing = [t for t in need if t not in have]
    print(f"\n{coh}: {A.shape} | 活性が得られたTF {len(have)}/{len(need)}" + (f" | 欠測 {missing}" if missing else ""))
    D = pd.DataFrame(index=A.obs_names)
    D["lesion"] = A.obs["lesion"].astype(str).values
    D["state"] = A.obs[statecol].astype(str).values if statecol in A.obs else "NA"
    for t in have: D[t] = ACT[t].values
    for k, v in SETS.items():
        cols = [t for t in v if t in have]
        assert cols, (k, v)
        D[k] = D[cols].mean(1)
    for iset in ("infl_new_NFkB", "infl_old"):
        for tset in ("tiss_new5_incl_SMAD3", "tiss_new4_excl_SMAD3"):
            for grp in ("lesion", "state"):
                for lab, g in list(D.groupby(grp, observed=True)) + [("POOLED", D)]:
                    if len(g) < 30: continue
                    sens.append(dict(cohort=coh, inflammatory_set=iset, tissueforming_set=tset,
                                     group=grp, label=str(lab), **shape2d(g[iset].values, g[tset].values)))
    M = pd.DataFrame(spearmanr(D[have].values)[0], index=have, columns=have)
    M.index.name = "TF"; M = M.reset_index(); M.insert(0, "cohort", coh); mats.append(M)
    for iset in ("infl_new_NFkB", "infl_old"):
        for a in [t for t in SETS[iset] if t in have]:
            for b in [t for t in TISS_NEW5 if t in have]:
                rs = [spearmanr(g[a], g[b]).statistic for _, g in D.groupby("lesion", observed=True) if len(g) >= 30]
                pairs.append(dict(cohort=coh, inflammatory_set=iset, infl_TF=a, tiss_TF=b,
                                  pooled_rho=round(float(spearmanr(D[a], D[b]).statistic), 4),
                                  n_lesions=len(rs), n_positive=int(sum(r > 0 for r in rs)),
                                  rho_min=round(min(rs), 4) if rs else np.nan,
                                  rho_max=round(max(rs), 4) if rs else np.nan,
                                  all_same_sign=bool(rs and (all(r > 0 for r in rs) or all(r < 0 for r in rs)))))

S = pd.DataFrame(sens).drop_duplicates()
S.to_csv(f"{TAB}/32_coactivity_ruleset_sensitivity.csv", index=False)
pd.concat(mats).to_csv(f"{TAB}/32_tf_activity_spearman_matrix.csv", index=False)
PR = pd.DataFrame(pairs)
PR.to_csv(f"{TAB}/32_crossside_pairs_sign_consistency.csv", index=False)

print("\n=== 病変別ρ・PC1分散比・dip（POOLED 行のみ） ===")
print(S[S.label == "POOLED"][["cohort", "inflammatory_set", "tissueforming_set", "group", "n",
                              "spearman_rho", "pc1_var_frac", "dip_p_infl", "dip_p_tiss"]].to_string(index=False))
print("\n=== 病変別ρ ===")
print(S[(S.group == "lesion") & (S.label != "POOLED")][["cohort", "inflammatory_set", "tissueforming_set",
                                                        "label", "n", "spearman_rho", "pc1_var_frac"]].to_string(index=False))
print("\n=== 群間ペアの符号一致 ===")
print(PR.groupby(["cohort", "inflammatory_set"]).agg(n_pairs=("all_same_sign", "size"),
                                                     all_same_sign=("all_same_sign", "sum"),
                                                     n_positive_pooled=("pooled_rho", lambda s: int((s > 0).sum()))).to_string())

# ---------- Comparison of old vs new conclusions ----------
old11 = pd.read_csv(f"{TAB}/11_coactivity_smad3_sensitivity.csv")
cmp_rows = []
for coh in S.cohort.unique():
    for tset, old_variant in [("tiss_new5_incl_SMAD3", "tissueforming_5TF_incl_SMAD3"),
                              ("tiss_new4_excl_SMAD3", "tissueforming_4TF_excl_SMAD3")]:
        o = old11[(old11.cohort == coh) & (old11.variant == old_variant) & (old11.group == "lesion") & (old11.label == "POOLED")]
        for iset in ("infl_old", "infl_new_NFkB"):
            n_ = S[(S.cohort == coh) & (S.inflammatory_set == iset) & (S.tissueforming_set == tset) &
                   (S.group == "lesion") & (S.label == "POOLED")]
            if not len(n_): continue
            cmp_rows.append(dict(cohort=coh, tissueforming_set=tset, inflammatory_set=iset,
                                 pooled_rho=float(n_.spearman_rho.iloc[0]),
                                 pc1_var_frac=float(n_.pc1_var_frac.iloc[0]),
                                 rho_in_script_11=float(o.spearman_rho.iloc[0]) if len(o) else np.nan,
                                 reproduces_script_11=bool(len(o) and abs(float(n_.spearman_rho.iloc[0]) - float(o.spearman_rho.iloc[0])) < 1e-3)
                                 if iset == "infl_old" else None))
C = pd.DataFrame(cmp_rows)
C.to_csv(f"{TAB}/32_old_vs_new_comparison.csv", index=False)
print("\n=== 旧セット（11 の値）と新セットの対照 ===")
print(C.to_string(index=False))
# ---------- Record that PC1 variance ratio is a deterministic function of |Pearson r| ----------
# Incorporated here verbatim from 32_pc1_is_function_of_r.csv, which the round-32 dialogue cell had written (identical).
# In PCA on two standardized variables, PC1 variance ratio = (1+|Pearson r|)/2, so
# the PC1 variance ratio is not evidence independent of correlation (explicitly write out |r| = 2*pc1_var_frac - 1).
pool = S[(S.label == "POOLED") & (S.group == "lesion")].set_index(
    ["cohort", "inflammatory_set", "tissueforming_set"])
PC = pd.DataFrame([(coh, i, tt, float(r.pc1_var_frac), 2 * float(r.pc1_var_frac) - 1)
                   for (coh, i, tt), r in pool.iterrows()],
                  columns=["cohort", "inflammatory_set", "tissueforming_set",
                           "pc1_var_frac", "implied_abs_pearson_r"])
PC["note"] = ("Standardized bivariate PCA has PC1 variance ratio = (1+|Pearson r|)/2, "
              "so it is not evidence independent of correlation")
PC.to_csv(f"{TAB}/32_pc1_is_function_of_r.csv", index=False)
print("\n=== PC1 分散比から逆算した |Pearson r| ===")
print(PC.to_string(index=False))

require_heartbeat("32 完了時")
print("\n32 DONE")
