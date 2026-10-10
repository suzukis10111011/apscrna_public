import os
#!/usr/bin/env python
"""Purpose: Define validation-cohort inflammation and stromal states within Leiden 0.5 clusters using
  validation-cohort medians; compare against ingest-projected states.
Input: {OUT}/all_cells_validation.h5ad, {OUT}/stromal_validation.h5ad
Output: {OUT}/stromal_validation.h5ad, {TAB}/25_state_crosstab_within_vs_ingest.csv,
  {TAB}/25_state_definition_validation_within.csv, {TAB}/25_validation_states_within_vs_ingest.csv

25: 検証側の状態ラベルを検証コホート内で定義する（ingest 投影の代替）。
規則は探索（04b）と同一 — Leiden 0.5 クラスタごとの炎症・組織形成スコア中央値を、
コホート全体の中央値と比べて4状態に命名する。閾値だけを検証コホート内で取る。
投影ラベル（state_projected）は残し、両者の差を表に出す。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("25")
import scanpy as sc, pandas as pd, numpy as np

OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
INFL, TISS = "p_A_cytokine_chemokine", "p_D_tissueforming"
St = sc.read_h5ad(f"{OUT}/stromal_validation.h5ad")
St.obs_names_make_unique()
if INFL not in St.obs or TISS not in St.obs:
    # Align with the same order as discovery: "score on the all-cells object → remove stroma"
    import json
    A = sc.read_h5ad(f"{OUT}/all_cells_validation.h5ad"); A.obs_names_make_unique()
    AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
    for k in ("A_cytokine_chemokine", "D_tissueforming"):
        gs = [g for g in sorted(set(AXJ["axes"][k])) if g in set(A.var_names)]
        sc.tl.score_genes(A, gs, ctrl_size=50, n_bins=25, random_state=0, score_name=f"p_{k}")
    for c in (INFL, TISS):
        St.obs[c] = A.obs[c].reindex(St.obs_names).values
    print(f"p_ スコアを全細胞オブジェクトで算出して転記（欠損 {int(St.obs[INFL].isna().sum())}）")
CL = "leiden_0.5"
assert CL in St.obs, f"{CL} が無い"
m_i = float(np.nanmedian(St.obs[INFL])); m_t = float(np.nanmedian(St.obs[TISS]))
print(f"検証コホート内の中央値: 炎症 {m_i:.4f} / 組織形成 {m_t:.4f}")
g = St.obs.groupby(CL, observed=True)
tab = pd.DataFrame({"n_cells": g.size(), "median_inflammatory": g[INFL].median(),
                    "median_tissueforming": g[TISS].median()})
tab["inflammatory_high"] = tab.median_inflammatory > m_i
tab["tissueforming_high"] = tab.median_tissueforming > m_t
def name(r):
    if r.inflammatory_high and not r.tissueforming_high: return "inflammatory"
    if r.tissueforming_high and not r.inflammatory_high: return "tissue-forming"
    if r.inflammatory_high and r.tissueforming_high:     return "intermediate"
    return "other"
tab["state"] = tab.apply(name, axis=1)
tab["reference_median_inflammatory"] = m_i; tab["reference_median_tissueforming"] = m_t
tab["threshold_source"] = "validation cohort internal median (no discovery value used)"
tab.to_csv(f"{TAB}/25_state_definition_validation_within.csv")
print("\nクラスタ -> 状態:")
print(tab[["n_cells","median_inflammatory","median_tissueforming","state"]].round(4).to_string())
St.obs["state_within"] = St.obs[CL].map(tab.state).astype(str)

rows = []
pj = St.obs["state_projected"].astype(str) if "state_projected" in St.obs else pd.Series("NA", index=St.obs_names)
for les in sorted(St.obs.lesion.astype(str).unique()) + ["POOLED"]:
    m = np.ones(len(St), bool) if les == "POOLED" else (St.obs.lesion.astype(str) == les).values
    for lab, col in [("within_cohort", St.obs["state_within"].values), ("ingest_projected", pj.values)]:
        vc = pd.Series(col[m]).value_counts()
        r = dict(lesion=les, labelling=lab, n_cells=int(m.sum()))
        for s in ("inflammatory","intermediate","tissue-forming","other"):
            r[f"n_{s}"] = int(vc.get(s, 0)); r[f"pct_{s}"] = round(100*vc.get(s, 0)/m.sum(), 2)
        rows.append(r)
    rows[-1]["agreement_pct"] = round(100*float((St.obs["state_within"].values[m] == pj.values[m]).mean()), 2)
R = pd.DataFrame(rows); R.to_csv(f"{TAB}/25_validation_states_within_vs_ingest.csv", index=False)
print("\n=== 状態比率（%）===")
print(R[["lesion","labelling","n_cells","pct_inflammatory","pct_intermediate","pct_tissue-forming","pct_other","agreement_pct"]].to_string(index=False))
ct = pd.crosstab(St.obs["state_within"], pj)
ct.to_csv(f"{TAB}/25_state_crosstab_within_vs_ingest.csv")
print("\n=== 混同表（行=コホート内定義, 列=投影）==="); print(ct.to_string())
St.write(f"{OUT}/stromal_validation.h5ad")
require_heartbeat("25 完了時")
print("\n25 DONE")
