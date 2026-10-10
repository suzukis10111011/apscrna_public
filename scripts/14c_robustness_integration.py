#!/usr/bin/env python
"""Purpose: Quantify stromal cell clustering and TF screening rank shifts under uncorrected PCA versus Harmony
  integration
Input: collectri_regulons.csv, all_cells_discovery.h5ad
Output: 14_robustness_integration.csv, 14_tf_rank_by_integration.csv

14c: How much the stromal set and TF screening rank shift between integration methods (uncorrected PCA vs Harmony).
scVI is omitted since it is not available in the environment (self-judgment).
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("14c")
import scanpy as sc, pandas as pd, numpy as np
from scipy.stats import spearmanr
from _utils import annotate_compartment, ulm_on_symbols, harmony_run

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables"); SEED = 0
A = sc.read_h5ad(f"{W}/data/processed/all_cells_discovery.h5ad")
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
RELA_T = set(NET.loc[NET.source == "RELA", "target"])
s2i = {}
for gid, s in zip(A.var_names, A.var["symbol"].astype(str)): s2i.setdefault(s, gid)
for k in ("A_cytokine_chemokine", "D_tissueforming"):
    ids = [s2i[g] for g in sorted(set(AXJ["axes"][k]) - RELA_T) if g in s2i]
    sc.tl.score_genes(A, ids, ctrl_size=50, n_bins=25, random_state=SEED, score_name=f"lf_{k}")

sets, labels = {}, {}
sets["harmony"] = set(A.obs_names[A.obs.compartment.astype(str) == "stromal"])
labels["harmony"] = A.obs.compartment.astype(str).copy()
B = A.copy()
sc.pp.neighbors(B, n_neighbors=15, use_rep="X_pca", random_state=SEED)
sc.tl.leiden(B, resolution=0.5, key_added="leiden_unc", random_state=SEED)
B, _v, pct = annotate_compartment(B, "leiden_unc", "uncorrected", TAB)
sets["uncorrected"] = set(B.obs_names[B.obs.compartment.astype(str) == "stromal"])
print("uncorrected の区画比率 (%):"); print(pct.to_string())

inter = len(sets["harmony"] & sets["uncorrected"]); uni = len(sets["harmony"] | sets["uncorrected"])
jac = inter / uni if uni else np.nan
print(f"\n間質集合: harmony {len(sets['harmony'])} / uncorrected {len(sets['uncorrected'])} "
      f"| intersection {inter} | Jaccard {jac:.4f}")

ACT = ulm_on_symbols(A, NET[["source","target","weight"]], tmin=5)
rank = {}
for meth, cells in sets.items():
    idx = A.obs_names.isin(cells)
    sub = ACT[idx]; ax = A.obs.loc[idx, "lf_D_tissueforming"].values
    r = {tf: spearmanr(sub[tf].values, ax).statistic for tf in ACT.columns}
    rank[meth] = pd.Series(r)
RK = pd.DataFrame(rank).dropna()
rho = spearmanr(RK["harmony"], RK["uncorrected"]).statistic
rr = RK.rank(ascending=False)
top20 = set(rr.nsmallest(20, "harmony").index) & set(rr.nsmallest(20, "uncorrected").index)
print(f"TF スクリーン（n={len(RK)}）: pooled ρ の Spearman 一致 = {rho:.4f} | 上位20の共通 {len(top20)}/20")
out = pd.DataFrame([dict(comparison="harmony_vs_uncorrected", n_tf=len(RK),
                         stromal_n_harmony=len(sets["harmony"]), stromal_n_uncorrected=len(sets["uncorrected"]),
                         stromal_intersection=inter, stromal_jaccard=round(jac, 4),
                         tf_rank_spearman=round(float(rho), 4), top20_overlap=len(top20),
                         note="scVI not available in this environment, omitted; axis = lf_D_tissueforming")])
out.to_csv(f"{TAB}/14_robustness_integration.csv", index=False)
RK.reset_index().rename(columns={"index":"TF"}).to_csv(f"{TAB}/14_tf_rank_by_integration.csv", index=False)
require_heartbeat("14c 完了時")
print("\n14c DONE")
