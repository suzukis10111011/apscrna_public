#!/usr/bin/env python
"""Purpose: Define stromal cells by marker genes (COL1A1/DCN/LUM high, immune/endothelial/epithelial markers
  low) and assess TF activity in discovery vs validation cohorts
Input: {TAB}/04_tf_screen_discovery_candidates.csv, {W}/frozen/net/collectri_regulons.csv
Output: {TAB}/18_stromal_markerdefined_vs_primary.csv, {TAB}/18_tf_sign_by_cellset.csv

18: Reproduction in a stromal set independent of embedding/clustering.
Rule: within the lesion, score_genes for COL1A1/DCN/LUM is high (>= Q percentile within the lesion) and
      each of the PTPRC / PECAM1 / RGS5 / KRT14 scores is below the within-lesion median.
Judgment call (14): Q is decided during discovery and applied as-is to validation. We adopt Q=80
      as a pre-specified round default value, and after confirming in discovery that the stromal
      proportion is comparable to the main analysis (13.99%), apply it to validation.
      Choosing Q to match the main analysis's set would eliminate independence, so no fitting is done.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("18")
import scanpy as sc, pandas as pd, numpy as np
from scipy.stats import spearmanr, rankdata
from _utils import ulm_on_symbols

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables"); SEED = 0; Q = 80
POS = ["COL1A1", "DCN", "LUM"]
NEG = {"immune": ["PTPRC"], "endothelial": ["PECAM1"], "pericyte": ["RGS5"], "epithelial": ["KRT14"]}
AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
RELA_T = set(NET.loc[NET.source == "RELA", "target"])
INF, TIS = "lf_A_cytokine_chemokine", "lf_D_tissueforming"
CAND = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
if "is_candidate" in CAND.columns: CAND = CAND[CAND.is_candidate]

def dctrl(x, y, dep, ng):
    M = np.column_stack([np.ones(len(x)), dep, ng])
    res = lambda v: rankdata(v) - M @ np.linalg.lstsq(M, rankdata(v), rcond=None)[0]
    return spearmanr(res(x), res(y)).statistic

def matched(df, col, min_n=25):
    nb = 10 if len(df) < 500 else 20
    d = df.dropna(subset=[col, "RELA"]).copy()
    if len(d) < 60: return dict(defA_n_bins=0, defA_frac_retained=np.nan, defA_bins_used=nb)
    d["bin"] = pd.qcut(d.RELA, nb, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin", observed=True) if len(g) >= min_n]
    if not B: return dict(defA_n_bins=0, defA_frac_retained=np.nan, defA_bins_used=nb)
    n = np.array([b[0] for b in B]); s = np.array([b[1] for b in B])
    return dict(defA_n_bins=len(B), defA_frac_retained=round(float(np.sqrt(np.average(s**2, weights=n))/tot), 4),
                defA_bins_used=nb)

rows, sign_rows = [], []
for path, coh in [(f"{W}/data/processed/all_cells_discovery.h5ad", "discovery_n5"),
                  (f"{W}/data/processed/all_cells_validation.h5ad", "validation_n3")]:
    if not os.path.exists(path): print(f"{coh}: {path} 無し"); continue
    A = sc.read_h5ad(path)
    s2i = {}
    for gid, s in zip(A.var_names, A.var["symbol"].astype(str)): s2i.setdefault(s, gid)
    def sco(name, gs):
        ids = [s2i[g] for g in gs if g in s2i]
        if not ids:
            print(f"  {name}: 遺伝子が行列に無い（{gs}）— このフィルタは適用しない")
            return 0
        sc.tl.score_genes(A, ids, ctrl_size=50, n_bins=25, random_state=SEED, score_name=name)
        return len(ids)
    sco("m_stromal", POS)
    neg_ok = [k for k, gs in NEG.items() if sco(f"m_{k}", gs) > 0]
    neg_missing = [k for k in NEG if k not in neg_ok]
    for k in AXJ["axes"]:                      # Since the candidate table has 6 axes, construct all axes
        sco(f"lf_{k}", sorted(set(AXJ["axes"][k]) - RELA_T))
    O = A.obs.copy(); O["lesion"] = O.lesion.astype(str)
    O.index = pd.RangeIndex(len(O))     # Fix to positional index (does not break even if there are duplicate barcodes)
    sel = np.zeros(len(O), bool)
    for les, g in O.groupby("lesion"):
        thr = np.percentile(g.m_stromal, Q)
        m = (g.m_stromal >= thr)
        for k in neg_ok: m &= (g[f"m_{k}"] < g[f"m_{k}"].median())
        sel[O.index.get_indexer(g.index[m])] = True
    O["marker_stromal"] = sel
    prim = (O.compartment.astype(str) == "stromal").values
    inter = int((sel & prim).sum()); uni = int((sel | prim).sum())
    print(f"\n{coh}: マーカー定義 {int(sel.sum())} / 主解析 {int(prim.sum())} | 共通 {inter} | Jaccard {inter/uni:.4f}")
    ACT = ulm_on_symbols(A, NET[["source","target","weight"]], tmin=5)
    O["RELA"] = ACT["RELA"].values if "RELA" in ACT.columns else np.nan
    O["logdepth"] = np.log10(O.total_counts.astype(float) + 1); O["n_genes"] = O.n_genes_by_counts.astype(float)
    for setlab, mask in [("marker_defined", sel), ("primary_cluster_based", prim)]:
        sub = O[mask]
        for les in sorted(O.lesion.unique()) + ["POOLED"]:
            g = sub if les == "POOLED" else sub[sub.lesion == les]
            allg = O if les == "POOLED" else O[O.lesion == les]
            if len(g) < 60: continue
            thr = np.quantile(allg[TIS], 0.90); top = allg[allg[TIS] >= thr]
            base = float(mask[O.index.get_indexer(allg.index)].mean())
            topf = float(mask[O.index.get_indexer(top.index)].mean())
            r = dict(cohort=coh, cellset=setlab, lesion=les, n_cells=len(g),
                     stromal_pct_overall=round(100*base, 2), stromal_pct_top_decile=round(100*topf, 2),
                     enrichment_fold=round(topf/base, 3) if base > 0 else np.nan,
                     raw_rho=round(float(spearmanr(g[INF], g[TIS]).statistic), 4),
                     depth_controlled_rho=round(float(dctrl(g[INF].values, g[TIS].values,
                                                            g.logdepth.values, g.n_genes.values)), 4),
                     jaccard_vs_primary=round(inter/uni, 4) if setlab == "marker_defined" else np.nan,
                     threshold=(f"within-lesion Q{Q} or above and {len(neg_ok)} negative markers below the within-lesion median"
                                + (f"({','.join(neg_missing)} absent from the matrix, not applied)" if neg_missing else ""))
                                if setlab == "marker_defined" else "cluster-level assignment")
            r.update(matched(g, TIS)); rows.append(r)
        idx = O.index.get_indexer(sub.index)
        for axis in CAND.axis.unique():
            for tf in CAND[CAND.axis == axis].TF:
                if tf not in ACT.columns: continue
                per = []
                for les in sorted(sub.lesion.unique()):
                    m2 = (sub.lesion == les).values
                    if m2.sum() < 30: continue
                    per.append(spearmanr(ACT[tf].values[idx][m2], sub[axis].values[m2]).statistic)
                pooled = spearmanr(ACT[tf].values[idx], sub[axis].values).statistic
                sign_rows.append(dict(cohort=coh, cellset=setlab, axis=axis, TF=tf,
                                      pooled_rho=round(float(pooled), 4), n_lesions=len(per),
                                      all_same_sign=bool(per and (all(p > 0 for p in per) or all(p < 0 for p in per)))))
R = pd.DataFrame(rows); R.to_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv", index=False)
SG = pd.DataFrame(sign_rows)
piv = SG.pivot_table(index=["cohort","axis"], columns="cellset", values="pooled_rho", aggfunc="mean")
cmp_ = SG.pivot_table(index=["cohort","axis","TF"], columns="cellset", values="pooled_rho")
cmp_ = cmp_.dropna()
agree = (np.sign(cmp_["marker_defined"]) == np.sign(cmp_["primary_cluster_based"])).groupby(level=[0,1]).agg(["sum","size"])
SG.to_csv(f"{TAB}/18_tf_sign_by_cellset.csv", index=False)
print("\n=== POOLED ===")
print(R[R.lesion == "POOLED"][["cohort","cellset","n_cells","stromal_pct_overall","enrichment_fold",
      "depth_controlled_rho","defA_frac_retained","jaccard_vs_primary"]].to_string(index=False))
print("\n=== 候補TFの符号が2つの集合で一致した数 ===")
print(agree.to_string())
require_heartbeat("18 完了時")
print("\n18 DONE")
