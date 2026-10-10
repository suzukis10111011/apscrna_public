import os
#!/usr/bin/env python
"""Purpose: Recompute leave-one-lesion-out cross-validation using 10 transcription factors defined by rule,
  including RELA and NFKB1 for Figure 3E
Input: stromal_discovery.h5ad, stromal_discovery_tf_activity.parquet, 04_tf_screen_discovery_candidates.csv
Output: 34_lolo_rulebased.csv

34: leave-one-lesion-out を「規則で定義した10因子」で再計算する。
28 の LOLO は旧セット（ATF1/MTA1 を含む）で走っており、RELA・NFKB1 を含まないため
論文図3E に使えなかった（表に無い因子が 0 本として描かれ、表題も誤った）。
推定量・候補規則は 28 と逐語で同じ。
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("34")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr, rankdata

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
INFL5 = ["RELB", "NFKB2", "RELA", "NFKB1", "REL"]          # Rule: inflammation-axis candidate x NF-kB family
TISS5 = ["SP7", "RUNX2", "DLX5", "MYRF", "SMAD3"]          # Rule: tissue-forming axis candidate x osteoblast lineage


def dctrl(x, y, dep, ng):
    M = np.column_stack([np.ones(len(x)), dep, ng])
    def res(v):
        b, *_ = np.linalg.lstsq(M, rankdata(v), rcond=None); return rankdata(v) - M @ b
    return spearmanr(res(x), res(y)).statistic


St = sc.read_h5ad(f"{OUT}/stromal_discovery.h5ad")
ACT = pd.read_parquet(f"{OUT}/stromal_discovery_tf_activity.parquet")
obs = St.obs
assert len(ACT) == len(obs), (ACT.shape, obs.shape)
S04 = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
chk = S04[(S04.axis == "lf_D_tissueforming") & (S04.TF == "SP7")].iloc[0]
got = dctrl(ACT["SP7"].values, obs["lf_D_tissueforming"].values, obs.logdepth.values, obs.n_genes.values)
assert abs(got - chk.pooled_rho) < 1e-6, f"dctrl 再現失敗: {got} vs {chk.pooled_rho}"
print(f"dctrl 再現確認 OK (SP7: {got:.6f}) | 因子 {ACT.shape[1]}")

les_all = sorted(obs.lesion.astype(str).unique())
lolo = []
for held in les_all:
    keep = (obs.lesion.astype(str) != held).values
    o = obs[keep]; A = ACT[keep]
    inc = sorted(o.lesion.astype(str).unique())
    lesarr = o.lesion.astype(str).values
    masks = {l: (lesarr == l) for l in inc}
    for axis in ("lf_A_cytokine_chemokine", "lf_D_tissueforming"):
        y = o[axis].values
        recs = []
        for tf in A.columns:
            x = A[tf].values
            per = np.array([dctrl(x[m], y[m], o.logdepth.values[m], o.n_genes.values[m])
                            for m in masks.values()])
            pooled = dctrl(x, y, o.logdepth.values, o.n_genes.values)
            recs.append(dict(TF=tf, pooled_rho=pooled, abs_pooled=abs(pooled),
                             same_sign=bool((per > 0).all() or (per < 0).all())))
        R = pd.DataFrame(recs)
        elig = R[R.same_sign].sort_values("abs_pooled", ascending=False)
        cand = set(elig.head(20).TF)
        rank = {t: i + 1 for i, t in enumerate(elig.TF)}
        for tf in (INFL5 if axis.startswith("lf_A") else TISS5):
            r = R[R.TF == tf]
            lolo.append(dict(held_out_lesion=held, n_lesions_used=len(inc), axis=axis, TF=tf,
                             pooled_rho=round(float(r.pooled_rho.iloc[0]), 4) if len(r) else np.nan,
                             sign_consistent=bool(r.same_sign.iloc[0]) if len(r) else False,
                             rank_within_axis=rank.get(tf, np.nan), is_candidate=tf in cand))
    print(f"  LOLO {held} 完了")
L = pd.DataFrame(lolo)
L.to_csv(f"{TAB}/34_lolo_rulebased.csv", index=False)
per = L.groupby(["axis", "TF"]).is_candidate.agg(["sum", "size"]).reset_index()
per["all_runs"] = per["sum"] == per["size"]
print(f"\n全 {len(les_all)} 回とも候補: {int(per.all_runs.sum())}/{len(per)}")
print(per.to_string(index=False))
bad = per[~per.all_runs]
if len(bad):
    print("\n外れた因子の各回順位:")
    print(L[L.TF.isin(bad.TF)][["TF", "held_out_lesion", "rank_within_axis", "is_candidate"]].to_string(index=False))
require_heartbeat("34 完了時")
print("\n34 DONE")
