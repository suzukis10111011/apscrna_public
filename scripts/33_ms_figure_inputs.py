import os
#!/usr/bin/env python
"""Purpose: Generate CSV inputs for manuscript figures from discovery/validation cohorts and computed
  permutation nulls, dip tests, and ingested UMAP coordinates
Input: {OUT}/all_cells_discovery.h5ad, all_cells_discovery_obs_scored.parquet, all_cells_validation.h5ad,
  stromal_discovery_labelled.h5ad, stromal_discovery_tf_activity.parquet, stromal_validation.h5ad
Output: {TAB}/33_compartment_marker_dotplot.csv, 33_defA_permutation_null.csv, 33_dip_by_lesion.csv,
  33_ingest_umap_validation.csv, 33_state_pct_8lesions.csv, 33_tf_validation_all730.csv, 33_umap_{tag}.csv

33: 論文図（results/figures_ms/）が読む CSV を作る。図は数値を必ず CSV から読む。
新規計算は4件:
  (1) 検証コホートでの 730 因子すべての pooled / 病変別 ρ（従来は候補76因子のみ）
  (2) 定義A（retained dispersion）の置換帰無帯
  (3) 病変ごとの dip 検定（2つのプログラムスコア）
  (4) 検証間質を探索の UMAP に投影した座標（05 の ingest と同一手順。05 は座標を破棄していた）
他は既存 h5ad / 表からの整形。
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("33")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr, rankdata
import diptest
from _utils import ulm_on_symbols

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")[["source", "target", "weight"]]
AX2 = ["lf_A_cytokine_chemokine", "lf_D_tissueforming"]
SEED = 0


def dctrl(x, y, dep, ng):
    M = np.column_stack([np.ones(len(x)), dep, ng])
    def res(v):
        b, *_ = np.linalg.lstsq(M, rankdata(v), rcond=None); return rankdata(v) - M @ b
    return spearmanr(res(x), res(y)).statistic


def matched(df, col, match, nbins=20, min_n=25):
    if len(df) < 500: nbins = 10
    q = pd.qcut(df[match].rank(method="first"), nbins, labels=False, duplicates="drop")
    tot = float(df[col].std(ddof=1)); num, den = 0.0, 0
    for _, g in df.groupby(q, observed=True):
        if len(g) < min_n: continue
        num += (len(g) - 1) * float(g[col].std(ddof=1)) ** 2; den += len(g) - 1
    return float(np.sqrt(num / den)) / tot if den and tot else np.nan


# ---------- (1) Validation: ρ for all 730 factors ----------
Sv = sc.read_h5ad(f"{OUT}/stromal_validation.h5ad")
Sd = sc.read_h5ad(f"{OUT}/stromal_discovery_labelled.h5ad")
ACTv = ulm_on_symbols(Sv, NET, tmin=5)
ACTd = pd.read_parquet(f"{OUT}/stromal_discovery_tf_activity.parquet")
S04 = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
tf730 = sorted(S04.TF.unique())
ov = Sv.obs; od = Sd.obs
rows = []
for axis in AX2:
    for tf in tf730:
        if tf not in ACTv.columns:
            rows.append(dict(axis=axis, TF=tf, in_validation=False)); continue
        x = ACTv[tf].values; y = ov[axis].values
        per = {}
        for les in sorted(ov.lesion.astype(str).unique()):
            m = (ov.lesion.astype(str) == les).values
            per[f"val_rho_{les}"] = round(dctrl(x[m], y[m], ov.logdepth.values[m], ov.n_genes.values[m]), 4)
        d0 = S04[(S04.axis == axis) & (S04.TF == tf)].iloc[0]
        pooled = dctrl(x, y, ov.logdepth.values, ov.n_genes.values)
        sgn = [np.sign(v) == np.sign(d0.pooled_rho) for v in per.values()]
        rows.append(dict(axis=axis, TF=tf, in_validation=True,
                         discovery_pooled_rho=round(float(d0.pooled_rho), 4),
                         validation_pooled_rho=round(float(pooled), 4),
                         discovery_is_candidate=bool(d0.is_candidate),
                         discovery_rank=float(d0.rank_within_axis) if pd.notna(d0.rank_within_axis) else np.nan,
                         same_sign_pooled=bool(np.sign(pooled) == np.sign(d0.pooled_rho)),
                         n_val_lesions_same_sign=int(sum(sgn)), **per))
A730 = pd.DataFrame(rows)
A730.to_csv(f"{TAB}/33_tf_validation_all730.csv", index=False)
ok = A730[A730.in_validation == True]
print(f"(1) 検証で活性が取れた因子 {ok.TF.nunique()}/{len(tf730)} | pooled 同符号 "
      + " / ".join(f"{a.split('_')[1]} {int(g.same_sign_pooled.sum())}/{len(g)}" for a, g in ok.groupby('axis')))

# ---------- (2) Permutation null for Definition A ----------
RELA_d = ACTd["RELA"].values; RELA_v = ACTv["RELA"].values
rng = np.random.default_rng(SEED)
NP = 200
rows = []
for coh, obs, rela in [("discovery_n5", od, RELA_d), ("validation_n3", ov, RELA_v)]:
    for les in sorted(obs.lesion.astype(str).unique()):
        m = (obs.lesion.astype(str) == les).values
        d = pd.DataFrame(dict(y=obs["lf_D_tissueforming"].values[m], RELA=rela[m]))
        obs_v = matched(d, "y", "RELA")
        null = np.array([matched(pd.DataFrame(dict(y=d.y.values, RELA=rng.permutation(d.RELA.values))), "y", "RELA")
                         for _ in range(NP)])
        rows.append(dict(cohort=coh, lesion=les, n=int(m.sum()), frac_retained=round(obs_v, 6),
                         null_mean=round(float(np.nanmean(null)), 6),
                         null_lo=round(float(np.nanquantile(null, 0.025)), 6),
                         null_hi=round(float(np.nanquantile(null, 0.975)), 6),
                         p_permutation=(1 + int((null <= obs_v).sum())) / (1 + NP),
                         n_perm=NP, one_sided="less", seed=SEED))
PN = pd.DataFrame(rows)
PN.to_csv(f"{TAB}/33_defA_permutation_null.csv", index=False)
print("(2) 定義A 置換帰無:"); print(PN[["cohort","lesion","frac_retained","null_mean","null_lo","null_hi","p_permutation"]].to_string(index=False))

# ---------- (3) Per-lesion dip test ----------
rows = []
for coh, obs in [("discovery_n5", od), ("validation_n3", ov)]:
    for les in sorted(obs.lesion.astype(str).unique()) + ["POOLED"]:
        m = np.ones(len(obs), bool) if les == "POOLED" else (obs.lesion.astype(str) == les).values
        for axis, lab in [(AX2[0], "Inflammatory program"), (AX2[1], "Tissue-forming program")]:
            v = obs[axis].values[m]
            rows.append(dict(cohort=coh, lesion=les, program=lab, n=int(m.sum()),
                             dip_statistic=round(float(diptest.diptest(v)[0]), 5),
                             dip_p=round(float(diptest.diptest(v)[1]), 4)))
DP = pd.DataFrame(rows)
DP.to_csv(f"{TAB}/33_dip_by_lesion.csv", index=False)
print(f"(3) dip: p の範囲 {DP.dip_p.min():.3f}–{DP.dip_p.max():.3f}、p<0.05 は {int((DP.dip_p<0.05).sum())}/{len(DP)}")

# ---------- (4) Project validation stroma onto the discovery UMAP (same procedure as in 05) ----------
_ds = Sd.var["symbol"].astype(str); _keep = ~_ds.duplicated(keep="first")
Dst = Sd[:, _keep.values].copy(); Dst.var_names = _ds[_keep].values
shared = [g for g in Dst.var_names if g in set(Sv.var_names)]
assert len(shared) > 5000, len(shared)
A = Dst[:, shared].copy(); B = Sv[:, shared].copy()
A.obsm["X_pca"] = A.obsm["X_pca_harmony"] if "X_pca_harmony" in A.obsm else A.obsm["X_pca"]
sc.tl.ingest(B, A, obs="state", embedding_method=("umap", "pca"))
pd.DataFrame(dict(cell=B.obs_names, lesion=B.obs.lesion.astype(str).values,
                  state_within=Sv.obs["state_within"].astype(str).values,
                  umap1=B.obsm["X_umap"][:, 0], umap2=B.obsm["X_umap"][:, 1])
             ).to_csv(f"{TAB}/33_ingest_umap_validation.csv", index=False)
print(f"(4) 投影座標 {B.n_obs} 細胞（共有遺伝子 {len(shared)}）")

# ---------- Formatting: UMAP, scores, markers, state proportions ----------
Od = pd.read_parquet(f"{OUT}/all_cells_discovery_obs_scored.parquet")
Ad = sc.read_h5ad(f"{OUT}/all_cells_discovery.h5ad", backed="r")
Av = sc.read_h5ad(f"{OUT}/all_cells_validation.h5ad", backed="r")
for tag, obj, sco in [("discovery", Ad, Od), ("validation", Av, None)]:
    u = obj.obsm["X_umap"]
    df = pd.DataFrame(dict(cohort=tag, cell=[str(x) for x in obj.obs_names],
                           lesion=obj.obs.lesion.astype(str).values,
                           compartment=obj.obs.compartment.astype(str).values,
                           umap1=u[:, 0], umap2=u[:, 1]))
    if sco is not None:
        df["inflammatory"] = sco["p_A_cytokine_chemokine"].values
        df["tissueforming"] = sco["p_D_tissueforming"].values
    else:
        df["inflammatory"] = obj.obs["lf_A_cytokine_chemokine"].values
        df["tissueforming"] = obj.obs["lf_D_tissueforming"].values
    df.to_csv(f"{TAB}/33_umap_{tag}.csv", index=False)
    print(f"    UMAP {tag}: {len(df)} 細胞")

MARK = {"Stromal cells": ["COL1A1", "DCN", "LUM"], "Immune (pan)": ["PTPRC"],
        "Myeloid cells": ["LYZ", "CD68"], "T cells": ["CD3D", "IL7R"],
        # IGKC is excluded here because it is excluded from the dot plot panel in Methods
        "Plasma cells": ["MZB1"], "Endothelium": ["PECAM1", "VWF"],
        "Pericytes/SMC": ["RGS5", "ACTA2"], "Epithelium": ["KRT14", "KRT5"]}
rows = []
for tag, obj, symcol in [("discovery", Ad, "symbol"), ("validation", Av, None)]:
    sym = obj.var[symcol].astype(str) if symcol and symcol in obj.var else pd.Series(obj.var_names, index=obj.var_names).astype(str)
    s2i = {}
    for gid, s in zip(obj.var_names, sym): s2i.setdefault(s, gid)
    comp = obj.obs.compartment.astype(str).values
    for lab, genes in MARK.items():
        for g in genes:
            if g not in s2i: continue
            v = np.asarray(obj[:, s2i[g]].X.todense()).ravel()
            for c in sorted(set(comp)):
                m = comp == c
                rows.append(dict(cohort=tag, marker_group=lab, gene=g, compartment=c, n_cells=int(m.sum()),
                                 mean_lognorm=round(float(v[m].mean()), 4),
                                 detection_pct=round(float(100 * (v[m] > 0).mean()), 2)))
pd.DataFrame(rows).to_csv(f"{TAB}/33_compartment_marker_dotplot.csv", index=False)
Ad.file.close(); Av.file.close()

sd = pd.read_csv(f"{TAB}/04b_state_pct_discovery.csv")
sv = pd.read_csv(f"{TAB}/25_validation_states_within_vs_ingest.csv")
sv = sv[sv.labelling == "within_cohort"]
# Both 04b and 25 have a POOLED row after the lesion rows. If this is left in, the 8-lesion table
# Row 9 is included, causing a 9th bar labeled "nan" to appear in Fig4c (this occurred once). Take only the lesion rows.
sd = sd[sd.lesion != "POOLED"]
sv = sv[sv.lesion != "POOLED"]
STATES = ["inflammatory", "intermediate", "tissue-forming", "other"]
rows = []
for _, r in sd.iterrows():
    rows.append(dict(cohort="discovery_n5", lesion=r["lesion"],
                     n_cells=int(sum(r[f"n_{s}"] for s in STATES if f"n_{s}" in r)),
                     **{s: float(r[f"pct_{s}"]) for s in STATES}))
for _, r in sv.iterrows():
    rows.append(dict(cohort="validation_n3", lesion=r["lesion"], n_cells=int(r["n_cells"]),
                     **{s: float(r[f"pct_{s}"]) for s in STATES}))
pd.DataFrame(rows).to_csv(f"{TAB}/33_state_pct_8lesions.csv", index=False)
print("\n33 DONE")
