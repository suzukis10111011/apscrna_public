#!/usr/bin/env python
"""Purpose: Compute axis scores (6-axis stage1 framework), TF activity (RELA regulon), and primary metrics for
  5 lesion types in discovery cohort; define stromal subset and generate sensitivity analyses.
Input: data/processed/all_cells_discovery.h5ad, frozen/genesets/stage1_geneset_definitions.json (6 axes,
  read-only), frozen/genesets/program_geneset_definitions.json (composite score sensitivity),
  frozen/net/collectri_regulons.csv (RELA 560 targets, TF activity)
Output: data/processed/stromal_discovery.h5ad, results/tables/04_*.csv (definition, enrichment, coverage,
  grant factors, correlation, TF screen)

04: Axis scores, TF activity, and primary metrics for the 5 discovery lesion types (1)(2)(3①)(4 discovery portion).
Input : data/processed/all_cells_discovery.h5ad
       frozen/genesets/stage1_geneset_definitions.json (6 axes, read-only)
       frozen/genesets/program_geneset_definitions.json (for composite score sensitivity analysis)
       frozen/net/collectri_regulons.csv (RELA 560 targets, TF activity)
Output : data/processed/stromal_discovery.h5ad
       results/tables/04_*.csv
Judgment call: There are two systems of definitions for the "inflammation score" and "tissue-forming score."
          We take as primary the 6 axes of stage1 named in the instructions, with inflammation = A_cytokine_chemokine
          and tissue-forming = D_tissueforming.
          The composite score from program_geneset_definitions (equal-weighted average after dividing by the
          frozen global_robust_sd) is additionally reported with the "_composite" suffix and treated as a
          sensitivity analysis.
Validation-side values are not used at all (this script reads only discovery data).
"""
import os, sys, json, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp
from scipy.stats import spearmanr, rankdata
warnings.filterwarnings("ignore")

import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
FZ = f"{W}/frozen"
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
# ===== 09: Data source switch (GEO matrix / kb requantification) =====
# APSC_SOURCE=geo (default) | kb (identical QC for all barcodes) | kbgeo (restricted to GEO barcodes)
# Keep the settings completely identical, only switching the input and output suffixes.
SRC = os.environ.get("APSC_SOURCE", "geo")
assert SRC in ("geo", "kb", "kbgeo", "kbtest"), SRC   # kbtest: for trying the pathway only with samples that are aligned
SFX = "" if SRC == "geo" else f"_{SRC}"
SEED = 0; np.random.seed(SEED)
sc.settings.verbosity = 1

S1 = json.load(open(f"{FZ}/genesets/stage1_geneset_definitions.json"))
PG = json.load(open(f"{FZ}/genesets/program_geneset_definitions.json"))
NET = pd.read_csv(f"{FZ}/net/collectri_regulons.csv")
RELA_TARGETS = set(NET.loc[NET.source == "RELA", "target"])
print("RELA targets:", len(RELA_TARGETS))

AX = {k: set(v) for k, v in S1["axes"].items()}
LF = {k: sorted(v - RELA_TARGETS) for k, v in AX.items()}          # leak-free variant
COMP = {k: v["genes"] for k, v in PG["components"].items()}
RSD = PG["global_robust_sd"]
INFL_C = PG["composite_inflammatory"]["components"]
TISS_C = PG["composite_tissueforming"]["components"]
print("軸:", {k: len(v) for k, v in AX.items()})
print("leak-free 後:", {k: len(v) for k, v in LF.items()})

D = sc.read_h5ad(f"{OUT}/all_cells_discovery{SFX}.h5ad")
print("all_cells_discovery:", D.shape, D.obs.lesion.value_counts().to_dict())
sym2id = {}
for gid, s in zip(D.var_names, D.var.symbol.astype(str)): sym2id.setdefault(s, gid)

def score(A, name, symbols):
    ids = [sym2id[s] for s in symbols if s in sym2id]
    if len(ids) < 5:
        print(f"  {name}: 遺伝子不足 ({len(ids)}) — skip"); return 0
    sc.tl.score_genes(A, ids, ctrl_size=50, n_bins=25, random_state=SEED, score_name=name)
    return len(ids)

cov = []
for k, v in AX.items():  cov.append(dict(set=f"p_{k}",  n_def=len(v), n_found=score(D, f"p_{k}",  sorted(v))))
for k, v in LF.items():  cov.append(dict(set=f"lf_{k}", n_def=len(v), n_found=score(D, f"lf_{k}", v)))
for k, v in COMP.items(): cov.append(dict(set=f"c_{k}", n_def=len(v), n_found=score(D, f"c_{k}", v)))
pd.DataFrame(cov).to_csv(f"{TAB}/04_geneset_coverage_discovery{SFX}.csv", index=False)

# Composite score: equal-weight average after dividing by frozen global robust SD (per frozen convention)
for lab, comps in [("inflammatory_composite", INFL_C), ("tissueforming_composite", TISS_C)]:
    D.obs[lab] = np.mean([D.obs[f"c_{c}"].values / RSD[c] for c in comps], axis=0)

D.obs["logdepth"] = np.log10(D.obs.total_counts.values)
D.obs["n_genes"] = D.obs.n_genes_by_counts.values
INFL, TISS = "p_A_cytokine_chemokine", "p_D_tissueforming"

# ---------- (1) Enrichment: stromal fraction in top 10% of tissue formation score / stromal fraction in entire lesion ----------
rows = []
for les, g in D.obs.groupby("lesion", observed=True):
    for lab, col in [("axis_D_tissueforming", TISS), ("composite_tissueforming", "tissueforming_composite")]:
        thr = np.quantile(g[col].values, 0.90)
        top = g[g[col].values >= thr]
        base = float((g.compartment == "stromal").mean())
        topf = float((top.compartment == "stromal").mean())
        rows.append(dict(lesion=les, cohort="discovery_n5", score=lab, n_cells=len(g), n_top_decile=len(top),
                         stromal_pct_overall=round(100*base, 2), stromal_pct_top_decile=round(100*topf, 2),
                         enrichment_fold=round(topf/base, 3) if base > 0 else np.nan))
# Since the 14.8% / 95.2% in the application document are pooled values from the discovery cohort, a pooled row is also output for comparison purposes only.
# (No combining across validation is performed.)
for lab, col in [("axis_D_tissueforming", TISS), ("composite_tissueforming", "tissueforming_composite")]:
    g = D.obs
    thr = np.quantile(g[col].values, 0.90); top = g[g[col].values >= thr]
    base = float((g.compartment == "stromal").mean()); topf = float((top.compartment == "stromal").mean())
    rows.append(dict(lesion="POOLED_discovery", cohort="discovery_n5", score=lab, n_cells=len(g), n_top_decile=len(top),
                     stromal_pct_overall=round(100*base, 2), stromal_pct_top_decile=round(100*topf, 2),
                     enrichment_fold=round(topf/base, 3) if base > 0 else np.nan))
E1 = pd.DataFrame(rows)
E1.to_csv(f"{TAB}/04_enrichment_by_lesion_discovery{SFX}.csv", index=False)
D.obs.to_parquet(f"{OUT}/all_cells_discovery_obs_scored{SFX}.parquet")
print("\n=== (1) 濃縮（病変別、コホート合算なし） ===")
print(E1[E1.score == "axis_D_tissueforming"][["lesion","n_cells","stromal_pct_overall","stromal_pct_top_decile","enrichment_fold"]].to_string(index=False))
print("全5病変で1倍超:", bool((E1[E1.score=="axis_D_tissueforming"].enrichment_fold > 1).all()))

# ---------- (2) Inflammation vs tissue formation ρ ----------
def dctrl(x, y, dep, ng):
    M = np.column_stack([np.ones(len(x)), dep, ng])
    def res(v):
        b, *_ = np.linalg.lstsq(M, rankdata(v), rcond=None); return rankdata(v) - M @ b
    return spearmanr(res(x), res(y)).statistic

rows = []
def add_rho(tag, sub, lesion):
    for lab, (a, b) in {"axis": (INFL, TISS), "composite": ("inflammatory_composite","tissueforming_composite")}.items():
        x, y = sub[a].values, sub[b].values
        rows.append(dict(cohort="discovery_n5", lesion=lesion, cellset=tag, definition=lab, n=len(sub),
                         raw_rho=round(spearmanr(x, y).statistic, 4),
                         depth_controlled_rho=round(dctrl(x, y, sub.logdepth.values, sub.n_genes.values), 4)))
for les, g in D.obs.groupby("lesion", observed=True): add_rho("all_cells", g, les)
add_rho("all_cells", D.obs, "POOLED_discovery")
strom_mask = (D.obs.compartment == "stromal").values
for les, g in D.obs[strom_mask].groupby("lesion", observed=True): add_rho("stromal", g, les)
add_rho("stromal", D.obs[strom_mask], "POOLED_discovery")
E2 = pd.DataFrame(rows); E2.to_csv(f"{TAB}/04_rho_inflam_vs_tissueforming_discovery{SFX}.csv", index=False)
print("\n=== (2) ρ（探索、8病変 pooled は出さない） ===")
print(E2[(E2.cellset=="stromal") & (E2.definition=="axis")][["lesion","n","raw_rho","depth_controlled_rho"]].to_string(index=False))

# ---------- Stromal subset: reintegration + Leiden 0.2-1.5 ----------
St = D[strom_mask].copy()
print("\n間質細胞:", St.shape, St.obs.lesion.value_counts().to_dict())
sc.pp.highly_variable_genes(St, n_top_genes=2000, flavor="seurat_v3", batch_key="lesion", layer="counts")
sc.pp.pca(St, n_comps=50, svd_solver="arpack", random_state=SEED, mask_var="highly_variable")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import harmony_run, load_grant_factors
harmony_run(St, "lesion", seed=SEED)
sc.pp.neighbors(St, n_neighbors=15, use_rep="X_pca_harmony", random_state=SEED)
sc.tl.umap(St, random_state=SEED)
for r in np.round(np.arange(0.2, 1.51, 0.1), 1):
    sc.tl.leiden(St, resolution=float(r), key_added=f"leiden_{r}", random_state=SEED, flavor="igraph", n_iterations=2)
print("leiden クラスタ数:", {f"{r}": St.obs[f"leiden_{r}"].nunique() for r in np.round(np.arange(0.2,1.51,0.1),1)})

# ---------- TF activity (ULM) ----------
from _utils import ulm_on_symbols
ACT = ulm_on_symbols(St, NET[["source","target","weight"]], tmin=5)
print("TF活性:", ACT.shape)
ACT.to_parquet(f"{OUT}/stromal_discovery_tf_activity{SFX}.parquet")
St.write(f"{OUT}/stromal_discovery{SFX}.h5ad")

# ---------- (3)① Discovery candidates (one-directional screen) ----------
LFCOLS = [f"lf_{k}" for k in AX]
obs = St.obs
res = []
for axis in LFCOLS:
    y_all = obs[axis].values
    for tf in ACT.columns:
        x_all = ACT[tf].values
        per = []
        for les in sorted(obs.lesion.unique()):
            m = (obs.lesion == les).values
            per.append(dctrl(x_all[m], y_all[m], obs.logdepth.values[m], obs.n_genes.values[m]))
        per = np.array(per)
        pooled = dctrl(x_all, y_all, obs.logdepth.values, obs.n_genes.values)
        res.append(dict(axis=axis, TF=tf, pooled_rho=pooled, abs_pooled=abs(pooled),
                        n_pos=int((per > 0).sum()), all5_same_sign=bool((per > 0).all() or (per < 0).all()),
                        **{f"rho_{l}": round(v, 4) for l, v in zip(sorted(obs.lesion.unique()), per)}))
S = pd.DataFrame(res)
S["is_candidate"] = False
for axis in LFCOLS:
    m = (S.axis == axis) & S.all5_same_sign
    top = S[m].sort_values("abs_pooled", ascending=False).head(20).index
    S.loc[top, "is_candidate"] = True
    S.loc[m, "rank_within_axis"] = S[m].abs_pooled.rank(ascending=False)
S.to_csv(f"{TAB}/04_tf_screen_discovery_candidates{SFX}.csv", index=False)
print("\n=== (3)① 探索候補（軸ごと上位20、5病変同符号のみ） ===")
print(S.groupby("axis").agg(n_tested=("TF","size"), n_sign_consistent=("all5_same_sign","sum"),
                            n_candidates=("is_candidate","sum")).to_string())

# The factor list for the application form is private (docs/private/grant_factors.txt). If absent, skip this block.
_gf = load_grant_factors(W)
if _gf:
    GRANT10 = [g for k in ("inflammatory", "tissueforming", "extra") for g in _gf.get(k, [])]
    g10 = S[S.TF.isin(GRANT10)].copy()
    _priv = os.environ.get("APSC_PRIVATE", f"{W}/docs/private/tables")
    os.makedirs(_priv, exist_ok=True)
    g10.to_csv(f"{_priv}/04_grant_factors_in_screen_discovery{SFX}.csv", index=False)
    print("\n=== (3)③ 申請書10因子+MYC が D_tissueforming 軸で候補入りするか ===")
    sub = g10[g10.axis == "lf_D_tissueforming"][["TF","pooled_rho","n_pos","all5_same_sign","rank_within_axis","is_candidate"]]
    print(sub.sort_values("pooled_rho").to_string(index=False))
else:
    print("\n(3)③ skipped: docs/private/grant_factors.txt is absent (grant material is not published)")

# ---------- (4) Definition A (discovery portion) ----------
def matched_residual(df, col, rela, nbins=20, min_n=25):
    d = df.dropna(subset=[col, rela]).copy()
    d["bin"] = pd.qcut(d[rela], nbins, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin") if len(g) >= min_n]
    n = np.array([b[0] for b in B]); s = np.array([b[1] for b in B])
    within = float(np.sqrt(np.average(s**2, weights=n)))
    return dict(n_cells=len(d), n_bins=len(B), median_bin_n=int(np.median(n)),
                total_sd=tot, matched_sd=within, frac_retained=within/tot)
obs2 = obs.copy(); obs2["RELA"] = ACT["RELA"].values
rows = []
for les in sorted(obs2.lesion.unique()):
    r = matched_residual(obs2[obs2.lesion == les], "lf_D_tissueforming", "RELA")
    r.update(lesion=les, cohort="discovery_n5",
             sd_reduction_pct=100*(1-r["frac_retained"]),
             variance_retained=r["frac_retained"]**2,
             variance_reduction_pct=100*(1-r["frac_retained"]**2))
    rows.append(r)
E4 = pd.DataFrame(rows); E4.to_csv(f"{TAB}/04_definitionA_discovery_v2{SFX}.csv", index=False)
print(f"\n=== (4) 定義A（探索、入力 = {'GEO 行列' if SRC == 'geo' else 'kb 再定量 (' + SRC + ')'}） ===")
print(E4[["lesion","n_cells","n_bins","total_sd","matched_sd","frac_retained","sd_reduction_pct","variance_reduction_pct"]].round(4).to_string(index=False))
print("\n04 DONE")
