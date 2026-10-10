#!/usr/bin/env python
"""Purpose: Preprocess, integrate, and annotate compartments for 5 discovery lesions from GEO (CellRanger
  matrices with mismatched gene references).
Input: data/geo/raw/GSM*_Sample{A..E}_{barcodes,features,matrix}.mtx.gz (SampleA: 33,694 genes; B–E: 33,538;
  aligned to 32,991 common Ensembl IDs).
Output: data/processed/all_cells_discovery.h5ad; results/tables/qc_cells_by_lesion_discovery.csv;
  results/tables/compartment_by_lesion_discovery.csv;
  results/tables/compartment_marker_zscores_discovery.csv.

03: Preprocessing, integration, and compartment annotation for the 5 discovery lesions (GEO processed matrices).
Input : data/geo/raw/GSM*_Sample{A..E}_{barcodes,features,matrix}
Output: data/processed/all_cells_discovery.h5ad
       results/tables/qc_cells_by_lesion_discovery.csv
       results/tables/compartment_by_lesion_discovery.csv
Judgment call: Only SampleA uses the old CellRanger reference (33,694 genes, 2-column features). B-E have 33,538.
          Integrate by aligning to the common set of 32,991 Ensembl gene IDs. Adopt the symbol annotations from the B-E side.
Defaults  : genes>=200, UMI>=500, mito<=20%, genes expressed in >=3 cells, scrublet 0.06, random_state=0
"""
import os, sys, glob, gzip, json
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.io import mmread
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # Make _utils readable before load

import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
RAW = f"{W}/data/geo/raw"
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
# ===== 09: Data source switch (GEO matrix / kb requantification) =====
# APSC_SOURCE=geo (default) | kb (identical QC for all barcodes) | kbgeo (restricted to GEO barcodes)
# Keep the settings completely identical, only switching the input and output suffixes.
SRC = os.environ.get("APSC_SOURCE", "geo")
assert SRC in ("geo", "kb", "kbgeo", "kbtest"), SRC   # kbtest: for trying the pathway only with samples that are aligned
SFX = "" if SRC == "geo" else f"_{SRC}"
os.makedirs(OUT, exist_ok=True); os.makedirs(TAB, exist_ok=True)
sc.settings.verbosity = 1
SEED = 0
np.random.seed(SEED)

LES = {"SampleA": "A", "SampleB": "B", "SampleC": "C", "SampleD": "D", "SampleE": "E"}
if SRC == "kbtest":
    # To allow verifying the wiring even partway through 08, restrict to samples that have verify.json.
    import glob as _g
    _done = {os.path.basename(os.path.dirname(x)) for x in _g.glob(f"{W}/data/counts_discovery/Sample*/verify.json")}
    LES = {k: v for k, v in LES.items() if k in _done}
    print(f"kbtest: 検体 {sorted(LES)} のみ使用")

def load_geo(sample):
    bc = glob.glob(f"{RAW}/*_{sample}_barcodes.tsv.gz")[0]
    ft = glob.glob(f"{RAW}/*_{sample}_features.tsv.gz")[0]
    mx = glob.glob(f"{RAW}/*_{sample}_matrix.mtx.gz")[0]
    with gzip.open(bc, "rt") as f: bcs = [l.strip() for l in f]
    with gzip.open(ft, "rt") as f:
        rows = [l.rstrip("\n").split("\t") for l in f]
    gid = [r[0] for r in rows]; sym = [r[1] if len(r) > 1 else r[0] for r in rows]
    with gzip.open(mx, "rt") as f: M = mmread(f)
    A = ad.AnnData(sp.csr_matrix(M.T.astype(np.float32)))            # cells x genes
    A.var_names = gid; A.var["symbol"] = sym
    A.obs_names = [f"{LES[sample]}_{b}" for b in bcs]
    A.obs["lesion"] = LES[sample]; A.obs["cohort"] = "discovery"
    A.obs["gsm"] = os.path.basename(bc).split("_")[0]
    A.var_names_make_unique()
    return A


def load_kb(sample):
    """08 の kb 出力。var_names は symbol（同一 symbol の gene id は合算済み）。
    GEO 行列は作成者側で cell calling 済みだが kb は counts_unfiltered なので、
    cell calling は 03 の QC 閾値（>=500 UMI, >=200 遺伝子, MT<=20%）が兼ねる。
    kbgeo はさらに GEO 行列に含まれるバーコードだけに絞り、定量の差だけを見る。"""
    from _utils import load_kb_counts
    A = load_kb_counts(f"{W}/data/counts_discovery/{sample}/counts_unfiltered", root=W)
    A.obs_names = [f"{LES[sample]}_{b}" for b in A.obs_names]
    A.obs["lesion"] = LES[sample]; A.obs["cohort"] = "discovery"; A.obs["gsm"] = sample
    A.var["symbol"] = list(A.var_names)
    if SRC == "kbgeo":
        bc = glob.glob(f"{RAW}/*_{sample}_barcodes.tsv.gz")[0]
        # GEO barcodes have the CellRanger "-1" suffix, kb barcodes are the raw 16mer. Strip the suffix to match them
        # (leaving it attached results in 0 matches. Confirmed empirically in SampleA: 8,107/8,107 = 100% match).
        with gzip.open(bc, "rt") as f:
            keep = {l.strip().split("-")[0] for l in f}
        m = [o.split("_", 1)[1] in keep for o in A.obs_names]
        assert sum(m) > 0.9 * len(keep), f"{sample}: GEO バーコードとの一致が少なすぎる {sum(m)}/{len(keep)}"
        print(f"  {sample}: GEO バーコードに限定 {sum(m)} / {A.n_obs}（GEO 側 {len(keep)}）")
        A = A[m].copy()
    else:
        # Filter by cell-side QC thresholds (>=500 UMI, >=200 genes) before concatenation. Since QC is per-cell,
        # the result is identical to the original 03 order (concatenate -> QC -> filter_genes), without having to hold 491k x 5.
        tot = np.asarray(A.X.sum(1)).ravel(); ng = np.asarray((A.X > 0).sum(1)).ravel()
        mm = (tot >= 500) & (ng >= 200)
        print(f"  {sample}: 連結前 QC {int(mm.sum())} / {A.n_obs} バーコード")
        A = A[mm].copy()
    return A


load = load_geo if SRC == "geo" else load_kb

print(f"== load (source={SRC}) ==")
ads = {s: load(s) for s in LES}
for s, A in ads.items(): print(f"  {s}: {A.shape}")

common = sorted(set.intersection(*[set(A.var_names) for A in ads.values()]))
# In the kb series, var_names are symbols, so output the ID-space name according to the source
print(("共通 symbol:" if SRC != "geo" else "共通 Ensembl ID:"), len(common))
sym_ref = ads["SampleB"].var["symbol"].to_dict()      # Adopt the newer annotation
for s in ads: ads[s] = ads[s][:, common].copy()

D = ad.concat(list(ads.values()), join="outer", label=None, index_unique=None)
D.var["symbol"] = [sym_ref.get(g, g) for g in D.var_names]
D.var["ensembl_id"] = D.var_names
print("concat:", D.shape, "| 病変別:", D.obs.lesion.value_counts().to_dict())

# ---- QC ----
mito = D.var.symbol.astype(str).str.startswith("MT-").values
D.var["mt"] = mito
sc.pp.calculate_qc_metrics(D, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
rec = []
for les, g in D.obs.groupby("lesion", observed=True):
    rec.append(dict(lesion=les, stage="0_raw", n_cells=len(g),
                    median_genes=float(g.n_genes_by_counts.median()),
                    median_umi=float(g.total_counts.median()),
                    median_pct_mt=float(g.pct_counts_mt.median())))
keep = (D.obs.n_genes_by_counts >= 200) & (D.obs.total_counts >= 500) & (D.obs.pct_counts_mt <= 20)
for les, g in D.obs[keep.values].groupby("lesion", observed=True):
    rec.append(dict(lesion=les, stage="1_qc_pass", n_cells=len(g),
                    median_genes=float(g.n_genes_by_counts.median()),
                    median_umi=float(g.total_counts.median()),
                    median_pct_mt=float(g.pct_counts_mt.median())))
D = D[keep.values].copy()
sc.pp.filter_genes(D, min_cells=3)
print("QC後:", D.shape)

# ---- doublet (per lesion) ----
D.obs["doublet_score"] = np.nan; D.obs["predicted_doublet"] = False
for les in sorted(D.obs.lesion.unique()):
    m = (D.obs.lesion == les).values
    sub = D[m].copy()
    try:
        sc.pp.scrublet(sub, expected_doublet_rate=0.06, random_state=SEED)
        D.obs.loc[m, "doublet_score"] = sub.obs["doublet_score"].values
        D.obs.loc[m, "predicted_doublet"] = sub.obs["predicted_doublet"].values
        print(f"  scrublet {les}: {int(sub.obs.predicted_doublet.sum())}/{m.sum()} doublet")
    except Exception as e:
        print(f"  scrublet {les} FAILED: {type(e).__name__}: {e}")
for les, g in D.obs[~D.obs.predicted_doublet.astype(bool)].groupby("lesion", observed=True):
    rec.append(dict(lesion=les, stage="2_singlet", n_cells=len(g),
                    median_genes=float(g.n_genes_by_counts.median()),
                    median_umi=float(g.total_counts.median()),
                    median_pct_mt=float(g.pct_counts_mt.median())))
D = D[~D.obs.predicted_doublet.astype(bool)].copy()
print("doublet除去後:", D.shape)
pd.DataFrame(rec).to_csv(f"{TAB}/qc_cells_by_lesion_discovery{SFX}.csv", index=False)

# ---- normalization ----
D.layers["counts"] = D.X.copy()
sc.pp.normalize_total(D, target_sum=1e4)
sc.pp.log1p(D)
D.raw = D

# ---- integration (discovery only) ----
sc.pp.highly_variable_genes(D, n_top_genes=2000, flavor="seurat_v3", batch_key="lesion", layer="counts")
sc.pp.pca(D, n_comps=50, svd_solver="arpack", random_state=SEED, mask_var="highly_variable")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import harmony_run
harmony_run(D, "lesion", seed=SEED)
sc.pp.neighbors(D, n_neighbors=15, use_rep="X_pca_harmony", random_state=SEED)
sc.tl.umap(D, random_state=SEED)
for r in (0.3, 0.5, 1.0):
    sc.tl.leiden(D, resolution=r, key_added=f"leiden_{r}", random_state=SEED, flavor="igraph", n_iterations=2)
    print(f"  leiden {r}: {D.obs[f'leiden_{r}'].nunique()} clusters")

# ---- compartment annotation ----
MARK = {
 "stromal":        ["COL1A1","COL1A2","DCN","LUM","PDGFRB"],
 "endothelial":    ["PECAM1","VWF","CDH5"],
 "pericyte_smc":   ["RGS5","ACTA2","MYH11"],
 "immune":         ["PTPRC","CD3E","MS4A1","CD79A","JCHAIN","LYZ","CD68","CD14","S100A8","KIT"],
 "epithelial":     ["KRT14","KRT5","EPCAM"],
 "schwann_neural": ["S100B","PLP1"],
 "erythrocyte":    ["HBB"],
}
sym2id = {}
for gid, s in zip(D.var_names, D.var.symbol.astype(str)): sym2id.setdefault(s, gid)
CLKEY = "leiden_0.5"
Z = sc.get.obs_df(D, keys=[sym2id[g] for gs in MARK.values() for g in gs if g in sym2id] + [CLKEY])
Z = Z.rename(columns={sym2id[g]: g for gs in MARK.values() for g in gs if g in sym2id})
means = Z.groupby(CLKEY, observed=True).mean()
zs = (means - means.mean()) / means.std(ddof=0).replace(0, np.nan)
score = pd.DataFrame({k: zs[[g for g in gs if g in zs.columns]].mean(1) for k, gs in MARK.items()})
assign = score.idxmax(1)
D.obs["compartment"] = D.obs[CLKEY].map(assign).astype(str)
score.assign(assigned=assign).to_csv(f"{TAB}/compartment_marker_zscores_discovery{SFX}.csv")
ct = pd.crosstab(D.obs.lesion, D.obs.compartment)
frac = (ct.T / ct.sum(1)).T
out = ct.copy(); out.columns = [f"n_{c}" for c in ct.columns]
for c in frac.columns: out[f"pct_{c}"] = (100 * frac[c]).round(2)
out["n_total"] = ct.sum(1)
out.to_csv(f"{TAB}/compartment_by_lesion_discovery{SFX}.csv")
print("\n区画×病変 (%):"); print((100 * frac).round(1).to_string())

D.write(f"{OUT}/all_cells_discovery{SFX}.h5ad")
print("\nsaved:", f"{OUT}/all_cells_discovery{SFX}.h5ad", D.shape)
json.dump(dict(n_cells=int(D.n_obs), n_genes=int(D.n_vars),
               common_ensembl=len(common),
               stromal_pct_by_lesion={k: float(v) for k, v in (100*frac.get("stromal", pd.Series(dtype=float))).items()}),
          open(f"{TAB}/discovery_summary{SFX}.json", "w"), indent=1)
