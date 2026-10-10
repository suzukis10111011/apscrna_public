#!/usr/bin/env python
"""Purpose: Reproduce TF regulatory signs and retained dispersion in GSE164241 gingival samples; skip cross-
  dataset level comparison due to gene-space divergence.
Input: 04_tf_screen_discovery_candidates.csv, collectri_regulons.csv
Output: 13_geneset_coverage_gse164241.csv, 13_gse164241_summary.csv, 13_gse164241_tf_sign_reproduction.csv

13 (Fig7): Apply the same procedure to gingival specimens from the periodontitis dataset GSE164241 (oral mucosa atlas).
Since the gene space differs, no level comparison (state ratios) is performed. What is examined is sign concordance and ρ / retained dispersion.
"""
import sys, os, re, io, tarfile, gzip, json, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("13")
import scanpy as sc, anndata as ad, pandas as pd, numpy as np, scipy.io as sio, scipy.sparse as sp
from scipy.stats import spearmanr, rankdata
from _utils import annotate_compartment, ulm_on_symbols

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables"); D13 = f"{W}/data/gse164241"; SEED = 0
os.makedirs(D13, exist_ok=True)
TAR = f"{D13}/GSE164241_RAW.tar"
if not os.path.exists(TAR):
    print("RAW.tar 取得中 (796 MB)")
    subprocess.run(["curl","-sS","-C","-","-o",TAR,
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164241/suppl/GSE164241_RAW.tar"], check=True)
print("tar:", os.path.getsize(TAR), "B")

# --- GSM metadata (tissue / disease state) ---
META = f"{D13}/gsm_meta.csv"
if not os.path.exists(META):
    s = subprocess.run(["curl","-s","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241&targ=self&form=text&view=brief"],
                       capture_output=True, text=True).stdout
    gsms = [l.split("= ")[1].strip() for l in s.split("\n") if l.startswith("!Series_sample_id")]
    rows = []
    for g in gsms:
        t = subprocess.run(["curl","-s",f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=self&form=text&view=brief"],
                           capture_output=True, text=True).stdout
        d = dict(gsm=g)
        for l in t.split("\n"):
            if l.startswith("!Sample_title"): d["title"] = l.split("= ",1)[1].strip()
            elif l.startswith("!Sample_characteristics_ch1"):
                kv = l.split("= ",1)[1].split(": ",1)
                if len(kv) == 2: d[kv[0].strip().replace(" ","_")] = kv[1].strip()
        rows.append(d)
    pd.DataFrame(rows).to_csv(META, index=False)
M = pd.read_csv(META)
print("GSM 数:", len(M), "| 列:", list(M.columns))
tcol = next((c for c in M.columns if c.lower() == "tissue"), None)
dcol = next((c for c in M.columns if "disease" in c.lower()), None)
print(M[[c for c in ("gsm","title",tcol,dcol) if c]].to_string(index=False))
if tcol is None:
    print("tissue 列が無い — 13 を中止"); raise SystemExit(0)
G = M[M[tcol].astype(str).str.contains("gingiv", case=False, na=False)]
print(f"\n歯肉検体: {len(G)}")
if len(G) == 0:
    print("歯肉検体が無い — 13 を中止"); raise SystemExit(0)

# --- Matrix loading ---
H5 = f"{D13}/gingiva.h5ad"
if not os.path.exists(H5):
    with tarfile.open(TAR) as tf:
        names = [n for n in tf.getnames() if not os.path.basename(n).startswith("._")]
        parts = []
        for _, r in G.iterrows():
            pre = [n for n in names if n.startswith(r.gsm)]
            f_bc = next((n for n in pre if "barcodes" in n), None)
            f_ft = next((n for n in pre if "features" in n or "genes" in n), None)
            f_mx = next((n for n in pre if "matrix" in n), None)
            if not (f_bc and f_ft and f_mx): print(f"  {r.gsm}: 3点揃わず — 除外"); continue
            bc = gzip.open(tf.extractfile(f_bc), "rt").read().split()
            ft = [l.split("\t") for l in gzip.open(tf.extractfile(f_ft), "rt").read().rstrip("\n").split("\n")]
            X = sio.mmread(gzip.open(tf.extractfile(f_mx), "rb")).T.tocsr()
            v = pd.DataFrame({"gene_id": [f[0] for f in ft],
                              "symbol": [f[1] if len(f) > 1 else f[0] for f in ft]})
            a = ad.AnnData(X, obs=pd.DataFrame(index=[f"{r.gsm}_{b}" for b in bc]), var=v.set_index("gene_id"))
            a.obs["donor"] = r.title; a.obs["gsm"] = r.gsm
            a.obs["disease"] = str(r[dcol]) if dcol else "NA"
            parts.append(a); print(f"  {r.gsm} {r.title}: {a.shape}")
        A = ad.concat(parts, join="inner", index_unique=None)
        A.var["symbol"] = parts[0].var["symbol"].reindex(A.var_names).values
    A.var_names_make_unique()
    A.var["mt"] = A.var["symbol"].astype(str).str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    A = A[(A.obs.n_genes_by_counts >= 200) & (A.obs.total_counts >= 500) & (A.obs.pct_counts_mt <= 20)].copy()
    sc.pp.filter_genes(A, min_cells=3)
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat_v3", layer="counts", batch_key="donor")
    sc.pp.pca(A, n_comps=50, random_state=SEED)
    sc.pp.neighbors(A, n_neighbors=15, random_state=SEED); sc.tl.leiden(A, resolution=0.5, key_added="leiden_0.5", random_state=SEED)
    A.obs["lesion"] = A.obs["donor"].astype(str)
    A, _v, pct = annotate_compartment(A, "leiden_0.5", "gse164241", TAB)
    print("\n区画比率 (%):"); print(pct.to_string())
    A.write(H5)
else:
    A = sc.read_h5ad(H5)
print("\n統合後:", A.shape)

AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
RELA_T = set(NET.loc[NET.source == "RELA", "target"])
LF = {k: sorted(set(v) - RELA_T) for k, v in AXJ["axes"].items()}
s2i = {}
for gid, s in zip(A.var_names, A.var["symbol"].astype(str)): s2i.setdefault(s, gid)
cov = []
for k, v in LF.items():
    ids = [s2i[g] for g in v if g in s2i]
    sc.tl.score_genes(A, ids, ctrl_size=50, n_bins=25, random_state=SEED, score_name=f"lf_{k}")
    cov.append(dict(set=f"lf_{k}", n_def=len(v), n_found=len(ids)))
pd.DataFrame(cov).to_csv(f"{TAB}/13_geneset_coverage_gse164241.csv", index=False)
print(pd.DataFrame(cov).to_string(index=False))

S = A[A.obs.compartment == "stromal"].copy()
print("間質:", S.shape, "| ドナー別:", S.obs.donor.value_counts().to_dict())
ACT = ulm_on_symbols(S, NET[["source","target","weight"]], tmin=5)
O = S.obs.copy(); O["RELA"] = ACT["RELA"].values if "RELA" in ACT.columns else np.nan
O["logdepth"] = np.log10(O.total_counts.astype(float) + 1); O["n_genes"] = O.n_genes_by_counts.astype(float)
INF, TIS = "lf_A_cytokine_chemokine", "lf_D_tissueforming"

def dctrl(x, y, dep, ng):
    Mx = np.column_stack([np.ones(len(x)), dep, ng])
    def res(v):
        b, *_ = np.linalg.lstsq(Mx, rankdata(v), rcond=None); return rankdata(v) - Mx @ b
    return spearmanr(res(x), res(y)).statistic

def matched(df, col, rela, min_n=25):
    nb = 10 if len(df) < 500 else 20
    d = df.dropna(subset=[col, rela]).copy()
    if len(d) < 60: return dict(n_bins=0, frac_retained=np.nan, bins_used=nb)
    d["bin"] = pd.qcut(d[rela], nb, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin", observed=True) if len(g) >= min_n]
    if not B: return dict(n_bins=0, frac_retained=np.nan, bins_used=nb)
    n = np.array([b[0] for b in B]); s = np.array([b[1] for b in B])
    w = float(np.sqrt(np.average(s**2, weights=n)))
    return dict(n_bins=len(B), total_sd=round(tot,4), matched_sd=round(w,4),
                frac_retained=round(w/tot,4), bins_used=nb)

rows = []
for grp in ["ALL"] + sorted(O.disease.astype(str).unique()):
    sub = O if grp == "ALL" else O[O.disease.astype(str) == grp]
    for don in sorted(sub.donor.astype(str).unique()) + ["POOLED"]:
        g = sub if don == "POOLED" else sub[sub.donor.astype(str) == don]
        if len(g) < 60: continue
        r = dict(group=grp, donor=don, n_cells=len(g),
                 raw_rho=round(float(spearmanr(g[INF], g[TIS]).statistic), 4),
                 depth_controlled_rho=round(float(dctrl(g[INF].values, g[TIS].values,
                                                        g.logdepth.values, g.n_genes.values)), 4))
        r.update({f"defA_{k}": v for k, v in matched(g, TIS, "RELA").items()})
        rows.append(r)
R = pd.DataFrame(rows).drop_duplicates(subset=["group","donor"])
R.to_csv(f"{TAB}/13_gse164241_summary.csv", index=False)
print("\n=== ドナー別 ==="); print(R.to_string(index=False))

CAND = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
CAND = CAND[CAND.is_candidate] if "is_candidate" in CAND.columns else CAND
srows = []
for grp in ["ALL"] + sorted(O.disease.astype(str).unique()):
    sub = O if grp == "ALL" else O[O.disease.astype(str) == grp]
    dons = sorted(sub.donor.astype(str).unique())
    for axis in CAND.axis.unique():
        tfs = CAND[CAND.axis == axis].TF.tolist()
        disc_sign = CAND.set_index(["axis","TF"]).pooled_rho.to_dict()
        for tf in tfs:
            if tf not in ACT.columns: 
                srows.append(dict(group=grp, axis=axis, TF=tf, present=False)); continue
            a = ACT[tf].values
            rs = []
            for d in dons:
                m = (sub.donor.astype(str) == d).values
                idx = O.index.get_indexer(sub.index[m])
                if m.sum() < 30: continue
                rs.append(spearmanr(a[idx], sub[axis].values[m]).statistic)
            ds = disc_sign.get((axis, tf), np.nan)
            srows.append(dict(group=grp, axis=axis, TF=tf, present=True, n_donors=len(rs),
                              pooled_rho_gse164241=round(float(spearmanr(
                                  a[O.index.get_indexer(sub.index)], sub[axis].values).statistic), 4),
                              discovery_pooled_rho=round(float(ds), 4) if ds == ds else np.nan,
                              all_donors_same_sign=bool(rs and (all(r > 0 for r in rs) or all(r < 0 for r in rs))),
                              sign_matches_discovery=bool(ds == ds and rs and
                                  np.sign(spearmanr(a[O.index.get_indexer(sub.index)], sub[axis].values).statistic) == np.sign(ds))))
SG = pd.DataFrame(srows); SG.to_csv(f"{TAB}/13_gse164241_tf_sign_reproduction.csv", index=False)
print("\n=== 探索候補TFの符号再現 ===")
print(SG[SG.present == True].groupby(["group","axis"]).agg(
    n=("TF","size"), same_sign_all_donors=("all_donors_same_sign","sum"),
    matches_discovery=("sign_matches_discovery","sum")).to_string())
require_heartbeat("13 完了時")
print("\n13 DONE")
