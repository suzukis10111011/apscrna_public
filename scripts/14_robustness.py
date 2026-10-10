#!/usr/bin/env python
"""Purpose: Evaluate robustness across gene sets, scoring methods, and integration approaches using Enrichr-
  sourced annotations
Input: collectri_regulons.csv, stromal_discovery_labelled.h5ad
Output: 14_geneset_provenance.csv, 14_robustness_genesets_scoring.csv

14: Robustness. (a) gene sets 3x3, (b) 3 scoring methods, (c) 2 integration methods.
Gene sets are obtained via Enrichr rather than gseapy's Msigdb (gsea-msigdb.org is
outside the sandbox allowlist. Self-judgment ③). Record the library name and retrieval date.
"""
import sys, os, json, datetime, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("14")
import scanpy as sc, pandas as pd, numpy as np
from scipy.stats import spearmanr, rankdata
from _utils import ulm_on_symbols, annotate_compartment

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables"); REF = f"{W}/ref"; SEED = 0
S = sc.read_h5ad(f"{W}/data/processed/stromal_discovery_labelled.h5ad")
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
RELA_T = set(NET.loc[NET.source == "RELA", "target"])

def enrichr_gmt(lib):
    p = f"{REF}/{lib}.gmt"
    if not os.path.exists(p):
        urllib.request.urlretrieve(
            f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={lib}", p)
    d = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t"); d[f[0]] = [g.split(",")[0].strip().upper() for g in f[2:] if g.strip()]
    return d

WANT = {"HALLMARK_TNFA_SIGNALING_VIA_NFKB": ("MSigDB_Hallmark_2020", "tnf-alpha signaling via nf-kb"),
        "GOBP_INFLAMMATORY_RESPONSE": ("GO_Biological_Process_2025", "inflammatory response (go:0006954)"),
        "GOBP_OSSIFICATION": ("GO_Biological_Process_2025", "ossification (go:0001503)"),
        "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION": ("Reactome_Pathways_2024", "extracellular matrix organization")}
EXT, prov = {}, []
for name, (lib, needle) in WANT.items():
    try:
        g = enrichr_gmt(lib)
    except Exception as e:
        prov.append(dict(set=name, library=lib, status=f"取得失敗 {type(e).__name__}", n_genes=0)); continue
    hit = [t for t in g if needle in t.lower()]
    if not hit:
        prov.append(dict(set=name, library=lib, status="該当語なし", n_genes=0)); continue
    EXT[name] = sorted(set(g[hit[0]]))
    prov.append(dict(set=name, library=lib, matched_term=hit[0], n_genes=len(EXT[name]), status="OK",
                     source="Enrichr geneSetLibrary API（gsea-msigdb.org は許可リスト外のため代替）",
                     downloaded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")))
P = pd.DataFrame(prov); P.to_csv(f"{TAB}/14_geneset_provenance.csv", index=False)
print(P.to_string(index=False))

INF_SETS = {"current_lf_A_cytokine_chemokine": sorted(set(AXJ["axes"]["A_cytokine_chemokine"]) - RELA_T)}
TIS_SETS = {"current_lf_D_tissueforming": sorted(set(AXJ["axes"]["D_tissueforming"]) - RELA_T)}
for k in ("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "GOBP_INFLAMMATORY_RESPONSE"):
    if k in EXT: INF_SETS[k] = EXT[k]
for k in ("GOBP_OSSIFICATION", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"):
    if k in EXT: TIS_SETS[k] = EXT[k]
print("\n炎症セット:", {k: len(v) for k, v in INF_SETS.items()})
print("組織形成セット:", {k: len(v) for k, v in TIS_SETS.items()})

s2i = {}
for gid, s in zip(S.var_names, S.var["symbol"].astype(str)): s2i.setdefault(s, gid)
ACT = ulm_on_symbols(S, NET[["source","target","weight"]], tmin=5)
O = S.obs.copy(); O["RELA"] = ACT["RELA"].values
O["logdepth"] = np.log10(O.total_counts.astype(float) + 1); O["n_genes"] = O.n_genes_by_counts.astype(float)

def dctrl(x, y, dep, ng):
    Mx = np.column_stack([np.ones(len(x)), dep, ng])
    def res(v):
        b, *_ = np.linalg.lstsq(Mx, rankdata(v), rcond=None); return rankdata(v) - Mx @ b
    return spearmanr(res(x), res(y)).statistic

def matched(df, col, min_n=25):
    nb = 10 if len(df) < 500 else 20
    d = df.dropna(subset=[col, "RELA"]).copy()
    d["bin"] = pd.qcut(d.RELA, nb, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin", observed=True) if len(g) >= min_n]
    if not B: return np.nan
    n = np.array([b[0] for b in B]); s = np.array([b[1] for b in B])
    return round(float(np.sqrt(np.average(s**2, weights=n)) / tot), 4)

# ---- 3 scoring methods ----
import decoupler as dc
def build_scores(sets, prefix):
    out = {}
    for nm, genes in sets.items():
        ids = [s2i[g] for g in genes if g in s2i]
        col = f"{prefix}__{nm}__score_genes"
        sc.tl.score_genes(S, ids, ctrl_size=50, n_bins=25, random_state=SEED, score_name=col)
        out[(nm, "score_genes")] = S.obs[col].values
    net = pd.DataFrame([dict(source=nm, target=g, weight=1.0) for nm, gs in sets.items() for g in gs])
    import anndata as ad_
    v = pd.DataFrame(index=pd.Index(S.var["symbol"].astype(str)))
    B = ad_.AnnData(S.X, obs=S.obs[[]].copy(), var=v); B.var_names_make_unique()
    for meth, fn in [("ulm", dc.mt.ulm), ("aucell", getattr(dc.mt, "aucell", None))]:
        if fn is None: continue
        try:
            fn(data=B, net=net, tmin=5, verbose=False)
            key = [k for k in B.obsm if k.startswith("score_")][-1]
            R = B.obsm[key]
            for nm in sets:
                if nm in R.columns: out[(nm, meth)] = R[nm].values
            del B.obsm[key]
        except Exception as e:
            print(f"  {meth} 失敗: {type(e).__name__} {e}")
    return out
IS = build_scores(INF_SETS, "INF"); TS = build_scores(TIS_SETS, "TIS")

rows = []
for (inm, imeth), iv in IS.items():
    for (tnm, tmeth), tv in TS.items():
        if imeth != tmeth: continue
        d = O.copy(); d["_i"] = iv; d["_t"] = tv
        for les in sorted(d.lesion.astype(str).unique()) + ["POOLED"]:
            g = d if les == "POOLED" else d[d.lesion.astype(str) == les]
            rows.append(dict(inflammatory_set=inm, tissueforming_set=tnm, scoring=imeth, lesion=les,
                             n=len(g), raw_rho=round(float(spearmanr(g._i, g._t).statistic), 4),
                             depth_controlled_rho=round(float(dctrl(g._i.values, g._t.values,
                                 g.logdepth.values, g.n_genes.values)), 4),
                             defA_frac_retained=matched(g, "_t")))
RB = pd.DataFrame(rows); RB.to_csv(f"{TAB}/14_robustness_genesets_scoring.csv", index=False)
print("\n=== 3x3 x スコア法（POOLED）===")
print(RB[RB.lesion == "POOLED"][["inflammatory_set","tissueforming_set","scoring",
      "depth_controlled_rho","defA_frac_retained"]].to_string(index=False))
require_heartbeat("14 完了時")
print("\n14a/b DONE")
