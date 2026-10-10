#!/usr/bin/env python
"""Purpose: Test inflammation–tissue formation axis coupling across tissue compartments, controlling for depth;
  assess myeloid-specific M1/M2 axes.
Input: all_cells.h5ad, collectri_regulons.csv (frozen regulatory network)
Output: 12_coupling_by_compartment.csv, 12_geneset_coverage_{COH}.csv

12 (Figure 6): Examine whether the inflammation axis and the tissue formation axis are coupled in compartments other than stroma.
For each compartment × lesion, depth-controlled ρ and retained dispersion of RELA matching (bins are 10 when n<500).
Compartments with <100 cells in a lesion are excluded. For myeloid, also construct dedicated M1/M2-repair axes and compute the same ρ.
Usage: python 12_coupling_by_compartment.py [all_cells.h5ad] [cohort_tag]
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("12")
import scanpy as sc, pandas as pd, numpy as np
from scipy.stats import spearmanr, rankdata
from _utils import ulm_on_symbols

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables"); SEED = 0
PATH = sys.argv[1] if len(sys.argv) > 1 else f"{W}/data/processed/all_cells_discovery.h5ad"
COH  = sys.argv[2] if len(sys.argv) > 2 else "discovery_n5"
AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
RELA_T = set(NET.loc[NET.source == "RELA", "target"])
AX = {k: set(v) for k, v in AXJ["axes"].items()}
LF = {k: sorted(v - RELA_T) for k, v in AX.items()}
M1 = ["IL1B","TNF","IL6","CXCL8","CCL3","CCL4","CXCL10","PTGS2","NFKBIA"]
M2 = ["MRC1","CD163","MSR1","IL10","TGFB1","CCL18","SPP1","VEGFA","MMP9","FOLR2"]

A = sc.read_h5ad(PATH)
s2i = {}
for gid, s in zip(A.var_names, A.var["symbol"].astype(str)): s2i.setdefault(s, gid)
def score(name, genes):
    ids = [s2i[g] for g in genes if g in s2i]
    sc.tl.score_genes(A, ids, ctrl_size=50, n_bins=25, random_state=SEED, score_name=name)
    return len(ids)
cov = [dict(set=f"lf_{k}", n_def=len(v), n_found=score(f"lf_{k}", v)) for k, v in LF.items()]
cov.append(dict(set="myeloid_M1", n_def=len(M1), n_found=score("myeloid_M1", M1)))
cov.append(dict(set="myeloid_M2repair", n_def=len(M2), n_found=score("myeloid_M2repair", M2)))
pd.DataFrame(cov).to_csv(f"{TAB}/12_geneset_coverage_{COH}.csv", index=False)
print(pd.DataFrame(cov).to_string(index=False))

ACT = ulm_on_symbols(A, NET[["source","target","weight"]], tmin=5)
O = A.obs.copy()
O["RELA"] = ACT["RELA"].values if "RELA" in ACT.columns else np.nan
O["logdepth"] = np.log10(O.total_counts.astype(float) + 1)
O["n_genes"] = O.n_genes_by_counts.astype(float)
INF, TIS = "lf_A_cytokine_chemokine", "lf_D_tissueforming"

def dctrl(x, y, dep, ng):
    M = np.column_stack([np.ones(len(x)), dep, ng])
    def res(v):
        b, *_ = np.linalg.lstsq(M, rankdata(v), rcond=None); return rankdata(v) - M @ b
    return spearmanr(res(x), res(y)).statistic

def matched_residual(df, col, rela, min_n=25):
    nbins = 10 if len(df) < 500 else 20
    d = df.dropna(subset=[col, rela]).copy()
    if len(d) < 60: return dict(n_bins=0, total_sd=np.nan, matched_sd=np.nan, frac_retained=np.nan, bins_used=nbins)
    d["bin"] = pd.qcut(d[rela], nbins, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin", observed=True) if len(g) >= min_n]
    if not B: return dict(n_bins=0, total_sd=tot, matched_sd=np.nan, frac_retained=np.nan, bins_used=nbins)
    n = np.array([b[0] for b in B]); s = np.array([b[1] for b in B])
    within = float(np.sqrt(np.average(s**2, weights=n)))
    return dict(n_bins=len(B), total_sd=round(tot,4), matched_sd=round(within,4),
                frac_retained=round(within/tot,4), bins_used=nbins)

COMPS = ["stromal","immune_myeloid","immune_plasma","immune_T","endothelial","pericyte_smc"]
key = "celltype_fine" if "celltype_fine" in O.columns else "compartment"
rows = []
for comp in COMPS:
    sub = O[O[key].astype(str) == comp]
    if len(sub) == 0: print(f"{comp}: no cells — excluded"); continue
    for les in sorted(sub.lesion.astype(str).unique()) + ["POOLED"]:
        g = sub if les == "POOLED" else sub[sub.lesion.astype(str) == les]
        if les != "POOLED" and len(g) < 100:
            rows.append(dict(cohort=COH, compartment=comp, lesion=les, n_cells=len(g),
                             excluded="fewer than 100 cells")); continue
        r = dict(cohort=COH, compartment=comp, lesion=les, n_cells=len(g), excluded="",
                 raw_rho=round(float(spearmanr(g[INF], g[TIS]).statistic), 4),
                 depth_controlled_rho=round(float(dctrl(g[INF].values, g[TIS].values,
                                                        g.logdepth.values, g.n_genes.values)), 4))
        r.update({f"defA_{k}": v for k, v in matched_residual(g, TIS, "RELA").items()})
        if comp == "immune_myeloid":
            r["myeloid_M1_vs_M2_raw_rho"] = round(float(spearmanr(g.myeloid_M1, g.myeloid_M2repair).statistic), 4)
            r["myeloid_M1_vs_M2_depth_controlled_rho"] = round(float(dctrl(
                g.myeloid_M1.values, g.myeloid_M2repair.values, g.logdepth.values, g.n_genes.values)), 4)
        rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(f"{TAB}/12_coupling_by_compartment.csv", index=False)
print("\n=== 区画別（POOLED 行）===")
c = ["compartment","n_cells","raw_rho","depth_controlled_rho","defA_frac_retained","defA_bins_used"]
print(R[(R.lesion == "POOLED")][c].to_string(index=False))
print("\n=== 病変別 depth 統制ρ ===")
print(R[R.excluded == ""].pivot_table(index="compartment", columns="lesion",
      values="depth_controlled_rho").round(3).to_string())
if "myeloid_M1_vs_M2_depth_controlled_rho" in R.columns:
    print("\n=== myeloid 固有軸 M1 × M2修復 ===")
    print(R[R.compartment == "immune_myeloid"][["lesion","n_cells","myeloid_M1_vs_M2_raw_rho",
          "myeloid_M1_vs_M2_depth_controlled_rho"]].to_string(index=False))
require_heartbeat("12 完了時")
print("\n12 DONE")
