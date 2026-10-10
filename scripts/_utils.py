"""Purpose: Common utility functions for single-cell analysis; direct harmonypy 2.0.0 integration to work
  around scanpy wrapper incompatibility
Input: none (module imported by analysis scripts)
Output: {tabdir}/03c_celltype_fine_by_lesion_{tag}.csv, {tabdir}/03c_compartment_by_lesion_{tag}.csv,
  {tabdir}/03c_compartment_votes_{tag}.csv

Common helpers.
harmony: scanpy.external.pp.harmony_integrate is incompatible with harmonypy 2.0.0
(the old API returns PCs x cells and the wrapper applies .T, but 2.0.0 returns cells x PCs, so
 the shape of obsm breaks: ValueError "Value had shape (50,)"). Call harmonypy directly.
"""
import os
import numpy as np


def harmony_run(adata, key, seed=0, basis="X_pca", out="X_pca_harmony", verbose=False):
    import harmonypy as hm
    X = np.asarray(adata.obsm[basis])
    ho = hm.run_harmony(X, adata.obs, [key], random_state=seed, verbose=verbose)
    Z = ho.result() if callable(getattr(ho, "result", None)) else getattr(ho, "result", None)
    if Z is None:
        Z = ho.Z_corr
    Z = np.asarray(Z)
    if Z.shape != X.shape:                      # Correct the orientation just in case
        if Z.T.shape == X.shape:
            Z = Z.T
        else:
            raise RuntimeError(f"harmony 出力形状が不正: {Z.shape}, 期待 {X.shape}")
    adata.obsm[out] = Z
    return adata


def ulm_on_symbols(A, net, tmin=5):
    """CollecTRI の target は gene symbol、こちらの var_names は Ensembl ID なので
    symbol を index にした軽量 AnnData を作って ULM を回す（04 初回失敗の原因）。"""
    import anndata as _ad, pandas as _pd, numpy as _np, decoupler as _dc
    v = _pd.DataFrame(index=_pd.Index(A.var["symbol"].astype(str), name=None))
    B = _ad.AnnData(A.X, obs=A.obs[[]].copy(), var=v)
    B.var_names_make_unique()
    _dc.mt.ulm(data=B, net=net, tmin=tmin, verbose=False)
    act = B.obsm["score_ulm"].copy(); act.index = A.obs_names
    return act


# ===== 03c: Compartment annotation (the single implementation shared by discovery and validation) =====
# Changes from 03b (z of cluster means, max within set):
#   Compute sc.tl.score_genes on a per-cell basis, take the per-cell argmax, and then
#   Determine compartments by majority vote per cluster. Since 03b is cross-cluster z,
#   when 1 cell type was scattered across many clusters they canceled each other out and unassigned increased
#   (33.5% in discovery lesion B). Cell-level scores are not affected by this.
MARK_COMP = {
 "stromal":        ["COL1A1","COL1A2","DCN","LUM","PDGFRB"],
 "endothelial":    ["PECAM1","VWF","CDH5"],
 "pericyte_smc":   ["RGS5","ACTA2","MYH11"],
 "immune":         ["PTPRC","CD3E","MS4A1","CD79A","JCHAIN","LYZ","CD68","CD14","S100A8","KIT",
                    "MZB1","DERL3","TNFRSF17"],
 "epithelial":     ["KRT14","KRT5","EPCAM"],
 "schwann_neural": ["S100B","PLP1"],
 "erythrocyte":    ["HBB"],
}
MINZ_COMP = 1.0
FINE_IMMUNE = {"plasma": ["MZB1","DERL3","TNFRSF17","JCHAIN"], "myeloid": ["LYZ","CD68","CD14","S100A8"],
               "T": ["CD3E"], "B": ["MS4A1","CD79A"], "mast": ["KIT"]}


def _symbol_view(A):
    """var_names が Ensembl ID の AnnData から symbol 索引の軽量 AnnData を作る。"""
    import anndata as ad, pandas as pd
    v = pd.DataFrame(index=pd.Index(A.var["symbol"].astype(str)))
    B = ad.AnnData(A.X, obs=A.obs[[]].copy(), var=v)
    B.var_names_make_unique()
    return B


def annotate_compartment(A, clkey="leiden_0.5", tag="discovery", tabdir=None, seed=0):
    """細胞単位 score_genes -> 細胞ごと argmax -> クラスタ多数決で区画を決める。
    戻り値: (A, per_cluster_table, pct_by_lesion)"""
    import scanpy as sc, pandas as pd, numpy as np
    B = _symbol_view(A)
    present = {k: [g for g in gs if g in B.var_names] for k, gs in MARK_COMP.items()}
    for k, gs in present.items():
        if gs:
            sc.tl.score_genes(B, gs, score_name=f"cs_{k}", ctrl_size=50, random_state=seed)
        else:
            B.obs[f"cs_{k}"] = 0.0
    S = B.obs[[f"cs_{k}" for k in MARK_COMP]].copy()
    S.columns = list(MARK_COMP)
    S.index = A.obs_names
    cell_best = S.idxmax(1).where(S.max(1) > 0, "unassigned")      # score<=0 is unclassified
    A.obs["compartment_cell"] = pd.Categorical(cell_best.values)
    for k in MARK_COMP: A.obs[f"cs_{k}"] = S[k].values

    votes = pd.crosstab(A.obs[clkey], A.obs["compartment_cell"])
    majority = votes.idxmax(1)
    frac_win = (votes.max(1) / votes.sum(1)).rename("majority_frac")

    # Cross-cluster z (03b method). This is used as the primary judgment since it reproduces the discovery stromal 6,215 and matches exactly
    # the definition that fully agrees with the original fibroblast object human_fibro.h5ad (6,215 cells).
    # If only cell-level majority vote is used, pericyte_smc 3,517 is absorbed into stromal,
    # and the stromal ratio shifts from 13.99% -> 22.20%, breaking comparability (confirmed in 03c first version).
    sym2id = {}
    for gid, s in zip(A.var_names, A.var["symbol"].astype(str)): sym2id.setdefault(s, gid)
    genes = sorted({g for gs in MARK_COMP.values() for g in gs if g in sym2id})
    Zdf = sc.get.obs_df(A, keys=[sym2id[g] for g in genes] + [clkey]).rename(
        columns={sym2id[g]: g for g in genes})
    means = Zdf.groupby(clkey, observed=True).mean()
    zs = (means - means.mean()) / means.std(ddof=0).replace(0, np.nan)
    zscore = pd.DataFrame({k: zs[[g for g in gs if g in zs.columns]].max(1) for k, gs in MARK_COMP.items()})
    zbest, zbestv = zscore.idxmax(1), zscore.max(1)
    # Compromise: for clusters where z reaches the threshold, use the z-based call; only fill in clusters that don't reach it by per-cell majority vote.
    # This resolves the unclassified cells from 03b (6,962 cells explored = 2 plasma cell clusters).
    assign = zbest.where(zbestv >= MINZ_COMP, majority.reindex(zbest.index))
    assign = assign.fillna("unassigned")

    Bi = _symbol_view(A)
    fine = {}
    fin_present = {k: [g for g in gs if g in Bi.var_names] for k, gs in FINE_IMMUNE.items()}
    for k, gs in fin_present.items():
        if gs: sc.tl.score_genes(Bi, gs, score_name=f"fs_{k}", ctrl_size=50, random_state=seed)
        else:  Bi.obs[f"fs_{k}"] = 0.0
    F = Bi.obs[[f"fs_{k}" for k in FINE_IMMUNE]]; F.columns = list(FINE_IMMUNE); F.index = A.obs_names
    fine_cell = F.idxmax(1)
    for cl in votes.index:
        if assign[cl] != "immune": fine[cl] = assign[cl]; continue
        m = (A.obs[clkey] == cl).values
        fine[cl] = "immune_" + fine_cell[m].value_counts().idxmax()

    A.obs["compartment"] = A.obs[clkey].map(assign).astype(str)
    A.obs["celltype_fine"] = A.obs[clkey].map(fine).astype(str)

    out = votes.copy()
    out["z_best"] = zbest; out["z_best_value"] = zbestv.round(3)
    out["cell_majority"] = majority; out["majority_frac"] = frac_win.round(3)
    out["source"] = np.where(zbestv >= MINZ_COMP, "cluster_z(03b)", "cell_majority(03c)")
    out["assigned"] = assign
    out["celltype_fine"] = pd.Series(fine)
    out["n_cells"] = A.obs[clkey].value_counts().reindex(out.index)
    ct = pd.crosstab(A.obs.lesion, A.obs.compartment); frac = (ct.T / ct.sum(1)).T
    o = ct.copy(); o.columns = [f"n_{c}" for c in ct.columns]
    for c in frac.columns: o[f"pct_{c}"] = (100 * frac[c]).round(2)
    o["n_total"] = ct.sum(1)
    if tabdir:
        out.to_csv(f"{tabdir}/03c_compartment_votes_{tag}.csv")
        o.to_csv(f"{tabdir}/03c_compartment_by_lesion_{tag}.csv")
        pd.crosstab(A.obs.lesion, A.obs.celltype_fine).to_csv(f"{tabdir}/03c_celltype_fine_by_lesion_{tag}.csv")
    return A, out, (100 * frac).round(2)


# ===== Single source of truth for the run-to-donor mapping =====
# Only ref/run_donor_map_v2.csv is authoritative. requant_run_donor_map.csv (derived from frozen)
# had CAP1 and CAP3 swapped, and trusting it had caused the data/counts labels to be swapped as well.
# The basis is the concordance rate with the previous stromal cell barcodes (94.6 / 96.1 / 96.9%, off-diagonal 1-7 cells).
RUN_DONOR_MAP = os.path.join(
    os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "ref", "run_donor_map_v2.csv")


def load_run_donor_map(path=RUN_DONOR_MAP):
    """→ (DataFrame, run->donor, donor->run)"""
    import pandas as pd
    m = pd.read_csv(path)
    return m, dict(zip(m.run_accession, m.donor)), dict(zip(m.donor, m.run_accession))


def assert_counts_mapping(cnt_dir, path=RUN_DONOR_MAP, strict=True):
    """data/counts/<donor>/run_info.json の call に、対応表が指すランが入っているか検証する。
    ディレクトリ名とランの取り違えを二度と通さないためのガード。→ {donor: accession}"""
    import json, os, re
    _, _, d2r = load_run_donor_map(path)
    found = {}
    for donor, acc in d2r.items():
        ri = os.path.join(cnt_dir, donor, "run_info.json")
        if not os.path.exists(ri):
            continue
        call = json.load(open(ri)).get("call", "")
        m = re.search(r"SRR\d+", call)
        got = m.group(0) if m else None
        found[donor] = got
        if got != acc:
            msg = f"{cnt_dir}/{donor} のランは {got}、対応表は {acc}（ref/run_donor_map_v2.csv）"
            if strict: raise RuntimeError("ラン→ドナー対応の不一致: " + msg)
            print("WARNING:", msg)
    return found


# ===== 20: Loading kb output (unifying symbol resolution) =====
# When kb's cells_x_genes.genes.names.txt has multiple gene ids sharing the same symbol,
# Drop all of that naming and retain gene id as is (because Ensembl 112 cDNA includes alt locus,
# 7,057 genes are affected, and PTPRC / MYH11 / DERL3 were also anonymized by this).
# Therefore, do not use names.txt; look it up again with ref/gene_id_to_symbol_ens112.tsv,
# and furthermore sum gene ids with the same symbol (MYH11 has 32.6% on a secondary copy, HLA-DQB1 has 47.8%).
SYMMAP = "ref/gene_id_to_symbol_ens112.tsv"


def load_kb_counts(d, symmap=None, aggregate=True, root=None):
    """Return the kb counts_unfiltered directory as AnnData.

    d        : .../counts_unfiltered
    aggregate: If True, sum gene ids sharing the same symbol and set var_names to symbol.
               Unresolved gene ids are left as gene id.
    var of the return value: symbol / gene_ids (comma-separated source ids) / n_ids / symbol_source
    """
    import os
    import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
    from scipy.io import mmread
    root = root or os.environ.get("APSC_W", ".")
    symmap = symmap or os.path.join(root, SYMMAP)
    X = sp.csr_matrix(mmread(os.path.join(d, "cells_x_genes.mtx")))
    bc = [l.strip() for l in open(os.path.join(d, "cells_x_genes.barcodes.txt"))]
    gid = [l.strip() for l in open(os.path.join(d, "cells_x_genes.genes.txt"))]
    A = ad.AnnData(X)
    A.obs_names = bc
    A.var_names = gid
    A.var["ensembl_id"] = gid
    M = pd.read_csv(symmap, sep="\t").set_index("gene_id")
    A.var["symbol"] = [M.symbol.get(g, g) for g in gid]
    A.var["symbol_source"] = [M.source.get(g, "absent_from_map") for g in gid]
    if not aggregate:
        return A
    key = pd.Index(A.var["symbol"].values)
    order = pd.unique(key)
    idx = pd.Series(np.arange(len(order)), index=order)
    cols = idx[key].values
    Agg = sp.csr_matrix((np.ones(len(cols)), (np.arange(len(cols)), cols)),
                        shape=(len(cols), len(order)))
    B = ad.AnnData((A.X @ Agg).tocsr())
    B.obs_names = A.obs_names
    B.var_names = [str(o) for o in order]
    grp = pd.DataFrame(dict(sym=key, gid=gid, src=A.var["symbol_source"].values)).groupby("sym", sort=False)
    B.var["symbol"] = B.var_names
    B.var["gene_ids"] = grp.gid.apply(lambda s: ",".join(s)).reindex(order).values
    B.var["n_ids"] = grp.size().reindex(order).values
    B.var["symbol_source"] = grp.src.first().reindex(order).values
    return B


# ===== 09: Normalization of barcode format =====
# GEO(CellRanger) is "AAAC...-1", kb is the raw 16mer, and after concatenation it is "<lesion>_<bc>".
# Always pass this through before set comparison.
def bare_barcode(x):
    s = str(x)
    if "_" in s:
        s = s.split("_", 1)[1]
    return s.split("-")[0]


def load_grant_factors(root=None):
    """Grant-proposal factor list, read from docs/private/grant_factors.txt.

    The list is grant material and is not published, so the file is absent from a public
    checkout; callers skip the grant-comparison block when this returns None.
    """
    import os
    root = root or os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = os.path.join(root, "docs", "private", "grant_factors.txt")
    if not os.path.exists(p):
        return None
    out = {}
    for line in open(p):
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        k, v = line.split(":", 1)
        out[k.strip()] = v.split()
    return out
