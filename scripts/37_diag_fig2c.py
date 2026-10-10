#!/usr/bin/env python
"""Purpose: Diagnostic plot for Fig2c embeddings colored by 6 covariates; tabulates 2-lobe composition
Input: 33_ingest_umap_validation.csv, 33_umap_discovery.csv, discovery & validation h5ad
Output: diag_fig2c.pdf, diag_fig2c.png, 34_fig2c_lobes.csv

37: Diagnostic plots coloring the Fig2c embedding in 6 ways, and a breakdown table of the 2 lobes.

Fig2c is the embedding of 6,215 discovery stromal cells, with the 1,779 validation cells projected onto it via ingest.
To determine what explains the 2 lobes, the same coordinates are colored by (a) cohort, (b) lesion,
(c) log10 total UMI, (d) number of detected genes, (e) mitochondrial fraction, and (f) state label.
Lobes are defined by merging leiden (resolution 0.02) clusters on the 2D embedding into 2 groups via single-linkage of centroids.
Data is not recomputed: coordinates are only read from 33_*.csv, and covariates and expression from the existing h5ad.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".numba_cache"))
from _msstyle import *                      # rcParams (Arial 7-8pt) and color definitions
import scanpy as sc, anndata as ad, pandas as pd, numpy as np, scipy.sparse as sp
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABD = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
FIGD = f"{W}/results/figures_ms"
GENES = ["DCN", "LUM", "POSTN", "IBSP", "SP7", "RUNX2", "PDGFRA", "PDGFRB", "ACTA2", "MKI67"]
LES = ["A", "B", "C", "D", "E", "CAP1", "CAP2", "CAP3"]
STATES = ["inflammatory", "intermediate", "tissue-forming", "other"]

Ud = pd.read_csv(f"{TABD}/33_umap_discovery.csv")
IN = pd.read_csv(f"{TABD}/33_ingest_umap_validation.csv")
sd = Ud[Ud.compartment == "stromal"]


def grab(fp, cells):
    a = sc.read_h5ad(fp)
    sym = a.var.symbol.astype(str) if "symbol" in a.var else pd.Series(a.var_names.astype(str), index=a.var_names)
    ids = {g: sym[sym == g].index[0] for g in GENES if (sym == g).any()}
    X = a[:, list(ids.values())].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    E = pd.DataFrame(X, index=a.obs_names.astype(str), columns=list(ids))
    O = a.obs.copy(); O.index = a.obs_names.astype(str)
    keep = [c for c in ["n_genes_by_counts", "total_counts", "pct_counts_mt", "state"] if c in O]
    return O[keep].join(E).reindex(list(cells))


Dd = grab(f"{W}/data/processed/stromal_discovery_labelled.h5ad", sd.cell.astype(str))
Dv = grab(f"{W}/data/processed/stromal_validation.h5ad", IN.cell.astype(str))
F = pd.concat([
    pd.DataFrame(dict(cell=sd.cell.astype(str).values, umap1=sd.umap1.values, umap2=sd.umap2.values,
                      cohort="discovery_n5", lesion=sd.lesion.astype(str).values,
                      state=Dd["state"].astype(str).values)),
    pd.DataFrame(dict(cell=IN.cell.astype(str).values, umap1=IN.umap1.values, umap2=IN.umap2.values,
                      cohort="validation_n3", lesion=IN.lesion.astype(str).values,
                      state=IN.state_within.astype(str).values))], ignore_index=True)
C = pd.concat([Dd, Dv], ignore_index=True)
for c in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] + GENES:
    F[c] = C[c].values if c in C else np.nan
F["log10_total_counts"] = np.log10(F.total_counts.astype(float))
assert len(F) == len(sd) + len(IN) and F.total_counts.notna().all()

# ---- Definition of the cluster ----
A = ad.AnnData(np.zeros((len(F), 1), dtype=np.float32), obs=F.set_index("cell"))
A.obsm["X_emb"] = F[["umap1", "umap2"]].values.astype(np.float32)
sc.pp.neighbors(A, n_neighbors=15, use_rep="X_emb", random_state=0)
sc.tl.leiden(A, resolution=0.02, key_added="leiden", flavor="igraph", n_iterations=2,
             directed=False, random_state=0)
F["leiden"] = A.obs["leiden"].astype(str).values
cen = F.groupby("leiden")[["umap1", "umap2"]].mean()
grp = pd.Series(fcluster(linkage(cen.values, method="single"), t=2, criterion="maxclust"), index=cen.index)
F["lobe"] = ["L" + str(g) for g in grp.reindex(F.leiden).values]
gap = pdist(cen.values)
print(f"leiden 0.02 -> {F.leiden.nunique()} クラスタ、単連結で 2 塊に統合 "
      f"(centroid distance {gap.min():.2f}-{gap.max():.2f})")
print(pd.crosstab(F.lobe, F.cohort).to_string())

# ---- Figure ----
f = plt.figure(figsize=(180 * MM, 118 * MM))
gs = f.add_gridspec(2, 3, hspace=0.30, wspace=0.40)
COH = {"discovery_n5": DISC, "validation_n3": VAL}
xy = dict(s=1.6, lw=0, rasterized=True)


def frame(ax, title, letter):
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1", fontsize=6.4); ax.set_ylabel("UMAP 2", fontsize=6.4)
    ax.set_title(title, loc="left", fontsize=7.2)
    panel(ax, letter, dx=-0.10)


a1 = f.add_subplot(gs[0, 0])
for k, col in COH.items():
    g = F[F.cohort == k]
    a1.scatter(g.umap1, g.umap2, c=col, label=f"{k.split('_')[0]} ({len(g):,})", **xy)
a1.legend(fontsize=5.6, markerscale=4, loc="lower left", handletextpad=0.2)
frame(a1, "Cohort", "a")

a2 = f.add_subplot(gs[0, 1])
for i, l in enumerate(LES):
    g = F[F.lesion == l]
    if len(g):
        a2.scatter(g.umap1, g.umap2, c=OKABE[i % len(OKABE)], label=l, **xy)
a2.legend(fontsize=5.0, markerscale=4, ncol=2, loc="lower left", handletextpad=0.2, columnspacing=0.6)
frame(a2, "Lesion", "b")

for (ax_pos, col, lab, letter, cm) in [((0, 2), "log10_total_counts", "log$_{10}$ total UMI", "c", "viridis"),
                                       ((1, 0), "n_genes_by_counts", "genes detected", "d", "viridis"),
                                       ((1, 1), "pct_counts_mt", "mitochondrial reads (%)", "e", "magma")]:
    ax = f.add_subplot(gs[ax_pos])
    o = np.argsort(F[col].values)
    sm = ax.scatter(F.umap1.values[o], F.umap2.values[o], c=F[col].values[o], cmap=cm, **xy)
    cb = f.colorbar(sm, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=6.2)   # The label is carried by the panel title (it overlaps the adjacent axis label)
    frame(ax, lab, letter)

a6 = f.add_subplot(gs[1, 2])
for s in STATES:
    g = F[F.state == s]
    if len(g):
        a6.scatter(g.umap1, g.umap2, c=STATE[s], label=f"{s} ({len(g):,})", **xy)
a6.legend(fontsize=5.0, markerscale=4, loc="lower left", handletextpad=0.2)
frame(a6, "State label", "f")

f.savefig(f"{FIGD}/diag_fig2c.png", dpi=600, bbox_inches="tight")
f.savefig(f"{FIGD}/diag_fig2c.pdf", bbox_inches="tight")
print("saved diag_fig2c.png / .pdf")

# ---- Lump table ----
rows = []
for lb, g in F.groupby("lobe"):
    r = dict(lobe=lb, n_cells=len(g), pct_of_all=round(100 * len(g) / len(F), 2),
             n_discovery=int((g.cohort == "discovery_n5").sum()),
             n_validation=int((g.cohort == "validation_n3").sum()))
    r["pct_discovery_within_lobe"] = round(100 * r["n_discovery"] / len(g), 2)
    r["pct_validation_within_lobe"] = round(100 * r["n_validation"] / len(g), 2)
    r["pct_of_all_discovery_6215"] = round(100 * r["n_discovery"] / len(sd), 2)
    r["pct_of_all_validation_1779"] = round(100 * r["n_validation"] / len(IN), 2)
    vc = g.lesion.value_counts()
    for l in LES:
        r[f"n_{l}"] = int(vc.get(l, 0))
    r["lesion_breakdown"] = "; ".join(f"{l} {int(vc.get(l, 0))}" for l in LES if vc.get(l, 0))
    sv = g.state.value_counts()
    for s in STATES:
        r[f"n_state_{s.replace('-', '_')}"] = int(sv.get(s, 0))
    r["state_breakdown"] = "; ".join(f"{s} {int(sv.get(s, 0))} ({100 * sv.get(s, 0) / len(g):.1f}%)"
                                     for s in STATES if sv.get(s, 0))
    r["median_n_genes"] = int(np.median(g.n_genes_by_counts.astype(float)))
    r["median_total_counts"] = int(np.median(g.total_counts.astype(float)))
    r["median_pct_mt"] = round(float(np.median(g.pct_counts_mt.astype(float))), 2)
    for gn in GENES:
        v = g[gn].astype(float); ok = v.notna()
        r[f"n_quantified_{gn}"] = int(ok.sum())
        r[f"mean_{gn}"] = round(float(v[ok].mean()), 3) if ok.sum() else np.nan
        r[f"pct_expr_{gn}"] = round(100 * float((v[ok] > 0).mean()), 1) if ok.sum() else np.nan
    rows.append(r)
L = pd.DataFrame(rows)
L["note"] = [("SP7 is absent from the validation gene space; in this lobe it is quantified only in the "
              f"{int(L.loc[i, 'n_quantified_SP7'])} discovery-derived cells")
             if L.loc[i, "n_quantified_SP7"] < L.loc[i, "n_cells"] else "all 10 genes quantified in all cells"
             for i in L.index]
L["lobe_definition"] = (f"leiden on the 2-D embedding (n_neighbors=15, resolution=0.02, seed 0): "
                        f"{F.leiden.nunique()} clusters merged into 2 groups by single linkage on centroids")
L.to_csv(f"{TABD}/34_fig2c_lobes.csv", index=False)
assert int(L.n_cells.sum()) == len(F) and int(L.n_validation.sum()) == len(IN)

# Control for whether depth can be an explanatory factor
thr = float(np.median(F.loc[F.cohort == "validation_n3", "n_genes_by_counts"].astype(float)))
hi = F[(F.cohort == "discovery_n5") & (F.n_genes_by_counts.astype(float) >= thr)]
print(f"\n検証の n_genes 中位数 {thr:.0f} 以上の探索細胞 {len(hi)} 個のうち L2 に落ちたのは "
      f"{int((hi.lobe == 'L2').sum())} 個 ({100 * (hi.lobe == 'L2').mean():.2f}%)")
print(L[["lobe", "n_cells", "n_discovery", "n_validation", "median_n_genes", "median_total_counts"]].to_string(index=False))
print("37 DONE")
