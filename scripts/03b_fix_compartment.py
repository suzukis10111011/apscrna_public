#!/usr/bin/env python
"""Purpose: Fix compartment annotations by using max marker z-scores instead of mean, add plasma cell markers
  to immune set, threshold unassigned clusters, and assign fine immune subtypes
Input: {OUT}/all_cells_discovery.h5ad
Output: {OUT}/all_cells_discovery.h5ad, {TAB}/03b_celltype_fine_by_lesion_discovery.csv,
  {TAB}/03b_compartment_by_lesion_{tag}.csv, {TAB}/03b_compartment_marker_maxz_{tag}.csv

03b: Correction of compartment annotations (provides functions usable for both discovery and validation, applied to discovery).

The initial implementation in 03 was incorrect. Evidence and corrections:
  (a) 5 clusters, 15,973 cells (35.9%), were assigned as erythrocyte, but **mean HBB expression
      is 0.00-0.01 across all clusters**, so erythrocytes do not actually exist. These show
      MZB1 2.47-2.98 / DERL3 1.56-1.80 / TNFRSF17 0.48-0.59 / XBP1 2.34-2.81 / IGHG1 6.0-6.4,
      characteristic of **plasma cells**, with PTPRC at 0.06 (plasma cells have low CD45 transcription).
      High UMI (7.3-8.7k), low gene count (1.0-1.2k), and extremely low mito (0.9-1.1%) are also consistent with plasma cells.
      It is known that plasma cells predominate in granulation tissue of chronic apical periodontitis.
  (b) Cluster 3 is LYZ 3.02 / CD68 1.08, i.e., myeloid lineage, but was assigned as schwann_neural.
      The cause was taking the **mean** of the 10 genes in the immune set. Since markers for
      T/B/plasma/myeloid/mast are mutually exclusive, averaging them buries every subtype.
  (c) Clusters where no marker lit up at all were also forcibly assigned via argmax.

Corrections:
  1. Changed each compartment's score to the **maximum** z within the set (instead of the mean).
  2. Added plasma cell markers MZB1, DERL3, TNFRSF17 to the immune set (supplementing the lineage
     for which JCHAIN is already present in the instructed immune subtypes. IGKC/IGHG1 are high across
     all clusters due to surrounding RNA contamination and are not used for discrimination).
  3. Clusters with max z below 1.0 are labeled "unassigned".
  4. Immune clusters are also assigned celltype_fine (plasma/myeloid/T/B/mast/other_immune).
Input : data/processed/all_cells_discovery.h5ad
Output : overwrite update of the same file (compartment, celltype_fine replaced)
       results/tables/03b_compartment_by_lesion_discovery.csv
       results/tables/03b_compartment_marker_maxz_discovery.csv
"""
import os, sys, json
import numpy as np, pandas as pd, scanpy as sc, warnings
warnings.filterwarnings("ignore")

import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")

MARK = {
 "stromal":        ["COL1A1","COL1A2","DCN","LUM","PDGFRB"],
 "endothelial":    ["PECAM1","VWF","CDH5"],
 "pericyte_smc":   ["RGS5","ACTA2","MYH11"],
 "immune":         ["PTPRC","CD3E","MS4A1","CD79A","JCHAIN","LYZ","CD68","CD14","S100A8","KIT",
                    "MZB1","DERL3","TNFRSF17"],          # The last 3 are additions (plasma cells)
 "epithelial":     ["KRT14","KRT5","EPCAM"],
 "schwann_neural": ["S100B","PLP1"],
 "erythrocyte":    ["HBB"],
}
FINE = {"plasma": ["MZB1","DERL3","TNFRSF17","JCHAIN"], "myeloid": ["LYZ","CD68","CD14","S100A8"],
        "T": ["CD3E"], "B": ["MS4A1","CD79A"], "mast": ["KIT"]}
MINZ = 1.0


def annotate(A, clkey="leiden_0.5", tag="discovery"):
    sym2id = {}
    for gid, s in zip(A.var_names, A.var.symbol.astype(str)): sym2id.setdefault(s, gid)
    genes = sorted({g for gs in MARK.values() for g in gs if g in sym2id})
    Z = sc.get.obs_df(A, keys=[sym2id[g] for g in genes] + [clkey]).rename(columns={sym2id[g]: g for g in genes})
    means = Z.groupby(clkey, observed=True).mean()
    zs = (means - means.mean()) / means.std(ddof=0).replace(0, np.nan)
    # Compartment score = maximum z within the set
    score = pd.DataFrame({k: zs[[g for g in gs if g in zs.columns]].max(1) for k, gs in MARK.items()})
    best = score.idxmax(1); bestv = score.max(1)
    assign = best.where(bestv >= MINZ, "unassigned")
    fine = {}
    for cl in score.index:
        if assign[cl] != "immune": fine[cl] = assign[cl]; continue
        sub = {k: zs.loc[cl, [g for g in gs if g in zs.columns]].max() for k, gs in FINE.items()}
        fine[cl] = f"immune_{max(sub, key=sub.get)}"
    A.obs["compartment"] = A.obs[clkey].map(assign).astype(str)
    A.obs["celltype_fine"] = A.obs[clkey].map(fine).astype(str)
    out = score.copy(); out["best"] = best; out["best_z"] = bestv
    out["assigned"] = assign; out["celltype_fine"] = pd.Series(fine)
    out["n_cells"] = A.obs[clkey].value_counts().reindex(out.index)
    out.to_csv(f"{TAB}/03b_compartment_marker_maxz_{tag}.csv")
    ct = pd.crosstab(A.obs.lesion, A.obs.compartment); frac = (ct.T / ct.sum(1)).T
    o = ct.copy(); o.columns = [f"n_{c}" for c in ct.columns]
    for c in frac.columns: o[f"pct_{c}"] = (100*frac[c]).round(2)
    o["n_total"] = ct.sum(1)
    o.to_csv(f"{TAB}/03b_compartment_by_lesion_{tag}.csv")
    return A, out, (100*frac).round(2)


if __name__ == "__main__":
    D = sc.read_h5ad(f"{OUT}/all_cells_discovery.h5ad")
    print("読み込み:", D.shape)
    D, tab, pct = annotate(D, "leiden_0.5", "discovery")
    print("\nクラスタ -> 区画（最大z方式）:")
    print(tab[["best","best_z","assigned","celltype_fine","n_cells"]].round(2).to_string())
    print("\n区画×病変 (%):"); print(pct.to_string())
    print("\n細分×病変 (%):")
    ctf = pd.crosstab(D.obs.lesion, D.obs.celltype_fine)
    print((100*(ctf.T/ctf.sum(1)).T).round(1).to_string())
    ctf.to_csv(f"{TAB}/03b_celltype_fine_by_lesion_discovery.csv")
    stro = pct["stromal"] if "stromal" in pct.columns else None
    pooled = 100*float((D.obs.compartment == "stromal").mean())
    print(f"\n探索 pooled の間質比率: {pooled:.2f}%  （申請書 14.8%、差 {pooled-14.8:+.2f} ポイント）")
    json.dump(dict(stromal_pct_pooled=round(pooled,2), grant=14.8, diff=round(pooled-14.8,2),
                   n_cells=int(D.n_obs)), open(f"{TAB}/03b_stromal_pct_discovery.json","w"), indent=1)
    D.write(f"{OUT}/all_cells_discovery.h5ad")
    print("\n03b DONE（all_cells_discovery.h5ad を更新）")
