#!/usr/bin/env python
"""Purpose: Extract the two per-cell figure source tables from the processed stromal objects.
Input: $APSC_OUTDIR/stromal_discovery_labelled.h5ad, $APSC_OUTDIR/stromal_validation.h5ad
  (obs columns and obsm["X_umap"] only; nothing is recomputed).
Output: $APSC_TABDIR/33_stromal_program_scores.csv (Fig4 D),
  $APSC_TABDIR/33_stromal_umap_by_cohort.csv (Fig2 C).

This is a pure extraction: the programme scores and the cohort-internal UMAP coordinates were
computed by 04/04b (discovery) and 05/25 (validation) and are read here as they stand. Both tables
were originally written by interactive cells; this script is the scripted form of exactly that code,
so re-running it must reproduce them byte-for-byte.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("33_extract_figure_sources")
import scanpy as sc, pandas as pd, numpy as np

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
os.makedirs(TAB, exist_ok=True)

# ---------- 33_stromal_program_scores.csv (Fig4 D) ----------
rows = []
for f, coh in [("stromal_discovery_labelled", "discovery_n5"), ("stromal_validation", "validation_n3")]:
    a = sc.read_h5ad(f"{OUT}/{f}.h5ad", backed="r")
    o = a.obs
    stcol = "state" if "state" in o else "state_within"
    rows.append(pd.DataFrame(dict(cohort=coh, lesion=o.lesion.astype(str).values,
                                  state=o[stcol].astype(str).values,
                                  inflammatory=o["lf_A_cytokine_chemokine"].values,
                                  tissueforming=o["lf_D_tissueforming"].values)))
    a.file.close()
S = pd.concat(rows, ignore_index=True)
S.to_csv(f"{TAB}/33_stromal_program_scores.csv", index=False)
print(S.groupby(["cohort", "lesion"]).size().to_string())
print("\n計", len(S), "細胞 | 列:", list(S.columns))

# ---------- 33_stromal_umap_by_cohort.csv (Fig2 C) ----------
a = sc.read_h5ad(f"{OUT}/stromal_discovery_labelled.h5ad", backed="r")
o = a.obs.copy(); o.index = a.obs_names.astype(str)
U = pd.DataFrame(np.asarray(a.obsm["X_umap"]), index=a.obs_names.astype(str), columns=["umap1", "umap2"]); a.file.close()
b = sc.read_h5ad(f"{OUT}/stromal_validation.h5ad", backed="r")
ob = b.obs.copy(); ob.index = b.obs_names.astype(str)
Ub = pd.DataFrame(np.asarray(b.obsm["X_umap"]), index=b.obs_names.astype(str), columns=["umap1", "umap2"]); b.file.close()
cols = ["lesion", "p_D_tissueforming", "lf_D_tissueforming"]
SU = pd.concat([
    U.join(o[cols + ["state"]]).assign(cohort="discovery_n5").reset_index().rename(columns={"index": "cell"}),
    Ub.join(ob[cols + ["state_within"]].rename(columns={"state_within": "state"})).assign(
        cohort="validation_n3").reset_index().rename(columns={"index": "cell"})],
    ignore_index=True)[["cohort", "cell", "lesion", "state", "umap1", "umap2",
                        "p_D_tissueforming", "lf_D_tissueforming"]]
assert len(SU) == len(S), (SU.shape, S.shape)
assert SU.notna().all().all()
SU.to_csv(f"{TAB}/33_stromal_umap_by_cohort.csv", index=False)
print("\n素材:", SU.shape, "| 列:", list(SU.columns))
print(SU.groupby("cohort").p_D_tissueforming.agg(["min", "max", "median"]).round(4).to_string())
print("座標範囲:"); print(SU.groupby("cohort")[["umap1", "umap2"]].agg(["min", "max"]).round(2).to_string())
require_heartbeat("33_extract_figure_sources 完了時")
print("\n33_extract_figure_sources DONE")
