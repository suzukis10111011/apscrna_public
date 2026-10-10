#!/usr/bin/env python
"""Purpose: Answer the manuscript-versus-tables questions raised by an independent reader
  of clean_v16, and record the quantities that are not already in the tables.
Input: results/tables/S11c_dispersion_baselines.csv, S11c_null_factor_selection.csv,
  04_rho_inflam_vs_tissueforming_discovery.csv, 14_robustness_genesets_scoring.csv,
  10_receptors_pseudobulk_log2cpm.csv, 17_receptors_mixedmodel.csv,
  data/processed/stromal_{discovery,validation}.h5ad and the matching TF-activity
  parquets, data/processed/all_cells_discovery_obs_scored.parquet,
  data/processed/all_cells_validation.h5ad
Output: E30_retained_dispersion_by_lesion.csv, E31_retained_dispersion_residualised.csv,
  E32_compartment_composite_rho.csv, E33_top_decile_composition.csv,
  E34_receptor_pseudobulk_vs_model.csv, E35_single_axis_pooled_values.csv

Nothing frozen is touched: the gene sets, the candidate factors, the thresholds and the
published tables are read only. The two quantities that are new - the retained dispersion
computed on a depth-residualised score, and the composite depth-controlled rho per
compartment - are written to this E series and nowhere else.
"""
import os
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/nbc")
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr, rankdata

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
T = f"{W}/results/tables"
P = f"{W}/data/processed"
OUT = f"{W}/08_external_review/tables"
ORDER = ["A", "B", "C", "D", "E", "CAP1", "CAP2", "CAP3"]
NOTE = ("post hoc, in response to a reading of clean_v16; frozen gene sets, candidate "
        "factors and thresholds are unchanged")


def matched_residual(df, col, by, nbins=20, min_n=25):
    """the published estimator: bin the cells on `by`, weighted within-bin SD of `col`"""
    d = df.dropna(subset=[col, by]).copy()
    d["bin"] = pd.qcut(d[by], nbins, duplicates="drop", labels=False)
    tot = d[col].std(ddof=1)
    B = [(len(g), g[col].std(ddof=1)) for _, g in d.groupby("bin") if len(g) >= min_n]
    n = np.array([b[0] for b in B])
    s = np.array([b[1] for b in B])
    return float(np.sqrt(np.average(s ** 2, weights=n))) / tot


def resid(v, dep, ng):
    M = np.column_stack([np.ones(len(v)), dep, ng])
    b, *_ = np.linalg.lstsq(M, v, rcond=None)
    return v - M @ b


def dctrl(x, y, dep, ng):
    """depth-controlled Spearman, as in 04/05/12: rank, regress out depth, correlate"""
    M = np.column_stack([np.ones(len(x)), dep, ng])

    def r(v):
        rv = rankdata(v)
        b, *_ = np.linalg.lstsq(M, rv, rcond=None)
        return rv - M @ b
    return float(spearmanr(r(x), r(y)).statistic)


# ---------------------------------------------------------------- (1) Table S11 as read
DB = pd.read_csv(f"{T}/S11c_dispersion_baselines.csv")
g = DB[DB.conditioning_group.isin(["RELA_activity", "null_tf_mean", "floor_self",
                                   "inflammation_score"])]
E30 = g.pivot_table(index=["cohort", "lesion", "n_cells"], columns="conditioning_group",
                    values=["frac_retained", "variance_retained"]).reset_index()
E30.columns = ["_".join([a for a in c if a]) for c in E30.columns.values]
E30["null_minus_RELA_variance"] = (E30.variance_retained_null_tf_mean
                                   - E30.variance_retained_RELA_activity).round(4)
E30["RELA_removed_more_than_the_null_reference"] = np.where(
    E30.null_minus_RELA_variance > 0, "yes", "no")
E30["score_conditioned"] = DB.score_conditioned.iloc[0]
E30["note"] = "read from Table S11 (dispersion baselines); " + NOTE
E30 = E30.set_index("lesion").loc[ORDER].reset_index()
E30.to_csv(f"{OUT}/E30_retained_dispersion_by_lesion.csv", index=False)

# ------------------------------------------- (2) the same statistic on residualised scores
SD = ad.read_h5ad(f"{P}/stromal_discovery.h5ad", backed="r").obs.copy()
SV = ad.read_h5ad(f"{P}/stromal_validation.h5ad", backed="r").obs.copy()
AD = pd.read_parquet(f"{P}/stromal_discovery_tf_activity.parquet")
AV = pd.read_parquet(f"{P}/stromal_validation_tf_activity.parquet")
NULLS = {c: sorted(DB[(DB.conditioning_group == "null_tf_individual")
                      & (DB.cohort == c)].conditioning.str.replace("null_tf_", "").unique())
         for c in ("discovery_n5", "validation_n3")}
NB = DB[DB.conditioning == "RELA_activity"].set_index("lesion").n_bins.to_dict()
rows = []
for coh, O, ACT in (("discovery_n5", SD, AD), ("validation_n3", SV, AV)):
    O["RELA"] = ACT["RELA"].values
    for tf in NULLS[coh]:
        O[f"nt_{tf}"] = ACT[tf].values
    O["logdepth"] = np.log10(O.total_counts.astype(float))
    O["n_genes"] = O.n_genes_by_counts.astype(float)
    for les, d in O.groupby("lesion", observed=True):
        d = d.copy()
        nb = int(NB[les])
        d["lf_D_residualised"] = resid(d.lf_D_tissueforming.values.astype(float),
                                       d.logdepth.values, d.n_genes.values)
        rec = dict(cohort=coh, lesion=str(les), n_cells=len(d), n_bins=nb)
        d["inflammation_score"] = d.lf_A_cytokine_chemokine.values.astype(float)
        for col, tag in (("lf_D_tissueforming", "raw_score"),
                         ("lf_D_residualised", "residualised_score")):
            rela = matched_residual(d, col, "RELA", nbins=nb)
            infl = matched_residual(d, col, "inflammation_score", nbins=nb)
            rec[f"{tag}_inflammation_score_frac_retained"] = round(infl, 4)
            rec[f"{tag}_inflammation_score_variance_retained"] = round(infl ** 2, 4)
            nul = [matched_residual(d, col, f"nt_{tf}", nbins=nb) for tf in NULLS[coh]]
            flo = matched_residual(d, col, col, nbins=nb)
            rec[f"{tag}_RELA_frac_retained"] = round(rela, 4)
            rec[f"{tag}_RELA_variance_retained"] = round(rela ** 2, 4)
            rec[f"{tag}_null_mean_frac_retained"] = round(float(np.mean(nul)), 4)
            rec[f"{tag}_null_mean_variance_retained"] = round(float(np.mean(nul)) ** 2, 4)
            rec[f"{tag}_null_min_frac_retained"] = round(float(np.min(nul)), 4)
            rec[f"{tag}_null_max_frac_retained"] = round(float(np.max(nul)), 4)
            rec[f"{tag}_floor_frac_retained"] = round(flo, 4)
            rec[f"{tag}_floor_variance_retained"] = round(flo ** 2, 4)
            rec[f"{tag}_null_minus_RELA_variance"] = round(
                float(np.mean(nul)) ** 2 - rela ** 2, 4)
        rows.append(rec)
E31 = pd.DataFrame(rows).set_index("lesion").loc[ORDER].reset_index()
E31["residualisation"] = ("score residualised on log10(total_counts) and "
                          "n_genes_by_counts before binning; the binning variable is the "
                          "raw ULM activity, as in the published statistic")
E31["note"] = NOTE
E31.to_csv(f"{OUT}/E31_retained_dispersion_residualised.csv", index=False)

# ------------------------------------- (6) composite depth-controlled rho per compartment
AO = pd.read_parquet(f"{P}/all_cells_discovery_obs_scored.parquet")
VO = ad.read_h5ad(f"{P}/all_cells_validation.h5ad", backed="r").obs.copy()
COMPS = ["stromal", "pericyte_smc", "endothelial", "immune_myeloid", "immune_plasma",
         "immune_T"]
rows = []
for coh, O in (("discovery_n5", AO), ("validation_n3", VO)):
    O = O.copy()
    O["logdepth"] = np.log10(O.total_counts.astype(float) + 1)
    O["n_genes"] = O.n_genes_by_counts.astype(float)
    key = "celltype_fine" if "celltype_fine" in O.columns else "compartment"
    for comp in COMPS:
        sub = O[O[key].astype(str) == comp]
        if not len(sub):
            continue
        for les in sorted(sub.lesion.astype(str).unique()) + ["POOLED"]:
            d = sub if les == "POOLED" else sub[sub.lesion.astype(str) == les]
            if les != "POOLED" and len(d) < 100:
                rows.append(dict(cohort=coh, compartment=comp, lesion=les, n_cells=len(d),
                                 excluded="fewer than 100 cells"))
                continue
            rows.append(dict(
                cohort=coh, compartment=comp, lesion=les, n_cells=len(d), excluded="",
                composite_raw_rho=round(float(spearmanr(d.inflammatory_composite,
                                                        d.tissueforming_composite).statistic), 4),
                composite_depth_controlled_rho=round(dctrl(
                    d.inflammatory_composite.values, d.tissueforming_composite.values,
                    d.logdepth.values, d.n_genes.values), 4)))
E32 = pd.DataFrame(rows)
E32["score_columns"] = "inflammatory_composite x tissueforming_composite"
E32["note"] = ("Table S8 carries the single-axis leak-free values only (lf_A x lf_D); "
               "these are the composite values for the same cells. ") + NOTE
E32.to_csv(f"{OUT}/E32_compartment_composite_rho.csv", index=False)

# ------------------------------------------- (7) what the discovery top decile is made of
rows = []
for les, d in list(AO.groupby("lesion", observed=True)) + [("POOLED_discovery", AO)]:
    thr = np.quantile(d.p_D_tissueforming.values, 0.90)
    top = d[d.p_D_tissueforming.values >= thr]
    for comp in sorted(AO.compartment.astype(str).unique()):
        rows.append(dict(lesion=str(les), compartment=comp, n_cells=len(d),
                         n_top_decile=len(top),
                         n_of_compartment_in_top_decile=int((top.compartment == comp).sum()),
                         pct_of_top_decile=round(100 * float((top.compartment == comp).mean()), 3),
                         pct_overall=round(100 * float((d.compartment == comp).mean()), 2)))
E33 = pd.DataFrame(rows)
E33["score"] = "p_D_tissueforming (as in 28_manuscript_numbers.py, Table S8)"
E33["note"] = NOTE
E33.to_csv(f"{OUT}/E33_top_decile_composition.csv", index=False)

# --------------------------------- (8) receptor pseudobulk differences and model estimates
PB = pd.read_csv(f"{T}/10_receptors_pseudobulk_log2cpm.csv")
MM = pd.read_csv(f"{T}/17_receptors_mixedmodel.csv")
genes = MM.gene.tolist()
piv = PB.pivot_table(index=["cohort", "lesion"], columns="state", values=genes)
dif = pd.DataFrame({g: piv[(g, "inflammatory")] - piv[(g, "tissue-forming")] for g in genes})
E34 = pd.DataFrame({
    "gene": genes,
    "pseudobulk_mean_diff_log2CPM": [round(dif[g].mean(), 3) for g in genes],
    "pseudobulk_min_diff": [round(dif[g].min(), 3) for g in genes],
    "pseudobulk_max_diff": [round(dif[g].max(), 3) for g in genes],
    "n_lesions_with_positive_diff": [int((dif[g] > 0).sum()) for g in genes],
    "n_lesions_evaluable": [int(dif[g].notna().sum()) for g in genes],
}).merge(MM[["gene", "coef_inflammatory_minus_tissueforming", "ci_low", "ci_high",
             "pvalue", "padj_BH"]], on="gene")
E34["units"] = ("pseudobulk difference in log2 CPM per lesion; model coefficient on "
                "log-normalised per-cell expression (natural log), so the two are not on "
                "the same scale")
E34["note"] = NOTE
E34.to_csv(f"{OUT}/E34_receptor_pseudobulk_vs_model.csv", index=False)

# ------------------------------------- (10) where each single-axis pooled value comes from
R04 = pd.read_csv(f"{T}/04_rho_inflam_vs_tissueforming_discovery.csv")
RB = pd.read_csv(f"{T}/14_robustness_genesets_scoring.csv")
p = R04[(R04.lesion == "POOLED_discovery") & (R04.cellset == "stromal")]
rb = RB[(RB.lesion == "POOLED") & (RB.inflammatory_set == "current_lf_A_cytokine_chemokine")
        & (RB.tissueforming_set == "current_lf_D_tissueforming")
        & (RB.scoring == "score_genes")]
E35 = pd.DataFrame([
    dict(value_in_the_manuscript="-0.065", where="Robustness",
         quantity=float(p[p.definition == "axis"].depth_controlled_rho.iloc[0]),
         score_convention="plain (p_) single-axis score, RELA targets retained",
         cell_set="discovery stromal, pooled over the five lesions",
         n_cells=int(p[p.definition == "axis"].n.iloc[0]),
         scored_on="all 44,433 discovery cells, then subset to stromal",
         source_file="04_rho_inflam_vs_tissueforming_discovery.csv"),
    dict(value_in_the_manuscript="-0.070", where="Figure S3a, first cell",
         quantity=float(rb.depth_controlled_rho.iloc[0]),
         score_convention="leak-free (lf_) single-axis score, RELA targets removed",
         cell_set="discovery stromal, pooled over the five lesions",
         n_cells=int(rb.n.iloc[0]),
         scored_on="the stromal object alone (control bins drawn within stromal cells)",
         source_file="14_robustness_genesets_scoring.csv"),
    dict(value_in_the_manuscript="-0.03", where="Results",
         quantity=float(p[p.definition == "leakfree"].depth_controlled_rho.iloc[0]),
         score_convention="leak-free (lf_) single-axis score, RELA targets removed",
         cell_set="discovery stromal, pooled over the five lesions",
         n_cells=int(p[p.definition == "leakfree"].n.iloc[0]),
         scored_on="all 44,433 discovery cells, then subset to stromal",
         source_file="04_rho_inflam_vs_tissueforming_discovery.csv"),
])
E35["note"] = ("all three are depth-controlled Spearman rho on the same cells; they "
               "differ in the score convention and in whether the scores were computed "
               "before or after taking the stromal subset. ") + NOTE
E35.to_csv(f"{OUT}/E35_single_axis_pooled_values.csv", index=False)

# ----------------------- (A3) provenance of every programme gene set, gene by gene
import json
S1J = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
PGJ = json.load(open(f"{W}/frozen/genesets/program_geneset_definitions.json"))
GMT = {}
for ln in open(f"{W}/ref/Reactome_Pathways_2024.gmt", errors="replace"):
    p = ln.rstrip("\n").split("\t")
    GMT[p[0]] = set(g for g in p[2:] if g)
DV = set(ad.read_h5ad(f"{P}/all_cells_discovery.h5ad", backed="r").var.symbol.astype(str))
VV = set(ad.read_h5ad(f"{P}/all_cells_validation.h5ad", backed="r").var.symbol.astype(str))


def covered(pathways):
    import difflib
    hit, miss = set(), []
    for nm in [x.strip() for x in pathways.split(";")] if pathways else []:
        key = nm if nm in GMT else (difflib.get_close_matches(nm, list(GMT), n=1,
                                                              cutoff=0.72) or [None])[0]
        if key:
            hit |= GMT[key]
        else:
            miss.append(nm)
    return hit, miss


rows = []
for k, v in S1J["axes"].items():
    src = S1J["axis_sources"][k]
    for g in sorted(v):
        rows.append(dict(set_name=k, set_type="programme axis (Table S2 screen axis)",
                         family="inflammatory" if k.startswith("A") else
                                ("tissue-forming" if k.startswith("D") else "other"),
                         gene_symbol=g, source_statement=src, reactome_pathways="",
                         n_genes_in_set=len(v), n_curated_added_declared="",
                         in_current_reactome_gmt="", in_discovery_gene_space=g in DV,
                         in_validation_gene_space=g in VV))
for m in PGJ["component_meta"]:
    comp = m["component"]
    genes = PGJ["components"][comp]["genes"]
    hit, miss = covered(m["reactome_sets"])
    for g in sorted(genes):
        rows.append(dict(set_name=comp, set_type="composite component",
                         family=m["family"].replace("_", "-"), gene_symbol=g,
                         source_statement=m["source"],
                         reactome_pathways=m["reactome_sets"]
                         + (f" [not found in the 2024 GMT: {'; '.join(miss)}]" if miss else ""),
                         n_genes_in_set=len(genes),
                         n_curated_added_declared=m["n_curated_added"],
                         in_current_reactome_gmt=g in hit,
                         in_discovery_gene_space=g in DV, in_validation_gene_space=g in VV))
E37 = pd.DataFrame(rows)
E37["note"] = ("the frozen record states how many genes were added by hand per component "
               "but not which ones, so in_current_reactome_gmt is a check against the "
               "2024 Reactome release and not the provenance of record. ") + NOTE
E37.to_csv(f"{OUT}/E37_programme_set_provenance.csv", index=False)

RSD = PGJ["global_robust_sd"]
E38 = pd.DataFrame([
    dict(composite=fam, component=c, weight=f"1/{len(PGJ[key]['components'])}",
         global_robust_sd=RSD[c],
         n_genes=PGJ["components"][c]["n"] if c in PGJ["components"] else "",
         divisor_definition="1.4826 * MAD over all cells of the annotated object",
         n_cells_used_for_the_divisor=42613,
         object_used="data/human_annot.h5ad, the pooled discovery object as annotated "
                     "when the sets were frozen (42,613 cells x 17,944 genes); the "
                     "re-processed discovery object of this paper has 44,433 cells, and "
                     "the divisors were deliberately not recomputed - they are constants "
                     "applied identically to every cell, donor and compartment",
         scoring="scanpy score_genes, 25 expression bins, 50 control genes per gene")
    for fam, key in (("inflammatory", "composite_inflammatory"),
                     ("tissue-forming", "composite_tissueforming"))
    for c in PGJ[key]["components"]]
    + [dict(composite="kept separate", component="TP2a_core44_FROZEN",
            weight="not averaged", global_robust_sd=RSD["TP2a_core44_FROZEN"], n_genes="",
            divisor_definition="1.4826 * MAD over all cells of the annotated object",
            n_cells_used_for_the_divisor=42613, object_used="as above",
            scoring="as above")])
E38["note"] = NOTE
E38.to_csv(f"{OUT}/E38_composite_constants.csv", index=False)

for f, d in [("E30", E30), ("E31", E31), ("E37", E37), ("E38", E38), ("E32", E32), ("E33", E33), ("E34", E34),
             ("E35", E35)]:
    print(f, d.shape)
print("E30 DONE")
