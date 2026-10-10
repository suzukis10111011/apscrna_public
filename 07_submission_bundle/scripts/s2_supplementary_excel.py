#!/usr/bin/env python
"""Submission bundle: Supplementary Tables S1-S12 in one workbook.

Reads the current results/tables (post-swap) and the S11/S12 sheets prepared in
06_submission_prep. results/ is read only; everything is written to
07_submission_bundle/tables.

Layout per sheet: row 1 column names, row 2 a definition of each column (grey
italic), row 3 onwards the data. Tables assembled from more than one file carry a
table_part column and a source_file column. An index sheet lists the twelve tables,
their parts and their row counts.

S2 carries the full-precision screen (the 4-dp rounding was removed when the files
were swapped into results/). S3 is the wide form whose column names state the score
convention (composite, plain p_, leak-free lf_). S5 states the leak-free convention
per row and carries the candidate sign counts on the depth-controlled estimator.
S11 and S12 are the post-audit tables.

Output: 07_submission_bundle/tables/Supplementary_Tables_S1-S12.xlsx
        07_submission_bundle/tables/supplementary_table_index.csv
"""
import os
import re
import sys
import time

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill

W = "/path/to/Documents/apscrna_v2"
TAB = f"{W}/results/tables"
SUB = f"{W}/06_submission_prep/tables"
EXT = f"{W}/08_external_review/tables"
OUT = f"{W}/07_submission_bundle/tables"
os.makedirs(OUT, exist_ok=True)

JP2EN = {
    "クラスタ判定": "cluster-based assignment",
    "候補は 5病変同符号かつ |pooled ρ| 上位20。旧セットから ATF1・MTA1 が外れ、RELA・NFKB1 が入る":
        "candidates are TFs with a consistent sign across all 5 lesions and |pooled rho| "
        "in the top 20; relative to the previous set, ATF1 and MTA1 drop out and RELA "
        "and NFKB1 enter",
    "標準化2変数の PCA では PC1 分散比 = (1+|Pearson r|)/2 なので、相関と独立な証拠ではない":
        "for a PCA of two standardised variables the PC1 variance ratio equals "
        "(1+|Pearson r|)/2, so it is a deterministic function of the correlation and "
        "not independent evidence",
    "n<20 のため算出せず": "not computed (n<20)",
    "細胞数<100": "n cells < 100",
    "scVI は環境に無いため省略。axis = lf_D_tissueforming":
        "scVI not available in this environment and was omitted; axis = lf_D_tissueforming",
    "病変内 Q80 以上 かつ 4陰性マーカーが病変内中央値未満":
        "within-lesion 80th percentile or above for the positive markers, and all 4 "
        "negative markers below the within-lesion median",
    "病変×状態で counts 合計 -> CPM -> log2(CPM+1)":
        "counts summed per lesion x state, then CPM, then log2(CPM+1)",
    "MixedLM: log正規化発現 ~ state, ランダム切片 = lesion（探索+検証を合算）":
        "MixedLM: log-normalised expression ~ state, random intercept = lesion "
        "(discovery and validation pooled)",
    "今回の区画判定による間質集合": "stromal set from the compartment assignment of this study",
    "前回の間質バーコード集合（今回の行列・正規化・スコアで再計算）":
        "previous stromal barcode set, recomputed with this study's matrix, "
        "normalisation and scores",
    "composite (inflammatory vs tissueforming)": "composite (inflammatory vs tissue-forming)",
}
JP_FALLBACK = [("病変内", "within-lesion"), ("以上 かつ", " or above and "),
               ("陰性マーカーが病変内中央値未満",
                "negative markers below the within-lesion median"),
               ("は行列に不在で未適用", " absent from the matrix, filter not applied"),
               ("細胞数", "n cells"), ("探索", "discovery"), ("検証", "validation"),
               ("病変", "lesion"), ("状態", "state"), ("区画", "compartment"),
               ("（", " ("), ("）", ")"), ("、", ", "), ("。", ". ")]
JP = re.compile(r"[\u3040-\u30ff\u4e00-\u9fff]")


def clean(d):
    def conv(v):
        if not isinstance(v, str) or not JP.search(v):
            return v
        if v in JP2EN:
            return JP2EN[v]
        out = v
        for a, b in JP_FALLBACK:
            out = out.replace(a, b)
        return out
    d = d.copy()
    for c in d.columns:
        if d[c].dtype == object:
            d[c] = d[c].map(conv)
    return d


SHEETS = {
    "S1": ("Cell counts at each QC stage, per lesion, both cohorts, and the coverage "
           "of the gene sets and of the marker panel in the quantified gene space",
           [("28_qc_all_stages.csv", TAB, "main"),
            ("20b_geneset_coverage_before_after.csv", TAB,
             "gene-set and marker-panel coverage")]),
    "S1b_validation_chemistry": (
        "Chemistry of the three validation runs checked from read structure: read-1 "
        "length, base composition around the barcode/UMI boundary, and barcode "
        "agreement with the 10x v2 and v3 inclusion lists",
        [("C6_validation_chemistry_check.csv", SUB, "per run"),
         ("C6_base_composition_by_position.csv", SUB, "base composition by position")]),
    "S2": ("Full 730-factor activity screen on all six axes, discovery ranking and "
           "validation sign reproduction, at full precision",
           [("04_tf_screen_discovery_candidates.csv", TAB, "discovery_n5"),
            ("05_tf_screen_validation_reproduction.csv", TAB, "validation_n3"),
            ("E21_figS8A_counts.csv", EXT,
             "candidate reproduction per programme, counted (as in Figure S8A)"),
            ("E24_geneset_members.csv", EXT, "programme gene sets, members"),
            ("E25_validation_all_factors.csv", EXT,
             "validation values of every factor of the two programmes "
             "(as in Figure 3A, B)")]),
    "S3": ("Raw and depth-controlled Spearman rho between the inflammatory and "
           "tissue-forming scores under all three score conventions, and the "
           "unimodality of the two scores per lesion",
           [("S3_rho_by_score_convention.csv", TAB, "main"),
            ("E9_programme_score_dip.csv", EXT,
             "unimodality of the programme scores per lesion (Hartigan dip test, as "
             "in Figure 4E)")]),
    "S4": ("Robustness across gene-set definitions and scoring methods",
           [("14_robustness_genesets_scoring.csv", TAB, "main")]),
    "S5": ("Robustness to the integration method and to an embedding-free, "
           "marker-defined stromal set, with the decile ranked on the LEAK-FREE score "
           "(lf_D_tissueforming) throughout the sensitivity parts; the enrichment "
           "series of the main analysis and of Figure 2 is the separate part "
           "'primary enrichment per lesion', whose decile is ranked on the full "
           "programme score p_D_tissueforming",
           [("14_robustness_integration.csv", TAB, "integration method"),
            ("18_stromal_markerdefined_vs_primary.csv", TAB,
             "stromal cell-set definition, decile ranked on the leak-free score "
             "(sensitivity analysis, not the Figure 2 series)"),
            ("18_tf_sign_by_cellset.csv", TAB, "candidate sign by cell set"),
            ("E8_projection_with_not_adopted.csv", EXT,
             "projected state labels vs the within-cohort definition"),
            ("E9_primary_enrichment.csv", EXT,
             "primary enrichment per lesion (as in Figure 2)"),
            ("E15_sign_agreement_part.csv", EXT,
             "candidate sign agreement between the two stromal sets, counted"),
            ("E23_projection_confusion.csv", EXT,
             "within-cohort against projected state labels, cell counts "
             "(as in Figure S7B)")]),
    "S6": ("GEO author matrices versus the independent kallisto|bustools "
           "re-quantification",
           [("09_geo_vs_kb.csv", TAB, "key quantities"),
            ("09_stromal_jaccard.csv", TAB, "stromal barcode agreement"),
            ("09_tf_agreement_geo_vs_kb.csv", TAB, "TF screen agreement")]),
    "S7": ("Cytokine-receptor expression between the inflammatory and tissue-forming "
           "states",
           [("10_receptors_by_state.csv", TAB, "pseudobulk, paired Wilcoxon"),
            ("17_receptors_mixedmodel.csv", TAB, "cell-level linear mixed model"),
            ("17_receptors_detection_fraction.csv", TAB, "detection fraction")]),
    "S8": ("Axis coupling per compartment, and the permutation test for top-decile "
           "compartment enrichment",
           [("12_coupling_by_compartment.csv", TAB, "main"),
            ("28_enrichment_permutation.csv", TAB, "top-decile permutation test"),
            ("E3_compartment_composition.csv", EXT,
             "compartment composition per lesion, both annotation versions")]),
    "S9": ("Independent periodontitis cohort GSE164241 (gingiva)",
           [("13_gse164241_summary.csv", TAB, "per donor"),
            ("13_gse164241_tf_sign_reproduction.csv", TAB, "TF sign reproduction"),
            ("E22_figS8B_counts.csv", EXT,
             "candidate reproduction per programme, counted (as in Figure S8B)")]),
    "S10": ("Co-activity of the two transcription-factor programmes, with the TF sets "
            "defined by rule from this study's own screen",
            [("32_set_definition.csv", TAB, "TF set definition by rule"),
             ("32_coactivity_ruleset_sensitivity.csv", TAB, "per-lesion and per-state statistics"),
             ("32_crossside_pairs_sign_consistency.csv", TAB, "cross-side TF pairs"),
             ("32_regulon_overlap.csv", TAB, "regulon overlap"),
             ("32_old_vs_new_comparison.csv", TAB, "rule-based vs previous TF set"),
             ("32_pc1_is_function_of_r.csv", TAB, "note on the PC1 variance ratio"),
             ("33_state_pct_8lesions.csv", TAB, "state proportions per lesion"),
             ("E6_myc_vs_set_means.csv", EXT,
              "MYC activity against the two TF sets"),
             ("E26_tf_activity_matrix.csv", EXT,
              "factor-by-factor activity correlation matrix (as in Figure 5D)")]),
    "S11": ("Random factor-set null, overlap-removed activities, and the null and "
            "floor references for retained dispersion (post hoc, added after "
            "independent audit)",
            [("S11a_random_factor_null.csv", SUB, "random-factor null"),
             ("S11b_overlap_removed_candidates.csv", SUB, "overlap-removed activities"),
             ("S11c_dispersion_baselines.csv", SUB, "dispersion baselines"),
             ("S11c_null_factor_selection.csv", SUB, "unrelated-factor selection"),
             ("33_defA_permutation_null.csv", TAB,
              "permutation null for the retained dispersion")]),
    "S12": ("Sex-marker check of single-donor origin, and the provenance of every "
            "reported value (script, output file, cell set, score definition)",
            [("S12a_single_donor_check.csv", SUB, "single-donor check"),
             ("S12b_number_provenance.csv", SUB, "provenance of reported values")]),
}

DEF = {
    "cohort": "discovery_n5 (lesions A-E) or validation_n3 (CAP1-CAP3)",
    "lesion": "lesion identifier; POOLED rows combine the lesions of that cohort",
    "donor": "donor identifier in the external cohort",
    "group": "disease group in the external cohort",
    "stage": "quality-control stage",
    "n_cells": "number of cells in this row's group",
    "n": "number of cells in this row's group",
    "n_cells_qc_passed": "cells passing quality control in this lesion",
    "cellset": "cell set the statistic was computed on",
    "cell_set": "cell set the statistic was computed on",
    "compartment": "compartment assigned by the marker-score majority rule",
    "axis": "programme axis; lf_ prefix means RELA target genes were removed",
    "TF": "transcription factor (CollecTRI regulon)",
    "raw_rho": "Spearman rho between the two programme scores, no covariates",
    "depth_controlled_rho": "Spearman rho of the residuals of the ranks on "
                            "[1, log10(total UMIs), number of detected genes]",
    "raw_rho_composite": "raw Spearman rho, composite scores",
    "raw_rho_plain_p": "raw Spearman rho, whole-axis p_ scores (RELA targets retained)",
    "raw_rho_leakfree_lf": "raw Spearman rho, leak-free lf_ scores (RELA targets removed)",
    "depth_controlled_rho_composite": "depth-controlled rho, composite scores",
    "depth_controlled_rho_plain_p": "depth-controlled rho, whole-axis p_ scores",
    "depth_controlled_rho_leakfree_lf": "depth-controlled rho, leak-free lf_ scores",
    "depth_controlled_rho_as_in_04_05": "depth-controlled rho with the depth covariate "
                                        "log10(total UMIs), as in the main scripts",
    "score_columns": "the two score columns correlated in this row",
    "score_columns_composite": "score columns behind the composite values",
    "score_columns_plain_p": "score columns behind the plain p_ values",
    "score_columns_leakfree_lf": "score columns behind the leak-free lf_ values",
    "depth_covariate": "covariates used for the depth-controlled estimate",
    "definition": "score convention: axis (plain p_), composite, or leakfree (lf_)",
    "pooled_rho": "depth-controlled rho over all cells of the cohort",
    "abs_pooled": "absolute value of pooled_rho, the ranking key of the screen",
    "n_pos": "number of lesions with a positive rho",
    "all5_same_sign": "the rho had the same sign in all five discovery lesions",
    "all_same_sign": "the rho had the same sign in every lesion of the cohort",
    "all_same_sign_published_plain_spearman": "the same flag on the published plain "
                                              "Spearman estimator",
    "pooled_rho_published_plain_spearman": "the published plain Spearman value",
    "estimator": "estimator behind pooled_rho in this row",
    "rank_within_axis": "rank by |pooled rho| among the sign-consistent factors",
    "is_candidate": "selected as a candidate (sign-consistent and in the top 20)",
    "precision": "whether the file was rounded",
    "in_validation": "the factor had an activity estimate in the validation cohort",
    "discovery_pooled_rho": "pooled depth-controlled rho in the discovery cohort",
    "validation_pooled_rho": "pooled depth-controlled rho in the validation cohort",
    "n_val_same_sign_as_discovery": "validation lesions matching the discovery sign",
    "reproduced": "the discovery sign held in all three validation lesions",
    "frac_retained": "within-bin SD of the tissue-forming score divided by its total SD",
    "variance_retained": "square of frac_retained",
    "sd_reduction_pct": "100 x (1 - frac_retained)",
    "variance_reduction_pct": "100 x (1 - variance_retained)",
    "total_sd": "SD of the tissue-forming score over all cells of the group",
    "matched_sd": "cell-count-weighted within-bin SD",
    "n_bins": "bins retained after dropping those with fewer than 25 cells",
    "bins_used": "bins requested (10 when the lesion has fewer than 500 cells)",
    "median_bin_n": "median number of cells per retained bin",
    "conditioning": "variable the cells were matched on",
    "conditioning_group": "family the conditioning variable belongs to",
    "defA_frac_retained": "retained dispersion, RELA activity from a ULM on the "
                          "reported cell set",
    "defA_frac_retained_published_allcell_ulm": "the published value, from a ULM on "
                                                "the whole all-cell object",
    "defA_RELA_input": "matrix the RELA activity was estimated on",
    "defA_n_bins": "bins retained for the definition A estimate",
    "defA_bins_used": "bins requested for the definition A estimate",
    "jaccard_vs_primary": "Jaccard index of the marker-defined and cluster-defined "
                          "stromal sets",
    "threshold": "rule defining the cell set",
    "stromal_pct_overall": "stromal cells as a percentage of the group",
    "stromal_pct_top_decile": "stromal cells as a percentage of the top decile of the "
                              "tissue-forming score",
    "enrichment_fold": "stromal_pct_top_decile divided by stromal_pct_overall",
    "n_top_decile": "cells in the top decile of the tissue-forming score",
    "p_permutation": "one-sided empirical p against the permutation null",
    "n_perm": "number of permutations",
    "one_sided": "direction of the one-sided test",
    "seed": "random seed",
    "pool": "set the random factors were drawn from",
    "pool_size": "number of factors available to draw from",
    "version": "activity version: as screened, or regulon overlap removed",
    "n_seeds": "number of random draws",
    "n_candidates_frozen": "candidates in the frozen list for this axis",
    "n_candidates_in_denominator": "candidates still evaluable in this activity version",
    "candidates_dropped": "candidates removed by the tmin=5 filter",
    "candidates_sign_agree_all3": "candidates keeping the discovery sign in all three "
                                  "validation lesions",
    "candidates_sign_agree_majority": "the same count under a two-of-three rule",
    "null_all3_mean": "mean of that count over the random draws",
    "null_all3_sd": "SD over the random draws",
    "null_all3_q025": "2.5th percentile of the null",
    "null_all3_q975": "97.5th percentile of the null",
    "null_all3_max": "largest count reached by any random draw",
    "empirical_p_sign_all3": "(draws at least as extreme + 1) / (draws + 1)",
    "empirical_p_sign_majority": "the same p under the two-of-three rule",
    "candidates_median_abs_validation_rho": "median |validation pooled rho| over the "
                                            "candidates",
    "null_median_abs_rho_mean": "the same statistic averaged over the random draws",
    "empirical_p_effect": "empirical p for the effect-size statistic",
    "observed_jaccard_disc_top20_vs_val_top20": "Jaccard index of the discovery and "
                                                "validation top-20 lists",
    "null_jaccard_mean": "the same index for the random draws",
    "empirical_p_jaccard": "empirical p for the overlap statistic",
    "discovery_top20": "the discovery top-20 factors",
    "validation_top20": "the validation top-20 factors",
    "in_denominator": "the candidate was evaluable in this activity version",
    "dropped_reason": "why the candidate was not evaluable",
    "regulon_size_full": "targets in the factor's regulon before removal",
    "regulon_size_after_removal": "targets left after removing the axis gene set",
    "reproduced_reported_rule": "the pre-registered all-three-lesion rule",
    "n_validation_lesions_same_sign": "validation lesions matching the discovery sign",
    "pooled_depth_controlled_rho_with_tissueforming": "correlation with the "
                                                      "tissue-forming axis used to "
                                                      "select the unrelated factors",
    "selection": "rule by which the factor was selected",
    "n_XIST_only": "cells with XIST and no Y-linked transcript",
    "n_Y_only": "cells with a Y-linked transcript and no XIST",
    "n_both_positive": "cells positive for both",
    "minority_ratio": "smaller of the two single-marker counts over their sum",
    "frac_any_Y_positive": "fraction of cells with any of RPS4Y1, DDX3Y, UTY, KDM5D",
    "judged_on": "markers available in that matrix",
    "declared_sex": "sex stated by the source record",
    "declared_donor": "donor attributes stated by the source",
    "marker_based_call": "call from the transcript counts",
    "agrees_with_declared": "the marker call matches the source record",
    "pooling_statement": "what the source publication says about pooling",
    "publication": "source publication",
    "accession": "repository accession for the sample",
    "category": "class of reported quantity",
    "quantity": "the reported quantity",
    "manuscript_location": "where the value appears in the manuscript",
    "source_script": "script that produced the value",
    "output_file": "file the value lives in",
    "score_definition": "score convention behind the value",
    "ulm_input": "matrix the TF activities were estimated on",
    "table_part": "which input file this block of rows came from",
    "source_file": "file this block of rows came from",
    "note": "provenance or caveat recorded with the row",
}
DEF.update({
    'median_genes': 'median detected genes per cell at this stage',
    'median_umi': 'median UMIs per cell at this stage',
    'median_pct_mt': 'median percentage of UMIs from mitochondrial genes',
    'source_file_original': 'source_file column of the input file, renamed to avoid a clash',
    'table_part_original': 'table_part column of the input file, renamed to avoid a clash',
    'inflammatory_set': 'gene set used for the inflammatory programme',
    'tissueforming_set': 'gene set used for the tissue-forming programme',
    'scoring': 'scoring method used for the per-cell programme score',
    'geo': 'value from the author-deposited GEO matrices',
    'kb_allbarcodes': 'value from the kallisto|bustools re-quantification, all barcodes',
    'kb_geobarcodes': 'value from the re-quantification restricted to the GEO barcodes',
    'diff_kb_allbarcodes_minus_geo': 're-quantification (all barcodes) minus GEO',
    'diff_kb_geobarcodes_minus_geo': 're-quantification (GEO barcodes) minus GEO',
    'comparison': 'the two things being compared in this row',
    'n_geo': 'size of the GEO-derived set',
    'n_other': 'size of the comparison set',
    'intersection': 'cells or factors in both sets',
    'union': 'cells or factors in either set',
    'jaccard': 'intersection over union',
    'pct_of_geo_recovered': 'percentage of the GEO set recovered',
    'source': 'input the row refers to',
    'n_candidates_geo': 'candidates from the GEO matrices',
    'n_candidates_kb': 'candidates from the re-quantification',
    'overlap': 'candidates common to both',
    'geo_only': 'candidates only from GEO',
    'kb_only': 'candidates only from the re-quantification',
    'n_TF_common': 'factors evaluable in both inputs',
    'same_sign': 'factors with the same sign of pooled rho in both inputs',
    'same_sign_pct': 'same_sign as a percentage of n_TF_common',
    'spearman_rho_of_pooled_rho': 'rank correlation of the pooled rho between the two inputs',
    'gene': 'receptor gene',
    'n_lesions': 'lesions contributing to this row',
    'mean_log2cpm_inflammatory': 'mean log2 CPM in inflammatory-state cells',
    'mean_log2cpm_tissueforming': 'mean log2 CPM in tissue-forming-state cells',
    'delta_infl_minus_tiss': 'inflammatory minus tissue-forming mean log2 CPM',
    'wilcoxon_stat': 'paired Wilcoxon statistic across lesions',
    'pvalue': 'nominal p value',
    'padj_BH': 'Benjamini-Hochberg adjusted p value',
    'pseudobulk': 'how the pseudobulk was formed',
    'genes_absent_from_matrix': 'requested genes not present in the matrix',
    'status': 'whether the model could be fitted',
    'coef_inflammatory_minus_tissueforming': 'mixed-model coefficient for state, inflammatory minus tissue-forming',
    'ci_low': 'lower bound of the 95% confidence interval',
    'ci_high': 'upper bound of the 95% confidence interval',
    'group_var': 'random-effect grouping variable',
    'converged': 'the mixed model converged',
    'model': 'model specification',
    'state': 'descriptive stromal state',
    'detection_fraction': 'fraction of cells with a non-zero count',
    'mean_lognorm': 'mean log-normalised expression',
    'included_in_model': 'the row entered the mixed model',
    'defA_total_sd': 'SD of the tissue-forming score before matching',
    'defA_matched_sd': 'cell-count-weighted within-bin SD after matching',
    'present': 'the factor had an activity estimate in this cohort',
    'n_donors': 'donors contributing to this row',
    'pooled_rho_gse164241': 'pooled depth-controlled rho in the gingival cohort',
    'all_donors_same_sign': 'the rho had the same sign in every donor',
    'sign_matches_discovery': 'the pooled sign matches the discovery sign',
    'set_name': 'name of the transcription-factor set',
    'included_by_rule': 'the factor satisfied the inclusion rule',
    'label': 'label of the series shown',
    'spearman_rho': 'Spearman rho between the two set-mean activities',
    'spearman_p': 'p value for that correlation',
    'pc1_var_frac': 'variance fraction on the first principal component of the two set means',
    'pc1_angle_deg': 'angle of the first principal component, degrees',
    'dip_p_infl': 'Hartigan dip-test p for the inflammatory set mean',
    'dip_p_tiss': 'Hartigan dip-test p for the tissue-forming set mean',
    'infl_TF': 'factor from the inflammatory set',
    'tiss_TF': 'factor from the tissue-forming set',
    'n_positive': 'lesions with a positive rho',
    'rho_min': 'smallest per-lesion rho',
    'rho_max': 'largest per-lesion rho',
    'TF_a': 'first factor of the pair',
    'TF_b': 'second factor of the pair',
    'n_a': 'regulon size of TF_a',
    'n_b': 'regulon size of TF_b',
    'n_shared': 'targets shared by the two regulons',
    'frac_of_smaller': 'shared targets as a fraction of the smaller regulon',
    'rho_in_script_11': 'the value obtained with the previous factor set',
    'reproduces_script_11': 'the rule-based set reproduces that value',
    'implied_abs_pearson_r': 'absolute Pearson r implied by the PC1 variance fraction',
    'n_tf': 'factors compared between the two embeddings',
    'stromal_n_harmony': 'stromal cells under Harmony integration',
    'stromal_n_uncorrected': 'stromal cells without integration',
    'stromal_intersection': 'cells stromal under both',
    'stromal_jaccard': 'Jaccard index of the two stromal sets',
    'tf_rank_spearman': 'rank correlation of the factor ranking between the two embeddings',
    'top20_overlap': 'candidates common to the two embeddings',
    'excluded': 'why the row was excluded from the estimate',
    'pct_overall': 'percentage of all cells of the group',
    'pct_top_decile': 'percentage of the top decile of the tissue-forming score',
    'null_mean_pct': 'mean of the permutation null, percentage',
    'null_p975_pct': '97.5th percentile of the permutation null, percentage',
    'myeloid_M1_vs_M2_raw_rho': 'raw rho between the myeloid M1-like and M2-like scores',
    'myeloid_M1_vs_M2_depth_controlled_rho': 'the same rho, depth-controlled',
    'reported_reproduced_count': 'count published in the pre-registered reproduction table',
    'n_draws_ge_candidates_all3': 'random draws reaching the candidate count',
    'null_majority_mean': 'null mean under the two-of-three rule',
    'n_draws_ge_candidates_majority': 'random draws reaching the candidate count under that rule',
    'null_median_abs_rho_q975': '97.5th percentile of the null effect size',
    'candidate_percentile_in_null_effect': 'percentile of the candidate effect size within the null',
    'null_jaccard_q975': '97.5th percentile of the null overlap',
    'null_tf_frac_retained_sd': 'SD of frac_retained over the ten unrelated factors',
    'null_tf_frac_retained_min': 'smallest frac_retained over the ten unrelated factors',
    'null_tf_frac_retained_max': 'largest frac_retained over the ten unrelated factors',
    'score_conditioned': 'score whose dispersion is being measured',
    'female_markers_used': 'female-specific transcripts present in that matrix',
    'male_markers_used': 'Y-linked transcripts present in that matrix',
    'value': 'the reported value as printed',
    'script': 'script that produced the value',
})

DEF.update({
    "run_accession": "SRA run accession",
    "donor_label": "lesion label this run was assigned to",
    "read_source": "how the reads were obtained for this check",
    "n_reads_examined": "reads examined (the first 200,000 of the run)",
    "read1_length": "length of read 1 in bases",
    "read1_length_counts": "observed read-1 lengths and their counts",
    "n_distinct_barcodes": "distinct 16-base barcodes among the reads examined",
    "n_top_barcodes": "barcodes ranked by read count that were tested against the "
                      "inclusion lists",
    "top5000_in_v2_737K_pct": "percentage of the top barcodes in the 10x v2 "
                              "(737K-august-2016) inclusion list",
    "top5000_in_v3_3M_pct": "percentage of the top barcodes in the 10x v3 "
                            "(3M-february-2018) inclusion list",
    "all_reads_in_v2_737K_pct": "percentage of all reads whose barcode is in the v2 list",
    "all_reads_in_v3_3M_pct": "percentage of all reads whose barcode is in the v3 list",
    "pos27_T_pct": "percentage of reads with T at position 27 of read 1",
    "pos28_T_pct": "percentage of reads with T at position 28 of read 1",
    "expected_if_v2": "what the two statistics look like if the library is v2",
    "expected_if_v3": "what they look like if the library is v3",
    "indirect_evidence": "evidence used where the reads could not be retrieved",
    "position": "position in read 1, counting from 1",
    "A": "percentage of reads with A at this position",
    "C": "percentage of reads with C at this position",
    "G": "percentage of reads with G at this position",
    "T": "percentage of reads with T at this position",
    "N": "percentage of reads with an ambiguous base at this position",
    "n_reads": "reads contributing to this position",
})


DEF.update({
    # S1, coverage of the gene sets and of the marker panel (20b)
    "gene_set": "gene set or marker panel whose coverage is counted",
    "n_def": "genes in the definition of the set",
    "found_before": "genes of the set present in the quantified gene space before the "
                    "symbol correction",
    "found_after": "genes of the set present after the symbol correction",
    "pct_before": "found_before as a percentage of n_def",
    "pct_after": "found_after as a percentage of n_def",
    "recovered": "found_after minus found_before",
    # S5, projected state labels (25)
    "labelling": "how the state label was obtained: within_cohort (clustering inside "
                 "the validation cohort) or ingest_projected (projected from the "
                 "discovery embedding)",
    "n_inflammatory": "stromal cells labelled inflammatory",
    "pct_inflammatory": "n_inflammatory as a percentage of the stromal cells of the row",
    "n_intermediate": "stromal cells labelled intermediate",
    "pct_intermediate": "n_intermediate as a percentage of the stromal cells of the row",
    "n_tissue-forming": "stromal cells labelled tissue-forming",
    "pct_tissue-forming": "n_tissue-forming as a percentage of the stromal cells of the "
                          "row",
    "n_other": "stromal cells labelled other",
    "pct_other": "n_other as a percentage of the stromal cells of the row",
    "agreement_pct": "percentage of cells given the same state label by the projection "
                     "and by the within-cohort definition",
    # S8, compartment composition (E3)
    "annotation_version": "compartment annotation version: 03c assigns every cell and "
                          "is the version used in the figures and in the other parts of "
                          "this table; 03b leaves part of the cells unassigned",
    "level": "compartment (six compartments) or fine (immune cells split into B, T, "
             "mast, myeloid and plasma)",
    "label": "compartment or fine cell type of the row",
    "n_total": "all cells of the lesion at this level of the annotation",
    "pct_of_all_cells": "n_cells as a percentage of n_total",
    "source_file_of_part": "file in results/tables the row was taken from",
    "scheme_adopted": "yes for the scheme used in the figures and in the rest of the "
                      "paper, no for a scheme the Methods describe as evaluated and "
                      "not adopted",
    "value_source": "how the row was obtained: recomputed from the named file, or "
                    "recorded in the Methods when the inputs of that run were not "
                    "retained",
    "gene_space": "shared gene space the projection ran over",
    # S5, enrichment: the sensitivity series (leak-free score) and the main series
    "stromal_pct_top_decile_leakfree_sensitivity":
        "percentage of the top decile that is stromal when the decile is ranked on the "
        "LEAK-FREE score lf_D_tissueforming; sensitivity analysis, NOT the series in "
        "Figure 2 (that one is stromal_pct_top_decile_primary)",
    "enrichment_fold_leakfree_sensitivity":
        "stromal_pct_top_decile_leakfree_sensitivity divided by stromal_pct_overall. "
        "SENSITIVITY ANALYSIS: the decile is ranked on the leak-free score "
        "lf_D_tissueforming and, in the marker_defined rows, the stromal set is the "
        "marker-defined one. The enrichment of the main analysis and of Figure 2 is "
        "enrichment_fold_primary, whose decile is ranked on the full programme score "
        "p_D_tissueforming; for the cluster-based rows the ranking score is the only "
        "difference between the two (CAP1 8.682 here against 8.787 there)",
    "stromal_pct_overall": "percentage of all cells of the lesion that belong to the "
                           "cell set named in cellset; independent of which score "
                           "ranks the decile, so it is shared by both series",
    "score_ranking_the_decile": "score whose top decile is taken",
    "stromal_set": "which stromal set the percentages are computed on",
    "stromal_pct_overall_primary": "percentage of all cells of the lesion that are "
                                   "stromal (main analysis)",
    "stromal_pct_top_decile_primary": "percentage of the top decile of the FULL "
                                      "programme score that is stromal (main analysis)",
    "enrichment_fold_primary": "stromal_pct_top_decile_primary divided by "
                               "stromal_pct_overall_primary; the series in Figure 2 and "
                               "in the ranges quoted in the Abstract and the Results",
    "figure": "figure the part is plotted in",
    # the counted parts behind Figures S7B and S8
    "figure_panel": "panel the counts are plotted in",
    "programme_gene_set": "programme the gene belongs to",
    "n_genes_in_set": "genes prespecified for that programme",
    "gene_symbol": "gene as it is named in the frozen definition",
    "in_discovery_gene_space": "TRUE when the symbol is present in the discovery matrix",
    "in_validation_gene_space": "TRUE when the symbol is present in the validation "
                                "matrix; the counts reproduce the coverage part of S1",
    "discovery_is_candidate": "TRUE for the 20 candidates of that programme",
    "discovery_rank": "rank of the factor within the programme in the discovery screen",
    "same_sign_pooled": "TRUE when the validation pooled rho has the discovery sign",
    "n_val_lesions_same_sign": "validation lesions in which the sign agrees with the "
                               "discovery pooled rho (0 to 3)",
    "val_rho_CAP1": "depth-controlled rho in CAP1",
    "val_rho_CAP2": "depth-controlled rho in CAP2",
    "val_rho_CAP3": "depth-controlled rho in CAP3",
    "n_reproduced": "candidates of that programme that reproduced",
    "n_evaluable": "candidates evaluable in that cohort, i.e. the denominator",
    "reproduction": "n_reproduced / n_evaluable, the form the legend prints",
    "definition_of_reproduced": "what counts as a reproduction in this part",
    "main_programme": "TRUE for the two programmes analysed in the main text",
    "state_projected": "state label carried over from the discovery cohort",
    "n_cells_confusion": "cells with that pair of labels",
    "agrees": "TRUE on the diagonal, where the two labellings give the same state",
    "n_cells_agreeing_total": "cells on the diagonal, i.e. the agreeing cells",
    "n_cells_total_confusion": "cells in the confusion matrix",
    "agreement_pct_confusion": "n_cells_agreeing_total as a percentage of "
                               "n_cells_total_confusion",
    # S5, the counted sign agreement between the two stromal sets
    "estimator_of_the_pooled_rho": "which estimate of the pooled rho the signs were "
                                   "taken from; the published file was written with "
                                   "the plain Spearman one",
    "estimator_of_the_pooled_rho_note": "what that estimator is",
    "n_agree": "candidate TFs whose pooled rho has the same sign in the marker-defined "
               "and in the cluster-based stromal set",
    "n_evaluable": "candidate TFs evaluable in both sets, i.e. the denominator",
    "agreement": "n_agree / n_evaluable, the form the Methods print; the rows with "
                 "axis = ALL SIX AXES carry 117/120 (discovery) and 110/112 "
                 "(validation) under the published plain-Spearman estimator; the rows "
                 "with axis = THE TWO PROGRAMMES (A and D) carry 40/40 in both cohorts "
                 "under the depth-controlled estimator and 40/40 and 39/40 under the "
                 "plain-Spearman one",
    "where_the_counts_come_from": "the part of this sheet the counts were taken from",
    # S3, unimodality of the programme scores
    "programme_score": "which of the two programme scores the test was run on",
    "dip_statistic_programme_score": "Hartigan dip statistic of that score in that "
                                     "lesion",
    "dip_p_programme_score": "p value of the dip test; a large value is consistent with "
                             "unimodality",
    # S10, state proportions (33) and MYC (E3)
    "inflammatory": "percentage of the stromal cells of the lesion in the inflammatory "
                    "state",
    "intermediate": "percentage in the intermediate state (high on both programmes)",
    "tissue-forming": "percentage in the tissue-forming state",
    "other": "percentage in neither state",
    "set_name": "TF set the row refers to: the rule-defined NF-kB set, the "
                "rule-defined osteoblast-lineage set, or the alternative inflammatory "
                "set with ATF1 and MTA1",
    "factors_in_set": "members of the set, or the single factor when the row is one "
                      "factor against MYC",
    "n_factors": "members of the set the activity estimate was available for",
    "summary": "which quantity the row reports: the correlation of MYC activity with "
               "the mean activity of the set, the mean or the median of the "
               "correlations with the individual factors, or one factor against MYC",
    "spearman_rho_MYC_vs_set_mean": "Spearman correlation described by the summary "
                                    "column, over the stromal cells of the row",
    "spearman_p": "p value of that Spearman correlation (not computed for the mean and "
                  "median summaries)",
    # S11, permutation null for the retained dispersion (33)
    "null_mean": "mean retained dispersion over the permutations",
    "null_lo": "2.5th percentile of the permutation null",
    "null_hi": "97.5th percentile of the permutation null",
    "p_permutation": "one-sided permutation p for the retained dispersion being below "
                     "the null",
    "n_perm": "number of permutations",
    "one_sided": "direction of the one-sided test",
    "seed": "random seed of the permutation",
})

PATTERNS = [(re.compile(r"^rho_(A|B|C|D|E)$"),
             "depth-controlled rho in discovery lesion {}"),
            (re.compile(r"^rho_(CAP[123])$"),
             "depth-controlled rho in validation lesion {}"),
            (re.compile(r"^val_rho_(CAP[123])$"),
             "depth-controlled rho in validation lesion {}")]


def define(col):
    if col in DEF:
        return DEF[col]
    for rx, tmpl in PATTERNS:
        m = rx.match(col)
        if m:
            return tmpl.format(m.group(1))
    return ""


TABLE_OF = {"S1b_validation_chemistry": "S1"}

index_rows = []
frames = {}
for sheet, (desc, parts) in SHEETS.items():
    blocks = []
    for fn, d, tag in parts:
        p = f"{d}/{fn}"
        t = clean(pd.read_csv(p))
        if "source_file" in t.columns:
            t = t.rename(columns={"source_file": "source_file_original"})
        if fn == "18_stromal_markerdefined_vs_primary.csv":
            # these rank the top decile on the leak-free score and belong to the
            # cell-set sensitivity analysis; the main-analysis series is its own part
            t = t.rename(columns={
                "stromal_pct_top_decile": "stromal_pct_top_decile_leakfree_sensitivity",
                "enrichment_fold": "enrichment_fold_leakfree_sensitivity"})
        if "table_part" in t.columns:
            t = t.rename(columns={"table_part": "table_part_original"})
        t.insert(0, "table_part", tag)
        t.insert(1, "source_file", fn)
        blocks.append(t)
        index_rows.append(dict(table=TABLE_OF.get(sheet, sheet), sheet=sheet,
                               description=desc, part=tag, source_file=fn,
                               n_rows=len(t), n_columns=t.shape[1]))
    frames[sheet] = pd.concat(blocks, ignore_index=True)

IDX = pd.DataFrame(index_rows)
IDX.to_csv(f"{OUT}/supplementary_table_index.csv", index=False)

xlsx = f"{OUT}/Supplementary_Tables_S1-S12.xlsx"
LEG = pd.read_csv(f"{EXT}/E20_table_legends.csv")

LEG = pd.read_csv(f"{EXT}/E20_table_legends.csv")
LEG.columns = ["Table", "Legend as it appears in the manuscript"]

with pd.ExcelWriter(xlsx, engine="openpyxl") as xl:
    # first sheet; the data sheets are written after it and keep every cell address
    LEG.to_excel(xl, sheet_name="Legends", index=False)
    IDX.to_excel(xl, sheet_name="index", index=False)
    for sheet, df in frames.items():
        df.to_excel(xl, sheet_name=sheet, index=False, startrow=0)
    wb = xl.book
    grey = Font(italic=True, color="555555", size=9)
    fill = PatternFill("solid", fgColor="F2F2F2")
    for sheet, df in frames.items():
        ws = wb[sheet]
        ws.insert_rows(2)
        for j, col in enumerate(df.columns, start=1):
            c = ws.cell(row=2, column=j, value=define(col))
            c.font = grey
            c.fill = fill
            c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "A3"
    ws = wb["index"]
    ws.freeze_panes = "A2"
    ws = wb["Legends"]
    ws.column_dimensions["A"].width = 11
    ws.column_dimensions["B"].width = 150
    ws.freeze_panes = "A2"
    ws["A1"].font = Font(bold=True)
    ws["B1"].font = Font(bold=True)
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
        for c in row:
            c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[row[0].row].height = 95

n_def = {s: sum(bool(define(c)) for c in d.columns) for s, d in frames.items()}
print(pd.DataFrame([dict(sheet=s, rows=len(d), cols=d.shape[1], defined=n_def[s])
                    for s, d in frames.items()]).to_string(index=False))
print("written", xlsx)
