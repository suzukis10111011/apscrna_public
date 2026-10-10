#!/usr/bin/env python
"""Pre-submission tidy-up: number record and new supplementary tables (requests 2, 6).

(2) results/tables/28_manuscript_numbers.csv is rebuilt with three new columns
    (score_convention, cell_set, source_file) plus cohort and source_script. Every
    row whose label claims a leak-free score is checked against the column the value
    actually comes from; rows 38-41 claim "(lf_)" but 05:143 correlates
    p_A_cytokine_chemokine with p_D_tissueforming, so the label is corrected to the
    plain score and the leak-free value for the same cells is added as a new row.
    Rows are also added for the post-audit quantities the revised manuscript cites.

(6) Table S11 (random-factor null, overlap-removed candidate results, residual
    dispersion null and floor) and Table S12 (sex-marker single-donor check, full
    provenance map) are assembled from 05_post_audit/tables.

Nothing in results/ is modified.

Output (06_submission_prep/tables):
    28_manuscript_numbers_corrected.csv
    28_manuscript_numbers_corrections_log.csv
    S11a_random_factor_null.csv  S11b_overlap_removed_candidates.csv
    S11c_dispersion_baselines.csv  S11c_null_factor_selection.csv
    S12a_single_donor_check.csv  S12b_number_provenance.csv
    Supplementary_Tables_S11_S12.xlsx
"""
import os
import re
import sys
import time

import numpy as np
import pandas as pd

W = "/path/to/Documents/apscrna_v2"
TAB = f"{W}/results/tables"
AUD = f"{W}/05_post_audit/tables"
OUT = f"{W}/06_submission_prep/tables"
NOTE = ("pre-submission correction; frozen gene sets, candidate factors and "
        "thresholds unchanged")


def log(m):
    s = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {m}"
    print(s, flush=True)
    open(f"{W}/06_submission_prep/logs/run_submission_prep.log", "a").write(s + "\n")


# ----------------------------------------------------------------- (2) record
MSN = pd.read_csv(f"{TAB}/28_manuscript_numbers.csv")
S3 = pd.read_csv(f"{OUT}/S3_rho_by_score_convention.csv")
NQ = pd.read_csv(f"{AUD}/Q1b_null_summary.csv")
BL = pd.read_csv(f"{AUD}/Q2_definitionA_baselines.csv")
DM = pd.read_csv(f"{AUD}/Q3_donor_mixture.csv")
UNI = pd.read_csv(f"{OUT}/definitionA_RELA_input_unification.csv")

SRC_META = {
    "28_lolo_grant_factors.csv": ("leak-free (lf_) axis score", "stromal",
                                  "discovery_n5", "34_lolo_rulebased.py"),
    "28_enrichment_permutation.csv": ("leak-free (lf_) axis score", "all cells",
                                      "discovery_n5", "28 (permutation test)"),
    "25_validation_states_within_vs_ingest.csv": (
        "plain (p_) axis score", "stromal", "validation_n3",
        "25_validation_states_within.py"),
    "17_receptors_detection_fraction.csv": ("not a score (detection fraction)",
                                            "stromal", "both", "17_receptor_power.py"),
    "13_gse164241_summary.csv": ("leak-free (lf_) axis score", "stromal (gingiva)",
                                 "GSE164241", "13_gse164241.py"),
    "05_rho_inflam_vs_tissueforming_validation.csv": (
        "plain (p_) axis score", "stromal", "validation_n3",
        "05_validation_and_projection.py"),
    "28_qc_all_stages.csv": ("not a score (cell counts)", "all cells", "both",
                             "28 (QC consolidation)"),
    "32_set_definition.csv": ("leak-free (lf_) axis score", "stromal", "both",
                              "32_coactivity_ruleset.py"),
    "32_coactivity_ruleset_sensitivity.csv": ("TF activity sets, not programme scores",
                                              "stromal", "both",
                                              "32_coactivity_ruleset.py"),
    "32_pc1_is_function_of_r.csv": ("TF activity sets, not programme scores",
                                    "stromal", "both", "32_coactivity_ruleset.py"),
    "32_crossside_pairs_sign_consistency.csv": ("TF activity sets, not programme "
                                                "scores", "stromal", "both",
                                                "32_coactivity_ruleset.py"),
    "32_regulon_overlap.csv": ("not a score (regulon overlap)", "n.a.", "n.a.",
                               "32_coactivity_ruleset.py"),
}

rows, corrections = [], []
for i, r in MSN.iterrows():
    item, src = str(r["item"]), str(r["source"])
    conv, cells, coh, script = SRC_META.get(
        src, ("(not classified)", "(not classified)", "(not classified)", ""))
    claims_lf = bool(re.search(r"\(lf_\)|leak-?free", item, re.I))
    actually_lf = conv.startswith("leak-free")
    corrected_item, cnote = item, ""
    if claims_lf and not actually_lf:
        corrected_item = item.replace("(lf_)", "(plain p_ score)")
        cnote = ("label claimed a leak-free score; 05:143 correlates "
                 "p_A_cytokine_chemokine with p_D_tissueforming, so the value is the "
                 "plain whole-axis score. Label corrected; the leak-free value for "
                 "the same cells is added as a separate row.")
        corrections.append(dict(row_index=i, original_item=item,
                                corrected_item=corrected_item, value=r["value"],
                                action="label corrected", reason=cnote, note=NOTE))
    rows.append(dict(item=corrected_item, value=r["value"], source=src,
                     score_convention=conv, cell_set=cells, cohort=coh,
                     source_script=script, source_file=src,
                     label_claimed_leak_free=claims_lf,
                     value_is_leak_free=actually_lf,
                     corrected_in_this_pass=bool(cnote), correction_note=cnote,
                     origin="published record", note=str(r["note"])))

LES_MAP = {"CAP1": "CAP1", "CAP2": "CAP2", "CAP3": "CAP3",
           "POOLED_validation": "POOLED_validation"}
for les in LES_MAP:
    q = S3[(S3.cohort == "validation_n3") & (S3.cell_set == "stromal")
           & (S3.lesion == les)].iloc[0]
    rows.append(dict(
        item=f"Validation single-axis (leak-free lf_ score) depth-controlled rho ({les})",
        value=f"{q.depth_controlled_rho_leakfree_lf:.4f} (raw "
              f"{q.raw_rho_leakfree_lf:+.4f}, n={int(q.n_cells)})",
        source="06_submission_prep/tables/S3_rho_by_score_convention.csv",
        score_convention="leak-free (lf_) axis score", cell_set="stromal",
        cohort="validation_n3", source_script="06_submission_prep/scripts/r1_tables.py",
        source_file="S3_rho_by_score_convention.csv",
        label_claimed_leak_free=True, value_is_leak_free=True,
        corrected_in_this_pass=True,
        correction_note="added so that the record carries both score conventions for "
                        "the same cells",
        origin="added in this pass", note=NOTE))
    corrections.append(dict(row_index=np.nan, original_item="(none)",
                            corrected_item=f"Validation single-axis (leak-free lf_ "
                                           f"score) depth-controlled rho ({les})",
                            value=f"{q.depth_controlled_rho_leakfree_lf:.4f}",
                            action="row added", reason="the record previously held "
                            "only the plain-score value for this quantity", note=NOTE))

for axis, lab in [("lf_A_cytokine_chemokine", "inflammatory"),
                  ("lf_D_tissueforming", "tissue-forming")]:
    q = NQ[(NQ.axis == axis) & (NQ.version == "as_screened")
           & (NQ.pool == "sign_consistent_non_candidates")].iloc[0]
    rows.append(dict(
        item=f"Random-factor null for candidate sign reproduction ({lab} programme)",
        value=f"candidates {int(q.candidates_sign_agree_all3)}/"
              f"{int(q.n_candidates_in_denominator)} vs null mean "
              f"{q.null_all3_mean:.2f} (max {int(q.null_all3_max)}), "
              f"empirical p = {q.empirical_p_sign_all3}",
        source="05_post_audit/tables/Q1b_null_summary.csv",
        score_convention="leak-free (lf_) axis score", cell_set="stromal",
        cohort="validation_n3",
        source_script="05_post_audit/scripts/q1b_null_signtest.py",
        source_file="Q1b_null_summary.csv", label_claimed_leak_free=True,
        value_is_leak_free=True, corrected_in_this_pass=True,
        correction_note="post-audit quantity now cited in the manuscript",
        origin="added in this pass", note=NOTE))
    q2 = NQ[(NQ.axis == axis) & (NQ.version == "overlap_removed")
            & (NQ.pool == "sign_consistent_non_candidates")].iloc[0]
    rows.append(dict(
        item=f"Candidate sign reproduction with regulon overlap removed ({lab})",
        value=f"{int(q2.candidates_sign_agree_all3)}/"
              f"{int(q2.n_candidates_in_denominator)} vs null mean "
              f"{q2.null_all3_mean:.2f}, empirical p = {q2.empirical_p_sign_all3}"
              + (f"; {q2.candidates_dropped} removed by tmin=5"
                 if isinstance(q2.candidates_dropped, str) and q2.candidates_dropped
                 else ""),
        source="05_post_audit/tables/Q1b_null_summary.csv",
        score_convention="leak-free (lf_) axis score", cell_set="stromal",
        cohort="validation_n3",
        source_script="05_post_audit/scripts/q1b_null_signtest.py",
        source_file="Q1b_null_summary.csv", label_claimed_leak_free=True,
        value_is_leak_free=True, corrected_in_this_pass=True,
        correction_note="post-audit quantity", origin="added in this pass", note=NOTE))

B = BL[BL.conditioning_group.isin(["RELA_activity", "null_tf_mean", "floor_self"])]
for grp, lab in [("RELA_activity", "RELA activity (reported conditioning)"),
                 ("null_tf_mean", "unrelated factors (mean of 10)"),
                 ("floor_self", "floor, conditioned on the score itself")]:
    g = B[B.conditioning_group == grp]
    rows.append(dict(
        item=f"Definition A retained dispersion, all 8 lesions - {lab}",
        value=f"SD ratio {g.frac_retained.min():.4f}-{g.frac_retained.max():.4f}; "
              f"variance ratio {g.variance_retained.min():.4f}-"
              f"{g.variance_retained.max():.4f}",
        source="05_post_audit/tables/Q2_definitionA_baselines.csv",
        score_convention="leak-free (lf_) tissue-forming score",
        cell_set="stromal", cohort="both",
        source_script="05_post_audit/scripts/q2_defA_baselines.py",
        source_file="Q2_definitionA_baselines.csv", label_claimed_leak_free=True,
        value_is_leak_free=True, corrected_in_this_pass=True,
        correction_note="post-audit baseline now cited in the manuscript",
        origin="added in this pass", note=NOTE))

d = DM[DM.judged_on.str.startswith("XIST")]
rows.append(dict(
    item="Single-donor check, discovery lesions (sex-marker discordance)",
    value=f"minority ratio {d.minority_ratio.min():.4f}-{d.minority_ratio.max():.4f} "
          f"over {len(d)} lesions",
    source="05_post_audit/tables/Q3_donor_mixture.csv",
    score_convention="not a score (marker counts)", cell_set="all cells",
    cohort="discovery_n5", source_script="05_post_audit/scripts/q3_donor_mixture.py",
    source_file="Q3_donor_mixture.csv", label_claimed_leak_free=False,
    value_is_leak_free=False, corrected_in_this_pass=True,
    correction_note="post-audit quantity now cited in the manuscript",
    origin="added in this pass", note=NOTE))
v = DM[~DM.judged_on.str.startswith("XIST")]
rows.append(dict(
    item="Single-donor check, validation lesions (Y-linked transcript fraction)",
    value="; ".join(f"{r.lesion} {100 * r.frac_any_Y_positive:.1f}%"
                    for _, r in v.iterrows()) + " (XIST absent from the feature space)",
    source="05_post_audit/tables/Q3_donor_mixture.csv",
    score_convention="not a score (marker counts)", cell_set="all cells",
    cohort="validation_n3", source_script="05_post_audit/scripts/q3_donor_mixture.py",
    source_file="Q3_donor_mixture.csv", label_claimed_leak_free=False,
    value_is_leak_free=False, corrected_in_this_pass=True,
    correction_note="post-audit quantity now cited in the manuscript",
    origin="added in this pass", note=NOTE))

R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/28_manuscript_numbers_corrected.csv", index=False)
pd.DataFrame(corrections).to_csv(f"{OUT}/28_manuscript_numbers_corrections_log.csv",
                                 index=False)
log(f"record: {len(MSN)} published rows -> {len(R)} rows, "
    f"{int(R.corrected_in_this_pass.sum())} corrected or added")

# ------------------------------------------------------------- (6) S11 and S12
s11a = pd.read_csv(f"{AUD}/Q1b_null_summary.csv")
s11a.insert(0, "table_part", "random-factor null")
s11a.to_csv(f"{OUT}/S11a_random_factor_null.csv", index=False)

det = pd.read_csv(f"{AUD}/Q1b_candidate_detail.csv")
drop = pd.read_csv(f"{AUD}/Q1a_tmin_dropped_factors.csv")
s11b = det[det.version == "overlap_removed"].copy()
s11b.insert(0, "table_part", "regulon-overlap-removed activities")
s11b = s11b.merge(drop[drop.is_frozen_candidate][["axis", "TF", "regulon_size_full",
                                                  "regulon_size_after_removal"]]
                  .drop_duplicates(), on=["axis", "TF"], how="left")
s11b.to_csv(f"{OUT}/S11b_overlap_removed_candidates.csv", index=False)

s11c = pd.read_csv(f"{AUD}/Q2_definitionA_baselines.csv")
s11c.insert(0, "table_part", "retained dispersion, conditioning baselines")
s11c.to_csv(f"{OUT}/S11c_dispersion_baselines.csv", index=False)
sel = pd.read_csv(f"{AUD}/Q2_null_tf_selection.csv")
sel.insert(0, "table_part", "factors used for the unrelated-factor baseline")
sel.to_csv(f"{OUT}/S11c_null_factor_selection.csv", index=False)

mix = pd.read_csv(f"{AUD}/Q3_donor_mixture.csv")
pool = pd.read_csv(f"{AUD}/Q3b_source_pooling_and_sex.csv")
s12a = pool.merge(mix[["lesion", "female_markers_used", "male_markers_used"]],
                  on="lesion", how="left")
s12a.insert(0, "table_part", "single-donor check")
s12a.to_csv(f"{OUT}/S12a_single_donor_check.csv", index=False)
s12b = pd.read_csv(f"{AUD}/Q4b_number_provenance.csv")
s12b.insert(0, "table_part", "provenance of every reported number")
s12b.to_csv(f"{OUT}/S12b_number_provenance.csv", index=False)

SHEETS = [("S11a_random_null", s11a), ("S11b_overlap_removed", s11b),
          ("S11c_dispersion_baselines", s11c), ("S11c_null_factors", sel),
          ("S12a_single_donor", s12a), ("S12b_provenance", s12b)]
with pd.ExcelWriter(f"{OUT}/Supplementary_Tables_S11_S12.xlsx") as xl:
    pd.DataFrame([
        dict(sheet="S11a_random_null",
             description="Sign reproduction of the 20 frozen candidates against "
                         "1,000 random draws of 20 non-candidate factors of the same "
                         "axis, with two specificity statistics"),
        dict(sheet="S11b_overlap_removed",
             description="The same candidates after removing, from the regulon "
                         "network, targets that belong to the axis gene set; factors "
                         "left with fewer than five targets are named"),
        dict(sheet="S11c_dispersion_baselines",
             description="Retained dispersion of the tissue-forming score under four "
                         "conditioning variables, all eight lesions"),
        dict(sheet="S11c_null_factors",
             description="The ten factors used for the unrelated-factor baseline and "
                         "their correlation with the tissue-forming axis"),
        dict(sheet="S12a_single_donor",
             description="Sex-specific transcript counts per lesion and the pooling "
                         "statements of the source publications"),
        dict(sheet="S12b_provenance",
             description="Every reported correlation, residual dispersion ratio, "
                         "candidate agreement count and cell count with its script, "
                         "file, cell set and score definition"),
    ]).to_excel(xl, sheet_name="index", index=False)
    for name, df in SHEETS:
        df.to_excel(xl, sheet_name=name[:31], index=False)
log(f"S11/S12 written: {[ (n, len(d)) for n, d in SHEETS ]}")
