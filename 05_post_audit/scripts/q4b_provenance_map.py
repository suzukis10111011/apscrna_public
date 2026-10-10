#!/usr/bin/env python
"""POST HOC (audit item 4, the mapping part).

One row per manuscript-facing number in the four categories the audit names
(correlation, residual variance ratio, candidate agreement count, cell count),
carrying the script that produced it, the file it lives in, the cell set it was
computed on, and the rho / conditioning definition behind it. Values are read from
results/tables, never retyped. Figure and supplementary-table locations come from
results/tables/figure_panel_n.csv and scripts/30_supplementary_tables.py.

A second table lists the places where two numbers that describe the same cells and
the same quantity were produced under different definitions.

Output: 05_post_audit/tables/Q4b_number_provenance.csv
        05_post_audit/tables/Q4b_definition_conflicts.csv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q0_common import AUD, POSTHOC_NOTE, TAB, log  # noqa: E402

RD = pd.read_csv(f"{TAB}/04_rho_inflam_vs_tissueforming_discovery.csv")
RV = pd.read_csv(f"{TAB}/05_rho_inflam_vs_tissueforming_validation.csv")
C12 = pd.read_csv(f"{TAB}/12_coupling_by_compartment.csv")
T18 = pd.read_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv")
S18 = pd.read_csv(f"{TAB}/18_tf_sign_by_cellset.csv")
E05 = pd.read_csv(f"{TAB}/05e_rho_prev_cellset_vs_clean.csv")
G13 = pd.read_csv(f"{TAB}/13_gse164241_summary.csv")
G13S = pd.read_csv(f"{TAB}/13_gse164241_tf_sign_reproduction.csv")
SCR = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
REP = pd.read_csv(f"{TAB}/05_tf_screen_validation_reproduction.csv")
A730 = pd.read_csv(f"{TAB}/33_tf_validation_all730.csv")
LOLO = pd.read_csv(f"{TAB}/34_lolo_rulebased.csv")
KB = pd.read_csv(f"{TAB}/09_tf_agreement_geo_vs_kb.csv")
SFX = pd.read_csv(f"{TAB}/26_tf_sign_reproduction_after_symbolfix.csv")
D4 = pd.read_csv(f"{TAB}/04_definitionA_discovery_v2.csv")
D5 = pd.read_csv(f"{TAB}/05_definitionA_validation_v2.csv")
PN = pd.read_csv(f"{TAB}/33_defA_permutation_null.csv")
QC = pd.read_csv(f"{TAB}/28_qc_all_stages.csv")
MSN = pd.read_csv(f"{TAB}/28_manuscript_numbers.csv")

R = []


def add(**kw):
    kw.setdefault("note", POSTHOC_NOTE)
    R.append(kw)


# ------------------------------------------------------------------ correlations
for T, script, coh, loc in [
        (RD, "04_discovery_scores_and_screen.py", "discovery_n5",
         "Fig4A (stromal, composite); FigS5A/B (stromal, both definitions); Table S3"),
        (RV, "05_validation_and_projection.py", "validation_n3",
         "Fig4A (stromal, composite); FigS5A/B (stromal, both definitions); Table S3; "
         "28_manuscript_numbers rows 38-41 (stromal, axis)")]:
    for _, r in T.iterrows():
        add(category="correlation", cohort=coh, group=r.lesion,
            quantity=f"depth-controlled rho ({r.definition} definition)",
            value=r.depth_controlled_rho, n_cells=r.n, script=script,
            output_file=os.path.basename(
                f"0{4 if coh.startswith('disc') else 5}_rho_inflam_vs_tissueforming_"
                f"{'discovery' if coh.startswith('disc') else 'validation'}.csv"),
            cell_set=r.cellset,
            score_definition=("p_A_cytokine_chemokine x p_D_tissueforming (whole axis)"
                              if r.definition == "axis" else
                              "inflammatory_composite x tissueforming_composite"),
            depth_covariate="log10(total_counts), n_genes_by_counts",
            ulm_input="n.a.", manuscript_location=loc)

for _, r in C12[C12.compartment == "stromal"].iterrows():
    add(category="correlation", cohort=r.cohort, group=r.lesion,
        quantity="depth-controlled rho (single axis, leak-free)",
        value=r.depth_controlled_rho, n_cells=r.n_cells,
        script="12_coupling_by_compartment.py",
        output_file="12_coupling_by_compartment.csv", cell_set="stromal",
        score_definition="lf_A_cytokine_chemokine x lf_D_tissueforming "
                         "(RELA targets removed)",
        depth_covariate="log10(total_counts + 1), n_genes_by_counts",
        ulm_input="n.a.", manuscript_location="Fig6A; Table S8")

for _, r in T18.iterrows():
    add(category="correlation", cohort=r.cohort, group=r.lesion,
        quantity="depth-controlled rho (single axis, leak-free)",
        value=r.depth_controlled_rho, n_cells=r.n_cells,
        script="18_stromal_markerdefined.py",
        output_file="18_stromal_markerdefined_vs_primary.csv", cell_set=r.cellset,
        score_definition="lf_A_cytokine_chemokine x lf_D_tissueforming "
                         "(RELA targets removed)",
        depth_covariate="log10(total_counts + 1), n_genes_by_counts",
        ulm_input="n.a.", manuscript_location="Table S5")

for _, r in E05.iterrows():
    add(category="correlation", cohort="validation_n3", group=r.lesion,
        quantity="depth-controlled rho (composite)", value=r.depth_controlled_rho,
        n_cells=r.n, script="05e_rho_prev_cellset.py",
        output_file="05e_rho_prev_cellset_vs_clean.csv", cell_set=r.cellset,
        score_definition="inflammatory_composite x tissueforming_composite",
        depth_covariate="log10(total_counts), n_genes_by_counts", ulm_input="n.a.",
        manuscript_location="diagnostic only; frozen/RESTORED.json states 05d/05e "
                            "are not used for manuscript numbers")

for _, r in G13.iterrows():
    add(category="correlation", cohort="GSE164241", group=f"{r.group}/{r.donor}",
        quantity="depth-controlled rho (single axis, leak-free)",
        value=r.depth_controlled_rho, n_cells=r.n_cells, script="13_gse164241.py",
        output_file="13_gse164241_summary.csv", cell_set="stromal",
        score_definition="lf_A_cytokine_chemokine x lf_D_tissueforming",
        depth_covariate="log10(total_counts + 1), n_genes_by_counts",
        ulm_input="stromal subset", manuscript_location="Fig7A; Table S9")

# ------------------------------------------------ residual variance ratios (def A)
for T, script, fn, coh, loc in [
        (D4, "04_discovery_scores_and_screen.py", "04_definitionA_discovery_v2.csv",
         "discovery_n5", "feeds Fig4B through 33_defA_permutation_null.csv"),
        (D5, "05_validation_and_projection.py", "05_definitionA_validation_v2.csv",
         "validation_n3", "feeds Fig4B through 33_defA_permutation_null.csv")]:
    for _, r in T.iterrows():
        add(category="residual_variance_ratio", cohort=coh, group=r.lesion,
            quantity="definition A frac_retained (SD ratio)", value=r.frac_retained,
            n_cells=r.n_cells, script=script, output_file=fn, cell_set="stromal",
            score_definition="lf_D_tissueforming conditioned on RELA activity",
            depth_covariate="n.a. (20 quantile bins of RELA; 10 when n < 500 in 05)",
            ulm_input="stromal subset", manuscript_location=loc)

for _, r in PN.iterrows():
    add(category="residual_variance_ratio", cohort=r.cohort, group=r.lesion,
        quantity="definition A frac_retained, with permutation null",
        value=r.frac_retained, n_cells=r.n, script="33_ms_figure_inputs.py",
        output_file="33_defA_permutation_null.csv", cell_set="stromal",
        score_definition="lf_D_tissueforming conditioned on RELA activity; "
                         f"null = RELA permuted within lesion, null mean {r.null_mean}",
        depth_covariate="n.a.", ulm_input="stromal subset",
        manuscript_location="Fig4B")

for _, r in C12[C12.compartment == "stromal"].iterrows():
    add(category="residual_variance_ratio", cohort=r.cohort, group=r.lesion,
        quantity="definition A frac_retained (SD ratio)", value=r.defA_frac_retained,
        n_cells=r.n_cells, script="12_coupling_by_compartment.py",
        output_file="12_coupling_by_compartment.csv", cell_set="stromal",
        score_definition="lf_D_tissueforming conditioned on RELA activity",
        depth_covariate="n.a.", ulm_input="all-cell object",
        manuscript_location="Fig6A; Table S8")

for _, r in T18.iterrows():
    add(category="residual_variance_ratio", cohort=r.cohort, group=r.lesion,
        quantity="definition A frac_retained (SD ratio)", value=r.defA_frac_retained,
        n_cells=r.n_cells, script="18_stromal_markerdefined.py",
        output_file="18_stromal_markerdefined_vs_primary.csv", cell_set=r.cellset,
        score_definition="lf_D_tissueforming conditioned on RELA activity",
        depth_covariate="n.a.", ulm_input="all-cell object",
        manuscript_location="Table S5")

# ------------------------------------------------------ candidate agreement counts
g = SCR.groupby("axis").agg(n_tested=("TF", "size"),
                            n_sign_consistent=("all5_same_sign", "sum"),
                            n_candidates=("is_candidate", "sum"))
for axis, r in g.iterrows():
    add(category="candidate_agreement", cohort="discovery_n5", group=axis,
        quantity="candidates selected (sign-consistent across 5 lesions, top 20 by "
                 "|pooled rho|)",
        value=f"{int(r.n_candidates)} of {int(r.n_sign_consistent)} sign-consistent "
              f"of {int(r.n_tested)} tested",
        n_cells=6215, script="04_discovery_scores_and_screen.py",
        output_file="04_tf_screen_discovery_candidates.csv", cell_set="stromal",
        score_definition="depth-controlled rho of TF activity vs the leak-free axis",
        depth_covariate="log10(total_counts), n_genes_by_counts",
        ulm_input="stromal subset", manuscript_location="Fig3C/D; Table S2")

for axis, r in REP.groupby("axis").agg(n=("TF", "size"),
                                       rep=("reproduced", "sum")).iterrows():
    add(category="candidate_agreement", cohort="validation_n3", group=axis,
        quantity="candidates whose sign reproduces in all 3 validation lesions",
        value=f"{int(r.rep)}/{int(r.n)}", n_cells=1779,
        script="05_validation_and_projection.py",
        output_file="05_tf_screen_validation_reproduction.csv", cell_set="stromal",
        score_definition="depth-controlled rho of TF activity vs the leak-free axis; "
                         "sign compared with the discovery pooled rho",
        depth_covariate="log10(total_counts), n_genes_by_counts",
        ulm_input="stromal subset",
        manuscript_location=("Fig3E (lf_A, lf_D); FigS8A (the other four axes); "
                             "Table S2"))

for (coh, cs, axis), r in S18.groupby(["cohort", "cellset", "axis"]).agg(
        n=("TF", "size"), same=("all_same_sign", "sum")).iterrows():
    add(category="candidate_agreement", cohort=coh, group=axis,
        quantity="candidates sign-consistent across the lesions of the cohort",
        value=f"{int(r.same)}/{int(r.n)}", n_cells=np.nan,
        script="18_stromal_markerdefined.py", output_file="18_tf_sign_by_cellset.csv",
        cell_set=cs,
        score_definition="PLAIN Spearman of TF activity vs the leak-free axis "
                         "(no depth control, unlike 04/05/33)",
        depth_covariate="none", ulm_input="all-cell object",
        manuscript_location="Table S5")

for (grp, axis), r in G13S.groupby(["group", "axis"]).agg(
        n=("TF", "size"),
        same=("sign_matches_discovery", lambda s: (s.astype(str) == "True").sum())
        ).iterrows():
    add(category="candidate_agreement", cohort="GSE164241", group=f"{grp}/{axis}",
        quantity="candidates whose pooled sign matches discovery",
        value=f"{int(r.same)}/{int(r.n)}", n_cells=np.nan, script="13_gse164241.py",
        output_file="13_gse164241_tf_sign_reproduction.csv", cell_set="stromal",
        score_definition="depth-controlled rho of TF activity vs the leak-free axis",
        depth_covariate="log10(total_counts + 1), n_genes_by_counts",
        ulm_input="stromal subset", manuscript_location="Fig7B; FigS8B; Table S9")

for axis, r in LOLO.groupby("axis").agg(
        n_rounds=("held_out_lesion", "nunique"),
        n_candidate_rows=("is_candidate", "sum")).iterrows():
    add(category="candidate_agreement", cohort="discovery_n5", group=axis,
        quantity="leave-one-lesion-out candidate rows (5 rounds x 5 named factors)",
        value=f"{int(r.n_candidate_rows)} of {int(r.n_rounds) * 5}", n_cells=np.nan,
        script="34_lolo_rulebased.py", output_file="34_lolo_rulebased.csv",
        cell_set="stromal",
        score_definition="depth-controlled rho of TF activity vs the leak-free axis, "
                         "4 lesions per round",
        depth_covariate="log10(total_counts), n_genes_by_counts",
        ulm_input="stromal subset",
        manuscript_location="Fig3F; FigS6; 28_manuscript_numbers rows 0-10")

for _, r in KB.iterrows():
    add(category="candidate_agreement", cohort="discovery_n5",
        group=f"{r.source}/{r.axis}",
        quantity="candidate overlap, author matrix vs re-quantification",
        value=f"{int(r.overlap)}/20 (Jaccard {r.jaccard})", n_cells=np.nan,
        script="09_geo_vs_kb.py", output_file="09_tf_agreement_geo_vs_kb.csv",
        cell_set="stromal",
        score_definition="depth-controlled rho of TF activity vs the leak-free axis",
        depth_covariate="log10(total_counts), n_genes_by_counts",
        ulm_input="stromal subset", manuscript_location="FigS4; Table S6")

for _, r in SFX.iterrows():
    add(category="candidate_agreement", cohort="validation_n3", group=r.axis,
        quantity="candidate sign reproduction before / after the symbol fix",
        value=f"{r.before_symbolfix} -> {r.after_symbolfix}", n_cells=1779,
        script="26 (symbol-fix check)",
        output_file="26_tf_sign_reproduction_after_symbolfix.csv", cell_set="stromal",
        score_definition="depth-controlled rho of TF activity vs the leak-free axis",
        depth_covariate="log10(total_counts), n_genes_by_counts",
        ulm_input="stromal subset", manuscript_location="methods text")

# ------------------------------------------------------------------- cell counts
for _, r in QC.iterrows():
    add(category="cell_count", cohort=r.cohort, group=f"{r.lesion}/{r.stage}",
        quantity="cells at this QC stage", value=r.n_cells, n_cells=r.n_cells,
        script="28 (QC consolidation)", output_file="28_qc_all_stages.csv",
        cell_set="all cells",
        score_definition="n.a.", depth_covariate="n.a.", ulm_input="n.a.",
        manuscript_location="FigS1A/B; Table S1; 28_manuscript_numbers rows 42-43")
for coh, n_all, n_str in [("discovery_n5", 44433, 6215), ("validation_n3", 26090, 1779)]:
    add(category="cell_count", cohort=coh, group="cohort total",
        quantity="all cells after QC / stromal cells", value=f"{n_all} / {n_str}",
        n_cells=n_all, script="03c / 04 / 05",
        output_file="data/processed/all_cells_*.h5ad, stromal_*.h5ad",
        cell_set="all cells; compartment == stromal", score_definition="n.a.",
        depth_covariate="n.a.", ulm_input="n.a.",
        manuscript_location="Fig1B/C, Fig2, Fig3, Fig4 (n = 6,215 + 1,779 = 7,994)")

P = pd.DataFrame(R)
P.to_csv(f"{AUD}/tables/Q4b_number_provenance.csv", index=False)

# ------------------------------------------------------------- definition clashes
def val(T, **q):
    d = T
    for k, v in q.items():
        d = d[d[k] == v]
    return float(d.iloc[0].depth_controlled_rho) if len(d) else np.nan


def dval(T, col, **q):
    d = T
    for k, v in q.items():
        d = d[d[k] == v]
    return float(d.iloc[0][col]) if len(d) else np.nan


K = []
for les in ["CAP1", "CAP2", "CAP3", "POOLED_validation"]:
    l18 = "POOLED" if les == "POOLED_validation" else les
    a = val(RV, lesion=les, cellset="stromal", definition="axis")
    b = val(T18, cohort="validation_n3", cellset="primary_cluster_based", lesion=l18)
    c = val(C12, cohort="validation_n3", compartment="stromal", lesion=l18)
    d = val(RV, lesion=les, cellset="stromal", definition="composite")
    K.append(dict(
        conflict="C1 same cells, same quantity, three score conventions",
        cohort="validation_n3", group=les,
        side_a="FigS5A / Table S3 / text: 05, p_ scores, log10(depth)", value_a=a,
        side_b="Table S5: 18, lf_ scores, log10(depth+1)", value_b=b,
        side_c="Fig6A / Table S8: 12 does not cover the validation cohort "
               "(12_coupling_by_compartment.csv is discovery only)",
        value_c=c,
        side_d="Fig4A: 05, composite scores", value_d=d,
        max_abs_difference=round(float(np.nanmax([a, b, c, d]) - np.nanmin([a, b, c, d])), 4),
        note=POSTHOC_NOTE))
for les in ["A", "B", "C", "D", "E", "POOLED_discovery"]:
    l18 = "POOLED" if les == "POOLED_discovery" else les
    a = val(RD, lesion=les, cellset="stromal", definition="axis")
    b = val(T18, cohort="discovery_n5", cellset="primary_cluster_based", lesion=l18)
    c = val(C12, cohort="discovery_n5", compartment="stromal", lesion=l18)
    d = val(RD, lesion=les, cellset="stromal", definition="composite")
    K.append(dict(
        conflict="C1 same cells, same quantity, three score conventions",
        cohort="discovery_n5", group=les,
        side_a="FigS5A / Table S3: 04, p_ scores, log10(depth)", value_a=a,
        side_b="Table S5: 18, lf_ scores, log10(depth+1)", value_b=b,
        side_c="Fig6A / Table S8: 12, lf_ scores, log10(depth+1), ULM on all cells",
        value_c=c,
        side_d="Fig4A: 04, composite scores", value_d=d,
        max_abs_difference=round(float(np.nanmax([a, b, c, d]) - np.nanmin([a, b, c, d])), 4),
        note=POSTHOC_NOTE))

for coh, lesions, DT in [("discovery_n5", ["A", "B", "C", "D", "E"], D4),
                         ("validation_n3", ["CAP1", "CAP2", "CAP3"], D5)]:
    for les in lesions:
        a = dval(DT, "frac_retained", lesion=les)
        b = dval(T18, "defA_frac_retained", cohort=coh,
                 cellset="primary_cluster_based", lesion=les)
        c = dval(C12, "defA_frac_retained", cohort=coh, compartment="stromal",
                 lesion=les)
        K.append(dict(
            conflict="C5 definition A: RELA activity from two different ULM inputs",
            cohort=coh, group=les,
            side_a="04/05 (feeds Fig4B): ULM on the stromal subset", value_a=a,
            side_b="Table S5 (18): ULM on the all-cell object", value_b=b,
            side_c=("Fig6A / Table S8 (12): ULM on the all-cell object"
                    if coh == "discovery_n5" else
                    "Fig6A / Table S8 (12): discovery only, no validation row"),
            value_c=c,
            side_d="", value_d=np.nan,
            max_abs_difference=round(float(np.nanmax([a, b, c]) - np.nanmin([a, b, c])), 4),
            note=POSTHOC_NOTE))

msn = MSN[MSN.item.str.contains("single-axis", na=False)]
for _, r in msn.iterrows():
    les = str(r["item"]).split("(")[-1].rstrip(")")
    K.append(dict(
        conflict="C3 label says lf_, value is the p_ column",
        cohort="validation_n3", group=les,
        side_a=f"28_manuscript_numbers: '{r['item']}' = {r['value']}",
        value_a=val(RV, lesion=les, cellset="stromal", definition="axis"),
        side_b="the lf_ value for the same cells (18)",
        value_b=val(T18, cohort="validation_n3", cellset="primary_cluster_based",
                    lesion="POOLED" if "POOLED" in les else les),
        side_c="", value_c=np.nan, side_d="", value_d=np.nan,
        max_abs_difference=np.nan,
        note="05 computes the 'axis' definition from p_A_cytokine_chemokine x "
             "p_D_tissueforming (05:143), not from the lf_ columns. " + POSTHOC_NOTE))

K.append(dict(
    conflict="C6 candidate agreement counted on two different estimators",
    cohort="both", group="lf_A / lf_D",
    side_a="Table S2 (05): depth-controlled rho, all 3 lesions, 20/20 and 19/20",
    value_a=np.nan,
    side_b="Table S5 (18): plain Spearman, no depth control", value_b=np.nan,
    side_c="", value_c=np.nan, side_d="", value_d=np.nan, max_abs_difference=np.nan,
    note="18:111-112 calls spearmanr directly on activity and score. " + POSTHOC_NOTE))
K.append(dict(
    conflict="C7 the repository's own reconciliation is applied to discovery only",
    cohort="both", group="single-axis pooled rho",
    side_a="19_rho_reconciliation.csv names 12 (lf_) canonical and 04 (p_) "
           "not canonical, for the pooled discovery value",
    value_a=dval(C12.rename(columns={"depth_controlled_rho": "v"}), "v",
                 cohort="discovery_n5", compartment="stromal", lesion="POOLED"),
    side_b="the validation counterpart printed in 28_manuscript_numbers is the "
           "p_ value from 05",
    value_b=val(RV, lesion="POOLED_validation", cellset="stromal", definition="axis"),
    side_c="", value_c=np.nan, side_d="", value_d=np.nan, max_abs_difference=np.nan,
    note=POSTHOC_NOTE))

KK = pd.DataFrame(K)
KK.to_csv(f"{AUD}/tables/Q4b_definition_conflicts.csv", index=False)
log(f"Q4b done: provenance rows {len(P)}, conflict rows {len(KK)}")
print(P.category.value_counts().to_string())
print(KK[["conflict", "group", "value_a", "value_b", "value_c", "value_d",
          "max_abs_difference"]].to_string(index=False))
