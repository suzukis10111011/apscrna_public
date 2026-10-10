#!/usr/bin/env python
"""Purpose: Prepare the quantities the external check could not find in S1-S12 so they
  can be added to the existing tables as further parts, and list where each one goes.
Input: results/tables/{20b_geneset_coverage_before_after, 25_validation_states_within_vs_
  ingest, 33_state_pct_8lesions, 33_defA_permutation_null, 03b/03c composition,
  05_compartment_by_lesion_validation, 11_tf_activity_spearman_matrix}.csv
Output: 08_external_review/tables/E3_compartment_composition.csv,
  08_external_review/tables/E3_myc_vs_tf_sets.csv,
  08_external_review/tables/E3_supplementary_additions.csv

Two of the six additions are assembled here because the quantity is spread over several
files: the compartment composition (two annotation versions, two cohorts, compartment
and fine level) and the MYC correlations (one row per factor of each set). The other
four are existing files added unchanged. Table numbers do not change: each addition is
a further part of the sheet that already carries the subject.
"""
import os

import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
TAB = f"{W}/results/tables"
OUT = f"{W}/08_external_review/tables"
os.makedirs(OUT, exist_ok=True)

INFL5 = ["NFKB2", "RELB", "REL", "ATF1", "MTA1"]
TISS5 = ["SP7", "DLX5", "RUNX2", "SMAD3", "MYRF"]

# ------------------------------------------------- compartment composition (S8)
rows = []
for fn, cohort, ver, level in (
        ("03c_compartment_by_lesion_discovery.csv", "discovery_n5", "03c", "compartment"),
        ("03b_compartment_by_lesion_discovery.csv", "discovery_n5", "03b", "compartment"),
        ("05_compartment_by_lesion_validation.csv", "validation_n3", "03c", "compartment"),
        ("03c_celltype_fine_by_lesion_discovery.csv", "discovery_n5", "03c", "fine"),
        ("03b_celltype_fine_by_lesion_discovery.csv", "discovery_n5", "03b", "fine")):
    d = pd.read_csv(f"{TAB}/{fn}")
    if level == "compartment":
        labels = [c[2:] for c in d.columns if c.startswith("n_") and c != "n_total"]
        for _, r in d.iterrows():
            tot = int(r["n_total"])
            for lab in labels:
                rows.append(dict(cohort=cohort, lesion=r["lesion"],
                                 annotation_version=ver, level=level, label=lab,
                                 n_cells=int(r[f"n_{lab}"]), n_total=tot,
                                 pct_of_all_cells=round(100 * r[f"n_{lab}"] / tot, 4),
                                 source_file_of_part=fn))
    else:
        labels = [c for c in d.columns if c != "lesion"]
        for _, r in d.iterrows():
            tot = int(sum(r[c] for c in labels))
            for lab in labels:
                rows.append(dict(cohort=cohort, lesion=r["lesion"],
                                 annotation_version=ver, level=level, label=lab,
                                 n_cells=int(r[lab]), n_total=tot,
                                 pct_of_all_cells=round(100 * r[lab] / tot, 4),
                                 source_file_of_part=fn))
COMP = pd.DataFrame(rows)
COMP["scheme_adopted"] = "yes"
COMP["value_source"] = "recomputed from the file named in source_file_of_part"

# The scheme the Methods describe as evaluated and rejected: the compartment of a
# cluster taken from the majority of the cell-level marker argmax, without the
# cluster-level z step. It absorbs the pericytes into the stromal compartment, which is
# why the manuscript quotes it. Reproduced from the vote table, pooled over the
# discovery lesions (the table is per cluster, so no per-lesion split is possible).
VOTES = pd.read_csv(f"{TAB}/03c_compartment_votes_discovery.csv")
tot = int(VOTES.n_cells.sum())
by = VOTES.groupby("cell_majority").n_cells.sum()
not_adopted = []
for lab, n in by.items():
    not_adopted.append(dict(
        cohort="discovery_n5", lesion="POOLED_discovery",
        annotation_version="03c_cell_majority_only", level="compartment", label=lab,
        n_cells=int(n), n_total=tot, pct_of_all_cells=round(100 * n / tot, 4),
        source_file_of_part="03c_compartment_votes_discovery.csv",
        scheme_adopted="no",
        value_source="purely cell-level scheme, evaluated and rejected (Methods): the "
                     "cluster takes the majority of the cell-level marker argmax, with "
                     "no cluster-level z step; the pericytes are absorbed into the "
                     "stromal compartment"))
COMP = pd.concat([COMP, pd.DataFrame(not_adopted)], ignore_index=True)
COMP.to_csv(f"{OUT}/E3_compartment_composition.csv", index=False)

# ------------------------------------- projected state labels, with the earlier run (S5)
PROJ = pd.read_csv(f"{TAB}/25_validation_states_within_vs_ingest.csv")
PROJ["gene_space"] = "15,820 shared symbols after the symbol correction"
PROJ["scheme_adopted"] = "yes"
PROJ["value_source"] = "recomputed from 25_validation_states_within_vs_ingest.csv"
PROJ = pd.concat([PROJ, pd.DataFrame([dict(
    lesion="CAP1", labelling="ingest_projected", n_cells=860,
    pct_inflammatory=57.7,
    gene_space="shared gene space before the symbol correction",
    scheme_adopted="no",
    value_source="recorded in the Methods as the figure this fraction moved from; the "
                 "pre-correction gene space was not retained, so this row is not "
                 "recomputed from the current files")])], ignore_index=True)
PROJ.to_csv(f"{OUT}/E8_projection_with_not_adopted.csv", index=False)

# ------------------------------------------------------- MYC correlations (S10)
M = pd.read_csv(f"{TAB}/11_tf_activity_spearman_matrix.csv")
myc = []
for _, r in M.iterrows():
    if r.TF == "MYC":
        continue
    grp = ("inflammatory TF set" if r.TF in INFL5 else
           "tissue-forming TF set" if r.TF in TISS5 else "other")
    myc.append(dict(cohort=r.cohort, factor=r.TF, programme_set=grp,
                    spearman_rho_with_MYC_activity=round(float(r.MYC), 6)))
MYC = pd.DataFrame(myc)
MYC.to_csv(f"{OUT}/E3_myc_vs_tf_sets.csv", index=False)

# --------------------------------------------------------------- the list itself
ADD = [
    dict(sheet="S1", table="S1", part="gene-set and marker-panel coverage",
         source_file="20b_geneset_coverage_before_after.csv", source_dir="results/tables",
         quantity="coverage of the programme gene sets and of the compartment marker "
                  "panel in the quantified gene space, before and after the symbol "
                  "correction",
         manuscript_numbers="No.51 marker-panel coverage 100%; No.53 programme gene-set "
                            "coverage 90.6-97.8%",
         columns="gene_set, n_def, found_before, found_after, pct_before, pct_after, "
                 "recovered"),
    dict(sheet="S5", table="S5",
         part="projected state labels vs the within-cohort definition",
         source_file="E8_projection_with_not_adopted.csv",
         source_dir="08_external_review/tables",
         quantity="state proportions of the validation stromal cells under the "
                  "within-cohort definition and under the ingest projection, and the "
                  "agreement between them",
         manuscript_numbers="No.147 agreement 40.9%; No.148 agreement in CAP3 10.8%; "
                            "No.145 the 57.7% the CAP1 inflammatory fraction moved "
                            "from, as a row marked not adopted",
         columns="lesion, labelling, n_cells, n_inflammatory, pct_inflammatory, "
                 "n_intermediate, pct_intermediate, n_tissue-forming, "
                 "pct_tissue-forming, n_other, pct_other, agreement_pct, gene_space, "
                 "scheme_adopted, value_source"),
    dict(sheet="S8", table="S8",
         part="compartment composition per lesion, both annotation versions",
         source_file="E3_compartment_composition.csv",
         source_dir="08_external_review/tables",
         quantity="cells per compartment and per fine cell type, as a count and as a "
                  "percentage of all cells of the lesion, for the annotation version "
                  "used in the figures (03c) and for the superseded one (03b)",
         manuscript_numbers="No.222-223 validation immune share 69-87%; the evidence "
                            "for No.220-221 and No.224-225 (see E1); No.94 the 22.2% "
                            "stromal share of the rejected purely cell-level scheme",
         columns="cohort, lesion, annotation_version, level, label, n_cells, n_total, "
                 "pct_of_all_cells, source_file_of_part, scheme_adopted, "
                 "value_source"),
    dict(sheet="S10", table="S10", part="state proportions per lesion",
         source_file="33_state_pct_8lesions.csv", source_dir="results/tables",
         quantity="percentage of stromal cells in each of the four descriptive states, "
                  "per lesion, cohort-internal definition",
         manuscript_numbers="No.345-346 inflammatory 16.7-39.1%; No.347-348 "
                            "tissue-forming 12.3-24.9%; No.349-352 validation 1.7-16.4% "
                            "and 25.4-39.7%",
         columns="cohort, lesion, n_cells, inflammatory, intermediate, tissue-forming, "
                 "other"),
    dict(sheet="S10", table="S10", part="MYC activity against the two TF sets",
         source_file="E3_myc_vs_tf_sets.csv", source_dir="08_external_review/tables",
         quantity="Spearman correlation of MYC activity with each factor of the "
                  "inflammatory and tissue-forming TF sets, in stromal cells of each "
                  "cohort",
         manuscript_numbers="No.385-386 MYC +0.45 and -0.06 (per-factor values only; "
                            "the correlation with the mean activity of each set is not "
                            "in any saved table - listed in E4)",
         columns="cohort, factor, programme_set, spearman_rho_with_MYC_activity"),
    dict(sheet="S11", table="S11",
         part="permutation null for the retained dispersion",
         source_file="33_defA_permutation_null.csv", source_dir="results/tables",
         quantity="within-lesion permutation null of the retained dispersion: null "
                  "mean, 95% interval, one-sided p and the number of permutations",
         manuscript_numbers="No.330 p = 0.005; No.332 p = 0.025; No.334 p = 0.76; "
                            "No.331/448/450 200 permutations; No.447 95% interval",
         columns="cohort, lesion, n, frac_retained, null_mean, null_lo, null_hi, "
                 "p_permutation, n_perm, one_sided, seed"),
]
E3 = pd.DataFrame(ADD)
E3.to_csv(f"{OUT}/E3_supplementary_additions.csv", index=False)
print(f"composition {COMP.shape}  myc {MYC.shape}")
print(E3[["sheet", "part", "source_file", "manuscript_numbers"]].to_string(index=False))
