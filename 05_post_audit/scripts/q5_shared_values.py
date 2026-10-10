#!/usr/bin/env python
"""POST HOC (audit item 5).

The dental-pulp study reuses this repository's clean stromal values in its own
comparison table (pulp_scrna/05_manuscript_tables/tables/
T8b_independence_two_studies.csv, variant 'comparability-matched to this study').
Those rows are reproduced here as a supplementary table of the apical-lesion paper,
and each value is checked twice: against the pulp file, and against the apical-lesion
source table it was taken from.

  rho, validation           18_stromal_markerdefined_vs_primary.csv,
                            cellset = primary_cluster_based (leak-free scores,
                            log10(total_counts + 1))
  frac_retained, validation 05_definitionA_validation_v2.csv (RELA from the ULM on
                            the stromal subset)
  rho and frac_retained,    12_coupling_by_compartment.csv, compartment = stromal
  discovery

The three validation lesions (n = 860 / 566 / 353) are the rows the audit asks for;
the five discovery lesions and the discovery pooled row are carried along because the
pulp table holds them under the same variant.

Output: 05_post_audit/tables/Q5_shared_values_with_pulp_study.csv
        05_post_audit/tables/Q5_value_identity_check.csv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q0_common import AUD, POSTHOC_NOTE, TAB, log  # noqa: E402

PULP = ("/path/to/Documents/pulp_scrna/05_manuscript_tables/tables/"
        "T8b_independence_two_studies.csv")
P = pd.read_csv(PULP)
A = P[(P.study.str.startswith("apscrna"))
      & (P.variant.str.startswith("comparability"))].copy()

T18 = pd.read_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv")
D5 = pd.read_csv(f"{TAB}/05_definitionA_validation_v2.csv")
C12 = pd.read_csv(f"{TAB}/12_coupling_by_compartment.csv")

rows, checks = [], []
for _, r in A.iterrows():
    val = r.cohort == "validation_n3"
    if val:
        s18 = T18[(T18.cohort == "validation_n3")
                  & (T18.cellset == "primary_cluster_based")
                  & (T18.lesion == r["sample"])].iloc[0]
        s05 = D5[D5.lesion == r["sample"]].iloc[0]
        src_rho, src_fr = float(s18.depth_controlled_rho), float(s05.frac_retained)
        n_src = int(s18.n_cells)
        rho_file = "18_stromal_markerdefined_vs_primary.csv"
        fr_file = "05_definitionA_validation_v2.csv"
        fr_cond = "RELA activity, ULM on the stromal subset (04:150 route)"
    else:
        s12 = C12[(C12.compartment == "stromal") & (C12.lesion == r["sample"])].iloc[0]
        src_rho, src_fr = float(s12.depth_controlled_rho), float(s12.defA_frac_retained)
        n_src = int(s12.n_cells)
        rho_file = fr_file = "12_coupling_by_compartment.csv"
        fr_cond = "RELA activity, ULM on the all-cell object (12:44)"
    rows.append(dict(
        cohort=r.cohort, lesion=r["sample"], n_cells=n_src,
        depth_controlled_rho=src_rho, frac_retained=src_fr,
        frac_retained_conditioned_on=fr_cond, n_bins=r.n_bins,
        cell_set=r.cell_set_short, rho_source_file=rho_file,
        frac_retained_source_file=fr_file,
        score_definition="lf_A_cytokine_chemokine x lf_D_tissueforming "
                         "(RELA targets removed)",
        depth_covariate="log10(total_counts + 1), n_genes_by_counts",
        audit_request_scope=bool(val),
        shared_with="pulp_scrna/05_manuscript_tables/tables/"
                    "T8b_independence_two_studies.csv (variant "
                    "'comparability-matched to this study')",
        note="These are the apical-lesion values the dental-pulp study places beside "
             "its own; the convention was chosen for comparability with that study "
             "and is not a statement about which apical-lesion value is correct. "
             + POSTHOC_NOTE))
    for q, pulp_v, src_v in [("depth_controlled_rho", float(r.depth_controlled_rho),
                              src_rho),
                             ("frac_retained", float(r.frac_retained), src_fr),
                             ("n_cells", float(r.n_cells), float(n_src))]:
        checks.append(dict(cohort=r.cohort, lesion=r["sample"], quantity=q,
                           value_in_pulp_table=pulp_v,
                           value_in_apscrna_source=src_v,
                           abs_difference=abs(pulp_v - src_v),
                           identical=bool(pulp_v == src_v),
                           pulp_file=os.path.basename(PULP),
                           apscrna_file=(rho_file if q == "depth_controlled_rho"
                                         else fr_file),
                           note=POSTHOC_NOTE))

S = pd.DataFrame(rows).sort_values(["cohort", "lesion"])
S.to_csv(f"{AUD}/tables/Q5_shared_values_with_pulp_study.csv", index=False)
K = pd.DataFrame(checks)
K.to_csv(f"{AUD}/tables/Q5_value_identity_check.csv", index=False)
log(f"Q5: {len(S)} shared rows, identity checks {int(K.identical.sum())}/{len(K)} "
    f"identical, max |diff| {K.abs_difference.max():.2e}")
print(S[["cohort", "lesion", "n_cells", "depth_controlled_rho", "frac_retained",
         "audit_request_scope"]].to_string(index=False))
print(K[~K.identical].to_string(index=False) if (~K.identical).any()
      else "all values identical")
