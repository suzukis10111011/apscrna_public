#!/usr/bin/env python
"""Purpose: Adjudicate the seven items the external check flagged as needing a
  manuscript-side adjustment, by tracing each number to the file, cell set and score
  definition it came from, and list the manuscript locations to change.
Input: 07_submission_bundle/tables/manuscript_number_agreement.csv,
  07_submission_bundle/manuscript_text/manuscript_clean_v6.txt,
  05_post_audit/tables/{Q4b_number_provenance,Q2_definitionA_baselines}.csv,
  results/tables/** (rho, enrichment, composition, co-activity, state tables)
Output: 08_external_review/tables/E1_item_adjudication.csv,
  08_external_review/tables/E1_alternatives_checked.csv,
  08_external_review/tables/E4_manuscript_fix_list.csv

No interpretation of the biology is written here: each row states the manuscript
wording, the value it prints, the value the saved tables hold for that quantity, which
file and definition that value comes from, and which of the two findings applies -
the manuscript number is wrong, or the number is right and the table it was compared
against uses a different definition.
"""
import os
import re

import numpy as np
import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
TAB = f"{W}/results/tables"
AUD = f"{W}/05_post_audit/tables"
BUN = f"{W}/07_submission_bundle"
OUT = f"{W}/08_external_review/tables"
os.makedirs(OUT, exist_ok=True)

AG = pd.read_csv(f"{BUN}/tables/manuscript_number_agreement.csv")
A = AG[AG.version == "clean_v6"].reset_index(drop=True)
Q4B = pd.read_csv(f"{AUD}/Q4b_number_provenance.csv")


def ctx(no, words=10):
    """the manuscript wording around the number, as the worksheet shows it"""
    return str(A.iloc[no - 1].context)


def printed(no):
    return str(A.iloc[no - 1].printed)


def q4b(pattern):
    """provenance rows of the number-provenance map matching a quantity"""
    m = Q4B.apply(lambda r: r.astype(str).str.contains(pattern, case=False, na=False)
                  .any(), axis=1)
    return Q4B[m]


# the per-lesion files store the correlations rounded to 4 dp, which changes the second
# decimal for a value like -0.35503; the wide table carries them at full precision
S3F = pd.read_csv(f"{TAB}/S3_rho_by_score_convention.csv")
S3F = S3F[S3F.cell_set == "stromal"].set_index("lesion")
RV = pd.read_csv(f"{TAB}/05_rho_inflam_vs_tissueforming_validation.csv")
RD = pd.read_csv(f"{TAB}/04_rho_inflam_vs_tissueforming_discovery.csv")
PREV = pd.read_csv(f"{TAB}/05e_rho_prev_cellset_vs_clean.csv")
M18 = pd.read_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv")
EV = pd.read_csv(f"{TAB}/05_enrichment_by_lesion_validation.csv")
ED = pd.read_csv(f"{TAB}/04_enrichment_by_lesion_discovery.csv")
C3B = pd.read_csv(f"{TAB}/03b_compartment_by_lesion_discovery.csv")
C3C = pd.read_csv(f"{TAB}/03c_compartment_by_lesion_discovery.csv")
F3B = pd.read_csv(f"{TAB}/03b_celltype_fine_by_lesion_discovery.csv")
F3C = pd.read_csv(f"{TAB}/03c_celltype_fine_by_lesion_discovery.csv")
CVAL = pd.read_csv(f"{TAB}/05_compartment_by_lesion_validation.csv")
D32 = pd.read_csv(f"{TAB}/32_coactivity_ruleset_sensitivity.csv")
Q2 = pd.read_csv(f"{AUD}/Q2_definitionA_baselines.csv")

rows = []
alt = []


def add(no, quantity, ms_value, table_value, source, cellset, definition, finding,
        evidence, correct=""):
    rows.append(dict(No=no, quantity=quantity, manuscript_context=ctx(no),
                     manuscript_value=ms_value, table_value=table_value,
                     source_file=source, cell_set=cellset, score_definition=definition,
                     finding=finding, correct_value=correct, evidence=evidence))


def add_alt(no, candidate, value, why):
    alt.append(dict(No=no, alternative_checked=candidate, value=value, verdict=why))


# ---------------------------------------------------------------- items 1 and 2
# Per-lesion and pooled composite correlation in the validation cohort. The same
# sentence quotes the five discovery lesions and both pooled values; those agree with
# the table, so the comparison is between like quantities.
st = RV[(RV.cellset == "stromal") & (RV.definition == "composite")].set_index("lesion")
sd = RD[(RD.cellset == "stromal") & (RD.definition == "composite")].set_index("lesion")
for no, les in ((313, "CAP1"), (314, "CAP2"), (315, "CAP3")):
    tv = float(S3F.loc[les, "depth_controlled_rho_composite"])
    add(no, f"composite depth-controlled rho, {les} stromal cells", printed(no),
        f"{tv:.4f}", "results/tables/05_rho_inflam_vs_tissueforming_validation.csv "
        f"(= Table S3 row {les}/stromal, column depth_controlled_rho_composite)",
        "stromal (compartment assignment of this study)", "composite, depth-controlled",
        "manuscript value not reproducible from any saved table",
        "the five discovery values and both pooled values in the same sentence agree "
        "with the same file and definition "
        f"(discovery A {sd.loc['A','depth_controlled_rho']:.4f}, "
        f"pooled {sd.loc['POOLED_discovery','depth_controlled_rho']:.4f}; validation "
        f"pooled {st.loc['POOLED_validation','depth_controlled_rho']:.4f}), so the "
        "three per-lesion validation values are the only ones that do not",
        correct=f"{tv:.2f}")
    for nm, df, col, key in (
            ("all cells, composite", RV[(RV.cellset == "all_cells")
                                        & (RV.definition == "composite")],
             "depth_controlled_rho", les),
            ("stromal, single-axis (plain)", RV[(RV.cellset == "stromal")
                                                & (RV.definition == "axis")],
             "depth_controlled_rho", les),
            ("stromal, leak-free", RV[(RV.cellset == "stromal")
                                      & (RV.definition == "leakfree")],
             "depth_controlled_rho", les),
            ("marker-defined stromal set", M18[(M18.cohort == "validation_n3")
                                               & (M18.cellset == "marker_defined")],
             "depth_controlled_rho", les),
            ("pre-correction stromal set (before the plasma-cell exclusion)",
             PREV[PREV.cellset == "prev_stromal_set"], "depth_controlled_rho", les),
            ("raw (not depth-controlled), composite",
             RV[(RV.cellset == "stromal") & (RV.definition == "composite")],
             "raw_rho", les)):
        s = df.set_index("lesion")
        v = float(s.loc[key, col]) if key in s.index else float("nan")
        add_alt(no, nm, round(v, 4),
                "matches the manuscript" if round(v, 2) == round(
                    float(printed(no).replace("−", "-")), 2) else "does not match")

pool_raw = float(S3F.loc["POOLED_validation", "raw_rho_composite"])
allc_raw = float(RV[(RV.cellset == "all_cells") & (RV.definition == "composite")
                    & (RV.lesion == "POOLED_validation")].raw_rho.iloc[0])
disc_raw = float(S3F.loc["POOLED_discovery", "raw_rho_composite"])
add(321, "raw composite correlation, pooled validation stromal cells", printed(321),
    f"{pool_raw:.4f}", "results/tables/05_rho_inflam_vs_tissueforming_validation.csv "
    "(= Table S3 row POOLED_validation/stromal, column raw_rho_composite)",
    "stromal", "composite, raw (no depth control)",
    "manuscript value not reproducible from any saved table",
    f"the discovery figure in the same clause agrees ({disc_raw:.4f} -> -0.05); the "
    f"same magnitude with the opposite sign is the all-cells pooled raw value "
    f"({allc_raw:+.4f}), a different cell set",
    correct=f"{pool_raw:.2f}")
add_alt(321, "all cells, pooled, raw composite", round(allc_raw, 4),
        "same magnitude, opposite sign, different cell set")
add_alt(321, "pre-correction stromal set, pooled, raw composite",
        round(float(PREV[(PREV.cellset == "prev_stromal_set")
                         & (PREV.lesion == "POOLED_validation")].raw_rho.iloc[0]), 4),
        "does not match")

# ---------------------------------------------------------------- items 3 and 4
# Compartment composition of the discovery lesions. Two annotation versions exist:
# 03b, which leaves part of the cells unassigned, and 03c, which assigns every cell
# and is the version Figure 1E and Table S8 are built from.
imm3b = C3B.pct_immune
imm3c = C3C.pct_immune
pl3b = 100 * F3B.immune_plasma / F3B.drop(columns=["lesion"]).sum(axis=1)
pl3c = 100 * F3C.immune_plasma / F3C.drop(columns=["lesion"]).sum(axis=1)
for no, lo in ((220, True), (221, False)):
    add(no, "immune share of all cells, discovery lesions "
            f"({'lower' if lo else 'upper'} bound)", printed(no),
        f"{imm3c.min():.2f}-{imm3c.max():.2f}",
        "results/tables/03c_compartment_by_lesion_discovery.csv (the version behind "
        "Figure 1E and Table S8); the printed range is that of "
        "results/tables/03b_compartment_by_lesion_discovery.csv",
        "all cells per discovery lesion",
        "compartment annotation, 03c (every cell assigned) vs 03b (part unassigned)",
        "manuscript value comes from the superseded annotation version",
        f"03b gives {imm3b.min():.2f}-{imm3b.max():.2f} per lesion, which is the "
        f"printed 25-79%; 03c gives {imm3c.min():.2f}-{imm3c.max():.2f}; 03b leaves "
        f"{C3B.pct_unassigned.min():.2f}-{C3B.pct_unassigned.max():.2f}% of cells "
        "unassigned while 03c assigns all of them",
        correct=f"{imm3c.min():.1f}-{imm3c.max():.1f}%")
for no, lo in ((224, True), (225, False)):
    add(no, "plasma-cell share of all cells, discovery lesions "
            f"({'lower' if lo else 'upper'} bound)", printed(no),
        f"{pl3c.min():.2f}-{pl3c.max():.2f}",
        "results/tables/03c_celltype_fine_by_lesion_discovery.csv; the printed range "
        "is that of results/tables/03b_celltype_fine_by_lesion_discovery.csv",
        "all cells per discovery lesion",
        "fine cell-type annotation, 03c vs 03b",
        "manuscript value comes from the superseded annotation version",
        f"03b gives {pl3b.min():.2f}-{pl3b.max():.2f}% (the printed 14.5-31.6%), 03c "
        f"gives {pl3c.min():.2f}-{pl3c.max():.2f}%",
        correct=f"{pl3c.min():.1f}-{pl3c.max():.1f}%")
add(222, "immune share of all cells, validation lesions (lower bound)", printed(222),
    f"{CVAL.pct_immune.min():.2f}", "results/tables/05_compartment_by_lesion_"
    "validation.csv", "all cells per validation lesion",
    "compartment annotation, 03c", "manuscript value agrees with the table",
    f"per-lesion immune share {CVAL.pct_immune.min():.2f}-{CVAL.pct_immune.max():.2f}%"
    " matches the printed 69-87%; the table is not part of S1-S12 (added, see E3)",
    correct=printed(222))
add(223, "immune share of all cells, validation lesions (upper bound)", printed(223),
    f"{CVAL.pct_immune.max():.2f}", "results/tables/05_compartment_by_lesion_"
    "validation.csv", "all cells per validation lesion",
    "compartment annotation, 03c", "manuscript value agrees with the table",
    "same row set as No.222", correct=printed(223))

# ---------------------------------------------------------------- item 5
ax = EV[EV.score == "axis_D_tissueforming"]
co = EV[EV.score == "composite_tissueforming"]
axd = ED[ED.score == "axis_D_tissueforming"]
per = M18[(M18.cohort == "validation_n3") & (M18.cellset == "primary_cluster_based")
          & (M18.lesion != "POOLED")]
add(240, "stromal top-decile enrichment, validation lesions (lower bound)",
    printed(240), f"{ax.enrichment_fold.min():.3f}-{ax.enrichment_fold.max():.3f}",
    "results/tables/05_enrichment_by_lesion_validation.csv (score = "
    "axis_D_tissueforming)", "all cells per validation lesion, stromal share of the "
    "top decile", "single-axis tissue-forming score, as in the discovery range of the "
    "same sentence", "manuscript value not reproducible from any saved table",
    f"the discovery range in the same sentence ({axd.enrichment_fold.min():.3f}-"
    f"{axd.enrichment_fold.max():.3f} -> 3.4-9.5) is the single-axis score, so the "
    f"matching validation range is {ax.enrichment_fold.min():.3f}-"
    f"{ax.enrichment_fold.max():.3f}; the composite score gives "
    f"{co.enrichment_fold.min():.3f}-{co.enrichment_fold.max():.3f} and the "
    f"marker-defined set {per.enrichment_fold.min():.3f}-"
    f"{per.enrichment_fold.max():.3f}; the cell the check pointed at (9.332) is the "
    "pooled marker-defined value, not a per-lesion one",
    correct=f"{ax.enrichment_fold.min():.1f}-{ax.enrichment_fold.max():.1f}")

# ---------------------------------------------------------------- item 6
les = D32[(D32.group == "lesion") & (D32.label != "POOLED")]
mins = {}
for s in les.inflammatory_set.unique():
    sub = les[les.inflammatory_set == s]
    j = sub[["dip_p_infl", "dip_p_tiss"]].min(axis=1).idxmin()
    mins[s] = (float(sub[["dip_p_infl", "dip_p_tiss"]].min(axis=1).min()),
               f"{sub.loc[j, 'cohort']}/{sub.loc[j, 'label']}")
pooled_old = D32[(D32.inflammatory_set == "infl_old") & (D32.cohort == "validation")
                 & (D32.label == "POOLED")].dip_p_infl.iloc[0]
add(380, "Hartigan dip test of the two TF-set activities, per lesion", printed(380),
    f"{mins['infl_new_NFkB'][0]:.4f}", "results/tables/32_coactivity_ruleset_"
    "sensitivity.csv (= Table S10, part 'per-lesion and per-state statistics'), rows "
    "with group = lesion and label != POOLED",
    "stromal cells per lesion, 8 lesions", "mean activity of the rule-defined NF-kB "
    "set and of the tissue-forming set",
    "manuscript bound not reproducible for the rule-defined set",
    f"the per-lesion minimum is {mins['infl_new_NFkB'][0]:.4f} "
    f"({mins['infl_new_NFkB'][1]}) for the rule-defined NF-kB set and "
    f"{mins['infl_old'][0]:.4f} ({mins['infl_old'][1]}) for the alternative set with "
    f"ATF1 and MTA1; the printed 0.78 is the pooled validation value of the "
    f"alternative set ({pooled_old:.4f}), not a per-lesion value, and Figure 5C shows "
    "both sets",
    correct=f"{mins['infl_new_NFkB'][0]:.2f} (rule-defined set); "
            f"{mins['infl_old'][0]:.2f} (alternative set)")
for s, (v, where) in mins.items():
    add_alt(380, f"per-lesion minimum, {s}", round(v, 4), where)
add_alt(380, "pooled validation, infl_old", round(float(pooled_old), 4),
        "equals the printed 0.78 but is a pooled value")

# ---------------------------------------------------------------- item 7
rela = Q2[Q2.conditioning == "RELA_activity"]
vr_min, vr_max = rela.variance_retained.min(), rela.variance_retained.max()
red_max = rela.variance_reduction_pct.max()
add(328, "variance of the tissue-forming score retained after matching on RELA "
         "activity (upper bound)", printed(328),
    f"{100 * vr_min:.2f}-{100 * vr_max:.2f}",
    "05_post_audit/tables/Q2_definitionA_baselines.csv (= Table S11, part 'dispersion "
    "baselines'); the same quantity as Table S3", "stromal cells per lesion",
    "definition A, matched on RELA activity in 20 quantile bins",
    "manuscript value is the table value rounded; the stated bound is exceeded",
    f"the retained variance runs {100 * vr_min:.2f}% (lesion "
    f"{rela.loc[rela.variance_retained.idxmin(), 'lesion']}) to {100 * vr_max:.2f}% "
    f"(lesion {rela.loc[rela.variance_retained.idxmax(), 'lesion']}); 85-100% is that "
    f"range rounded, and the upper bound is above 100% by "
    f"{100 * vr_max - 100:.2f} percentage points",
    correct="85-101% (or 'about 100%')")
add(344, "variance of the tissue-forming score attributable to inflammatory activity "
         "(upper bound)", printed(344), f"{red_max:.4f}",
    "05_post_audit/tables/Q2_definitionA_baselines.csv (= Table S11, part 'dispersion "
    "baselines'), column variance_reduction_pct", "stromal cells per lesion",
    "definition A, matched on RELA activity",
    "manuscript value is the table value rounded; the stated bound is exceeded",
    f"the largest reduction is {red_max:.4f}% (lesion "
    f"{rela.loc[rela.variance_reduction_pct.idxmax(), 'lesion']}), so 'no more than "
    f"15%' is exceeded by {red_max - 15:.3f} percentage points; the Discussion "
    "sentence (No.414) already says 'about 15%'",
    correct="about 15% (or 'no more than 16%')")

E1 = pd.DataFrame(rows)
E1.to_csv(f"{OUT}/E1_item_adjudication.csv", index=False)
pd.DataFrame(alt).to_csv(f"{OUT}/E1_alternatives_checked.csv", index=False)

# ------------------------------------------------------------------------- E4
# the manuscript-side list: only the rows whose finding is that the printed number
# cannot be reproduced for the quantity the sentence states
FIX = E1[~E1.finding.str.startswith("manuscript value agrees with the table")]
E4 = FIX.assign(**{
    "manuscript sentence (10 words each side)": FIX.manuscript_context,
    "current value": FIX.manuscript_value,
    "correct value": FIX.correct_value,
    "evidence file": FIX.source_file,
})[["No", "quantity", "manuscript sentence (10 words each side)", "current value",
    "correct value", "evidence file", "finding"]]
E4.to_csv(f"{OUT}/E4_manuscript_fix_list.csv", index=False)
print(E1[["No", "quantity", "manuscript_value", "table_value", "finding",
          "correct_value"]].to_string(index=False))
print(f"\nE1 {E1.shape}  alternatives {len(alt)}  E4 {E4.shape}")
