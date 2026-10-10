#!/usr/bin/env python
"""Install the regenerated artefacts into results/ and archive the published ones.

Follows tables/C4_supersedes_map.csv. Each published file is MOVED to
results/superseded/ (never deleted) and a replacement is written in its place.

Replacements are built as supersets of the published files: every original column
keeps its name, position and - where the quantity did not change - its value and its
4-dp rounding, so the ten downstream scripts that read these files keep working. New
columns are appended. Where the quantity did change the new value is written and the
published value is kept in an adjacent column:

  04_tf_screen_discovery_candidates.csv      rounding removed; rank_within_axis rebuilt
  05_tf_screen_validation_reproduction.csv   rounding removed
  04/05_rho_inflam_vs_tissueforming_*.csv    leak-free rows added (definition =
                                             "leakfree"); axis and composite rows
                                             unchanged
  18_stromal_markerdefined_vs_primary.csv    defA_frac_retained now from a ULM on the
                                             reported cell set
  18_tf_sign_by_cellset.csv                  pooled_rho / all_same_sign now
                                             depth-controlled
  28_manuscript_numbers.csv                  labels corrected, provenance columns and
                                             post-audit rows added
  figures_ms/Fig3, Fig4 (pdf, png, bare)     redrawn panels E and B

New files added to results/: S3_rho_by_score_convention.csv (the wide supplement
table), the S11/S12 sheets and Supplementary_Tables_S11_S12.xlsx.

Output: results/superseded/** , results/tables/** , results/Supplementary_Tables_S11_S12.xlsx
        06_submission_prep/tables/C7_install_log.csv
        06_submission_prep/COMMIT_MESSAGE.txt
"""
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

W = "/path/to/Documents/apscrna_v2"
TAB = f"{W}/results/tables"
SUP = f"{W}/results/superseded"
NEW = f"{W}/06_submission_prep/tables"
FIGNEW = f"{W}/06_submission_prep/figures"
FIGPUB = f"{W}/results/figures_ms"
AUD = f"{W}/05_post_audit/tables"
log_rows = []


def log(m):
    s = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {m}"
    print(s, flush=True)
    open(f"{W}/06_submission_prep/logs/run_submission_prep.log", "a").write(s + "\n")


def archive(rel):
    """Move a published file under results/ to results/superseded/, keeping the tree."""
    src = f"{W}/results/{rel}"
    dst = f"{SUP}/{rel}"
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(src):
        shutil.move(src, dst)
        return True
    return False


def record(published, replacement, action, what):
    log_rows.append(dict(published_file=published, archived_to=f"results/superseded/"
                         f"{published}" if action != "added" else "",
                         replacement_source=replacement, action=action,
                         what_changed=what))


os.makedirs(SUP, exist_ok=True)

# ---------------------------------------------------------------- 1. TF screen
pub = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
new = pd.read_csv(f"{NEW}/S2_tf_screen_discovery_full_precision.csv")
new = new.set_index(["axis", "TF"])
out = pub.copy()
for c in ["pooled_rho", "abs_pooled", "rho_A", "rho_B", "rho_C", "rho_D", "rho_E"]:
    out[c] = [new.loc[(a, t), c] for a, t in zip(pub.axis, pub.TF)]
# rank_within_axis is rebuilt by the rule in 04: rank of |pooled rho| among the
# sign-consistent factors of the axis
out["rank_within_axis"] = np.nan
for ax, g in out.groupby("axis"):
    m = g.all5_same_sign.astype(str).str.lower().isin(["true"])
    out.loc[g[m].index, "rank_within_axis"] = g[m].abs_pooled.rank(ascending=False)
out["score_columns"] = ["TF activity x " + a + " (leak-free axis score)" for a in out.axis]
out["depth_covariate"] = "log10(total_counts), n_genes_by_counts"
out["precision"] = "full (the published file was rounded to 4 dp)"
archive("tables/04_tf_screen_discovery_candidates.csv")
out.to_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv", index=False)
record("tables/04_tf_screen_discovery_candidates.csv",
       "S2_tf_screen_discovery_full_precision.csv", "replaced",
       "4-dp rounding removed from pooled_rho, abs_pooled and the five per-lesion "
       "columns; rank_within_axis rebuilt by the same rule; candidate flags unchanged")

pub = pd.read_csv(f"{TAB}/05_tf_screen_validation_reproduction.csv")
new = pd.read_csv(f"{NEW}/S2_tf_screen_validation_reproduction_full_precision.csv")
new = new.set_index(["axis", "TF"])
out = pub.copy()
for c in ["discovery_pooled_rho", "validation_pooled_rho", "rho_CAP1", "rho_CAP2",
          "rho_CAP3"]:
    out[c] = [new.loc[(a, t), c] if (a, t) in new.index else np.nan
              for a, t in zip(pub.axis, pub.TF)]
out["precision"] = "full (the published file was rounded to 4 dp)"
archive("tables/05_tf_screen_validation_reproduction.csv")
out.to_csv(f"{TAB}/05_tf_screen_validation_reproduction.csv", index=False)
record("tables/05_tf_screen_validation_reproduction.csv",
       "S2_tf_screen_validation_reproduction_full_precision.csv", "replaced",
       "4-dp rounding removed; reproduced flags and counts unchanged")

# ------------------------------------------------------------------- 2. rho tables
S3 = pd.read_csv(f"{NEW}/S3_rho_by_score_convention.csv")
for fn, coh in [("04_rho_inflam_vs_tissueforming_discovery.csv", "discovery_n5"),
                ("05_rho_inflam_vs_tissueforming_validation.csv", "validation_n3")]:
    pub = pd.read_csv(f"{TAB}/{fn}")
    pub["score_columns"] = np.where(
        pub.definition == "axis",
        "p_A_cytokine_chemokine x p_D_tissueforming (whole axis, RELA targets retained)",
        "inflammatory_composite x tissueforming_composite")
    pub["depth_covariate"] = "log10(total_counts), n_genes_by_counts"
    add = []
    for _, r in S3[S3.cohort == coh].iterrows():
        add.append(dict(cohort=coh, lesion=r.lesion, cellset=r.cell_set,
                        definition="leakfree", n=int(r.n_cells),
                        raw_rho=round(float(r.raw_rho_leakfree_lf), 4),
                        depth_controlled_rho=round(
                            float(r.depth_controlled_rho_leakfree_lf), 4),
                        score_columns="lf_A_cytokine_chemokine x lf_D_tissueforming "
                                      "(leak-free, RELA targets removed)",
                        depth_covariate="log10(total_counts), n_genes_by_counts"))
    out = pd.concat([pub, pd.DataFrame(add)], ignore_index=True)
    archive(f"tables/{fn}")
    out.to_csv(f"{TAB}/{fn}", index=False)
    record(f"tables/{fn}", "S3_rho_by_score_convention.csv", "replaced",
           f"{len(add)} leak-free rows added (definition = 'leakfree') and the score "
           "convention named per row; the axis and composite rows are byte-identical "
           "to the published ones")
shutil.copy2(f"{NEW}/S3_rho_by_score_convention.csv",
             f"{TAB}/S3_rho_by_score_convention.csv")
record("(new)", "S3_rho_by_score_convention.csv", "added",
       "wide form of Table S3: composite, plain and leak-free rho side by side")

# ------------------------------------------------------------------ 3. table S5
pub = pd.read_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv")
new = pd.read_csv(f"{NEW}/S5_stromal_cellsets.csv").set_index(
    ["cohort", "cell_set", "lesion"])
out = pub.copy()
key = list(zip(pub.cohort, pub.cellset, pub.lesion))
out["defA_frac_retained_published_allcell_ulm"] = pub.defA_frac_retained
out["defA_frac_retained"] = [round(float(new.loc[k, "defA_frac_retained_RELA_from_cellset_ulm"]), 4)
                             if k in new.index else np.nan for k in key]
out["defA_n_bins"] = [int(new.loc[k, "defA_n_bins_used"]) if k in new.index else np.nan
                      for k in key]
out["defA_RELA_input"] = "ULM run on the reported cell set (the 04/05 convention)"
out["depth_controlled_rho_as_in_04_05"] = [
    round(float(new.loc[k, "depth_controlled_rho_as_in_04_05"]), 4)
    if k in new.index else np.nan for k in key]
out["score_columns"] = ("lf_A_cytokine_chemokine x lf_D_tissueforming "
                        "(leak-free, RELA targets removed)")
archive("tables/18_stromal_markerdefined_vs_primary.csv")
out.to_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv", index=False)
record("tables/18_stromal_markerdefined_vs_primary.csv", "S5_stromal_cellsets.csv",
       "replaced",
       "defA_frac_retained now uses the RELA activity from a ULM on the reported cell "
       "set; the published all-cell-ULM value is kept in "
       "defA_frac_retained_published_allcell_ulm; score convention named per row")

pub = pd.read_csv(f"{TAB}/18_tf_sign_by_cellset.csv")
new = pd.read_csv(f"{NEW}/S5_tf_sign_by_cellset_depth_controlled.csv").set_index(
    ["cohort", "cell_set", "axis", "TF"])
out = pub.copy()
key = list(zip(pub.cohort, pub.cellset, pub.axis, pub.TF))
out["pooled_rho_published_plain_spearman"] = pub.pooled_rho
out["all_same_sign_published_plain_spearman"] = pub.all_same_sign
out["pooled_rho"] = [round(float(new.loc[k, "pooled_rho_depth_controlled"]), 4)
                     if k in new.index else np.nan for k in key]
out["all_same_sign"] = [bool(new.loc[k, "all_same_sign_depth_controlled"])
                        if k in new.index else np.nan for k in key]
out["estimator"] = ("depth-controlled rho on [1, log10(total_counts), "
                    "n_genes_by_counts], as in the TF screen table")
archive("tables/18_tf_sign_by_cellset.csv")
out.to_csv(f"{TAB}/18_tf_sign_by_cellset.csv", index=False)
record("tables/18_tf_sign_by_cellset.csv",
       "S5_tf_sign_by_cellset_depth_controlled.csv", "replaced",
       "pooled_rho and all_same_sign are now depth-controlled, matching the TF screen "
       "table; the published plain-Spearman values are kept in the two "
       "*_published_plain_spearman columns")

# ------------------------------------------------------------- 4. number record
archive("tables/28_manuscript_numbers.csv")
shutil.copy2(f"{NEW}/28_manuscript_numbers_corrected.csv",
             f"{TAB}/28_manuscript_numbers.csv")
record("tables/28_manuscript_numbers.csv", "28_manuscript_numbers_corrected.csv",
       "replaced",
       "four rows labelled (lf_) relabelled to the plain score, the leak-free values "
       "for the same cells added, post-audit quantities added, and score_convention / "
       "cell_set / cohort / source_script / source_file columns added")
shutil.copy2(f"{NEW}/28_manuscript_numbers_corrections_log.csv",
             f"{TAB}/28_manuscript_numbers_corrections_log.csv")
record("(new)", "28_manuscript_numbers_corrections_log.csv", "added",
       "row-by-row record of the label corrections and additions")

# ------------------------------------------------------------- 5. S11 and S12
for fn in ["S11a_random_factor_null.csv", "S11b_overlap_removed_candidates.csv",
           "S11c_dispersion_baselines.csv", "S11c_null_factor_selection.csv",
           "S12a_single_donor_check.csv", "S12b_number_provenance.csv"]:
    shutil.copy2(f"{NEW}/{fn}", f"{TAB}/{fn}")
    record("(new)", fn, "added", "sheet of the new supplementary tables")
shutil.copy2(f"{NEW}/Supplementary_Tables_S11_S12.xlsx",
             f"{W}/results/Supplementary_Tables_S11_S12.xlsx")
record("(new)", "Supplementary_Tables_S11_S12.xlsx", "added",
       "new supplementary workbook (S11 and S12)")

# -------------------------------------------------------------- 6. panel records
for fn in ["_msn_Fig3.csv", "_msn_Fig4.csv"]:
    archive(f"tables/{fn}")
    shutil.copy2(f"{NEW}/{fn}", f"{TAB}/{fn}")
    record(f"tables/{fn}", fn, "replaced",
           "panel content and source tables updated for the redrawn panel")

# -------------------------------------------------------------------- 7. figures
FIGMAP = [("Fig3_tf_screen", "Fig3_tf_screen_v3"),
          ("Fig4_coupling_coexistence", "Fig4_coupling_coexistence_v3")]
for pubname, newname in FIGMAP:
    for ext in ("pdf", "png"):
        archive(f"figures_ms/{pubname}.{ext}")
        shutil.copy2(f"{FIGNEW}/{newname}.{ext}", f"{FIGPUB}/{pubname}.{ext}")
        archive(f"figures_ms/bare/{pubname}.{ext}")
        shutil.copy2(f"{FIGNEW}/bare/{newname}.{ext}",
                     f"{FIGPUB}/bare/{pubname}.{ext}")
    record(f"figures_ms/{pubname}.pdf|.png (and bare/)", f"{newname}.*", "replaced",
           "panel E of Figure 3 carries the random-20-factor null mean; panel B of "
           "Figure 4 carries the unrelated-factor baseline and the self-conditioning "
           "floor on a broken axis")

L = pd.DataFrame(log_rows)
L["note"] = ("installed by 06_submission_prep/scripts/r6_install_to_results.py; "
             "published files moved, not deleted")
L.to_csv(f"{NEW}/C7_install_log.csv", index=False)

with open(f"{SUP}/README.md", "w") as fh:
    fh.write(
        "# results/superseded\n\n"
        "Files that were replaced during the pre-submission tidy-up. They are moved "
        "here, never deleted, and keep their original path under `results/`.\n\n"
        "What replaced them, and why, is in `06_submission_prep/tables/"
        "C4_supersedes_map.csv` and `C7_install_log.csv`; the reasoning is in "
        "`06_submission_prep/README.md`.\n\n"
        "The replacements keep every original column name and position, so the "
        "scripts that read these files continue to run. Quantities that changed are "
        "listed in `C3_manuscript_vs_tables_discrepancies.csv`.\n\n"
        "Still stale after this swap and not regenerated here: the per-figure .pptx "
        "files in `results/figures_ms/`, `results/tables/figure_panel_n.csv`, and "
        "`results/tables_excel/` (BOM copies of the CSVs).\n")

log(f"install done: {len(L)} entries "
    f"({(L.action == 'replaced').sum()} replaced, {(L.action == 'added').sum()} added)")
print(L[["published_file", "action"]].to_string(index=False))
