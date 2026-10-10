# Inflammatory activity poorly predicts tissue-forming transcription in stromal cells of human apical periodontitis

Code, frozen inputs and result tables for the reanalysis of eight human apical
periodontitis lesions reported in the paper above. Everything needed to redo the
analysis from the public sequencing data is here; the processed objects are deposited
separately (see *Data*).

## Layout

| path | what it holds |
|---|---|
| `scripts/` | the analysis, figure and bundle scripts, numbered in running order |
| `frozen/` | the gene sets, the CollecTRI network and the per-cell covariates fixed before the analysis |
| `reference_records/` | the run-to-sample tables, the chemistry decision and the symbol map |
| `results/tables/` | every result table the scripts write |
| `results/figures_ms/` | the final figures (PDF and 600 dpi PNG) and the text-free versions used for the editable slides |
| `05_post_audit/` | the analyses added after the independent audit |
| `06_submission_prep/` | the final figure and table preparation |
| `08_external_review/` | the two external checks and the E-series records |
| `07_submission_bundle/` | the supplementary workbook, the figure slides and the checking record |
| `logs/` | the run logs of the quantification and the analysis steps |

## Running order

1. `scripts/01a_mkidx.sh`, `01b_validation_pipeline.sh`, `01c_fetch_sra_parallel.sh` -
   build the kallisto index and fetch the runs
2. `scripts/02_quant.sh`, `02b_quant.sh`, `08_discovery_requant.sh` - `kb count`
3. `scripts/03_preprocess_discovery.py`, `03b_fix_compartment.py`,
   `03c_recheck_discovery.py` - QC, doublet removal, compartment annotation
4. `scripts/03_preprocess_discovery.py` (Harmony) and
   `scripts/14c_robustness_integration.py` - integration and its sensitivity
5. `scripts/04_discovery_scores_and_screen.py` - programme scores
6. `scripts/_utils.py` (`ulm_on_symbols`) inside steps 5-7 - TF activities
7. `scripts/04_discovery_scores_and_screen.py` - the TF screen and the candidates
8. `scripts/05_validation_and_projection.py` - the validation cohort
9. `05_post_audit/scripts/q*.py` and `08_external_review/scripts/e*.py` - the analyses
   added after the audit and after the external checks
10. `06_submission_prep/scripts/r*.py` and `07_submission_bundle/scripts/s*.py` -
    the figures, the supplementary workbook and the number checks

## Input data

| cohort | accession |
|---|---|
| discovery lesions A-C | GEO GSE181688 |
| discovery lesions D, E | GEO GSE197680 |
| validation lesions CAP1-CAP3 | SRA PRJNA1213680 (SRR32119703, SRR32119704, SRR32119705) |
| independent gingival cohort | GEO GSE164241 |
| reference | Ensembl 112 GRCh38 cDNA; CollecTRI through OmniPath |

## Environment

Python 3.11.16; see `requirements.txt` and `environment.yml` for the package versions the
analysis was run with (anndata 0.12.19, decoupler 2.2.0, diptest 0.11.0, h5py 3.13.0, harmonypy 2.0.0, kb-python 0.30.2 and
the rest), and bustools 0.45.1, kallisto 0.52.0, sra-tools 3.4.1 for read processing.

## Where each supplementary table comes from

| table | source files | written by |
|---|---|---|
| S1 | 20b_geneset_coverage_before_after.csv; 28_qc_all_stages.csv; C6_base_composition_by_position.csv; C6_validation_chemistry_check.csv | 06_submission_prep/scripts/r5_validation_chemistry.py; scripts/28_manuscript_numbers.py |
| S2 | 04_tf_screen_discovery_candidates.csv; 05_tf_screen_validation_reproduction.csv; E21_figS8A_counts.csv; E24_geneset_members.csv; E25_validation_all_factors.csv | 08_external_review; 08_external_review/scripts/e24_new_parts.py; scripts/04_discovery_scores_and_screen.py; scripts/05_validation_and_projection.py |
| S3 | E9_programme_score_dip.csv; S3_rho_by_score_convention.csv | 06_submission_prep/scripts/r1_tables.py; 08_external_review/scripts/e9_references.py |
| S4 | 14_robustness_genesets_scoring.csv | scripts/14_robustness.py |
| S5 | 14_robustness_integration.csv; 18_stromal_markerdefined_vs_primary.csv; 18_tf_sign_by_cellset.csv; E15_sign_agreement_part.csv; E23_projection_confusion.csv; E8_projection_with_not_adopted.csv; E9_primary_enrichment.csv | 08_external_review; 08_external_review/scripts/e9_references.py; scripts/14c_robustness_integration.py; scripts/18_stromal_markerdefined.py |
| S6 | 09_geo_vs_kb.csv; 09_stromal_jaccard.csv; 09_tf_agreement_geo_vs_kb.csv | scripts/09_geo_vs_kb.py |
| S7 | 10_receptors_by_state.csv; 17_receptors_detection_fraction.csv; 17_receptors_mixedmodel.csv | scripts/10_receptors_by_state.py; scripts/17_receptor_power.py |
| S8 | 12_coupling_by_compartment.csv; 28_enrichment_permutation.csv; E3_compartment_composition.csv | 08_external_review/scripts/e3_table_additions.py; scripts/12_coupling_by_compartment.py; scripts/28_manuscript_numbers.py |
| S9 | 13_gse164241_summary.csv; 13_gse164241_tf_sign_reproduction.csv; E22_figS8B_counts.csv | 08_external_review; scripts/13_gse164241.py |
| S10 | 32_coactivity_ruleset_sensitivity.csv; 32_crossside_pairs_sign_consistency.csv; 32_old_vs_new_comparison.csv; 32_pc1_is_function_of_r.csv; 32_regulon_overlap.csv; 32_set_definition.csv; 33_state_pct_8lesions.csv; E26_tf_activity_matrix.csv; E6_myc_vs_set_means.csv | 08_external_review; 08_external_review/scripts/e6_myc_set_means.py |
| S11 | 33_defA_permutation_null.csv; S11a_random_factor_null.csv; S11b_overlap_removed_candidates.csv; S11c_dispersion_baselines.csv; S11c_null_factor_selection.csv | 06_submission_prep/scripts/r3_record_and_S11_S12.py |
| S12 | S12a_single_donor_check.csv; S12b_number_provenance.csv | 06_submission_prep/scripts/r3_record_and_S11_S12.py |

Full list, part by part, in `table_and_figure_provenance.csv`.

## Where each figure comes from

| figure | written by |
|---|---|
| Fig3_tf_screen (final version) | 06_submission_prep/scripts/r2_fig3.py |
| Fig4_coupling_coexistence (final version) | 06_submission_prep/scripts/r2_fig4.py |
| FigS4_matrix_vs_requantification (axis label corrected) | 08_external_review/scripts/e27_figS4_relabel.py |
| fig1_cohort_composition | scripts/fig1.py |
| fig2_tf_screen | scripts/fig2.py |
| fig3_definitionA | scripts/fig3.py |
| fig4_validation | scripts/fig4.py |
| fig5_receptors_coactivity | scripts/fig5.py |
| fig6_compartments | scripts/fig6.py |
| fig7_gse164241 | scripts/fig7.py |
| figS1_qc | scripts/figS1_qc.py |
| figS2_integration | scripts/figS2_integration.py |
| figS3_genesets | scripts/figS3_genesets.py |
| figS4_geo_vs_kb | scripts/figS4_geo_vs_kb.py |
| figS5_raw_vs_depth | scripts/figS5_raw_vs_depth.py |
| figS6_lolo | scripts/figS6_lolo.py |
| figS_robustness | scripts/figS_robustness.py |

## Data

The processed objects are deposited separately (4.29 GB, 8 files; see the
`manifest.csv` and `README.md` of the deposit). They are the AnnData objects and per-cell
tables the scripts read. The raw reads are not redistributed here.

## What is not in this repository

The manuscript files, the working notes kept in Japanese, the individual worksheet the
checker filled in, the raw and intermediate sequencing files and every object above 2 GB.
In the files that are here, the home directory of the machine the analysis ran on is
written as `/path/to`, and the person who checked the numbers is named by role.

The analysis was written and run in Japanese. For this release the comments and
docstrings of every `.py` and `.sh` file were translated into English, the code itself
left untouched (`translation_check.csv` records, file by file, that the syntax tree with
the docstrings removed is identical before and after); the speaker notes of the
PowerPoint figure files were removed, the slides themselves unchanged; and the Japanese
note columns of the record tables were translated. What is still in Japanese is the run
logs under `logs/` and the progress messages the scripts print as they run, which are
string literals in the code and were left alone so that the published code is the code
that produced the results. `release_scan.csv` lists, file by file, how much Japanese text
each file still contains.

## Licence

Code: MIT (`LICENSE`). Tables, figures and derived data: CC BY 4.0 (`LICENSE-DATA`).
The primary sequencing data stay under the terms of GEO and SRA.

## Citation

> The authors. *Inflammatory activity poorly predicts tissue-forming transcription in stromal cells of human apical periodontitis*. Zenodo, DOI: 10.5281/zenodo.XXXXXXX (to be filled in on
> deposition).

Please cite the paper as well once it is published.
