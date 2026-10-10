#!/usr/bin/env python
"""Purpose: Prepare the two quantities the second external check could not reach - the
  per-lesion enrichment of the main analysis, and the unimodality of the programme
  scores - as parts to be added to Tables S5 and S3.
Input: results/tables/{04_enrichment_by_lesion_discovery,05_enrichment_by_lesion_
  validation,33_dip_by_lesion}.csv
Output: 08_external_review/tables/E9_primary_enrichment.csv,
  08_external_review/tables/E9_programme_score_dip.csv

The enrichment of the main analysis ranks the top decile on the full programme score
(p_D_tissueforming) over all cells of the lesion, and the stromal set is the
cluster-based one; this is the series behind Figure 2 and the ranges quoted in the
Abstract and the Results. The enrichment columns already in Table S5 rank the decile on
the leak-free score and belong to the cell-set sensitivity analysis, which is why their
per-lesion values differ.
"""
import os

import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
TAB = f"{W}/results/tables"
OUT = f"{W}/08_external_review/tables"
os.makedirs(OUT, exist_ok=True)

rows = []
for fn, cohort in (("04_enrichment_by_lesion_discovery.csv", "discovery_n5"),
                   ("05_enrichment_by_lesion_validation.csv", "validation_n3")):
    d = pd.read_csv(f"{TAB}/{fn}")
    for _, r in d[d.score == "axis_D_tissueforming"].iterrows():
        rows.append(dict(
            cohort=cohort, lesion=str(r.lesion),
            score_ranking_the_decile="p_D_tissueforming (full programme score)",
            stromal_set="cluster-based (the set used throughout the paper)",
            n_cells=int(r.n_cells), n_top_decile=int(r.n_top_decile),
            stromal_pct_overall_primary=float(r.stromal_pct_overall),
            stromal_pct_top_decile_primary=float(r.stromal_pct_top_decile),
            enrichment_fold_primary=float(r.enrichment_fold),
            figure="Figure 2A, B", source_file_of_part=fn))
PRIM = pd.DataFrame(rows)
PRIM.to_csv(f"{OUT}/E9_primary_enrichment.csv", index=False)

DIP = pd.read_csv(f"{TAB}/33_dip_by_lesion.csv")
DIP = DIP.rename(columns={"program": "programme_score",
                          "dip_statistic": "dip_statistic_programme_score",
                          "dip_p": "dip_p_programme_score"})
DIP["programme_score"] = DIP.programme_score.str.replace("program", "programme",
                                                         regex=False)
DIP["figure"] = "Figure 4E"
DIP["source_file_of_part"] = "33_dip_by_lesion.csv"
DIP.to_csv(f"{OUT}/E9_programme_score_dip.csv", index=False)

print(PRIM[["cohort", "lesion", "stromal_pct_overall_primary",
            "stromal_pct_top_decile_primary", "enrichment_fold_primary"]]
      .to_string(index=False))
per = PRIM[~PRIM.lesion.str.upper().str.startswith("POOL")]
print(f"\nper-lesion range: discovery "
      f"{per[per.cohort == 'discovery_n5'].enrichment_fold_primary.min():.3f}-"
      f"{per[per.cohort == 'discovery_n5'].enrichment_fold_primary.max():.3f}, "
      f"validation {per[per.cohort == 'validation_n3'].enrichment_fold_primary.min():.3f}-"
      f"{per[per.cohort == 'validation_n3'].enrichment_fold_primary.max():.3f}, "
      f"all eight {per.enrichment_fold_primary.min():.3f}-"
      f"{per.enrichment_fold_primary.max():.3f}")
print(f"dip rows {DIP.shape}, p range "
      f"{DIP.dip_p_programme_score.min():.4f}-{DIP.dip_p_programme_score.max():.4f}")
