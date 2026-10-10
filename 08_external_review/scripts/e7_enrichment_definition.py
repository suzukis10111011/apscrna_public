#!/usr/bin/env python
"""Purpose: Record that the discovery and the validation enrichment ranges in the
  manuscript come from the same score definition and the same cell set, and why the
  per-lesion values in Table S5 differ.
Input: results/tables/{04_enrichment_by_lesion_discovery,05_enrichment_by_lesion_
  validation,18_stromal_markerdefined_vs_primary}.csv, scripts/{04_discovery_scores_and_
  screen,05_validation_and_projection,18_stromal_markerdefined}.py
Output: 08_external_review/tables/E7_enrichment_definition.csv

The enrichment is the stromal share of the top decile of the tissue-forming score
divided by the stromal share of all cells of the lesion. Three things can differ between
two such numbers: which score ranks the decile, which cells count as stromal, and which
cells form the denominator. Each row states all three, with the line of the script that
fixes the score.
"""
import os
import re
import subprocess

import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
TAB = f"{W}/results/tables"
OUT = f"{W}/08_external_review/tables"
os.makedirs(OUT, exist_ok=True)


def score_line(script, pattern):
    r = subprocess.run(["grep", "-nE", pattern, f"{W}/scripts/{script}"],
                       capture_output=True, text=True)
    first = r.stdout.splitlines()[0]
    no, text = first.split(":", 1)
    return f"scripts/{script}:{no}", text.strip()


ED = pd.read_csv(f"{TAB}/04_enrichment_by_lesion_discovery.csv")
EV = pd.read_csv(f"{TAB}/05_enrichment_by_lesion_validation.csv")
M18 = pd.read_csv(f"{TAB}/18_stromal_markerdefined_vs_primary.csv")

rows = []
for label, d, script, pattern, score_col, cellset in (
        ("discovery, manuscript range 3.4-9.5", ED[ED.score == "axis_D_tissueforming"],
         "04_discovery_scores_and_screen.py", r"^INFL, TISS =",
         "p_D_tissueforming", "primary cluster-based stromal set"),
        ("validation, manuscript range 8.8-10.0",
         EV[EV.score == "axis_D_tissueforming"], "05_validation_and_projection.py",
         r"^INFL, TISS =", "p_D_tissueforming", "primary cluster-based stromal set"),
        ("validation, Table S5 per lesion",
         M18[(M18.cohort == "validation_n3")
             & (M18.cellset == "primary_cluster_based") & (M18.lesion != "POOLED")],
         "18_stromal_markerdefined.py", r"^INF, TIS =",
         "lf_D_tissueforming", "primary cluster-based stromal set"),
        ("validation, Table S5 marker-defined set",
         M18[(M18.cohort == "validation_n3") & (M18.cellset == "marker_defined")
             & (M18.lesion != "POOLED")], "18_stromal_markerdefined.py",
         r"^INF, TIS =", "lf_D_tissueforming", "marker-defined stromal set"),
):
    loc, text = score_line(script, pattern)
    per = d.set_index("lesion")
    fold = per["enrichment_fold"]
    over = per["stromal_pct_overall"]
    top = per["stromal_pct_top_decile"]
    rows.append(dict(
        series=label,
        lesions=";".join(f"{k} {v:.3f}" for k, v in fold.items()),
        range_rounded=f"{fold.min():.1f}-{fold.max():.1f}",
        score_ranking_the_decile=score_col,
        score_fixed_at=f"{loc}  ({text})",
        stromal_set=cellset,
        denominator="all cells of the lesion",
        stromal_pct_overall=";".join(f"{k} {v:.2f}" for k, v in over.items()),
        stromal_pct_top_decile=";".join(f"{k} {v:.2f}" for k, v in top.items()),
        source_file=("04_enrichment_by_lesion_discovery.csv" if "discovery" in label
                     else "05_enrichment_by_lesion_validation.csv"
                     if "manuscript range" in label
                     else "18_stromal_markerdefined_vs_primary.csv")))
E7 = pd.DataFrame(rows)
E7.to_csv(f"{OUT}/E7_enrichment_definition.csv", index=False)
print(E7[["series", "range_rounded", "score_ranking_the_decile", "stromal_set",
          "stromal_pct_overall"]].to_string(index=False))
print("\nper lesion:")
print(E7[["series", "lesions", "stromal_pct_top_decile"]].to_string(index=False))
