#!/usr/bin/env python
"""Submission bundle: which file is which figure or table in the manuscript.

One row per deliverable file, with the manuscript number it carries, the panels it
contains and the tables the numbers come from. Figure rows take their panel list and
source tables from the regenerated panel record; supplementary-table rows take their
parts from the workbook index.

Output: 07_submission_bundle/tables/file_to_manuscript_map.csv
"""
import hashlib
import os

import pandas as pd

W = "/path/to/Documents/apscrna_v2"
B = f"{W}/07_submission_bundle"
FIG = f"{W}/results/figures_ms"

FIGURES = [("Fig1", "Figure 1", "Fig1_design_composition"),
           ("Fig2", "Figure 2", "Fig2_enrichment"),
           ("Fig3", "Figure 3", "Fig3_tf_screen"),
           ("Fig4", "Figure 4", "Fig4_coupling_coexistence"),
           ("Fig5", "Figure 5", "Fig5_receptors_coactivity"),
           ("Fig6", "Figure 6", "Fig6_compartments"),
           ("Fig7", "Figure 7", "Fig7_independent_cohort"),
           ("FigS1", "Figure S1", "FigS1_qc"),
           ("FigS2", "Figure S2", "FigS2_integration"),
           ("FigS3", "Figure S3", "FigS3_genesets_scoring"),
           ("FigS4", "Figure S4", "FigS4_matrix_vs_requantification"),
           ("FigS5", "Figure S5", "FigS5_raw_vs_depth"),
           ("FigS6", "Figure S6", "FigS6_leave_one_out"),
           ("FigS7", "Figure S7", "FigS7_state_labelling"),
           ("FigS8", "Figure S8", "FigS8_other_programs")]
CHANGED = {"Fig3": "panel E now carries the random-20-factor null mean",
           "Fig4": "panel B now carries the unrelated-factor baseline and the "
                   "self-conditioning floor on a broken axis"}

P = pd.read_csv(f"{B}/tables/figure_panel_n.csv")
IDX = pd.read_csv(f"{B}/tables/supplementary_table_index.csv")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16] if os.path.exists(p) else ""


rows = []
for key, number, stem in FIGURES:
    g = P[P.figure == key]
    for kind, path in [("print figure (vector)", f"{FIG}/{stem}.pdf"),
                       ("print figure (600 dpi)", f"{FIG}/{stem}.png"),
                       ("text-free base (vector)", f"{FIG}/bare/{stem}.pdf"),
                       ("text-free base (600 dpi)", f"{FIG}/bare/{stem}.png"),
                       ("editable layout", f"{B}/figures/pptx/{key}.pptx")]:
        rows.append(dict(
            manuscript_item=number, item_type="figure", file_kind=kind,
            file=os.path.relpath(path, W), exists=os.path.exists(path),
            sha256_16=sha(path),
            panels="".join(sorted(g.panel.astype(str).unique())),
            n_panels=len(g),
            source_tables="; ".join(sorted({t for s in g.source_table.dropna()
                                            for t in str(s).split("; ")})),
            changed_in_this_round=CHANGED.get(key, ""),
            note="figures read from results/figures_ms; the pptx is rebuilt from the "
                 "text-free base and its text-position record"))
rows.append(dict(manuscript_item="Figures 1-S8 combined", item_type="figure",
                 file_kind="editable layout, all figures",
                 file=os.path.relpath(f"{B}/figures/figures_editable.pptx", W),
                 exists=os.path.exists(f"{B}/figures/figures_editable.pptx"),
                 sha256_16=sha(f"{B}/figures/figures_editable.pptx"), panels="",
                 n_panels=len(P), source_tables="", changed_in_this_round="",
                 note="one slide per figure"))

for t, g in IDX.groupby("table", sort=False):
    num = "Table " + t
    sheets = list(dict.fromkeys(g.sheet))
    rows.append(dict(
        manuscript_item=num, item_type="supplementary table",
        file_kind="sheet " + " and ".join(f"'{x}'" for x in sheets) + " of the workbook",
        file=os.path.relpath(f"{B}/tables/Supplementary_Tables_S1-S12.xlsx", W),
        exists=True, sha256_16=sha(f"{B}/tables/Supplementary_Tables_S1-S12.xlsx"),
        panels="", n_panels=int(g.n_rows.sum()),
        source_tables="; ".join(g.source_file),
        changed_in_this_round=({"S1": "second sheet added: the chemistry check of the "
                                      "three validation runs",
                                "S2": "full precision, all six axes",
                                "S3": "column names state the score convention; "
                                      "leak-free columns added",
                                "S5": "leak-free stated per row; definition A on one "
                                      "RELA input; sign counts depth-controlled",
                                "S11": "new", "S12": "new"}).get(t, ""),
        note=g.description.iloc[0]))

M = pd.DataFrame(rows)
M.to_csv(f"{B}/tables/file_to_manuscript_map.csv", index=False)
print(M.groupby("item_type").agg(files=("file", "nunique"), rows=("file", "size")).to_string())
print(M[M.item_type == "supplementary table"][["manuscript_item", "n_panels",
                                               "changed_in_this_round"]].to_string(index=False))
missing = M[~M.exists]
print("missing files:", missing.file.tolist() if len(missing) else "none")
