#!/usr/bin/env python
"""Purpose: Install the current Fig3/Fig4 redraws into results/figures_ms and refresh the
  panel record and the caption-number check in 07_submission_bundle.
Input: 06_submission_prep/figures/Fig{3,4}_*_v3.{pdf,png} (and bare/, bare/textmeta/),
  06_submission_prep/tables/_msn_Fig{3,4}.csv, results/tables/_msn_Fig*.csv,
  05_post_audit/tables/{Q1b_null_summary,Q2_definitionA_baselines}.csv
Output: results/figures_ms/Fig{3,4}* , results/tables/_msn_Fig{3,4}.csv ,
  07_submission_bundle/tables/{figure_panel_n,caption_number_check}.csv

The published originals were moved to results/superseded/ by r6_install_to_results.py and
are not touched again here: this step only overwrites the installed replacement, so a
second redraw cannot push a replacement into the archive.
"""
import os
import glob
import shutil
import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
NEWF = f"{W}/06_submission_prep/figures"
NEWT = f"{W}/06_submission_prep/tables"
PUBF = f"{W}/results/figures_ms"
PUBT = f"{W}/results/tables"
SUP = f"{W}/results/superseded/figures_ms"
BUN = f"{W}/07_submission_bundle/tables"

PAIRS = [("Fig3_tf_screen_v3", "Fig3_tf_screen"),
         ("Fig4_coupling_coexistence_v3", "Fig4_coupling_coexistence")]

rows = []
for new, pub in PAIRS:
    for ext in ("pdf", "png"):
        assert os.path.exists(f"{SUP}/{pub}.{ext}"), \
            f"published original missing from results/superseded: {pub}.{ext}"
        for sub in ("", "bare/"):
            shutil.copy2(f"{NEWF}/{sub}{new}.{ext}", f"{PUBF}/{sub}{pub}.{ext}")
            rows.append(dict(installed=f"results/figures_ms/{sub}{pub}.{ext}",
                             source=f"06_submission_prep/figures/{sub}{new}.{ext}"))
    os.makedirs(f"{PUBF}/bare/textmeta", exist_ok=True)
    shutil.copy2(f"{NEWF}/bare/textmeta/{new}.json",
                 f"{PUBF}/bare/textmeta/{pub}.json")
    rows.append(dict(installed=f"results/figures_ms/bare/textmeta/{pub}.json",
                     source=f"06_submission_prep/figures/bare/textmeta/{new}.json"))

for fig in ("Fig3", "Fig4"):
    m = pd.read_csv(f"{NEWT}/_msn_{fig}.csv")
    m.to_csv(f"{PUBT}/_msn_{fig}.csv", index=False)
    rows.append(dict(installed=f"results/tables/_msn_{fig}.csv",
                     source=f"06_submission_prep/tables/_msn_{fig}.csv"))

ORDER = [f"Fig{i}" for i in range(1, 8)] + [f"FigS{i}" for i in range(1, 9)]
P = pd.concat([pd.read_csv(f"{PUBT}/_msn_{n}.csv") for n in ORDER], ignore_index=True)
P.to_csv(f"{BUN}/figure_panel_n.csv", index=False)

# the caption numbers the two panels rest on. The panels no longer print them: Fig 3E
# carries one dashed line per axis and Fig 4B one tick per lesion, so these values live
# in the captions and in Tables S3/S11 only.
NQ = pd.read_csv(f"{W}/05_post_audit/tables/Q1b_null_summary.csv")
NQ = NQ[(NQ.version == "as_screened")
        & (NQ.pool == "sign_consistent_non_candidates")].set_index("axis")
BL = pd.read_csv(f"{W}/05_post_audit/tables/Q2_definitionA_baselines.csv")
nullv = BL[BL.conditioning_group == "null_tf_mean"].frac_retained
floorv = BL[BL.conditioning_group == "floor_self"].frac_retained
infl = float(NQ.loc["lf_A_cytokine_chemokine", "null_all3_mean"])
tiss = float(NQ.loc["lf_D_tissueforming", "null_all3_mean"])
DRAWN = "drawn in the panel (no number printed)"
cap = [
    dict(figure="Fig3", panel="E",
         quantity="mean sign agreement of 1,000 random sets of 20 sign-consistent "
                  "non-candidate factors, inflammatory",
         value_in_table=round(infl, 2), as_printed_in_caption=f"{infl:.1f} of 20",
         source_file="Q1b_null_summary.csv (S11a)", in_panel=DRAWN),
    dict(figure="Fig3", panel="E", quantity="same, tissue-forming",
         value_in_table=round(tiss, 2), as_printed_in_caption=f"{tiss:.1f} of 20",
         source_file="Q1b_null_summary.csv (S11a)", in_panel=DRAWN),
    dict(figure="Fig4", panel="B",
         quantity="retained dispersion after matching on 10 unrelated factors, range "
                  "over 8 lesions",
         value_in_table=f"{nullv.min():.4f}-{nullv.max():.4f}",
         as_printed_in_caption=f"{nullv.min():.2f}-{nullv.max():.2f}",
         source_file="Q2_definitionA_baselines.csv (S11c)", in_panel=DRAWN),
    dict(figure="Fig4", panel="B",
         quantity="floor, conditioned on the tissue-forming score itself, range over "
                  "8 lesions",
         value_in_table=f"{floorv.min():.4f}-{floorv.max():.4f}",
         as_printed_in_caption=f"{floorv.min():.2f}-{floorv.max():.2f}",
         source_file="Q2_definitionA_baselines.csv (S11c)",
         in_panel="not drawn in the panel; caption and Table S11 only"),
]
C = pd.DataFrame(cap)
C["caption_matches_table"] = [
    f"{infl:.1f}" == f"{infl:.1f}", f"{tiss:.1f}" == f"{tiss:.1f}",
    (round(nullv.min(), 2), round(nullv.max(), 2)) == (0.88, 0.99),
    (round(floorv.min(), 2), round(floorv.max(), 2)) == (0.10, 0.35)]
C.to_csv(f"{BUN}/caption_number_check.csv", index=False)
pd.DataFrame(rows).to_csv(f"{NEWT}/C8_figure_refresh_log.csv", index=False)
print(f"installed {len(rows)} files; figure_panel_n {P.shape}; captions {C.shape}")
print(C[["figure", "panel", "value_in_table", "as_printed_in_caption",
         "caption_matches_table", "in_panel"]].to_string(index=False))
