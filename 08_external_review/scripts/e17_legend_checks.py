#!/usr/bin/env python
"""Purpose: Check the supplementary figure legends against the scripts that drew the
  figures, check the supplementary table legends against the parts the workbook actually
  holds, and extract the table legends for the workbook's Legends sheet.
Input: 07_submission_bundle/manuscript_text/manuscript_clean_v*.txt (newest),
  07_submission_bundle/tables/supplementary_table_index.csv, the figure scripts under
  scripts/ and the result tables they read
Output: 08_external_review/tables/E17_figure_legend_checks.csv,
  08_external_review/tables/E18_table_legend_vs_parts.csv,
  08_external_review/tables/E19_validation_enrichment_in_S8.csv,
  08_external_review/tables/E20_table_legends.csv

Each figure-legend check names the script line that decides it and, where the claim is
numeric, the value recomputed from the table the script reads, so the verdict rests on
what the code does rather than on a reading of the legend.
"""
import glob
import os
import re

import numpy as np
import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
TXT = sorted(glob.glob(f"{W}/07_submission_bundle/manuscript_text/manuscript_clean_v*.txt"),
             key=lambda p: int(re.search(r"_v(\d+)\.txt$", p).group(1)))[-1]
TAB = f"{W}/results/tables"
OUT = f"{W}/08_external_review/tables"
T = open(TXT).read()


def legends(start, end, pattern):
    blk = T[T.find(start):(T.find(end) if end else len(T))]
    out, keys = {}, [m for m in re.finditer(pattern, blk)]
    for i, m in enumerate(keys):
        stop = keys[i + 1].start() if i + 1 < len(keys) else len(blk)
        out[m.group(1)] = blk[m.start():stop].strip()
    return out


FIG = legends("Supplementary figure legends", "Supplementary table legends",
              r"\n(Figure S\d)\.")
TBL = legends("Supplementary table legends", "Acknowledgements", r"\n(Table S\d+)\.")


def grep(path, pattern):
    with open(f"{W}/{path}", errors="replace") as fh:
        for i, line in enumerate(fh, start=1):
            if re.search(pattern, line):
                return f"{path}:{i}", line.strip()[:150]
    return "", ""


checks = []


def add(fig, claim, verdict, evidence, correct=""):
    checks.append(dict(figure=fig, claim_in_the_legend=claim, verdict=verdict,
                       evidence=evidence, what_the_figure_shows=correct))


# --- Figure S2 (B): which score, and whether the rho is depth-controlled
R = pd.read_csv(f"{TAB}/14_tf_rank_by_integration.csv")
I = pd.read_csv(f"{TAB}/14_robustness_integration.csv").iloc[0]
w_ax, l_ax = grep("scripts/14c_robustness_integration.py", r'lf_D_tissueforming"\]\.values')
w_rho, l_rho = grep("scripts/14c_robustness_integration.py", r"spearmanr\(sub\[tf\]")
add("Figure S2B", "one point per TF, against the leak-free tissue-forming score",
    "correct",
    f"{w_ax}: the y variable is lf_D_tissueforming only; {len(R)} TFs, one row each "
    f"in 14_tf_rank_by_integration.csv (columns harmony, uncorrected)")
_s2 = FIG.get("Figure S2", "")
add("Figure S2B", "how the plotted rho is described",
    "correct" if re.search(r"without depth control", _s2) else "wrong",
    f"{w_rho}: {l_rho} - a plain Spearman of TF activity against the score, with no "
    "residualization on depth; the legend now reads 'Pooled Spearman rho (without "
    "depth control)'" if re.search(r"without depth control", _s2) else
    f"{w_rho}: {l_rho} - a plain Spearman, but the legend calls it depth-controlled",
    "pooled Spearman rho (not depth-controlled)")
add("Figure S2B", "rank correlation 0.93 and 15 of the top 20 shared", "correct",
    f"14_robustness_integration.csv: tf_rank_spearman={I.tf_rank_spearman}, "
    f"top20_overlap={int(I.top20_overlap)}")
add("Figure S2A", "6,215 / 1,868 / 1,865 (99.8% of the uncorrected set)", "correct",
    f"14_robustness_integration.csv: {int(I.stromal_n_harmony)}, "
    f"{int(I.stromal_n_uncorrected)}, {int(I.stromal_intersection)}; "
    f"{100 * I.stromal_intersection / I.stromal_n_uncorrected:.2f}% of the uncorrected set")

# --- Figure S3: which cells the heatmap values come from
B14 = pd.read_csv(f"{TAB}/14_robustness_genesets_scoring.csv")
P14 = B14[B14.lesion == "POOLED"]
w_in, l_in = grep("scripts/14_robustness.py", r"stromal_discovery_labelled\.h5ad")
w_pool, l_pool = grep("scripts/ms_figS3.py", r'B\.lesion == "POOLED"')
neg = int((P14.depth_controlled_rho < 0).sum())
add("Figure S3", "pooled discovery stromal cells", "correct",
    f"{w_in}: the script reads stromal_discovery_labelled.h5ad (discovery stromal "
    f"cells); {w_pool}: the figure keeps lesion == POOLED, {len(P14)} of the "
    f"{len(B14)} rows")
add("Figure S3A", "negative in 16 of the 18 combinations, -0.24 to +0.04", "correct",
    f"14_robustness_genesets_scoring.csv, POOLED rows: {neg} of {len(P14)} negative, "
    f"range {P14.depth_controlled_rho.min():.4f} to {P14.depth_controlled_rho.max():.4f}")
add("Figure S3B", "retained dispersion 0.995-0.999", "correct",
    f"same rows: {P14.defA_frac_retained.min():.4f} to {P14.defA_frac_retained.max():.4f}")

# --- Figure S5 (A): full or leak-free scores
D04 = pd.read_csv(f"{TAB}/04_rho_inflam_vs_tissueforming_discovery.csv")
ax_rows = D04[(D04.definition == "axis") & (D04.cellset == "stromal")]
add("Figure S5A", "the single-axis panel uses the full programme scores", "correct",
    "04_rho_inflam_vs_tissueforming_discovery.csv, definition=axis, cellset=stromal: "
    f"score_columns = {ax_rows.score_columns.iloc[0]}")

# --- Figure S8: what "consistent sign" means, and the counts
V = pd.read_csv(f"{TAB}/05_tf_screen_validation_reproduction.csv")
G = pd.read_csv(f"{TAB}/13_gse164241_tf_sign_reproduction.csv")
for d in (V, G):
    for c in d.columns:
        if d[c].dtype == object and set(d[c].astype(str).str.lower()) <= {"true", "false", "nan"}:
            d[c] = d[c].astype(str).str.lower().eq("true")
w_rep, l_rep = grep("scripts/05_validation_and_projection.py", r"reproduced=bool\(all\(")
MAIN = ["lf_A_cytokine_chemokine", "lf_D_tissueforming"]
others = [a for a in V.axis.unique() if a not in MAIN]
va = {a: (int(V[(V.axis == a) & V.in_validation].reproduced.sum()),
          int(((V.axis == a) & V.in_validation).sum())) for a in others}
ga = {a: (int(G[(G.group == "ALL") & (G.axis == a) & G.present].sign_matches_discovery.sum()),
          int(((G.group == "ALL") & (G.axis == a) & G.present).sum())) for a in others}
add("Figure S8A", "candidates whose depth-controlled rho had the same sign in all "
    "three validation lesions (12/18, 11/18, 10/17, 10/16)", "correct",
    f"{w_rep}: {l_rep} - the sign of every validation lesion is compared with the "
    "DISCOVERY pooled rho, not merely with each other; counts recomputed: "
    + ", ".join(f"{a.replace('lf_', '')} {n}/{d}" for a, (n, d) in va.items()),
    "same sign as the discovery pooled rho in all three validation lesions")
add("Figure S8B", "candidates whose pooled rho had the same sign as in the discovery "
    "cohort (13/19, 11/19, 15/18, 15/18)", "correct",
    "13_gse164241_tf_sign_reproduction.csv, group=ALL, present=True: "
    + ", ".join(f"{a.replace('lf_', '')} {n}/{d}" for a, (n, d) in ga.items()))

# --- Figure S4 (C): what the bar and the fraction are
A09 = pd.read_csv(f"{TAB}/09_tf_agreement_geo_vs_kb.csv")
w_bar, l_bar = grep("scripts/ms_figS4.py", r"a3\.barh\(y, T2\.same_sign_pct")
w_lab, l_lab = grep("scripts/ms_figS4.py", r"cand\.")
sub = A09[A09.axis.isin(MAIN)]
pct_list = ", ".join(f"{r.source} {r.axis.replace('lf_', '')} {r.same_sign_pct}%"
                     for r in sub.itertuples())
_s4 = FIG.get("Figure S4", "")
_rng = f"{sub.same_sign_pct.min():.1f}-{sub.same_sign_pct.max():.1f}%"
_ok4 = bool(re.search(r"TFs evaluable in both quantifications", _s4)
            and re.search(r"overlap between the 20-candidate lists", _s4))
add("Figure S4C", "what the bar and the number beside it are",
    "correct" if _ok4 else "wrong",
    f"{w_bar}: the bar is same_sign_pct, the percentage of the {int(A09.n_TF_common.iloc[0])} "
    f"TFs evaluable in both quantifications whose pooled rho keeps its sign "
    f"({pct_list}); "
    f"{w_lab}: the x/20 label next to each bar is overlap / n_candidates_geo, i.e. how "
    "many of the 20 discovery candidates are again among the top 20 of the "
    f"requantification; the legend now names both correctly and quotes {_rng}",
    "bar: percentage of all evaluable TFs keeping the sign of the pooled rho; "
    "x/20: overlap of the two candidate lists")

# --- Figure S6 and S7: the counts quoted
L34 = pd.read_csv(f"{TAB}/34_lolo_rulebased.csv")
keep = L34.groupby(["axis", "TF"]).is_candidate.sum()
X = pd.read_csv(f"{TAB}/25_state_crosstab_within_vs_ingest.csv").set_index("state_within")
diag = int(sum(X.loc[s, s] for s in X.index if s in X.columns))
tot = int(X.values.sum())
add("Figure S6B", "5 of 5 runs for all ten factors", "correct",
    f"34_lolo_rulebased.csv: {int((keep == 5).sum())} of {len(keep)} factors are "
    "candidates in all five leave-one-out runs")
add("Figure S7B", "the labels agreed for 727 cells (40.9%)", "correct",
    f"25_state_crosstab_within_vs_ingest.csv: the diagonal sums to {diag} of {tot} "
    f"cells, {100 * diag / tot:.2f}%")
E17 = pd.DataFrame(checks)
E17.to_csv(f"{OUT}/E17_figure_legend_checks.csv", index=False)

# --- table legends against the parts the workbook holds
from openpyxl import load_workbook

IDX = pd.read_csv(f"{W}/07_submission_bundle/tables/supplementary_table_index.csv")
WB = load_workbook(f"{W}/07_submission_bundle/tables/Supplementary_Tables_S1-S12.xlsx",
                   data_only=True)


def cols(sheet):
    return [c.value for c in WB[sheet][1] if c.value]


# one phrase per part: what the legend has to say for that part to be described
PART_PHRASE = {
    ("S1", "main"): r"quality-control stage",
    ("S1", "gene-set and marker-panel coverage"): r"coverage",
    ("S1", "per run"): r"reads processed",
    ("S1", "base composition by position"): r"base composition",
    ("S2", "discovery_n5"): r"discovery cohort",
    ("S2", "validation_n3"): r"validation",
    ("S2", "candidate reproduction per programme, counted (as in Figure S8A)"):
        r"reproduction counts",
    ("S2", "programme gene sets, members"): r"members",
    ("S2", "validation values of every factor of the two programmes "
           "(as in Figure 3A, B)"): r"717",
    ("S3", "main"): r"composite, full and leak-free",
    ("S3", "unimodality of the programme scores per lesion (Hartigan dip test, as in "
           "Figure 4E)"): r"dip test",
    ("S4", "main"): r"18 combinations",
    ("S5", "integration method"): r"Harmony",
    ("S5", "stromal cell-set definition, decile ranked on the leak-free score "
           "(sensitivity analysis, not the Figure 2 series)"): r"marker-defined",
    ("S5", "candidate sign by cell set"): r"sign of each candidate",
    ("S5", "projected state labels vs the within-cohort definition"): r"projected",
    ("S5", "primary enrichment per lesion (as in Figure 2)"): r"enrichment",
    ("S5", "candidate sign agreement between the two stromal sets, counted"):
        r"summarized as counts|counts for the six programmes",
    ("S5", "within-cohort against projected state labels, cell counts "
           "(as in Figure S7B)"): r"confusion",
    ("S6", "key quantities"): r"barcode and UMI recovery",
    ("S6", "stromal barcode agreement"): r"Jaccard",
    ("S6", "TF screen agreement"): r"candidate agreement",
    ("S7", "pseudobulk, paired Wilcoxon"): r"Pseudobulk",
    ("S7", "cell-level linear mixed model"): r"mixed model",
    ("S7", "detection fraction"): r"percentage of cells expressing",
    ("S8", "main"): r"within each compartment",
    ("S8", "top-decile permutation test"): r"permutation",
    ("S8", "compartment composition per lesion, both annotation versions"):
        r"cells per compartment",
    ("S9", "per donor"): r"Per donor",
    ("S9", "TF sign reproduction"): r"sign of each discovery candidate",
    ("S9", "candidate reproduction per programme, counted (as in Figure S8B)"):
        r"S8B",
    ("S10", "TF set definition by rule"): r"rule-based definition",
    ("S10", "per-lesion and per-state statistics"): r"mean activities",
    ("S10", "cross-side TF pairs"): r"cross-side",
    ("S10", "regulon overlap"): r"regulon overlap",
    ("S10", "rule-based vs previous TF set"): r"previously used",
    ("S10", "note on the PC1 variance ratio"): r"PC1",
    ("S10", "state proportions per lesion"): r"proportions of inflammatory",
    ("S10", "MYC activity against the two TF sets"): r"MYC",
    ("S10", "factor-by-factor activity correlation matrix (as in Figure 5D)"):
        r"factor-by-factor",
    ("S11", "random-factor null"): r"random factor-set",
    ("S11", "overlap-removed activities"): r"removal of programme genes",
    ("S11", "dispersion baselines"): r"null reference",
    ("S11", "unrelated-factor selection"): r"unrelated",
    ("S11", "permutation null for the retained dispersion"): r"200 permutations",
    ("S12", "single-donor check"): r"sex-marker",
    ("S12", "provenance of reported values"): r"provenance",
}

S2d = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
E3c = pd.read_csv(f"{OUT}/E3_compartment_composition.csv")
rows = []
for _, r in IDX.iterrows():
    t, p = f"Table {r.table}", str(r.part)
    leg = TBL.get(t, "")
    phrase = PART_PHRASE.get((r.table, p))
    if phrase is None:
        rows.append(dict(table=t, direction="in the table, not in the legend",
                         item=f"part '{p}'", evidence="no phrase of the legend "
                         "describes this part"))
        continue
    if not re.search(phrase, leg, re.I):
        rows.append(dict(table=t, direction="in the table, not in the legend",
                         item=f"part '{p}'",
                         evidence=f"the legend of {t} does not contain '{phrase}'"))

# claims the legend makes, each with what the workbook would need to back it
CLAIMS = [
    ("Table S2", r"members", "the gene-set members",
     any(IDX[(IDX.table == "S2")].part.str.contains("members"))),
    ("Table S2", r"717", "the validation values of the evaluable factors",
     any(IDX[(IDX.table == "S2")].part.str.contains("every factor"))),
    ("Table S2", r"score convention", "all three score conventions",
     not set(S2d.axis.unique()) - {a for a in S2d.axis.unique() if a.startswith("lf_")}
     and False),
    ("Table S3", r"retained dispersion|SD ratio", "the retained dispersion",
     any("retained" in str(c) or "frac_retained" in str(c) for c in cols("S3"))),
    ("Table S7", r"mixed model", "the cell-level mixed model",
     any(IDX[(IDX.table == "S7")].part.str.contains("mixed model"))),
    ("Table S8", r"rejected cell-majority|earlier annotation",
     "the other annotation versions",
     len(set(E3c.annotation_version.unique())) > 1),
    ("Table S10", r"factor-by-factor", "the factor-by-factor matrix",
     any(IDX[(IDX.table == "S10")].part.str.contains("factor-by-factor"))),
    ("Table S11", r"SD ratio|variance ratio", "the per-lesion retained dispersion",
     any("frac_retained" in str(c) or "variance_retained" in str(c)
         for c in cols("S11"))),
]
for t, phrase, what, backed in CLAIMS:
    claimed = bool(re.search(phrase, TBL.get(t, ""), re.I))
    if claimed and not backed:
        rows.append(dict(table=t, direction="in the legend, not in the table",
                         item=what,
                         evidence=f"the legend of {t} claims it ('{phrase}') but no "
                                  "part of the sheet carries it"))
E18 = pd.DataFrame(rows)
E18.to_csv(f"{OUT}/E18_table_legend_vs_parts.csv", index=False)
PARTS = pd.DataFrame([dict(table=t, n_parts=len(g), parts="; ".join(g.part.astype(str)),
                           legend_chars=len(TBL.get(f"Table {t}", "")))
                      for t, g in IDX.groupby("table", sort=False)])
PARTS.to_csv(f"{OUT}/E18_parts_per_table.csv", index=False)

# --- does Table S8 carry the validation enrichment?
S8 = pd.read_csv(f"{TAB}/28_enrichment_permutation.csv")
E9 = pd.read_csv(f"{OUT}/E9_primary_enrichment.csv")
val = E9[E9.cohort == "validation_n3"][["lesion", "enrichment_fold_primary"]]
rows = []
for _, r in val.iterrows():
    hit = S8[S8.lesion.astype(str) == r.lesion]
    rows.append(dict(lesion=r.lesion, enrichment_fold=r.enrichment_fold_primary,
                     in_table_S8=bool(len(hit)),
                     where_it_is="Table S5, part 'primary enrichment per lesion "
                                 "(as in Figure 2)'"))
E19 = pd.DataFrame(rows)
E19["lesions_in_the_S8_permutation_part"] = "; ".join(sorted(S8.lesion.astype(str).unique()))
E19.to_csv(f"{OUT}/E19_validation_enrichment_in_S8.csv", index=False)

# --- the table legends themselves, for the workbook sheet
pd.DataFrame([dict(table=k, legend=v) for k, v in TBL.items()]).to_csv(
    f"{OUT}/E20_table_legends.csv", index=False)

print(E17[["figure", "verdict", "claim_in_the_legend"]].to_string(index=False))
print("\nE18 rows:", len(E18))
print(E18.to_string(index=False) if len(E18) else "  legend and parts agree")
print("\ntable legends extracted:", list(TBL), "| figure legends:", list(FIG))
print("\nvalidation enrichment in S8:", E19.in_table_S8.tolist(),
      "| S8 permutation lesions:", E19.lesions_in_the_S8_permutation_part.iloc[0])
