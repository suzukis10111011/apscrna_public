#!/usr/bin/env python
"""Pre-submission tidy-up: manuscript reconciliation (requests 1, 2b, 7).

The manuscript is manusript/20261006_IEJ_Manuscript_post-audit_redline_v2.docx, read
through pandoc twice: --track-changes=accept (the revised text) and
--track-changes=reject (the text before the tracked revisions). Both are reconciled
against the tables.

C1  every location that reports a within-sample rho, with the score convention it
    uses, whether the current wording states it, and the wording to use if not
C2  every number in the two manuscript versions, matched against the table corpus
C3  the numbers that disagree with the tables, or that this tidy-up changes
C4  which regenerated file supersedes which published one
C5  the bracketed placeholders left in the text and what resolves them

Output: 06_submission_prep/tables/C1..C5_*.csv
"""
import os
import re
import sys
import time

import numpy as np
import pandas as pd

W = "/path/to/Documents/apscrna_v2"
TAB = f"{W}/results/tables"
AUD = f"{W}/05_post_audit/tables"
OUT = f"{W}/06_submission_prep/tables"
MST = f"{W}/06_submission_prep/manuscript_text"
NOTE = ("pre-submission reconciliation; frozen gene sets, candidate factors and "
        "thresholds unchanged")


def log(m):
    s = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {m}"
    print(s, flush=True)
    open(f"{W}/06_submission_prep/logs/run_submission_prep.log", "a").write(s + "\n")


S3 = pd.read_csv(f"{OUT}/S3_rho_by_score_convention.csv")
S5 = pd.read_csv(f"{OUT}/S5_stromal_cellsets.csv")
SG = pd.read_csv(f"{OUT}/S5_tf_sign_by_cellset_depth_controlled.csv")
NQ = pd.read_csv(f"{AUD}/Q1b_null_summary.csv")
UNI = pd.read_csv(f"{OUT}/definitionA_RELA_input_unification.csv")

# ------------------------------------------------------------------ C1 wording
def s3v(coh, les, cellset, conv):
    q = S3[(S3.cohort == coh) & (S3.lesion == les) & (S3.cell_set == cellset)]
    return float(q.iloc[0][f"depth_controlled_rho_{conv}"]) if len(q) else np.nan


C1 = [
    dict(location="Abstract and Results, within-stromal coupling",
         quantity="pooled depth-controlled rho per cohort",
         score_used="composite (inflammatory_composite x tissueforming_composite)",
         value=f"discovery {s3v('discovery_n5','POOLED_discovery','stromal','composite'):.4f}, "
               f"validation {s3v('validation_n3','POOLED_validation','stromal','composite'):.4f}",
         current_wording_states_the_score=True,
         wording_to_use="unchanged: the text already says composite",
         source_file="06_submission_prep/tables/S3_rho_by_score_convention.csv"),
    dict(location="Figure 4A legend",
         quantity="per-lesion and pooled depth-controlled rho in stromal cells",
         score_used="composite", current_wording_states_the_score=True,
         value="per lesion, see S3_rho_by_score_convention.csv",
         wording_to_use="unchanged: 'between the composite inflammatory and "
                        "tissue-forming scores'",
         source_file="S3_rho_by_score_convention.csv"),
    dict(location="Figure S5A legend",
         quantity="raw versus depth-controlled rho per lesion, panel A",
         score_used="plain single-axis (p_A_cytokine_chemokine x p_D_tissueforming)",
         value="per lesion, see S3_rho_by_score_convention.csv",
         current_wording_states_the_score=False,
         wording_to_use="replace 'single-axis and composite definitions' with "
                        "'single-axis (whole-axis p_ scores, RELA targets retained) "
                        "and composite definitions'; panel A is the single-axis "
                        "panel and uses the plain score, not the leak-free score",
         source_file="S3_rho_by_score_convention.csv"),
    dict(location="Table S3",
         quantity="raw and depth-controlled correlations",
         score_used="published: composite and plain only; regenerated table carries "
                    "composite, plain and leak-free side by side",
         value="20 rows (2 cohorts x 2 cell sets x lesions and pooled)",
         current_wording_states_the_score=False,
         wording_to_use="column names now carry the convention "
                        "(raw_rho_composite / _plain_p / _leakfree_lf and the "
                        "matching depth_controlled_rho_* columns); legend to read "
                        "'raw and depth-controlled correlations under all three "
                        "score conventions'",
         source_file="S3_rho_by_score_convention.csv"),
    dict(location="Table S5",
         quantity="marker-defined versus cluster-defined stromal set",
         score_used="leak-free (lf_A x lf_D)",
         value="20 rows; score_columns column states the convention in every row",
         current_wording_states_the_score=False,
         wording_to_use="legend to read 'leak-free (lf_) programme scores "
                        "throughout'; the regenerated table states it per row and "
                        "gives the depth-controlled rho under both depth covariate "
                        "forms",
         source_file="S5_stromal_cellsets.csv"),
    dict(location="Figure 6A legend (not named in the request, same class of issue)",
         quantity="depth-controlled rho per compartment",
         score_used="leak-free (lf_A x lf_D), 12_coupling_by_compartment.py:49",
         value="see results/tables/12_coupling_by_compartment.csv",
         current_wording_states_the_score=False,
         wording_to_use="'single-axis inflammatory and tissue-forming scores' is "
                        "ambiguous; 12 uses the leak-free scores, so read "
                        "'leak-free single-axis scores'",
         source_file="results/tables/12_coupling_by_compartment.csv"),
]
for r in C1:
    r["note"] = NOTE
pd.DataFrame(C1).to_csv(f"{OUT}/C1_score_convention_annotations.csv", index=False)

# ------------------------------------------------------- numeric corpus
CORPUS_FILES = [
    (f"{OUT}/28_manuscript_numbers_corrected.csv", "corrected number record"),
    (f"{OUT}/S3_rho_by_score_convention.csv", "S3 regenerated"),
    (f"{OUT}/S5_stromal_cellsets.csv", "S5 regenerated"),
    (f"{OUT}/S5_tf_sign_by_cellset_depth_controlled.csv", "S5 sign regenerated"),
    (f"{OUT}/S2_tf_screen_discovery_full_precision.csv", "S2 discovery regenerated"),
    (f"{OUT}/S2_tf_screen_validation_reproduction_full_precision.csv",
     "S2 validation regenerated"),
    (f"{OUT}/definitionA_RELA_input_unification.csv", "definition A unified"),
    (f"{AUD}/Q1b_null_summary.csv", "S11a random-factor null"),
    (f"{AUD}/Q2_definitionA_baselines.csv", "S11c dispersion baselines"),
    (f"{AUD}/Q3_donor_mixture.csv", "S12a single-donor check"),
    (f"{TAB}/28_qc_all_stages.csv", "published QC"),
    (f"{TAB}/28_enrichment_permutation.csv", "published enrichment permutation"),
    (f"{TAB}/33_defA_permutation_null.csv", "published definition A null"),
    (f"{TAB}/33_state_pct_8lesions.csv", "published state proportions"),
    (f"{TAB}/33_dip_by_lesion.csv", "published dip test"),
    (f"{TAB}/12_coupling_by_compartment.csv", "published compartment coupling"),
    (f"{TAB}/13_gse164241_summary.csv", "published gingival cohort"),
    (f"{TAB}/13_gse164241_tf_sign_reproduction.csv", "published gingival signs"),
    (f"{TAB}/34_lolo_rulebased.csv", "published leave-one-lesion-out"),
    (f"{TAB}/32_coactivity_ruleset_sensitivity.csv", "published co-activity"),
    (f"{TAB}/32_regulon_overlap.csv", "published regulon overlap"),
    (f"{TAB}/17_receptors_detection_fraction.csv", "published receptor detection"),
    (f"{TAB}/10_receptors_by_state.csv", "published receptor pseudobulk"),
    (f"{TAB}/09_geo_vs_kb.csv", "published GEO vs kb"),
    (f"{TAB}/09_tf_agreement_geo_vs_kb.csv", "published GEO vs kb TF agreement"),
    (f"{TAB}/14_robustness_genesets_scoring.csv", "published gene-set robustness"),
    (f"{TAB}/14_robustness_integration.csv", "published integration robustness"),
    (f"{TAB}/03c_compartment_by_lesion_discovery.csv", "published compartments"),
    (f"{TAB}/05_compartment_by_lesion_validation.csv", "published compartments"),
    (f"{TAB}/03c_compartment_by_lesion_validation.csv", "published compartments"),
    (f"{TAB}/17_receptors_mixedmodel.csv", "published receptor mixed model"),
]
corpus = {}


def register(v, label):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return
    try:
        x = float(v)
    except (TypeError, ValueError):
        return
    for dp in range(0, 5):
        corpus.setdefault(round(abs(x), dp), set()).add(label)
    if abs(x) < 1:
        for dp in (1, 2, 3):
            corpus.setdefault(round(abs(x) * 100, dp), set()).add(label + " (as %)")


for path, label in CORPUS_FILES:
    if not os.path.exists(path):
        log(f"corpus: missing {os.path.basename(path)}")
        continue
    T = pd.read_csv(path)
    for c in T.columns:
        s = pd.to_numeric(T[c], errors="coerce")
        for v in s.dropna().unique():
            register(v, f"{os.path.basename(path)}:{c}")
        if T[c].dtype == object:
            for cell in T[c].dropna().astype(str).unique():
                for tok in re.findall(r"-?\d+(?:\.\d+)?", cell):
                    register(tok, f"{os.path.basename(path)}:{c}")
# derived quantities that the methods text quotes but no table stores as a cell
import json as _json
_ax = _json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))["axes"]
for k, v in _ax.items():
    register(len(v), f"frozen/genesets/stage1_geneset_definitions.json:axes[{k}]")
_net = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
register(len(_net), "frozen/net/collectri_regulons.csv:edges")
register(_net.source.nunique(), "frozen/net/collectri_regulons.csv:sources")
register(int((_net.source == "RELA").sum()), "frozen/net/collectri_regulons.csv:RELA targets")
_qc = pd.read_csv(f"{TAB}/28_qc_all_stages.csv")
for (coh, st), g in _qc.groupby(["cohort", "stage"]):
    register(int(g.n_cells.sum()), f"28_qc_all_stages.csv:{coh}/{st} (cohort sum)")
register(23993, "data/processed/all_cells_validation.h5ad:n_features")
register(22669, "data/processed/all_cells_discovery.h5ad:n_features")
log(f"corpus: {len(corpus)} distinct rounded values")

# ------------------------------------------------------- C2 manuscript numbers
HEAD = re.compile(r"^(Abstract|Introduction|Materials and methods|Results|Discussion|"
                  r"Conclusion[s]?|References|Keywords|Acknowledgements|"
                  r"Ethical approval|Data availability statement|"
                  r"Supplementary figures|Supplementary tables|Figure \d|Figure S\d)")
NUM = re.compile(r"(?<![\w.])([−–-]?\d{1,3}(?:[, ]\d{3})*(?:\.\d+)?|[−–-]?\d*\.\d+)"
                 r"(?![\d])\s*(%|/\s*\d+)?")
SKIP_CTX = re.compile(r"\(\d{4}\)|et al|Nature|Journal|doi|pp\.|ISBN")

rows = []
for version in ["accept", "reject"]:
    txt = open(f"{MST}/manuscript_{version}.txt").read()
    section = "front matter"
    for ln, line in enumerate(txt.split("\n"), 1):
        if not line.strip():
            continue
        m = HEAD.match(line.strip())
        if m:
            section = m.group(1)
        if section == "References":
            continue
        for mo in NUM.finditer(line):
            raw, suf = mo.group(1), (mo.group(2) or "").replace(" ", "")
            ctx = line[max(0, mo.start() - 90):mo.end() + 90]
            if SKIP_CTX.search(line[max(0, mo.start() - 30):mo.end() + 30]):
                continue
            s = raw.replace("−", "-").replace("–", "-").replace(",", "").replace(" ", "")
            try:
                val = abs(float(s))
            except ValueError:
                continue
            printed = raw + suf
            dp = len(s.split(".")[1]) if "." in s else 0
            key = round(val, dp)
            hit = corpus.get(key, set())
            if suf == "%" and not hit:
                hit = corpus.get(round(val, dp), set())
            rows.append(dict(
                version=version, section=section, line=ln, printed=printed,
                value=val, decimals=dp, matched=bool(hit),
                n_sources=len(hit),
                sources="; ".join(sorted(hit)[:4]),
                context=re.sub(r"\s+", " ", ctx).strip(), note=NOTE))
C2 = pd.DataFrame(rows)
C2.to_csv(f"{OUT}/C2_manuscript_numbers_extracted.csv", index=False)
log(f"C2: {len(C2)} numbers ({C2.groupby('version').matched.mean().round(3).to_dict()} "
    f"matched fraction)")

# --------------------------------------------------- C3 actionable differences
acc = open(f"{MST}/manuscript_accept.txt").read()
rej = open(f"{MST}/manuscript_reject.txt").read()
C3 = []


def find_ctx(txt, needle, width=150):
    i = txt.find(needle)
    return re.sub(r"\s+", " ", txt[max(0, i - width):i + width]).strip() if i >= 0 else ""


nq = NQ[(NQ.version == "as_screened")
        & (NQ.pool == "sign_consistent_non_candidates")].set_index("axis")
nq_all = NQ[(NQ.version == "as_screened")
            & (NQ.pool == "all_non_candidates")].set_index("axis")
for axis, lab in [("lf_A_cytokine_chemokine", "inflammatory"),
                  ("lf_D_tissueforming", "tissue-forming")]:
    old = f"{nq_all.loc[axis, 'null_all3_mean']:.1f}"
    new = f"{nq.loc[axis, 'null_all3_mean']:.1f}"
    C3.append(dict(
        item=f"Figure 3E legend, random-set mean sign agreement ({lab})",
        in_manuscript=old, after_this_pass=new,
        changes=old != new,
        reason="Figure 3E legend only. The Results sentence already gives both "
               "pools (7.5 / 10.3 and 10.4 / 13.2); the legend quotes the "
               "all-non-candidate pool, while the redrawn panel marks the "
               "sign-consistent pool the request names",
        where=find_ctx(acc, "random sets of 20 non-candidate factors"),
        source_file="Q1b_null_summary.csv", note=NOTE))

for _, r in UNI[UNI.replacement_needed].iterrows():
    if r.cell_set != "primary_cluster_based" or r.lesion == "POOLED":
        continue
    C3.append(dict(
        item=f"Definition A, {r.lesion}: value carried in Table S5 and Table S8",
        in_manuscript=f"{r.published_18_allcell_ulm:.4f} (RELA from the all-cell ULM)",
        after_this_pass=f"{r.unified_value_RELA_from_cellset_ulm:.4f} "
                        f"(RELA from the ULM on the reported cell set)",
        changes=True,
        reason="the main-text and Figure 4B values already use the ULM on the stromal "
               "subset; S5 and S8 used a ULM on the whole all-cell object for the "
               "same lesions",
        where="Table S5 / Table S8 only; the main text and Figure 4B are unchanged",
        source_file="definitionA_RELA_input_unification.csv", note=NOTE))

old_cnt = pd.read_csv(f"{TAB}/18_tf_sign_by_cellset.csv")
for (coh, cs, axis), g in SG.groupby(["cohort", "cell_set", "axis"]):
    if axis not in ("lf_A_cytokine_chemokine", "lf_D_tissueforming"):
        continue
    o = old_cnt[(old_cnt.cohort == coh) & (old_cnt.cellset == cs)
                & (old_cnt.axis == axis)]
    ov = int(o.all_same_sign.astype(str).str.lower().isin(["true"]).sum()) if len(o) else np.nan
    nv = int(g.all_same_sign_depth_controlled.sum())
    if ov != nv:
        C3.append(dict(
            item=f"Table S5 sign agreement, {coh} / {cs} / {axis}",
            in_manuscript=f"{ov}/{len(g)} (plain Spearman)",
            after_this_pass=f"{nv}/{len(g)} (depth-controlled, as in Table S2)",
            changes=True,
            reason="request 4: the two supplementary tables counted agreement on "
                   "different estimators",
            where="Table S5 only", source_file="S5_tf_sign_by_cellset_depth_controlled.csv",
            note=NOTE))

for q, lab in [("[figure to be regenerated]", "figure placeholder"),
               ("[verify for CAP1", "methods placeholder")]:
    n_acc = acc.count(q)
    if n_acc:
        C3.append(dict(
            item=f"{lab}: '{q}'", in_manuscript=f"{n_acc} occurrence(s)",
            after_this_pass="resolved",
            changes=True,
            reason=("Figure 3E and Figure 4B are regenerated in "
                    "06_submission_prep/figures" if "figure" in lab else
                    "the validation source publication states three named donors and "
                    "describes loading the suspension on a 10x Chromium controller; "
                    "no pooling, multiplexing or hashing is reported "
                    "(Q3b_source_pooling_and_sex.csv)"),
            where=find_ctx(acc, q), source_file=("06_submission_prep/figures"
                                                 if "figure" in lab
                                                 else "Q3b_source_pooling_and_sex.csv"),
            note=NOTE))

A730 = pd.read_csv(f"{TAB}/33_tf_validation_all730.csv")
CV = pd.read_csv(f"{TAB}/03c_compartment_by_lesion_validation.csv")
n717 = int(A730[A730.axis == A730.axis.iloc[0]].in_validation.astype(str)
           .str.lower().isin(["true"]).sum())
sch = [c for c in CV.columns if "schwann" in c.lower() and c.startswith("n_")]
n365 = int(CV[sch[0]].sum()) if sch else np.nan
for item, printed, got, src in [
        ("Results and Figure 3 legend: factors evaluable in both cohorts",
         717, n717, "33_tf_validation_all730.csv (in_validation)"),
        ("Results: Schwann/neural cells, validation cohort",
         365, n365, "03c_compartment_by_lesion_validation.csv")]:
    C3.append(dict(item=item, in_manuscript=str(printed),
                   after_this_pass=str(got), changes=bool(printed != got),
                   reason=("derived count not stored as a single table cell; "
                           "recomputed here from the source table"),
                   where="", source_file=src, note=NOTE))
C3.append(dict(
    item="Methods placeholder: '[verify selection rule against 05_post_audit README]'",
    in_manuscript="1 occurrence", after_this_pass="resolved", changes=True,
    reason="the ten factors are those with the smallest absolute pooled "
           "depth-controlled rho against lf_D_tissueforming within that cohort, "
           "selected post hoc and listed per cohort in Q2_null_tf_selection.csv; the "
           "discovery and validation sets are therefore different factors",
    where="Materials and methods, definition A reference values",
    source_file="Q2_null_tf_selection.csv", note=NOTE))
_e = NQ[(NQ.axis == "lf_A_cytokine_chemokine") & (NQ.version == "overlap_removed")
        & (NQ.pool == "all_non_candidates")].iloc[0]
C3.append(dict(
    item="Results placeholder: '[verify effect-size p]' (overlap-removed, inflammatory)",
    in_manuscript="1 occurrence", after_this_pass=f"p = {_e.empirical_p_effect} "
                  f"(median |rho| {_e.candidates_median_abs_validation_rho:.2f} "
                  f"against {_e.null_median_abs_rho_mean:.2f})",
    changes=True, reason="value available in the null summary",
    where="Results, overlap-removed activities",
    source_file="Q1b_null_summary.csv", note=NOTE))
_s1 = NQ[(NQ.axis == "lf_D_tissueforming") & (NQ.version == "overlap_removed")]
C3.append(dict(
    item="Consistency: which draw pool is quoted for the overlap-removed results",
    in_manuscript="p = 0.67 (tissue-forming sign), from the all-non-candidate pool",
    after_this_pass="the sign-consistent pool gives p = "
                    f"{float(_s1[_s1.pool == 'sign_consistent_non_candidates'].empirical_p_sign_all3.iloc[0])}",
    changes=False,
    reason="the as-screened results quote both pools; the overlap-removed sentence "
           "quotes only the all-non-candidate pool. Both are in Table S11",
    where="Results, overlap-removed activities", source_file="Q1b_null_summary.csv",
    note=NOTE))

unmatched = C2[(C2.version == "accept") & (~C2.matched)
               & (C2.section.isin(["Results", "Abstract"]) | C2.section.str.startswith("Figure"))]
for _, r in unmatched.iterrows():
    C3.append(dict(item=f"unmatched number in {r.section}: {r.printed}",
                   in_manuscript=r.printed, after_this_pass="(no table value found)",
                   changes=False,
                   reason="no value in the table corpus rounds to this number at the "
                          "printed precision; check whether it is a derived or "
                          "narrative figure",
                   where=r.context, source_file="", note=NOTE))
pd.DataFrame(C3).to_csv(f"{OUT}/C3_manuscript_vs_tables_discrepancies.csv", index=False)
log(f"C3: {len(C3)} rows")

# ------------------------------------------------------------- C4 supersedes
C4 = [
    ("results/tables/04_rho_inflam_vs_tissueforming_discovery.csv; "
     "05_rho_inflam_vs_tissueforming_validation.csv",
     "06_submission_prep/tables/S3_rho_by_score_convention.csv",
     "Table S3", "adds the leak-free columns and names the convention in every "
     "column; the composite and plain values are reproduced exactly "
     "(max |difference| 0)"),
    ("results/tables/04_tf_screen_discovery_candidates.csv",
     "06_submission_prep/tables/S2_tf_screen_discovery_full_precision.csv",
     "Table S2", "same values without the 4-dp rounding, all six axes"),
    ("results/tables/05_tf_screen_validation_reproduction.csv",
     "06_submission_prep/tables/"
     "S2_tf_screen_validation_reproduction_full_precision.csv",
     "Table S2", "same values without the 4-dp rounding"),
    ("results/tables/18_stromal_markerdefined_vs_primary.csv",
     "06_submission_prep/tables/S5_stromal_cellsets.csv",
     "Table S5", "definition A now uses the RELA activity from a ULM on the reported "
     "cell set; rho given under both depth covariate forms; score convention stated "
     "per row"),
    ("results/tables/18_tf_sign_by_cellset.csv",
     "06_submission_prep/tables/S5_tf_sign_by_cellset_depth_controlled.csv",
     "Table S5", "sign agreement counted with the depth-controlled estimator used in "
     "Table S2; the plain Spearman column is retained for comparison"),
    ("results/tables/28_manuscript_numbers.csv",
     "06_submission_prep/tables/28_manuscript_numbers_corrected.csv",
     "number record", "score convention, cell set and source file added per row; the "
     "four rows labelled (lf_) relabelled to the plain score and the leak-free values "
     "added; post-audit quantities added"),
    ("results/figures_ms/Fig3_tf_screen.pdf",
     "06_submission_prep/figures/Fig3_tf_screen_v3.pdf",
     "Figure 3", "panel E carries the random-20-factor null mean"),
    ("results/figures_ms/Fig4_coupling_coexistence.pdf",
     "06_submission_prep/figures/Fig4_coupling_coexistence_v3.pdf",
     "Figure 4", "panel B carries the unrelated-factor baseline and the "
     "self-conditioning floor on a broken axis; the permutation band is unchanged"),
    ("(new)", "06_submission_prep/tables/Supplementary_Tables_S11_S12.xlsx",
     "Tables S11 and S12", "new supplementary tables assembled from 05_post_audit"),
]
pd.DataFrame(C4, columns=["published_file", "replacement_file", "manuscript_location",
                          "what_changed"]).assign(note=NOTE).to_csv(
    f"{OUT}/C4_supersedes_map.csv", index=False)

# ------------------------------------------------------------- C5 placeholders
C5 = []
for version, txt in [("accept", acc), ("reject", rej)]:
    for mo in re.finditer(r"\[[^\]\n]{3,120}\]", txt):
        C5.append(dict(version=version, placeholder=mo.group(0),
                       context=re.sub(r"\s+", " ",
                                      txt[max(0, mo.start() - 140):mo.end() + 60]).strip(),
                       note=NOTE))
pd.DataFrame(C5).to_csv(f"{OUT}/C5_placeholders.csv", index=False)
log(f"C4/C5 written ({len(C4)} supersede rows, {len(C5)} placeholders)")
