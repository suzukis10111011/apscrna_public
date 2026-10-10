#!/usr/bin/env python
"""Purpose: Consolidate supplementary tables S1-S9 into one English-column xlsx workbook
Input: results/tables/ CSV files for S1-S9 (GEO matrix primary analysis, kb requant as S6 robustness check)
Output: results/Supplementary_Tables_S1-S9.xlsx

30: Consolidate supplementary tables S1-S9 into a single xlsx (one table per sheet, column names in English).
The primary analysis is the GEO matrix version. The kb requantification appears only as a supplementary robustness check in S6.
Output: results/Supplementary_Tables_S1-S9.xlsx (+ CSV of each sheet are the original files in results/tables/)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np

import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
T = os.environ.get("APSC_TABDIR", f"{W}/results/tables")   # Source from which the input (table) is read
OUTX = f"{W}/results/Supplementary_Tables_S1-S10.xlsx"

# English translation of description columns that may contain Japanese (translation of values; column names are already in English)
JP2EN = {
    "クラスタ判定": "cluster-based assignment",
    "候補は 5病変同符号かつ |pooled ρ| 上位20。旧セットから ATF1・MTA1 が外れ、RELA・NFKB1 が入る":
        "candidates are TFs with a consistent sign across all 5 lesions and |pooled rho| in the top 20; relative to the previous set, ATF1 and MTA1 drop out and RELA and NFKB1 enter",
    "標準化2変数の PCA では PC1 分散比 = (1+|Pearson r|)/2 なので、相関と独立な証拠ではない":
        "for a PCA of two standardised variables the PC1 variance ratio equals (1+|Pearson r|)/2, so it is a deterministic function of the correlation and not independent evidence",
    "n<20 のため算出せず": "not computed (n<20)",
    "細胞数<100": "n cells < 100",
    "scVI は環境に無いため省略。axis = lf_D_tissueforming":
        "scVI not available in this environment and was omitted; axis = lf_D_tissueforming",
    "病変内 Q80 以上 かつ 4陰性マーカーが病変内中央値未満":
        "within-lesion 80th percentile or above for the positive markers, and all 4 negative markers below the within-lesion median",
    "病変×状態で counts 合計 -> CPM -> log2(CPM+1)":
        "counts summed per lesion x state, then CPM, then log2(CPM+1)",
    "MixedLM: log正規化発現 ~ state, ランダム切片 = lesion（探索+検証を合算）":
        "MixedLM: log-normalised expression ~ state, random intercept = lesion (discovery and validation pooled)",
    "今回の区画判定による間質集合": "stromal set from the compartment assignment of this study",
    "前回の間質バーコード集合（今回の行列・正規化・スコアで再計算）":
        "previous stromal barcode set, recomputed with this study's matrix, normalisation and scores",
    "composite (inflammatory vs tissueforming)": "composite (inflammatory vs tissue-forming)",
}
JP_FALLBACK = [
    ("病変内", "within-lesion"), ("以上 かつ", " or above and "),
    ("陰性マーカーが病変内中央値未満", "negative markers below the within-lesion median"),
    ("は行列に不在で未適用", " absent from the matrix, filter not applied"),
    ("細胞数", "n cells"), ("探索", "discovery"), ("検証", "validation"),
    ("病変", "lesion"), ("状態", "state"), ("区画", "compartment"),
    ("（", " ("), ("）", ")"), ("、", ", "), ("。", ". "),
]


def clean(d):
    """説明列に残る日本語を英訳する（列名は元から英語）。辞書に無いものは語句置換で英語化し、
    それでも日本語が残る値は最後に検出して報告する（黙って出さない）。"""
    import re
    jp = re.compile(r"[\u3040-\u30ff\u4e00-\u9fff]")
    def conv(v):
        if not isinstance(v, str) or not jp.search(v):
            return v
        if v in JP2EN:
            return JP2EN[v]
        out = v
        for a, b in JP_FALLBACK:
            out = out.replace(a, b)
        return out
    d = d.copy()
    for c in d.columns:
        if d[c].dtype == object:
            d[c] = d[c].map(conv)
    return d


SHEETS = {
    "S1_QC": ([("28_qc_all_stages.csv", None)],
              "Cell counts at each QC stage, per sample, both cohorts"),
    "S2_TF_screen": ([("04_tf_screen_discovery_candidates.csv", "discovery_n5"),
                      ("05_tf_screen_validation_reproduction.csv", "validation_n3")],
                     "Full ranking of the 730-TF activity screen on both named axes; discovery screen and validation sign reproduction"),
    "S3_rho_raw_vs_depth": ([("04_rho_inflam_vs_tissueforming_discovery.csv", "discovery_n5"),
                             ("05_rho_inflam_vs_tissueforming_validation.csv", "validation_n3")],
                            "Raw and depth-controlled Spearman rho between the inflammatory and tissue-forming scores"),
    "S4_genesets_scoring": ([("14_robustness_genesets_scoring.csv", None)],
                            "Robustness across gene-set definitions and scoring methods"),
    "S5_integration_cellsets": ([("14_robustness_integration.csv", "integration method"),
                                 ("18_stromal_markerdefined_vs_primary.csv", "stromal cell-set definition")],
                                "Robustness to the integration method and to an embedding-free, marker-defined stromal set"),
    "S6_GEO_vs_kb": ([("09_geo_vs_kb.csv", "key quantities"),
                      ("09_stromal_jaccard.csv", "stromal barcode agreement"),
                      ("09_tf_agreement_geo_vs_kb.csv", "TF screen agreement")],
                     "Supplementary robustness check: GEO author matrices vs independent kallisto|bustools re-quantification"),
    "S7_receptors": ([("10_receptors_by_state.csv", "pseudobulk, paired Wilcoxon"),
                      ("17_receptors_mixedmodel.csv", "cell-level linear mixed model"),
                      ("17_receptors_detection_fraction.csv", "detection fraction")],
                     "Cytokine-receptor expression between the inflammatory and tissue-forming states"),
    "S8_compartment_coupling": ([("12_coupling_by_compartment.csv", None),
                                 ("28_enrichment_permutation.csv", "top-decile permutation test")],
                                "Axis coupling per compartment, and the permutation test for top-decile compartment enrichment"),
    "S10_coactivity": ([("32_set_definition.csv", "TF set definition by rule"),
                        ("32_coactivity_ruleset_sensitivity.csv", "per-lesion and per-state statistics"),
                        ("32_crossside_pairs_sign_consistency.csv", "cross-side TF pairs"),
                        ("32_regulon_overlap.csv", "regulon overlap"),
                        ("32_old_vs_new_comparison.csv", "rule-based vs previous TF set"),
                        ("32_pc1_is_function_of_r.csv", "note on the PC1 variance ratio")],
                       "Co-activity of the two transcription-factor programmes, with the TF sets defined by rule from this study's own screen"),
    "S9_GSE164241": ([("13_gse164241_summary.csv", "per donor"),
                      ("13_gse164241_tf_sign_reproduction.csv", "TF sign reproduction")],
                     "Independent periodontitis cohort GSE164241 (gingiva)"),
}

meta, made = [], {}
for sheet, (parts, desc) in SHEETS.items():
    frames = []
    for fn, tag in parts:
        p = f"{T}/{fn}"
        if not os.path.exists(p):
            print(f"  !! {sheet}: {fn} が無い — 飛ばす"); continue
        d = clean(pd.read_csv(p))
        if "source_file" in d.columns: d = d.rename(columns={"source_file": "source_file_original"})
        if "table_part" in d.columns: d = d.rename(columns={"table_part": "table_part_original"})
        d.insert(0, "table_part", tag if tag else "main")
        d.insert(1, "source_file", fn)
        frames.append(d)
    if not frames:
        continue
    D = pd.concat(frames, ignore_index=True, sort=False)
    made[sheet] = D
    meta.append(dict(sheet=sheet, description=desc, n_rows=len(D), n_columns=D.shape[1],
                     source_files="; ".join(f for f, _ in parts if os.path.exists(f"{T}/{f}"))))
M = pd.DataFrame(meta)
with pd.ExcelWriter(OUTX, engine="openpyxl") as xw:
    M.to_excel(xw, sheet_name="S0_index", index=False)
    for sheet, D in made.items():
        D.to_excel(xw, sheet_name=sheet[:31], index=False)
print(M.to_string(index=False))
print(f"\n書き出し: {OUTX} ({os.path.getsize(OUTX)/1e6:.2f} MB, {len(made)+1} シート)")
jp = {s: [c for c in D.columns if any(ord(ch) > 0x3000 for ch in str(c))] for s, D in made.items()}
print("日本語の列名:", {k: v for k, v in jp.items() if v} or "なし")
left = {s: int(sum(D[c].astype(str).str.contains(r"[぀-ヿ一-鿿]", regex=True).sum() for c in D.columns if D[c].dtype == object))
        for s, D in made.items()}
print("日本語を含む値の残存:", {k: v for k, v in left.items() if v} or "なし")
print("\n30 DONE")
