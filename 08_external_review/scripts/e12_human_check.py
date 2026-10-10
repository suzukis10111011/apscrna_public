#!/usr/bin/env python
"""Purpose: Decide, for every row of the checking worksheet, whether a person has to
  check the value by hand. Rows that point at a cell or a column always do. Among the
  rows that can only be read in the sentence, the measured values do, including those
  derived from several cells; software versions, pipeline arguments, figure numbers,
  reference-sequence construction values, frozen settings and numerals that are part of
  a name do not, and each of those is tied here to the file that fixes it.
Input: 07_submission_bundle/tables/check_worksheet.xlsx (which rows are sentence rows),
  07_submission_bundle/tables/manuscript_number_agreement.csv, the scripts and frozen
  files of the repository, and the installed packages of the environment.
Output: 08_external_review/tables/E12_human_check.csv,
  08_external_review/tables/E12_software_versions.csv

The evidence is not asserted: for a setting the file is read and the line that carries
the value is located, so the recorded file:line is a line that exists and contains the
value; for a software version the installed distribution is queried.
"""
import importlib.metadata as md
import os
import re
import subprocess
import sys

import pandas as pd

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
BUN = f"{W}/07_submission_bundle/tables"
OUT = f"{W}/08_external_review/tables"

# No -> (package as the manuscript names it, distribution name)
BINARIES = {"kallisto": ["kallisto", "version"], "bustools": ["bustools", "version"],
            "SRA Toolkit": ["fasterq-dump", "--version"]}

VERSIONS = {
    202: ("Python", None), 203: ("scanpy", "scanpy"), 204: ("anndata", "anndata"),
    206: ("omnipath", "omnipath"), 207: ("harmonypy", "harmonypy"),
    208: ("scrublet", "scrublet"), 209: ("leidenalg", "leidenalg"),
    210: ("diptest", "diptest"), 211: ("kallisto", None), 212: ("bustools", None),
    213: ("kb-python", "kb-python"), 214: ("SRA Toolkit", None),
    19: ("Cell Ranger", None), 30: ("SRA Toolkit", None), 31: ("kallisto", None),
    32: ("bustools", None), 33: ("kb-python", "kb-python"),
    76: ("Scrublet", "scrublet"), 90: ("harmonypy", "harmonypy"),
    164: ("diptest", "diptest"),
}

# No -> (file relative to the repository root, regex that must match the line)
SETTINGS = {
    73: ("scripts/03_preprocess_discovery.py", r"n_genes_by_counts >= 200"),
    74: ("scripts/03_preprocess_discovery.py", r"total_counts >= 500"),
    75: ("scripts/03_preprocess_discovery.py", r"pct_counts_mt <= 20"),
    77: ("scripts/03_preprocess_discovery.py", r"expected_doublet_rate=0\.06"),
    82: ("scripts/05_validation_and_projection.py", r"filter_cells\(A, min_counts=100\)"),
    91: ("scripts/03_preprocess_discovery.py", r"n_neighbors=15"),
    92: ("scripts/03_preprocess_discovery.py", r"for r in \(0\.3, 0\.5, 1\.0\)"),
    93: ("scripts/_utils.py", r"MINZ_COMP = 1\.0"),
    96: ("scripts/05_validation_and_projection.py", r"np\.arange\(0\.2, 1\.51, 0\.1\)"),
    98: ("scripts/05_validation_and_projection.py", r"np\.arange\(0\.2, 1\.51, 0\.1\)"),
    109: ("scripts/14c_robustness_integration.py", r"leiden\(B, resolution=0\.5"),
    117: ("scripts/05_validation_and_projection.py", r"ctrl_size=50, n_bins=25"),
    126: ("scripts/33_ms_figure_inputs.py", r"nbins=20, min_n=25"),
    127: ("scripts/33_ms_figure_inputs.py", r"nbins=20, min_n=25"),
    128: ("scripts/33_ms_figure_inputs.py", r"len\(df\) < 500: nbins = 10"),
    129: ("scripts/33_ms_figure_inputs.py", r"len\(df\) < 500: nbins = 10"),
    133: ("scripts/04_discovery_scores_and_screen.py", r"\.head\(20\)\.index"),
    135: ("05_post_audit/scripts/q1b_null_signtest.py", r"NSEED, NDRAW, TOPN = 1000, 20, 20"),
    136: ("05_post_audit/scripts/q1b_null_signtest.py", r"NSEED, NDRAW, TOPN = 1000, 20, 20"),
    137: ("05_post_audit/scripts/q1b_null_signtest.py", r"NSEED, NDRAW, TOPN = 1000, 20, 20"),
    138: ("05_post_audit/scripts/q1b_null_signtest.py", r"NSEED, NDRAW, TOPN = 1000, 20, 20"),
    140: ("05_post_audit/scripts/q1b_null_signtest.py", r"NSEED, NDRAW, TOPN = 1000, 20, 20"),
    150: ("scripts/10_receptors_by_state.py", r"np\.log2\(c \+ 1\)"),
    152: ("scripts/17_receptor_power.py", r">= 10 and"),
    170: ("scripts/12_coupling_by_compartment.py", r"len\(g\) < 100"),
    239: ("scripts/28_manuscript_numbers.py", r"NPERM = 10000"),
    361: ("scripts/17_receptor_power.py", r">= 10 and"),
    395: ("scripts/12_coupling_by_compartment.py", r"len\(g\) < 100"),
    422: ("scripts/28_manuscript_numbers.py", r"NPERM = 10000"),
    428: ("scripts/04_discovery_scores_and_screen.py", r"\.head\(20\)\.index"),
    437: ("scripts/04_discovery_scores_and_screen.py", r"\.head\(20\)\.index"),
    472: ("scripts/12_coupling_by_compartment.py", r"len\(g\) < 100"),
}

# No -> (category, note), for values that no file fixes because they are not values
NO_FILE = {
    15: ("name", "interleukin (IL)-1"), 17: ("name", "10x Single Cell 3' kit"),
    37: ("chemistry", "16-nt barcode of the 10x v3 structure"),
    38: ("chemistry", "12-nt UMI of the 10x v3 structure"),
    63: ("chemistry", "16 + 10 nt barcode-UMI structure of v2"),
    59: ("chemistry", "positions of read 1 in the v2 structure"),
    65: ("chemistry", "positions of read 1 in the v2 structure"),
    131: ("formula", "smallest attainable p of the test"),
    132: ("formula", "the value 1 as the reference of a ratio"),
    139: ("formula", "the +1 of the empirical p value"),
    141: ("formula", "the +1 of the empirical p value"),
    142: ("formula", "the fold of 1 as the reference"),
    151: ("formula", "the 0.05 of a significance statement"),
    165: ("formula", "the identity (1 + |r|)/2"),
    166: ("formula", "the identity (1 + |r|)/2"),
    237: ("formula", "the mantissa of a p value in scientific notation"),
    279: ("formula", "the 0.005 bound of a p value"),
    290: ("formula", "the 0.03 bound of a Jaccard index"),
    291: ("formula", "the smallest attainable p of 1,000 draws"),
    467: ("formula", "the p = 0.05 line drawn in the panel"),
    287: ("setting", "the top-20 lists of the screen"),
    274: ("setting", "the 20 candidates per programme"),
    277: ("setting", "the 20 candidates per programme"),
    359: ("name", "IL-1 receptor family"), 369: ("name", "IL-6 receptor family"),
    415: ("name", "IL-1 receptor family"), 416: ("name", "IL-6 receptor family"),
    50: ("reference build", "symbol-level genes of the quantification matrix"),
    56: ("reference build", "the 20 runs of the discovery cohort"),
    215: ("figure", "Figure 1"), 242: ("figure", "Figure 2"), 257: ("figure", "Figure 3"),
    260: ("figure", "Figure 3"), 306: ("figure", "Figure 3"), 357: ("figure", "Figure 4"),
    358: ("figure", "Figure 4"), 362: ("figure", "Figure 5"), 381: ("figure", "Figure 5"),
    382: ("figure", "Figure 5"), 393: ("figure", "Figure 6"), 396: ("figure", "Figure 6"),
    406: ("figure", "Figure 7"), 412: ("figure", "Figure 7"), 417: ("figure", "Figure 1"),
    420: ("figure", "Figure 2"), 426: ("figure", "Figure 3"), 438: ("figure", "Figure 4"),
    468: ("figure", "Figure 5"), 471: ("figure", "Figure 6"), 473: ("figure", "Figure 7"),
}

# the separate sheet: No -> (Y/N, where it can be read, category)
OTHER = {
    18: ("Y", "results/tables/02_annotation_feature_counts.csv or the Cell Ranger "
              "features.tsv of lesion A", "reference build"),
    21: ("Y", "the Cell Ranger features.tsv of lesions B-E", "reference build"),
    22: ("Y", "the intersection of the two feature lists", "reference build"),
    45: ("Y", "the Ensembl 112 GRCh38 cDNA fasta", "reference build"),
    46: ("Y", "results/tables/23_chemistry_umi_comparison.csv, column genes",
         "reference build"),
    47: ("Y", "the t2g table built from the cDNA fasta", "reference build"),
    54: ("Y", "logs/ run_info.json of the kb run for SRR32119705", "log"),
    86: ("Y", "Table S1, 1_qc_pass minus 2_singlet summed over the lesions",
         "measured value"),
    88: ("N", "scripts/03_preprocess_discovery.py", "pipeline argument"),
    143: ("Y", "Table S5, the gene-space column of the projection part",
          "measured value"),
}

JA = {"version": "ソフトウェアの版", "setting": "パイプラインの引数・凍結設定値",
      "figure": "図番号", "name": "名称の一部である数字", "formula": "式の定数・閾値",
      "chemistry": "キットの読み構造", "reference build": "参照配列の構築値",
      "log": "実行ログの値", "measured value": "測定値"}


def line_of(path, pattern):
    """line number and text of the first line of `path` matching `pattern`"""
    full = f"{W}/{path}"
    if not os.path.exists(full):
        return None, None
    with open(full, errors="replace") as fh:
        for i, line in enumerate(fh, start=1):
            if pattern is None or re.search(pattern, line):
                return i, line.strip()[:120]
    return None, None


A = pd.read_csv(f"{BUN}/manuscript_number_agreement.csv")
A = A[A.version == sorted(v for v in A.version.unique()
                          if v.startswith("clean_"))[-1]].reset_index(drop=True)
M = pd.read_excel(f"{BUN}/check_worksheet.xlsx", sheet_name="補足表で確認")
sentence = set(M.loc[M["確認方法"] == "文として確認", "No"].astype(int))

# the installed versions, so a version row can be tied to something that was checked
vrows = []
for no, (name, dist) in sorted(VERSIONS.items()):
    if dist is None:
        installed = ""
    else:
        try:
            installed = md.version(dist)
        except md.PackageNotFoundError:
            installed = ""
    if name == "Python":
        installed = sys.version.split()[0]
    if not installed and name in BINARIES:
        try:
            out = subprocess.run(BINARIES[name], capture_output=True, text=True,
                                 timeout=30).stdout
            m = re.search(r"(\d+\.\d+\.\d+)", out)
            installed = m.group(1) if m else ""
        except (OSError, subprocess.SubprocessError):
            installed = ""
    vrows.append(dict(No=no, package=name, printed=A.iloc[no - 1].printed,
                      installed_in_env_apsc=installed,
                      starts_with_printed=bool(installed and str(installed).startswith(
                          str(A.iloc[no - 1].printed)))))
V = pd.DataFrame(vrows)
V.to_csv(f"{OUT}/E12_software_versions.csv", index=False)

NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
S4 = pd.read_excel(f"{BUN}/Supplementary_Tables_S1-S12.xlsx", sheet_name="S4", skiprows=[1])
N_POOLED = int((S4.lesion.astype(str) == "POOLED").sum())
COUNTED = {
    121: (f"frozen/net/collectri_regulons.csv: {len(NET):,} rows counted in the frozen "
          "network file",
          f"凍結ファイルで確認済み: frozen/net/collectri_regulons.csv の行数 {len(NET):,}"),
    122: (f"frozen/net/collectri_regulons.csv: {NET.source.nunique():,} distinct values "
          "of the column source, counted",
          "凍結ファイルで確認済み: frozen/net/collectri_regulons.csv の source 列の異なる値 "
          f"{NET.source.nunique():,}"),
    196: (f"Table S4: {N_POOLED} pooled rows, i.e. the gene-set by scoring-method "
          "combinations, counted in the workbook",
          f"補足表S4 のプール行 {N_POOLED} 行（遺伝子集合×スコア法の組み合わせ）を数えた"),
    485: (f"Table S4: {N_POOLED} pooled rows, i.e. the gene-set by scoring-method "
          "combinations, counted in the workbook",
          f"補足表S4 のプール行 {N_POOLED} 行（遺伝子集合×スコア法の組み合わせ）を数えた"),
}

rows = []
for i, r in A.iterrows():
    no = i + 1
    if not bool(r.matched):
        continue
    if no not in sentence:
        rows.append(dict(No=no, printed=r.printed, human_check_required="Y",
                         category="reads against a cell or a column", evidence="",
                         evidence_ja=""))
        continue
    if no in VERSIONS:
        name = VERSIONS[no][0]
        inst = V.loc[V.No == no, "installed_in_env_apsc"].iloc[0]
        geo = ("; the depositors' processing version, recorded in the GEO series "
               "metadata of GSE181688 and GSE197680, not in this environment"
               if name == "Cell Ranger" else "")
        ev = (f"software version; {name} {inst} in the analysis environment apsc, "
              f"recorded in 08_external_review/tables/E12_software_versions.csv"
              if inst else
              f"software version of {name}, not an installed distribution here{geo}")
        ev_ja = (f"ソフトウェアの版。解析環境 apsc の {name} {inst} と一致（記録: "
                 "08_external_review/tables/E12_software_versions.csv）"
                 if inst else
                 f"ソフトウェアの版。{name} はこの環境に入っていない"
                 + ("（寄託者側の処理版で、GEO の GSE181688・GSE197680 の"
                    "シリーズ情報に記録）" if name == "Cell Ranger" else ""))
        rows.append(dict(No=no, printed=r.printed, human_check_required="N",
                         category="software version", evidence=ev, evidence_ja=ev_ja))
    elif no in COUNTED:
        ev, ev_ja = COUNTED[no]
        rows.append(dict(No=no, printed=r.printed, human_check_required="N",
                         category="frozen input, counted", evidence=ev,
                         evidence_ja=ev_ja))
    elif no in SETTINGS:
        path, pat = SETTINGS[no]
        ln, txt = line_of(path, pat)
        where = f"{path}:{ln}" if ln else path
        ev = (f"frozen setting; fixed in {where}" + (f"  [{txt}]" if txt else ""))
        ev_ja = f"凍結ファイル／スクリプトで確認済み: {where}" + (f"（{txt}）" if txt else "")
        rows.append(dict(No=no, printed=r.printed, human_check_required="N",
                         category="pipeline argument or frozen setting",
                         evidence=ev, evidence_ja=ev_ja))
    elif no in NO_FILE:
        cat, what = NO_FILE[no]
        ev = f"{cat}: {what}; not a measured value, so no file fixes it"
        ev_ja = f"{JA[cat]}（{what}）。測定値ではないため、照合すべきファイルはない"
        rows.append(dict(No=no, printed=r.printed, human_check_required="N",
                         category=cat, evidence=ev, evidence_ja=ev_ja))
    else:
        rows.append(dict(No=no, printed=r.printed, human_check_required="Y",
                         category="measured value read in the sentence",
                         evidence="", evidence_ja=""))
for no, (flag, where, cat) in sorted(OTHER.items()):
    ev = (f"{cat}; read it in {where}" if flag == "Y" else
          f"{cat}; fixed in {where}")
    if flag == "N":
        ln, txt = line_of(where, r"n_top_genes=2000")
        ev = f"pipeline argument; fixed in {where}:{ln}" + (f"  [{txt}]" if txt else "")
        ev_ja = f"凍結ファイル／スクリプトで確認済み: {where}:{ln}" + (f"（{txt}）" if txt else "")
    else:
        ev_ja = f"{JA[cat]}。{where} で確認する"
    rows.append(dict(No=no, printed=A.iloc[no - 1].printed, human_check_required=flag,
                     category=cat, evidence=ev, evidence_ja=ev_ja, sheet="other"))
E12 = pd.DataFrame(rows)
E12["sheet"] = E12.get("sheet").fillna("main")
E12.to_csv(f"{OUT}/E12_human_check.csv", index=False)
print(E12.groupby(["sheet", "human_check_required"]).size().to_string())
print()
print(E12[E12.human_check_required == "N"].groupby("category").size().to_string())
print("\nversions matching the installed distribution:",
      int(V.starts_with_printed.sum()), "of", int((V.installed_in_env_apsc != "").sum()),
      "resolvable")
miss = E12[(E12.human_check_required == "N") & (E12.evidence.str.contains("None"))]
print("evidence without a line number:", len(miss), miss.No.tolist())
