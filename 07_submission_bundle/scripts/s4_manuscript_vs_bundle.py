#!/usr/bin/env python
"""Submission bundle: every number in the manuscript against the regenerated tables.

Three texts are extracted with pandoc at run time: the clean submission manuscript
(the primary comparison target) and the tracked-changes file in both accept and
reject modes. The comparison corpus is the
twelve sheets of Supplementary_Tables_S1-S12.xlsx, the regenerated figure panel record
and caption check, the corrected number record in results/tables, and the derived
quantities that no table stores as a cell (gene-set sizes, regulon network totals,
QC cohort sums, feature counts).

Output: 07_submission_bundle/tables/manuscript_number_agreement.csv
        07_submission_bundle/tables/manuscript_unmatched_numbers.csv
        07_submission_bundle/manuscript_text/manuscript_{clean_v4,redline_v4_accept,redline_v4_reject}.txt
        07_submission_bundle/tables/manuscript_source_files.csv
"""
import decimal
import hashlib
import json
import os
import re
import subprocess
import time
from collections import Counter

import numpy as np
import pandas as pd

W = "/path/to/Documents/apscrna_v2"
B = f"{W}/07_submission_bundle"
TAB = f"{W}/results/tables"
MST = f"{B}/manuscript_text"
os.makedirs(MST, exist_ok=True)
NOTE = "submission bundle; compared against the regenerated Tables S1-S12"


def log(m):
    s = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {m}"
    print(s, flush=True)
    open(f"{B}/logs/run_bundle.log", "a").write(s + "\n")


# ---- the two v4 files: the clean submission text and the tracked-changes version
MDIR = f"{W}/manuscript"


def latest(kind):
    """Highest-numbered _v<N>.docx of that kind in manuscript/."""
    pat = re.compile(rf"{kind}_v(\d+)\.docx$")
    hits = [(int(m.group(1)), f"{MDIR}/{f}")
            for f in os.listdir(MDIR) if (m := pat.search(f))]
    if not hits:
        raise SystemExit(f"no {kind}_v<N>.docx in {MDIR}")
    return max(hits)


VER, CLEAN = latest("clean")
_, REDLINE = latest("redline")
VERSIONS = [(f"clean_v{VER}", CLEAN, "accept"),
            (f"redline_v{VER}_accept", REDLINE, "accept"),
            (f"redline_v{VER}_reject", REDLINE, "reject")]
prov = []
for label, path, mode in VERSIONS:
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}")
    subprocess.run(["pandoc", f"--track-changes={mode}", "-t", "plain", "--wrap=none",
                    path, "-o", f"{MST}/manuscript_{label}.txt"], check=True)
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]
    prov.append(dict(version=label, file=os.path.relpath(path, W),
                     track_changes=mode, bytes=os.path.getsize(path), sha256_16=h,
                     mtime_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                             time.gmtime(os.path.getmtime(path))),
                     n_chars=len(open(f"{MST}/manuscript_{label}.txt").read())))
    log(f"{label}: {os.path.basename(path)} {os.path.getsize(path)} B sha {h}")
pd.DataFrame(prov).to_csv(f"{B}/tables/manuscript_source_files.csv", index=False)

# ------------------------------------------------------------------- corpus
# A table value of 39.65 is printed as 39.7, which is round-half-up; Python's round() is
# half-to-even and would file it under 39.6, so the manuscript value would not match.
# The worksheet generator uses the same rule (s5_check_worksheet.py).
def rnd(v, dp):
    return float(decimal.Decimal(repr(float(v))).quantize(
        decimal.Decimal(1).scaleb(-dp), rounding=decimal.ROUND_HALF_UP))


corpus = {}


def register(v, label):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return
    if not np.isfinite(x):
        return
    for dp in range(0, 5):
        corpus.setdefault(rnd(abs(x), dp), set()).add(label)
    if abs(x) < 1:
        for dp in (1, 2, 3):
            corpus.setdefault(rnd(abs(x) * 100, dp), set()).add(label + " (as %)")


def scan(df, label):
    for c in df.columns:
        s = pd.to_numeric(df[c], errors="coerce")
        for v in s.dropna().unique():
            register(v, f"{label}:{c}")
        if df[c].dtype == object:
            for cell in df[c].dropna().astype(str).unique():
                for tok in re.findall(r"-?\d+(?:\.\d+)?", cell):
                    register(tok, f"{label}:{c}")


xl = pd.ExcelFile(f"{B}/tables/Supplementary_Tables_S1-S12.xlsx")
for sh in xl.sheet_names:
    if sh == "index":
        continue
    d = pd.read_excel(xl, sheet_name=sh, skiprows=[1])
    scan(d, f"Table {sh}")
for path, label in [(f"{B}/tables/figure_panel_n.csv", "figure_panel_n"),
                    (f"{B}/tables/caption_number_check.csv", "caption_number_check"),
                    (f"{TAB}/28_manuscript_numbers.csv", "number record"),
                    (f"{W}/06_submission_prep/tables/C6_validation_chemistry_check.csv",
                     "validation chemistry check"),
                    (f"{W}/06_submission_prep/tables/C6_base_composition_by_position.csv",
                     "validation chemistry base composition"),
                    (f"{TAB}/33_defA_permutation_null.csv", "33_defA_permutation_null"),
                    (f"{TAB}/33_state_pct_8lesions.csv", "33_state_pct_8lesions"),
                    (f"{TAB}/33_dip_by_lesion.csv", "33_dip_by_lesion"),
                    (f"{TAB}/33_tf_validation_all730.csv", "33_tf_validation_all730"),
                    (f"{TAB}/34_lolo_rulebased.csv", "34_lolo_rulebased"),
                    (f"{TAB}/25_validation_states_within_vs_ingest.csv", "25_states"),
                    (f"{TAB}/03c_compartment_by_lesion_discovery.csv", "03c_compartments"),
                    (f"{TAB}/05_compartment_by_lesion_validation.csv", "05_compartments"),
                    (f"{TAB}/03c_compartment_by_lesion_validation.csv", "03c_compartments_val")]:
    if os.path.exists(path):
        scan(pd.read_csv(path), label)

AX = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))["axes"]
for k, v in AX.items():
    register(len(v), f"frozen gene sets:axes[{k}]")
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
register(len(NET), "frozen network:edges")
register(NET.source.nunique(), "frozen network:sources")
register(int((NET.source == "RELA").sum()), "frozen network:RELA targets")
QC = pd.read_csv(f"{TAB}/28_qc_all_stages.csv")
for (coh, st), g in QC.groupby(["cohort", "stage"]):
    register(int(g.n_cells.sum()), f"28_qc_all_stages:{coh}/{st} (cohort sum)")
register(23993, "validation matrix:n_features")
register(22669, "discovery matrix:n_features")
A730 = pd.read_csv(f"{TAB}/33_tf_validation_all730.csv")
register(int(A730[A730.axis == A730.axis.iloc[0]].in_validation.astype(str)
             .str.lower().isin(["true"]).sum()),
         "33_tf_validation_all730:factors evaluable in both cohorts")
CV = pd.read_csv(f"{TAB}/03c_compartment_by_lesion_validation.csv")
sch = [c for c in CV.columns if "schwann" in c.lower() and c.startswith("n_")]
if sch:
    register(int(CV[sch[0]].sum()), "03c_compartment_by_lesion_validation:schwann total")
log(f"corpus: {len(corpus)} distinct rounded values")
# Extractions of earlier manuscript versions stay on disk when a version is replaced.
# They are listed, not removed: deleting files under the repository is the user's call.
_stale = sorted(f for f in os.listdir(MST)
                if f.endswith(".txt") and f"_v{VER}" not in f)
if _stale:
    log("superseded extractions still in manuscript_text/ (not used by this run): "
        + ", ".join(_stale))


# ------------------------------------------------------------- extraction
HEAD = re.compile(r"^(Abstract|Introduction|Materials and methods|Results|Discussion|"
                  r"Conclusion[s]?|References|Keywords|Acknowledgements|"
                  r"Ethical approval|Data availability statement|"
                  r"Supplementary figures|Supplementary tables|Figure \d|Figure S\d)")
NUM = re.compile(r"(?<![\w.])([−–-]?\d{1,3}(?:[, ]\d{3})*(?:\.\d+)?|[−–-]?\d*\.\d+)"
                 r"(?![\d])\s*(%|/\s*\d+)?")
SKIP = re.compile(r"\(\d{4}\)|et al|doi|pp\.|ISBN")

rows = []
for version, _, _ in VERSIONS:
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
            if SKIP.search(line[max(0, mo.start() - 30):mo.end() + 30]):
                continue
            s = raw.replace("−", "-").replace("–", "-").replace(",", "").replace(" ", "")
            try:
                val = abs(float(s))
            except ValueError:
                continue
            dp = len(s.split(".")[1]) if "." in s else 0
            hit = corpus.get(rnd(val, dp), set())
            rows.append(dict(
                version=version, section=section, line=ln, printed=raw + suf,
                value=val, decimals=dp, matched=bool(hit), n_sources=len(hit),
                # direct matches first, percent conversions after, so the four shown
                # are the ones that carry the quantity rather than a coincidence
                sources="; ".join(sorted(hit, key=lambda t: ("(as %)" in t, t))[:4]),
                context=re.sub(r"\s+", " ",
                               line[max(0, mo.start() - 90):mo.end() + 90]).strip(),
                note=NOTE))
C = pd.DataFrame(rows)
# the rows a reader flagged as absent from the tables carry, from here on, the file
# that does hold them (08_external_review/scripts/e14_setting_provenance.py)
_pv = f"{W}/08_external_review/tables/E14_setting_provenance.csv"
if os.path.exists(_pv):
    _p = pd.read_csv(_pv)[["No", "category", "value_in_the_file", "where_it_lives",
                           "status"]]
    _p = _p.rename(columns={"category": "provenance_category",
                            "value_in_the_file": "provenance_value",
                            "where_it_lives": "provenance_source",
                            "status": "provenance_status"})
    # the worksheet numbering is the row order of the clean version, so the
    # provenance is attached there only
    _clean = sorted(v for v in C.version.unique() if v.startswith("clean_"))[-1]
    C["No"] = C.groupby("version").cumcount() + 1
    C.loc[C.version != _clean, "No"] = -1
    C = C.merge(_p, on="No", how="left").drop(columns=["No"])

# values whose source of record is a particular part rather than whichever file the
# matcher happened to list first: the pointer is written into the provenance columns
# (08_external_review/tables/E39_manual_pointers.csv)
_mp = f"{W}/08_external_review/tables/E39_manual_pointers.csv"
if os.path.exists(_mp):
    _m = pd.read_csv(_mp)
    _clean = sorted(v for v in C.version.unique() if v.startswith("clean_"))[-1]
    _n = 0
    for _, _r in _m.iterrows():
        sel = ((C.version == _clean) & (C.section.astype(str) == str(_r.section))
               & (C.printed.astype(str) == str(_r.printed))
               & C.context.astype(str).str.contains(str(_r.context_contains), regex=False))
        if not sel.any():
            print(f"  [manual pointer] no row for {_r.section} / {_r.printed} / "
                  f"{_r.context_contains}")
            continue
        C.loc[sel, "provenance_category"] = "figure value, pointer recorded by hand"
        C.loc[sel, "provenance_value"] = str(_r.what)
        C.loc[sel, "provenance_source"] = (f"Supplementary_Tables_S1-S12.xlsx "
                                           f"{_r.sheet_cell} - {_r.part}")
        C.loc[sel, "provenance_status"] = str(_r.note)
        _n += int(sel.sum())
    print(f"manual pointers written into the provenance columns: {_n} rows")
C.to_csv(f"{B}/tables/manuscript_number_agreement.csv", index=False)
U = C[~C.matched].copy()

# where each unmatched number actually lives. None of these is a supplementary-table
# cell; the entry says which file or derivation carries it, so the check is still
# possible. Keyed on the string as printed in the manuscript.
_qc = pd.read_csv(f"{TAB}/28_qc_all_stages.csv").pivot_table(
    index=["cohort", "lesion"], columns="stage", values="n_cells")
_rm = (_qc["1_qc_pass"] - _qc["2_singlet"])
HOME = {
    "33 694": "参照配列の構築: 病変A の Cell Ranger アノテーションの feature 数。"
              "Materials and methods の記載値で、補足表には入っていない",
    "33 538": "参照配列の構築: 病変B-E の Cell Ranger 3.0.2/3.1.0 アノテーションの feature 数",
    "32 991": "参照配列の構築: 探索2系統のアノテーションで共通する Ensembl 遺伝子 ID 数",
    "206,723": "参照配列の構築: Ensembl 112 GRCh38 cDNA セットの転写産物数",
    "41,116": "参照配列の構築: t2g 表のバージョン付き遺伝子 ID 数。"
              "results/tables/23_chemistry_umi_comparison.csv の genes 列にも同じ値がある",
    "7,057": "参照配列の構築: Ensembl cDNA セットで同じ symbol を共有する ID 数",
    "90.6": "プログラム遺伝子集合の被覆率。コホート別の被覆率表にあり、S1-S12 には含めていない",
    "435,676,976": "SRR32119705 の処理リード数。logs/ の kb run_info から",
    "970": "ケミストリー検証で v2 として扱ったときの上位5,000バーコードの UMI 中央値。"
           "results/tables/23_chemistry_umi_comparison.csv",
    "406": "表S1 から導出できる: 病変ごとに 1_qc_pass - 2_singlet を計算すると "
           "{RANGE} の範囲になる",
    "2,000": "パイプラインの引数: コホートごとに選んだ highly variable genes の数",
    "15,820": "ingest 射影で共有された symbol 数。スクリプト 05 が logs/ に出力する",
}
HOME["406"] = HOME["406"].replace("{RANGE}", f"{int(_rm.min())}-{int(_rm.max())}")
HOME_EN = {
    "33 694": "reference build: features in the Cell Ranger annotation used for lesion A; "
              "a Materials and methods figure, not part of the supplementary tables",
    "33 538": "reference build: features in the Cell Ranger 3.0.2/3.1.0 annotation "
              "(lesions B-E)",
    "32 991": "reference build: Ensembl gene identifiers shared by the two discovery "
              "annotations",
    "206,723": "reference build: transcripts in the Ensembl 112 GRCh38 cDNA set",
    "41,116": "reference build: versioned gene identifiers in the t2g table; the same "
              "value is in the genes column of "
              "results/tables/23_chemistry_umi_comparison.csv",
    "7,057": "reference build: identifiers sharing a symbol in the Ensembl cDNA set",
    "90.6": "programme gene-set coverage, in the per-cohort coverage tables; not "
            "included in S1-S12",
    "435,676,976": "reads processed for SRR32119705, from the kb run_info in logs/",
    "970": "median UMIs of the top 5,000 barcodes when the chemistry test treated the "
           "run as v2; results/tables/23_chemistry_umi_comparison.csv",
    "406": "derivable from Table S1: 1_qc_pass minus 2_singlet per lesion gives the "
           "range {RANGE}",
    "2,000": "pipeline parameter: highly variable genes selected per cohort",
    "15,820": "symbols shared by the ingest projection, printed by script 05 to logs/",
}
HOME_EN["406"] = HOME_EN["406"].replace("{RANGE}", f"{int(_rm.min())}-{int(_rm.max())}")
U["where_the_value_lives"] = U.printed.map(lambda x: HOME.get(str(x).strip(), ""))
U["where_the_value_lives_en"] = U.printed.map(lambda x: HOME_EN.get(str(x).strip(), ""))
U["is_derivable_from_a_supplementary_table"] = U.printed.str.strip().eq("406")
_un = U[U.version == f"clean_v{VER}"]
if (_un.where_the_value_lives == "").any():
    log("unannotated unmatched numbers: "
        + ", ".join(_un[_un.where_the_value_lives == ""].printed))
U.to_csv(f"{B}/tables/manuscript_unmatched_numbers.csv", index=False)
g = C.groupby("version").agg(n=("printed", "size"), matched=("matched", "sum"))
g["pct"] = (100 * g.matched / g.n).round(1)
print(g.to_string())
print(U[U.version == f"clean_v{VER}"].groupby("section").size().to_string())
log("s4 done")
