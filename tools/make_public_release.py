#!/usr/bin/env python
"""Purpose: Assemble the public release of the apical periodontitis reanalysis - the code,
  the frozen inputs, the result tables and figures, the submission bundle and the records
  of the two checks - under ~/Documents/apscrna_public, with the processed objects kept
  separately for deposition, and scan the result for anything that should not be public.
Input: the working repository (APSC_ROOT, read only)
Output: ~/Documents/apscrna_public/ (+ zip), ~/Documents/apscrna_public_data/ (+ zip),
  and the scan report release_scan.csv inside the release

Nothing in the working repository is written to. Text files are copied with the home
directory of the machine replaced by a placeholder, and the name of the person who
checked the worksheet is replaced by their role; both substitutions are recorded in the
scan report and in the README of the release.
"""
import ast
import datetime
import hashlib
import io
import json
import tokenize
import importlib.metadata as md
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

import pandas as pd

SRC = os.environ.get("APSC_ROOT", os.path.expanduser("~/Documents/apscrna_v2"))
HOME = os.path.expanduser("~")
DEST = os.path.expanduser("~/Documents/apscrna_public")
DOI = "10.5281/zenodo.23277211"      # reserved on Zenodo for this release
DATA = os.path.expanduser("~/Documents/apscrna_public_data")
TITLE = ("Inflammatory activity poorly predicts tissue-forming transcription in "
         "stromal cells of human apical periodontitis")
JP = re.compile(r"[\u3040-\u30ff\u4e00-\u9fff]")
TEXT_EXT = {".py", ".sh", ".md", ".csv", ".tsv", ".json", ".txt", ".log", ".yml",
            ".cff", ".gmt"}
PERSON = re.compile(r"independent checker|independent checker|independent checker|independent checker")
SECRET = re.compile(r"(sk-[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{20,}|"
                    r"api[_-]?key\s*[:=]\s*['\"][^'\"]{12,})", re.I)
# a real home directory, not the bare string the scan itself looks for
ABSPATH = re.compile(r"/Users/[A-Za-z0-9._-]+/")

copied, notes = [], []
_subs = {"home": 0, "person": 0}


def sanitise(text, where=""):
    """the home directory and the checker's name do not travel into the release"""
    n = text.count(HOME)
    text = text.replace(HOME, "/path/to")
    text, m = PERSON.subn("independent checker", text)
    _subs["home"] += n
    _subs["person"] += m
    return text


TRANS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "comment_translations.json")
TRANS = (json.load(open(TRANS_PATH, encoding="utf8"))
         if os.path.exists(TRANS_PATH) else {})
JP = re.compile(r"[\u3040-\u30ff\u4e00-\u9fff]")
checks = []          # one row per file whose comments were translated


def _offsets(src):
    off, a = [0], 0
    for ln in src.splitlines(keepends=True):
        a += len(ln)
        off.append(a)
    return off


def _apply(src, spans):
    for a, b, new in sorted(spans, reverse=True):
        src = src[:a] + new + src[b:]
    return src


def translate_py(src):
    """replace the Japanese comments and docstrings, leaving the code untouched"""
    off, spans = _offsets(src), []
    for t in tokenize.generate_tokens(io.StringIO(src).readline):
        if t.type == tokenize.COMMENT and t.string in TRANS:
            a = off[t.start[0] - 1] + t.start[1]
            spans.append((a, a + len(t.string), TRANS[t.string]))
    n = len(spans)
    src = _apply(src, spans)

    off, spans = _offsets(src), []
    for node in ast.walk(ast.parse(src)):
        if not isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                                 ast.ClassDef)):
            continue
        doc = ast.get_docstring(node, clean=False)
        if doc is None or doc not in TRANS:
            continue
        e = node.body[0].value
        a = off[e.lineno - 1] + e.col_offset
        b = off[e.end_lineno - 1] + e.end_col_offset
        m = re.match(r'([A-Za-z]*)("""|\'\'\'|"|\')', src[a:b])
        q = m.group(2)
        if src[a + len(m.group(0)):b - len(q)] != doc:
            continue                      # escape sequences: leave the docstring alone
        new = TRANS[doc]
        if q in new:
            continue
        spans.append((a, b, m.group(0) + new + q))
    return _apply(src, spans), n + len(spans)


def translate_sh(src):
    out, n = [], 0
    for ln in src.split("\n"):
        s = ln.lstrip()
        if s.startswith("#") and s in TRANS:
            out.append(ln[:len(ln) - len(s)] + TRANS[s])
            n += 1
        elif "#" in ln and not JP.search(ln.split("#", 1)[0]) \
                and "#" + ln.split("#", 1)[1] in TRANS:
            out.append(ln.split("#", 1)[0] + TRANS["#" + ln.split("#", 1)[1]])
            n += 1
        else:
            out.append(ln)
    return "\n".join(out), n


def _ast_key(src):
    """the syntax tree with the docstrings emptied, so only the code is compared"""
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) \
                    and isinstance(first.value.value, str):
                first.value.value = ""
    return ast.dump(tree)


def _sh_check(orig, new):
    """the two scripts may differ only from a '#' onwards, and only where the original
    was Japanese: everything before that column must be identical"""
    a, c = orig.split("\n"), new.split("\n")
    if len(a) != len(c):
        return False
    for x, y in zip(a, c):
        if x == y:
            continue
        i = x.find("#")
        if i < 0 or x[:i] != y[:i] or not y[i:].startswith("#") \
                or not JP.search(x[i:]) or JP.search(y[i:]):
            return False
    return True


def _pptx_sig(f):
    from pptx import Presentation
    pr = Presentation(f)
    return ([[(s.shape_type, s.name) for s in sl.shapes] for sl in pr.slides],
            [[s.text_frame.text for s in sl.shapes if s.has_text_frame]
             for sl in pr.slides])


def _xlsx_sig(f):
    import openpyxl
    wb = openpyxl.load_workbook(f)
    return {ws.title: [[c.value for c in row] for row in ws.iter_rows()] for ws in wb}


PREV_DIR = tempfile.mkdtemp(prefix="apscrna_previous_build_")


def prev_bytes(rel):
    """the file of the previous build, before this one replaced it"""
    q = f"{PREV_DIR}/{rel}"
    return open(q, "rb").read() if os.path.exists(q) else None


def keep_if_same(dst, prev, sig):
    """a rewritten Office file is byte-different even when its content is identical;
    put the published bytes back so that rebuilding does not churn the repository"""
    if prev is None:
        return False
    try:
        if sig(io.BytesIO(prev)) == sig(dst):
            open(dst, "wb").write(prev)
            return True
    except Exception:
        pass
    return False


def strip_notes(src_path, dst_path):
    """drop the speaker notes, leaving the slides themselves as they are"""
    from pptx import Presentation
    before = Presentation(src_path)
    sig_before = [[sh.shape_type for sh in sl.shapes] for sl in before.slides]
    txt_before = [[sh.text_frame.text for sh in sl.shapes if sh.has_text_frame]
                  for sl in before.slides]
    dropped = 0
    for slide in before.slides:
        for rid, rel in list(slide.part.rels.items()):
            if rel.reltype.endswith("/notesSlide"):
                slide.part.drop_rel(rid)
                dropped += 1
    before.save(dst_path)
    after = Presentation(dst_path)
    same = ([[sh.shape_type for sh in sl.shapes] for sl in after.slides] == sig_before
            and [[sh.text_frame.text for sh in sl.shapes if sh.has_text_frame]
                 for sl in after.slides] == txt_before)
    left = sum(1 for sl in after.slides
               for rid, rel in sl.part.rels.items()
               if rel.reltype.endswith("/notesSlide"))
    return dropped, same, left


def put(rel_src, rel_dst=None, rename=None):
    """copy one file into the release, substituting the home path in text files"""
    s, d = f"{SRC}/{rel_src}", f"{DEST}/{rel_dst or rel_src}"
    if not os.path.exists(s):
        notes.append(f"missing in the working repository, not copied: {rel_src}")
        return
    os.makedirs(os.path.dirname(d), exist_ok=True)
    ext = os.path.splitext(s)[1].lower()
    if ext == ".csv" and os.path.getsize(s) > 1:
        try:
            t = pd.read_csv(s, dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            open(d, "w").write(sanitise(open(s, errors="replace").read(),
                                        rel_dst or rel_src))
            copied.append(rel_dst or rel_src)
            return
        drop = [c for c in t.columns
                if c.endswith("_ja") or f"{c}_en" in t.columns]
        if drop:
            t = t.drop(columns=drop)
            notes.append(f"Japanese-language columns dropped from "
                         f"{rel_dst or rel_src}: {', '.join(drop)}")
        t = t.map(lambda v: TRANS.get(v, v) if isinstance(v, str) else v)
        t = t.map(lambda v: sanitise(v, rel_dst or rel_src) if isinstance(v, str) else v)
        t.to_csv(d, index=False)
    elif ext in TEXT_EXT or ext == ".csv":
        text = sanitise(open(s, errors="replace").read(), rel_dst or rel_src)
        if ext in (".py", ".sh"):
            out, n = (translate_py(text) if ext == ".py" else translate_sh(text))
            if ext == ".py":
                ok = _ast_key(text) == _ast_key(out)
            else:
                ok = _sh_check(text, out)
            checks.append(dict(file=rel_dst or rel_src, kind=ext.lstrip("."),
                               units_translated=n, code_identical="Y" if ok else "N",
                               japanese_left=len(JP.findall(out))))
            assert ok, f"the code changed while translating {rel_src}"
            text = out
        elif ext == ".json":
            for k, v in TRANS.items():
                if k in text:
                    text = text.replace(k, v)
        open(d, "w").write(text)
    elif ext == ".pptx":
        prev = prev_bytes(rel_dst or rel_src)
        dropped, same, left = strip_notes(s, d)
        checks.append(dict(file=rel_dst or rel_src, kind="pptx",
                           units_translated=f"{dropped} notes pages removed",
                           code_identical="Y" if same and left == 0 else "N",
                           japanese_left=len(JP.findall(
                               subprocess.run(["unzip", "-p", d, "ppt/notesSlides/*"],
                                              capture_output=True,
                                              text=True).stdout))))
        assert same and left == 0, f"the slides changed while stripping notes: {rel_src}"
        keep_if_same(d, prev, _pptx_sig)
    else:
        shutil.copy2(s, d)
    copied.append(rel_dst or rel_src)


def put_tree(rel_dir, exts, rel_dst=None, skip=()):
    base = f"{SRC}/{rel_dir}"
    for root, dirs, files in os.walk(base):
        dirs[:] = [x for x in dirs if x not in ("__pycache__",) and
                   not any(s in os.path.join(root, x) for s in skip)]
        for f in sorted(files):
            if os.path.splitext(f)[1].lower() not in exts or f.endswith(".bak"):
                continue
            rel = os.path.relpath(os.path.join(root, f), base)
            put(f"{rel_dir}/{rel}", f"{(rel_dst or rel_dir)}/{rel}")


# ---------------------------------------------------------------- the release
if os.path.exists(DEST):
    for e in sorted(os.listdir(DEST)):
        if e.startswith("."):          # .git, .gitattributes: the user's repository
            continue
        shutil.move(f"{DEST}/{e}", f"{PREV_DIR}/{e}")   # kept for keep_if_same
os.makedirs(DEST, exist_ok=True)

put_tree("scripts", {".py", ".sh"})
put_tree("frozen", {".json", ".csv", ".sha256"})
for f in ("chemistry_decision.json", "discovery_runs.tsv", "run_donor_map_v2.csv",
          "ena_filereport.tsv", "gene_id_to_symbol_ens112.tsv"):
    put(f"ref/{f}", f"reference_records/{f}")
put_tree("results/tables", {".csv"})
put_tree("results/figures_ms", {".pdf", ".png", ".json"},
         skip=("figures_draft", "superseded"))
STAGE_README = {
 "05_post_audit": "# Analyses added after the independent audit\n\n"
   "`scripts/q*.py` and the tables they write. These answer the points the audit "
   "raised: the correlation matrices under every score convention, the random-factor "
   "null for candidate reproduction, the baselines for the retained dispersion, the "
   "donor-mixture check and the provenance of the numbers reported in the paper.\n",
 "06_submission_prep": "# Figure and table preparation for submission\n\n"
   "`scripts/r*.py` redraw Figures 3 and 4 in their final form, build the "
   "post-audit supplementary tables S11 and S12, reconcile the manuscript numbers "
   "with the tables, verify the chemistry of the validation runs and install the "
   "final figures into `results/figures_ms/`.\n",
 "08_external_review": "# External checks of the reported numbers\n\n"
   "`scripts/e*.py` and the E-series tables. Two independent checks were run over "
   "every number in the manuscript: an automated one and a cell-by-cell check by a "
   "person. The scripts adjudicate the points each check raised, add the parts the "
   "supplementary tables were missing, tie every parameter and software version to "
   "the file that fixes it, recover a base-composition measurement that had not been "
   "kept, verify the figure and table legends against the code, and assemble the "
   "final checking record in `07_submission_bundle/checks/`.\n",
}

for stage in ("05_post_audit", "06_submission_prep", "08_external_review"):
    put_tree(f"{stage}/tables", {".csv"})
    put_tree(f"{stage}/scripts", {".py"})
    put_tree(f"{stage}/logs", {".log"})
    # the working notes of each stage are kept in Japanese and are not published;
    # a short English note takes their place
    open(f"{DEST}/{stage}/README.md", "w").write(STAGE_README[stage])
    copied.append(f"{stage}/README.md")
BUNDLE = ["Supplementary_Tables_S1-S12.xlsx", "supplementary_table_index.csv",
          "manuscript_number_agreement.csv", "manuscript_unmatched_numbers.csv",
          "manuscript_source_files.csv", "file_to_manuscript_map.csv",
          "figure_panel_n.csv", "caption_number_check.csv", "pptx_layout_checks.csv"]
for f in BUNDLE:
    put(f"07_submission_bundle/tables/{f}")
put_tree("07_submission_bundle/figures/pptx", {".pptx"})
put("07_submission_bundle/figures/figures_editable.pptx")
put("07_submission_bundle/figures/colour_keys.csv")
put_tree("07_submission_bundle/scripts", {".py"})
open(f"{DEST}/07_submission_bundle/README.md", "w").write(
    "# Submission bundle\n\n"
    "`tables/Supplementary_Tables_S1-S12.xlsx` is the supplementary workbook: a "
    "`Legends` sheet, an `index` sheet and one sheet per table, each sheet holding its "
    "parts one after another with the column definitions in row 2. "
    "`supplementary_table_index.csv` lists every part and the file it comes from; "
    "`table_and_figure_provenance.csv` in the repository root adds the script that "
    "writes each one. `manuscript_number_agreement.csv` is the machine comparison of "
    "every number in the manuscript against the bundle, `manuscript_unmatched_numbers."
    "csv` the ones it cannot place and where they live instead. `figures/pptx/` holds "
    "one editable slide per figure, built from the text-free renderings and the "
    "position metadata. `checks/check_worksheet_final.xlsx` is the record of the "
    "cell-by-cell check.\n")
copied.append("07_submission_bundle/README.md")
put_tree("logs", {".log"})

# the checking record, with the checker named by role rather than by name
import openpyxl

src_x = f"{SRC}/07_submission_bundle/checks/check_worksheet_final.xlsx"
wb = openpyxl.load_workbook(src_x)
n_anon = 0
for sheet in wb:
    for row in sheet.iter_rows():
        for c in row:
            if isinstance(c.value, str) and PERSON.search(c.value):
                c.value = PERSON.sub("independent checker", c.value)
                n_anon += 1
os.makedirs(f"{DEST}/07_submission_bundle/checks", exist_ok=True)
dst_x = f"{DEST}/07_submission_bundle/checks/check_worksheet_final.xlsx"
prev_x = prev_bytes("07_submission_bundle/checks/check_worksheet_final.xlsx")
wb.save(dst_x)
keep_if_same(dst_x, prev_x, _xlsx_sig)
copied.append("07_submission_bundle/checks/check_worksheet_final.xlsx")
notes.append(f"the checker's name was replaced by their role in {n_anon} cells of "
             "check_worksheet_final.xlsx")

# ---------------------------------------------------------------- environment
PKGS = ["scanpy", "anndata", "decoupler", "omnipath", "harmonypy", "scrublet",
        "leidenalg", "diptest", "numpy", "scipy", "pandas", "matplotlib",
        "statsmodels", "openpyxl", "python-pptx", "h5py", "pyarrow", "kb-python"]
vers = {}
for p in PKGS:
    try:
        vers[p] = md.version(p)
    except md.PackageNotFoundError:
        pass
BIN = {}
for name, cmd in (("kallisto", ["kallisto", "version"]),
                  ("bustools", ["bustools", "version"]),
                  ("sra-tools", ["fasterq-dump", "--version"])):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=30).stdout
        m = re.search(r"(\d+\.\d+\.\d+)", out)
        if m:
            BIN[name] = m.group(1)
    except (OSError, subprocess.SubprocessError):
        pass
py = ".".join(sys.version.split()[0].split(".")[:3])
open(f"{DEST}/requirements.txt", "w").write(
    "\n".join(f"{k}=={v}" for k, v in sorted(vers.items())) + "\n")
open(f"{DEST}/environment.yml", "w").write(
    "name: apscrna\nchannels:\n  - conda-forge\n  - bioconda\ndependencies:\n"
    f"  - python={py}\n"
    + "".join(f"  - {k}={v}\n" for k, v in sorted(BIN.items()))
    + "  - pip\n  - pip:\n"
    + "".join(f"      - {k}=={v}\n" for k, v in sorted(vers.items())))
copied += ["requirements.txt", "environment.yml"]

YEAR = datetime.date.today().year
open(f"{DEST}/LICENSE", "w").write(f"""MIT License

Copyright (c) {YEAR} The authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
""")
open(f"{DEST}/LICENSE-DATA", "w").write(
    "The tables, figures and derived data files in this repository "
    "(results/, 05_post_audit/, 06_submission_prep/, 08_external_review/, "
    "07_submission_bundle/ and the deposited processed objects) are released under "
    "the Creative Commons Attribution 4.0 International licence (CC BY 4.0), "
    "https://creativecommons.org/licenses/by/4.0/.\n\n"
    "The code (scripts/ and the scripts/ folders of the staged analyses) is released "
    "under the MIT licence; see LICENSE.\n\n"
    "The primary sequencing data are not redistributed here; they remain under the "
    "terms of their original repositories (GEO and SRA, see README.md).\n")
copied += ["LICENSE", "LICENSE-DATA"]
# ------------------------------------------------- where each table and figure is made
IDX = pd.read_csv(f"{SRC}/07_submission_bundle/tables/supplementary_table_index.csv")
SI = pd.read_csv(f"{SRC}/results/tables/31_script_index.csv")
STAGED = {  # the files written by the staged analyses rather than by scripts/
    "S3_rho_by_score_convention.csv": "06_submission_prep/scripts/r1_tables.py",
    "S11a_random_factor_null.csv": "06_submission_prep/scripts/r3_record_and_S11_S12.py",
    "S11b_overlap_removed_candidates.csv":
        "06_submission_prep/scripts/r3_record_and_S11_S12.py",
    "S11c_dispersion_baselines.csv":
        "06_submission_prep/scripts/r3_record_and_S11_S12.py",
    "S11c_null_factor_selection.csv":
        "06_submission_prep/scripts/r3_record_and_S11_S12.py",
    "S12a_single_donor_check.csv": "06_submission_prep/scripts/r3_record_and_S11_S12.py",
    "S12b_number_provenance.csv": "06_submission_prep/scripts/r3_record_and_S11_S12.py",
    "C6_validation_chemistry_check.csv":
        "06_submission_prep/scripts/r5_validation_chemistry.py",
    "C6_base_composition_by_position.csv":
        "06_submission_prep/scripts/r5_validation_chemistry.py",
}


def maker(source_file):
    if source_file in STAGED:
        return STAGED[source_file]
    if source_file.startswith("E"):
        n = re.match(r"E(\d+)", source_file)
        hit = [f for f in os.listdir(f"{SRC}/08_external_review/scripts")
               if f.startswith(f"e{n.group(1)}_")] if n else []
        return f"08_external_review/scripts/{hit[0]}" if hit else "08_external_review"
    stem = os.path.splitext(source_file)[0]
    for _, r in SI.iterrows():
        if stem in str(r.tables):
            return f"scripts/{r.script}"
    return ""


PROV = IDX.assign(made_by=[maker(f) for f in IDX.source_file])
PROV.to_csv(f"{DEST}/table_and_figure_provenance.csv", index=False)
copied.append("table_and_figure_provenance.csv")
rows = []
for t, g in PROV.groupby("table", sort=False):
    rows.append(f"| {t} | {'; '.join(sorted(set(g.source_file)))} | "
                f"{'; '.join(sorted(set(x for x in g.made_by if x)))} |")
TABLE_MAP = "\n".join(rows)

FIGS = []
for _, r in SI.iterrows():
    if isinstance(r.figures, str) and r.figures.strip():
        for f in str(r.figures).split(";"):
            if f.strip():
                FIGS.append((f.strip(), f"scripts/{r.script}"))
FIGS += [("Fig3_tf_screen (final version)", "06_submission_prep/scripts/r2_fig3.py"),
         ("Fig4_coupling_coexistence (final version)",
          "06_submission_prep/scripts/r2_fig4.py"),
         ("FigS4_matrix_vs_requantification (axis label corrected)",
          "08_external_review/scripts/e27_figS4_relabel.py")]
FIG_MAP = "\n".join(f"| {a} | {b} |" for a, b in sorted(set(FIGS)))

# ---------------------------------------------------------------- the data folder
os.makedirs(DATA, exist_ok=True)
DATA_FILES = ["all_cells_discovery.h5ad", "all_cells_validation.h5ad",
              "stromal_discovery.h5ad", "stromal_discovery_labelled.h5ad",
              "stromal_validation.h5ad", "all_cells_discovery_obs_scored.parquet",
              "stromal_discovery_tf_activity.parquet",
              "stromal_validation_tf_activity.parquet"]
WHAT = {
    "all_cells_discovery.h5ad": "all 44,433 discovery cells after QC and doublet "
                                "removal, with compartment labels, Leiden clusterings "
                                "and QC covariates; the programme scores for these "
                                "cells are in all_cells_discovery_obs_scored.parquet",
    "all_cells_validation.h5ad": "all 26,090 validation cells, same processing",
    "stromal_discovery.h5ad": "the 6,215 discovery stromal cells",
    "stromal_discovery_labelled.h5ad": "the discovery stromal cells with the four state "
                                       "labels used in Figure 4",
    "stromal_validation.h5ad": "the 1,779 validation stromal cells",
    "all_cells_discovery_obs_scored.parquet": "the per-cell table of the discovery "
                                              "cohort (scores, labels, covariates)",
    "stromal_discovery_tf_activity.parquet": "CollecTRI activities of the discovery "
                                             "stromal cells",
    "stromal_validation_tf_activity.parquet": "CollecTRI activities of the validation "
                                              "stromal cells",
}
drows = []
for f in DATA_FILES:
    s = f"{SRC}/data/processed/{f}"
    if not os.path.exists(s):
        notes.append(f"processed object missing, not deposited: {f}")
        continue
    dst_data = f"{DATA}/{f}"
    if not (os.path.exists(dst_data)
            and os.path.getsize(dst_data) == os.path.getsize(s)):
        shutil.copy2(s, dst_data)
    h = hashlib.sha256()
    with open(s, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    drows.append(dict(file=f, size_bytes=os.path.getsize(s),
                      size_GB=round(os.path.getsize(s) / 1e9, 3),
                      sha256=h.hexdigest(), contents=WHAT.get(f, "")))
DD = pd.DataFrame(drows)
DD.to_csv(f"{DATA}/manifest.csv", index=False)
DATA_MAP = "\n".join(f"| {r.file} | {r.size_GB} | {r.contents} |"
                     for r in DD.itertuples())
total_gb = round(DD.size_bytes.sum() / 1e9, 2)
open(f"{DATA}/README.md", "w").write(
    f"# Processed objects for *{TITLE}*\n\n"
    f"{len(DD)} files, {total_gb} GB in total, with sha256 sums in `manifest.csv`. "
    "They are the processed objects the analysis scripts read; the raw reads are in GEO "
    "and SRA (see the code repository). The requantification variants "
    "(`*_kb.h5ad`, `*_kbgeo.h5ad`, `*_kbtest.h5ad`) are not deposited because they are "
    "rebuilt from the raw reads by `scripts/08_discovery_requant.sh` and "
    "`scripts/09_geo_vs_kb.py`.\n\n"
    f"| file | GB | contents |\n|---|---|---|\n{DATA_MAP}\n")

README = f"""# {TITLE}

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

Python {py}; see `requirements.txt` and `environment.yml` for the package versions the
analysis was run with ({', '.join(f'{k} {v}' for k, v in sorted(vers.items())[:6])} and
the rest), and {', '.join(f'{k} {v}' for k, v in sorted(BIN.items())) or 'kallisto, bustools and the SRA toolkit'} for read processing.

## Where each supplementary table comes from

| table | source files | written by |
|---|---|---|
{TABLE_MAP}

Full list, part by part, in `table_and_figure_provenance.csv`.

## Where each figure comes from

| figure | written by |
|---|---|
{FIG_MAP}

## Data

The processed objects are deposited separately on Zenodo under DOI {DOI}
({total_gb} GB, {len(DD)} files; see the `manifest.csv` and `README.md` of the deposit).
They are the AnnData objects and per-cell tables the scripts read. The raw reads are not redistributed here.

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

> The authors. *{TITLE}*. Zenodo, DOI: {DOI}.

Please cite the paper as well once it is published.
"""
open(f"{DEST}/README.md", "w").write(README)
copied.append("README.md")
os.makedirs(f"{DEST}/tools", exist_ok=True)
open(f"{DEST}/tools/make_public_release.py", "w").write(
    sanitise(open(os.path.abspath(__file__)).read(), "tools/make_public_release.py"))
copied.append("tools/make_public_release.py")
if os.path.exists(TRANS_PATH):
    shutil.copy2(TRANS_PATH, f"{DEST}/tools/comment_translations.json")
    copied.append("tools/comment_translations.json")
CK = pd.DataFrame(checks)
CK.to_csv(f"{DEST}/translation_check.csv", index=False)
copied.append("translation_check.csv")

# ---------------------------------------------------------------- the scan
scan = []
for root, dirs, files in os.walk(DEST):
    dirs[:] = [d for d in dirs
                   if d != "__pycache__" and not d.startswith(".")]
    for f in sorted(files):
        p = os.path.join(root, f)
        rel = os.path.relpath(p, DEST)
        ext = os.path.splitext(f)[1].lower()
        jp = home = person = secret = 0
        if ext in TEXT_EXT:
            t = open(p, errors="replace").read()
            jp = len(JP.findall(t))
            home = len(ABSPATH.findall(t))
            person = len(PERSON.findall(t))
            secret = len(SECRET.findall(t))
        elif ext in (".xlsx", ".pptx"):
            with zipfile.ZipFile(p) as z:
                t = " ".join(z.read(n).decode("utf8", "replace")
                             for n in z.namelist() if n.endswith((".xml", ".rels")))
            jp = len(JP.findall(t))
            home = len(ABSPATH.findall(t))
            person = len(PERSON.findall(t))
        scan.append(dict(file=rel, size_bytes=os.path.getsize(p), japanese_chars=jp,
                         absolute_home_paths=home, personal_names=person,
                         secret_like=secret,
                         filename_has_personal_name=bool(PERSON.search(f)),
                         filename_has_japanese=bool(JP.search(f))))
SC = pd.DataFrame(scan)
SC.to_csv(f"{DEST}/release_scan.csv", index=False)
pd.DataFrame([dict(note=n) for n in notes]).to_csv(f"{DEST}/release_notes.csv",
                                                   index=False)

# ---------------------------------------------------------------- the archives
repo_zip = os.path.expanduser("~/Documents/apscrna_public_v1.0.zip")
zips = [repo_zip]
if os.environ.get("SKIP_DATA_ZIP") != "1":   # the deposit archive is rebuilt below
    zips.append(os.path.expanduser("~/Documents/apscrna_public_data_v1.0.zip"))
for z in zips:
    if os.path.exists(z):
        os.remove(z)
with zipfile.ZipFile(repo_zip, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for root, dirs, files in os.walk(DEST):
        dirs[:] = [d for d in dirs
                   if d != "__pycache__" and not d.startswith(".")]
        files[:] = [f for f in files if not f.startswith(".")]
        for f in sorted(files):
            p = os.path.join(root, f)
            z.write(p, os.path.join("apscrna_public", os.path.relpath(p, DEST)))

data_zip = os.path.expanduser("~/Documents/apscrna_public_data_v1.0.zip")
if os.environ.get("SKIP_DATA_ZIP") != "1":
    with zipfile.ZipFile(data_zip, "w", zipfile.ZIP_DEFLATED, compresslevel=1,
                         allowZip64=True) as z:
        for f in sorted(os.listdir(DATA)):
            z.write(f"{DATA}/{f}", f"apscrna_public_data/{f}")
    print(f"deposit zip: {os.path.getsize(data_zip)/1e9:.2f} GB  {data_zip}")

print(f"files in the release: {len(SC)}  ({SC.size_bytes.sum()/1e6:.1f} MB)")
print(f"release zip: {os.path.getsize(repo_zip)/1e6:.1f} MB  {repo_zip}")
print(f"deposit: {len(DD)} files, {total_gb} GB in {DATA}")
print("\nscan:")
print(f"  files with Japanese text        : "
      f"{int((SC.japanese_chars > 0).sum())} of {len(SC)}")
print(f"  files with an absolute home path: {int((SC.absolute_home_paths > 0).sum())}")
print(f"  files mentioning a person       : {int((SC.personal_names > 0).sum())}")
print(f"  files with a secret-like string : {int((SC.secret_like > 0).sum())}")
print(f"  file names with a personal name : {int(SC.filename_has_personal_name.sum())}")
print(f"  file names with Japanese        : {int(SC.filename_has_japanese.sum())}")
for c in ("japanese_chars", "absolute_home_paths", "personal_names"):
    top = SC[SC[c] > 0].sort_values(c, ascending=False).head(8)
    if len(top):
        print(f"\n  {c}:")
        for r in top.itertuples():
            print(f"    {getattr(r, c):6d}  {r.file}")

shutil.rmtree(PREV_DIR, ignore_errors=True)
