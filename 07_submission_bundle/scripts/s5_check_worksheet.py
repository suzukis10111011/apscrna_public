#!/usr/bin/env python
"""Submission bundle: a check-by-hand worksheet for a second reader.

One row per number in the clean manuscript. Each row carries the sentence context,
the value as printed, and the address (sheet and cell) in
Supplementary_Tables_S1-S12.xlsx where that value can be read, plus blank columns for
the reader to write what they read, whether it agrees, and any remark.

Addresses come from a full scan of the workbook: every data cell is indexed by its
value rounded to 0-4 decimals, and numbers inside text cells are indexed too. A
manuscript number is matched at the precision it is printed with. Values that occur
in more than five cells (small integers such as 20 or 100) get no address - pointing
at one arbitrary cell would be worse than saying the check is a reading of the
sentence, not a cell lookup.

The eleven numbers that no supplementary-table cell carries go on their own sheet
with the file and row where each can be checked instead.

Output: 07_submission_bundle/tables/check_worksheet.xlsx
"""
import decimal
import os
import re
import sys

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

W = "/path/to/Documents/apscrna_v2"
B = f"{W}/07_submission_bundle"
XLSX = f"{B}/tables/Supplementary_Tables_S1-S12.xlsx"
MAX_CELLS = 5

AG = pd.read_csv(f"{B}/tables/manuscript_number_agreement.csv")
VER = sorted({v for v in AG.version if v.startswith("clean_")})[-1]
A = AG[AG.version == VER].reset_index(drop=True)
UN = pd.read_csv(f"{B}/tables/manuscript_unmatched_numbers.csv")
U = UN[UN.version == VER].reset_index(drop=True)
TXT = open(f"{B}/manuscript_text/manuscript_{VER}.txt").read()

# references the external check corrected; applied after the automatic assignment so
# the worksheet points at the cell that carries the quantity, not one that happens to
# hold the same number (08_external_review/scripts/e2_reference_corrections.py)
FIXES = {}
for _fx in (f"{W}/08_external_review/tables/E2_worksheet_reference_corrections.csv",
            f"{W}/08_external_review/tables/E9_reference_corrections.csv",
            f"{W}/08_external_review/tables/E11_reference_corrections.csv"):
    if not os.path.exists(_fx):
        continue
    _f = pd.read_csv(_fx)
    _f = _f[_f.action.str.startswith(("applied", "cleared"))]
    if "old_reference" not in _f.columns:
        _f = _f.assign(old_reference="")
    # a later file wins, so the second round's targets replace the first round's
    FIXES.update({int(r.No): r for _, r in _f.iterrows()})

# A table value of 39.65 is printed as 39.7 in the manuscript, which is round-half-up;
# Python's round() is half-to-even and would register it under 39.6, so the manuscript
# value would never find its cell. All rounding for the index and the lookup goes
# through here.
def rnd(v, dp):
    q = decimal.Decimal(1).scaleb(-dp)
    return float(decimal.Decimal(repr(float(v))).quantize(
        q, rounding=decimal.ROUND_HALF_UP))


# ------------------------------------------------------------------- wording
# Everything the reader sees is held here so the same construction can be written in
# Japanese or in English; only the manuscript context, the values and the cell
# addresses are language-independent and copied through unchanged.
LANG = sys.argv[1] if len(sys.argv) > 1 else "ja"
assert LANG in ("ja", "en"), LANG

STR = {
 "ja": {
  "sheet_howto": "手順", "sheet_main": "補足表で確認", "sheet_other": "補足表にない値",
  "how_cell": "セル照合", "how_col": "列で照合", "how_sent": "文として確認",
  "c_no": "No", "c_section": "原稿の節", "c_context": "原稿中の文脈（前後10語）",
  "c_hcr": "human_check_required",
  "c_printed": "原稿の値", "c_how": "確認方法", "c_file": "参照先ファイル",
  "c_sheet": "シート名", "c_cell": "セル番地", "c_col": "列名", "c_tabval": "表の値",
  "c_ncand": "候補セル数", "c_note": "補足", "c_read": "自分で読んだ値",
  "c_agree": "一致（Y/N）", "c_remark": "備考", "c_where": "どこで確認するか",
  "c_derivable": "補足表から導出できるか", "c_item": "項目", "c_value": "内容",
  "yes": "はい", "no": "いいえ",
  "k_lesion": "病変名", "k_cohort": "コホート", "k_TF": "TF名", "k_gene": "遺伝子名",
  "k_compartment": "区画名", "k_axis": "軸", "k_cell_set": "細胞集合",
  "k_version": "活性版", "k_pool": "母集団", "k_stage": "QC段階", "k_donor": "ドナー",
  "k_run_accession": "ラン", "k_comparison": "比較", "k_set_name": "遺伝子集合",
  "sep": "・",
  "m_more": " ほか{n}件",
  "m_ties": "文脈の語と同じだけ合う候補が他にもある: {list}",
  "m_best_of_n": "同じ値のセルはブック内に {n} 個あるが、列名が文脈に最も合うのはこの1つ",
  "m_loose": "原稿と表で印字精度が違う。丸めて一致するかを見る。",
  "m_col": "この列で値を探す{rows}。候補の列: {cols}（該当セル {n} 個）",
  "m_col_rows": "（該当する行: {hint}）",
  "m_pairs": "文脈に近い順の候補列: {cols}",
  "m_narrowed": "文脈の手がかり（{keys}）で {n} 個の候補から絞った。",
  "m_narrowed_done": "文脈の手がかり（{keys}）で絞込み済み。",
  "m_absent": "補足表のブックには同じ値のセルが無い。原稿照合表では次のファイルで一致している: {src}",
  "m_absent_head": "補足表のブックには同じ値のセルが無い",
  "m_nosrc": "(記録なし)",
  "m_pct": "原稿照合表での一致は相関係数をパーセントに読み替えたもので、単位が違う可能性がある。文脈に合う表を探して確認する。",
  "c_note2": "備考（自動）",
  "m_nohuman": "人の確認は不要（human_check_required = N）。{ev}。",
  "m_cleared": "外部照合で、ここに指定していたセル（{old}）は別の量を指していると分かった。{why}。{how}",
  "m_howelse": "表で確認するには: {how}。",
  "m_fixed": "外部照合（08_external_review/tables/E2）で参照先を修正した。修正前の指定: {old}。",
  "w_param": "閾値・クラスタリングの設定値・式の中の定数で、表が持つ量ではない。文の記述が方法の記載と合っているかを読んで確認する。",
  "w_range": "複数行にまたがる値（範囲・平均・合計など）で、1つのセルに対応しない。候補セル {n} 個、{p} 通りの列。",
  "w_nokeys": "前後10語に病変名・コホート・TF名・区画名が現れず、行を絞る手がかりが無い。候補セル {n} 個、{p} 通りの列。",
  "w_many": "文脈の手がかりで絞っても候補が多く、列名からも1つに決められない。候補セル {n} 個、{p} 通りの列。",
  "title": "原稿の数値を補足表で突き合わせるための確認シート",
  "subtitle": "1行が原稿中の数値1つに対応する。黄色の3列（自分で読んだ値・一致（Y/N）・備考）だけを記入する。それ以外の列は書き換えない。",
  "steps": [
   "1. このブックと同じフォルダにある Supplementary_Tables_S1-S12.xlsx を別ウィンドウで開く。",
   "2. まず「human_check_required」列を見る。Y の行だけを確認すればよい。N の行は"
   "ソフトウェアの版・パイプラインの引数・図番号・参照配列の構築値・凍結設定値・"
   "名称の一部である数字で、凍結ファイルまたはスクリプトの該当行で確認済みなので、"
   "手で確認する必要はない（確認済みの場所は「補足」に書いてある）。",
   "3. 「補足表で確認」シートの human_check_required が Y の行を、上から順に見る。",
   "4. 「シート名」と「セル番地」の指すセルを補足表のブックで開く。たとえば S3 と H7 なら、"
   "S3 シートの H7 セル。各シートは1行目が列名、2行目が列の定義、3行目からがデータ。",
   "5. そのセルに表示されている値を「自分で読んだ値」に書き写す。四捨五入せず見えたまま書く。",
   "6. 「原稿の値」と見比べて、同じなら「一致（Y/N）」に Y、違えば N を入れる。"
   "原稿が丸めている場合（表 0.9759 → 原稿 0.98 など）は、丸めて一致するなら Y とし、"
   "「備考」に「丸め」と書く。",
   "7. 「確認方法」が「列で照合」の行は、セルを1つに絞れない。指定のシートの指定の列を"
   "フィルタして該当行の値を読む。「補足」に「該当する行」として、セル番地と"
   "その行の病変名・コホート・TF名などを挙げてあるので、文脈に合う行を選ぶ。",
   "8. 「確認方法」が「文として確認」の行は、セルを指定できない。理由が「補足」の先頭に"
   "書いてある（範囲や平均などで複数行にまたがる／閾値や設定値で表の量ではない／"
   "文脈に行を絞る手がかりが無い／手がかりで絞っても候補が多い）。"
   "文の述べている内容が表と矛盾しないかを読んで判断し、「自分で読んだ値」には"
   "確認に使った表の値を書く。",
   "9. 「補足」が「補足表のブックには同じ値のセルが無い」で始まる行は、補足表の外の"
   "ファイルで一致している。そのファイル名が書いてあるので、見られない場合は"
   "「備考」に「要確認」と書いて残す。",
   "10. 「補足表にない値」シートは、補足表のセルでは確認できない{u_n}件。"
   "「どこで確認するか」に書いてあるファイルや導出方法で確認する。",
   "11. N を付けた行、判断に迷った行は「備考」に理由を書いて残す。空欄のままにしない。",
   "12. 記入が終わったら集計する。「補足表で確認」シートの1行目でオートフィルタを開き、"
   "「一致（Y/N）」で N だけを絞って行数を数える（画面左下の「○○レコード中○○個」で分かる）。"
   "同じ要領で Y の件数、空欄の件数も数え、N と空欄の行を一覧にして報告する。"
   "このブックには数式を入れていないので、件数は自動では出ない。",
  ],
  "meta_head": "（以下は作成時点の内訳。記入欄の集計ではない）",
  "meta_ms": "照合対象の原稿", "meta_ms_v": "manuscript/{v} の docx（清書版）",
  "meta_n": "原稿中の数値", "meta_matched": "補足表のセルで確認できる",
  "meta_other": "補足表にない", "meta_other_v": "{n} 件（別シート）",
  "meta_cell": "うちセル番地を指定できた（確認方法=セル照合）",
  "meta_col": "シートと列を指定できた（確認方法=列で照合）",
  "meta_absent": "うち補足表のブックに同じ値のセルが無い（別ファイルで一致）",
  "meta_upg": "うち文脈の手がかりで格上げできた（文→列・文→セル・列→セル）",
  "meta_hcr_y": "うち human_check_required = Y",
  "meta_hcr_n": "うち human_check_required = N（凍結ファイル／スクリプトで確認済み）",
  "meta_sent": "セル指定ができず文として確認",
  "meta_src": "作成元", "unit": " 件",
 },
 "en": {
  "sheet_howto": "How to use", "sheet_main": "Check against tables",
  "sheet_other": "Not in the tables",
  "how_cell": "cell lookup", "how_col": "column lookup", "how_sent": "read the sentence",
  "c_no": "No", "c_section": "Manuscript section",
  "c_context": "Context in the manuscript (10 words each side)",
  "c_hcr": "human_check_required",
  "c_printed": "Value in the manuscript", "c_how": "How to check",
  "c_file": "File to open", "c_sheet": "Sheet", "c_cell": "Cell", "c_col": "Column",
  "c_tabval": "Value in the table", "c_ncand": "Candidate cells", "c_note": "Notes",
  "c_read": "Value you read", "c_agree": "Agrees (Y/N)", "c_remark": "Your remarks",
  "c_where": "Where to check it", "c_derivable": "Derivable from a table",
  "c_item": "Item", "c_value": "Value", "yes": "yes", "no": "no",
  "k_lesion": "lesion", "k_cohort": "cohort", "k_TF": "TF", "k_gene": "gene",
  "k_compartment": "compartment", "k_axis": "axis", "k_cell_set": "cell set",
  "k_version": "activity version", "k_pool": "null pool", "k_stage": "QC stage",
  "k_donor": "donor", "k_run_accession": "run", "k_comparison": "comparison",
  "k_set_name": "gene set", "sep": ", ",
  "m_more": " and {n} more",
  "m_ties": "other candidates fit the context words equally well: {list}",
  "m_best_of_n": "the workbook holds {n} cells with this value; this is the one whose column name fits the context best",
  "m_loose": "the manuscript and the table print this at different precision; check that it agrees after rounding. ",
  "m_col": "look for the value in this column{rows}. Candidate columns: {cols} ({n} cells carry the value)",
  "m_col_rows": " (rows to look at: {hint})",
  "m_pairs": "candidate columns, closest to the context first: {cols}",
  "m_narrowed": "narrowed from {n} candidates using the identifiers named in the sentence ({keys}). ",
  "m_narrowed_done": "already narrowed using the identifiers named in the sentence ({keys}). ",
  "m_absent": "no cell in the supplementary workbook carries this value; the manuscript agreement table matches it in: {src}",
  "m_absent_head": "no cell in the supplementary workbook carries this value",
  "m_nosrc": "(not recorded)",
  "m_pct": "the agreement table matched this by reading a correlation coefficient as a percentage, so the units may differ; find the table that fits the context and check there. ",
  "c_note2": "Notes (automatic)",
  "m_nohuman": "No hand check needed (human_check_required = N). {ev}. ",
  "m_cleared": "the external check found that the cell named here ({old}) holds another quantity: {why}. {how}",
  "m_howelse": "To check it against the tables: {how}. ",
  "m_fixed": "reference corrected after the external check (08_external_review/tables/E2); it previously pointed at {old}. ",
  "w_param": "a threshold, a software setting or a constant inside a formula, not a quantity the tables carry; read it against the methods text. ",
  "w_range": "a value spanning several rows (a range, mean or total), so it has no single cell; {n} candidate cells in {p} columns. ",
  "w_nokeys": "the 10 words each side name no lesion, cohort, TF or compartment, so there is nothing to narrow the rows with; {n} candidate cells in {p} columns. ",
  "w_many": "even after narrowing on the identifiers in the sentence the candidates are many, and the column names do not single one out; {n} candidate cells in {p} columns. ",
  "title": "Worksheet for checking the manuscript numbers against the supplementary tables",
  "subtitle": "One row is one number in the manuscript. Fill in only the three yellow columns (Value you read, Agrees (Y/N), Your remarks) and leave the rest unchanged.",
  "steps": [
   "1. Open Supplementary_Tables_S1-S12.xlsx, in the same folder as this workbook, in a second window.",
   "2. Start from the 'human_check_required' column: only the Y rows need checking. "
   "An N row is a software version, a pipeline argument, a figure number, a "
   "reference-sequence construction value, a frozen setting or a numeral that is part "
   "of a name; each has already been checked against the line of the frozen file or "
   "script named in 'Notes', so it needs no hand check.",
   "3. Work down the Y rows of the 'Check against tables' sheet from the first one.",
   "4. Go to the cell named by 'Sheet' and 'Cell' in the supplementary workbook; S3 and H7 means "
   "cell H7 of sheet S3. In every sheet row 1 holds the column names, row 2 the column "
   "definitions, and the data start at row 3.",
   "5. Copy the value shown in that cell into 'Value you read', exactly as it appears, without rounding.",
   "6. Compare it with 'Value in the manuscript': write Y in 'Agrees (Y/N)' if they are the same, "
   "N if they differ. Where the manuscript has rounded (table 0.9759, manuscript 0.98), write Y if "
   "it agrees after rounding and put 'rounded' in 'Your remarks'.",
   "7. Rows whose 'How to check' is 'column lookup' cannot be narrowed to one cell. Filter the named "
   "column of the named sheet and read the row that fits. 'Notes' lists the rows to look at, each with "
   "its cell address and the lesion, cohort or TF that row belongs to.",
   "8. Rows whose 'How to check' is 'read the sentence' cannot be pointed at a cell; the reason opens "
   "the 'Notes' entry (spans several rows as a range or mean / is a threshold or setting rather than a "
   "tabulated quantity / the sentence names nothing to narrow the rows with / the candidates stay "
   "many). Decide by reading whether the sentence contradicts the tables, and put the table value you "
   "used in 'Value you read'.",
   "9. Rows whose 'Notes' begin 'no cell in the supplementary workbook carries this value' are matched "
   "in a file outside S1-S12. The file name is given; if you cannot open it, write 'to be confirmed' in "
   "'Your remarks'.",
   "10. The 'Not in the tables' sheet holds the {u_n} numbers that no supplementary-table cell carries. "
   "Check each one in the file or the derivation named in 'Where to check it'.",
   "11. For every N, and for anything you were unsure about, write the reason in 'Your remarks'. Do not leave it blank.",
   "12. When you have finished, count the results yourself: open the auto-filter on row 1 of the "
   "'Check against tables' sheet, filter 'Agrees (Y/N)' to N and read the row count at the bottom left "
   "of the window. Count the Y rows and the blank rows the same way, then list the N and blank rows in "
   "your report. This workbook holds no formulas, so no count appears on its own.",
  ],
  "meta_head": "(the figures below describe the worksheet as built; they do not count your answers)",
  "meta_ms": "Manuscript checked", "meta_ms_v": "manuscript/{v} docx (clean version)",
  "meta_n": "Numbers in the manuscript",
  "meta_matched": "Checkable in a supplementary-table cell",
  "meta_other": "Not in the supplementary tables", "meta_other_v": "{n} (separate sheet)",
  "meta_cell": "of which a cell address could be given (How to check = cell lookup)",
  "meta_col": "a sheet and column could be given (How to check = column lookup)",
  "meta_absent": "of which no supplementary-table cell carries the value (matched in another file)",
  "meta_upg": "of which the identifiers in the sentence allowed a more precise reference",
  "meta_hcr_y": "of these, human_check_required = Y",
  "meta_hcr_n": "of these, human_check_required = N (already checked against the "
                "frozen file or script)",
  "meta_sent": "no cell could be given, read the sentence instead",
  "meta_src": "Built by", "unit": "",
 },
}


def T(key, **kw):
    t = STR[LANG][key]
    return t.format(**kw) if kw else t


# ---------------------------------------------------------------- cell index
wb = load_workbook(XLSX, data_only=True)
index = {}


def put(key, entry):
    index.setdefault(key, []).append(entry)


for ws in wb.worksheets:
    if ws.title == "index":
        continue
    header = [c.value for c in ws[1]]
    for row in ws.iter_rows(min_row=3):
        for c in row:
            if c.value is None:
                continue
            col = header[c.column - 1] if c.column - 1 < len(header) else ""
            if isinstance(c.value, (int, float)) and not isinstance(c.value, bool):
                v = float(c.value)
                for dp in range(0, 5):
                    put((rnd(abs(v), dp), dp),
                        (ws.title, c.coordinate, col, v, "numeric"))
            elif isinstance(c.value, str):
                for tok in re.findall(r"-?\d+(?:\.\d+)?", c.value):
                    try:
                        v = float(tok)
                    except ValueError:
                        continue
                    dp = len(tok.split(".")[1]) if "." in tok else 0
                    put((rnd(abs(v), dp), dp),
                        (ws.title, c.coordinate, col, v, "in text"))
print(f"cell index: {len(index)} value+precision keys over {len(wb.sheetnames) - 1} sheets")

# ------------------------------------------------- word context, 10 each side
NUM = re.compile(r"(?<![\w.])([−–-]?\d{1,3}(?:[, ]\d{3})*(?:\.\d+)?|[−–-]?\d*\.\d+)"
                 r"(?![\d])\s*(%|/\s*\d+)?")
HEAD = re.compile(r"^(Abstract|Introduction|Materials and methods|Results|Discussion|"
                  r"Conclusion[s]?|References|Keywords|Acknowledgements|"
                  r"Ethical approval|Data availability statement|"
                  r"Supplementary figures|Supplementary tables|Figure \d|Figure S\d)")
SKIP = re.compile(r"\(\d{4}\)|et al|doi|pp\.|ISBN")


def word_window(line, start, end, n=10):
    before = line[:start].split()[-n:]
    after = line[end:].split()[:n]
    return " ".join(before) + "  [[ " + line[start:end].strip() + " ]]  " + " ".join(after)


ctx, section = [], "front matter"
for line in TXT.split("\n"):
    if not line.strip():
        continue
    m = HEAD.match(line.strip())
    if m:
        section = m.group(1)
    if section == "References":
        continue
    for mo in NUM.finditer(line):
        if SKIP.search(line[max(0, mo.start() - 30):mo.end() + 30]):
            continue
        ctx.append((section, word_window(line, mo.start(), mo.end())))
assert len(ctx) == len(A), (len(ctx), len(A))

# ---------------------------------------------------------------- main sheet
# A number is anchored to a cell only when few cells carry it. Most manuscript
# numbers are small integers that recur throughout the workbook, so three levels are
# used: a cell address, a sheet-and-column to filter, or a reading of the sentence.
SENT = re.compile(r"Table (S\d+)")

# A value alone does not identify a cell: 0.9759 and 94.2 occur in several columns.
# Candidates are therefore ranked by how well the column name matches the words around
# the number in the manuscript, with a bonus when the sentence cites that table.
ABBREV = {"tf": "transcription factor", "rho": "correlation", "sd": "dispersion",
          "pct": "percent", "n": "number", "defa": "definition a", "lf": "leak free",
          "qc": "quality", "mt": "mitochondrial", "hvg": "variable"}


# Excel reads a string that starts with "=", "+", "-" or "@" as a formula; openpyxl
# stores it as one, and Excel then strips it and reports a repaired file. Every string
# written here is plain text, so a leading formula character is prefixed with an
# ellipsis (the context windows are cut mid-sentence anyway).
def safe_text(x):
    if isinstance(x, str) and x[:1] in ("=", "+", "-", "@"):
        return "… " + x
    return x


def words(text):
    w = set(re.findall(r"[a-z]+", str(text).lower()))
    return w | {x[:-1] for x in w if x.endswith("s")} | {x + "s" for x in w}


def column_score(window, sheet, column, cited):
    ctx = words(window)
    cw = set()
    for w in re.split(r"[_\s]+", str(column).lower()):
        if w:
            cw.add(w)
            cw |= words(ABBREV.get(w, ""))
    bonus = 2 if cited and (sheet == cited or sheet.startswith(cited + "b")) else 0
    return len(cw & ctx) + bonus



# ------------------------------------------------- row keys and context entities
# A column name tells you which quantity a number is; it does not tell you which row.
# The manuscript usually names the row in words - a lesion, a cohort, a transcription
# factor, a compartment, an axis - so those words are extracted from the sentence and
# used to narrow the candidate cells to the rows that carry the same identifiers.
KEY_COLUMNS = ["lesion", "cohort", "TF", "compartment", "axis", "cell_set", "cellset",
               "stage", "version", "pool", "donor", "run_accession", "gene", "set_name",
               "comparison"]
rowkeys = {}
for ws in wb.worksheets:
    if ws.title == "index":
        continue
    header = [c.value for c in ws[1]]
    cols = {h: j for j, h in enumerate(header, start=1) if h in KEY_COLUMNS}
    if not cols:
        continue
    for rr in range(3, ws.max_row + 1):
        d = {}
        for h, j in cols.items():
            v = ws.cell(rr, j).value
            if isinstance(v, str) and v.strip():
                d[h] = v.strip()
        if d:
            rowkeys[(ws.title, rr)] = d

TF_NAMES = {v for d in rowkeys.values() for k, v in d.items() if k in ("TF", "gene")}
LESIONS = ["A", "B", "C", "D", "E", "CAP1", "CAP2", "CAP3"]
AXIS_WORDS = {
    "lf_A_cytokine_chemokine": ["inflammatory", "cytokine", "chemokine"],
    "lf_B_survival_antiapoptotic": ["survival", "apoptotic", "antiapoptotic"],
    "lf_C1_proliferation": ["proliferation"],
    "lf_C2_migration": ["migration"],
    "lf_C_prolif_migration": ["proliferation and migration"],
    "lf_D_tissueforming": ["tissue-forming", "tissueforming", "tissue forming",
                           "hard-tissue", "osteogenic"],
}
COMPARTMENT_WORDS = {"stromal": ["stromal"], "immune": ["immune"],
                     "endothelial": ["endothelial"], "pericyte_smc": ["pericyte", "smooth"],
                     "immune_myeloid": ["myeloid"], "immune_T": ["T cell", "T-cell"],
                     "immune_plasma": ["plasma cell"]}
CELLSET_WORDS = {"marker_defined": ["marker-defined", "marker defined"],
                 "primary_cluster_based": ["cluster-based", "cluster based", "primary"],
                 "all cells": ["all cells"], "all_cells": ["all cells"],
                 "stromal": ["stromal"]}


def context_entities(window):
    """identifiers named in the words around the number, by key-column name"""
    w = str(window)
    low = w.lower()
    out = {}
    les = set()
    for cap in re.findall(r"CAP\s?([123])", w):
        les.add("CAP" + cap)
    # "lesion A", "lesions A and B", "lesions A-E": single capitals in a lesion phrase
    for mo in re.finditer(r"lesion[s]?\s+((?:[A-E](?:\s*[-\u2013,/and ]+\s*)?)+)", w):
        seg = mo.group(1)
        letters = re.findall(r"[A-E]", seg)
        if re.search(r"[-\u2013]", seg) and len(letters) == 2:
            a, b = letters
            les |= {chr(c) for c in range(ord(a), ord(b) + 1)}
        else:
            les |= set(letters)
    if les:
        out["lesion"] = les
    coh = set()
    if "discovery" in low:
        coh |= {"discovery", "discovery_n5", "POOLED_discovery"}
    if "validation" in low:
        coh |= {"validation", "validation_n3", "POOLED_validation"}
    if "GSE164241" in w or "gingiv" in low:
        coh |= {"GSE164241"}
    if coh:
        out["cohort"] = coh
    tfs = {t for t in re.findall(r"\b[A-Z][A-Z0-9]{2,}\b", w) if t in TF_NAMES}
    if tfs:
        out["TF"] = tfs
        out["gene"] = tfs
    ax = {k for k, ws_ in AXIS_WORDS.items() if any(x in low for x in ws_)}
    if ax:
        out["axis"] = ax
    cp = {k for k, ws_ in COMPARTMENT_WORDS.items()
          if any(x.lower() in low for x in ws_)}
    if cp:
        out["compartment"] = cp
    cs = {k for k, ws_ in CELLSET_WORDS.items() if any(x in low for x in ws_)}
    if cs:
        out["cell_set"] = cs
        out["cellset"] = {x for x in cs if x in ("marker_defined", "primary_cluster_based")}
    if "overlap-removed" in low or "overlap removed" in low:
        out["version"] = {"overlap_removed"}
    elif "as screened" in low or "as-screened" in low:
        out["version"] = {"as_screened"}
    if "sign-consistent" in low or "sign consistency" in low:
        out["pool"] = {"sign_consistent_non_candidates"}
    return {k: v for k, v in out.items() if v}


def row_of(coord):
    return int(re.sub(r"[^0-9]", "", coord))


def entity_filter(cands, ents):
    """candidates whose row carries the identifiers named in the sentence"""
    if not ents:
        return cands, []
    kept, used = [], set()
    for h in cands:
        keys = rowkeys.get((h[0], row_of(h[1])), {})
        ok, hit = True, False
        for k, allowed in ents.items():
            if k in keys:
                if keys[k] in allowed:
                    hit = True
                    used.add(k)
                else:
                    ok = False
                    break
        if ok and hit:
            kept.append(h)
    return (kept, sorted(used)) if kept else (cands, [])


def cited_sheet(line, start, end):
    ls = line.rfind(". ", 0, start)
    le = line.find(". ", end)
    seg = line[(ls + 1 if ls >= 0 else 0):(le if le >= 0 else len(line))]
    h = SENT.findall(seg)
    return h[0] if h else ""


cites = []
section = "front matter"
for line in TXT.split("\n"):
    if not line.strip():
        continue
    m = HEAD.match(line.strip())
    if m:
        section = m.group(1)
    if section == "References":
        continue
    for mo in NUM.finditer(line):
        if SKIP.search(line[max(0, mo.start() - 30):mo.end() + 30]):
            continue
        cites.append(cited_sheet(line, mo.start(), mo.end()))
assert len(cites) == len(A)

JP_KEY = {k: T("k_" + ("cell_set" if k == "cellset" else k))
          for k in ("lesion", "cohort", "TF", "gene", "compartment", "axis",
                    "cell_set", "cellset", "version", "pool", "stage", "donor",
                    "run_accession", "comparison", "set_name")}


# Thresholds, clustering parameters and constants inside formulae are printed as
# numbers but are not quantities the tables carry, so a value match to some cell is
# a coincidence. These are kept at the sentence level whatever the candidates look like.
GUARD_NOUNS = ["threshold", "resolution", "leiden", "scrublet", "version", "harmony"]
GUARD_PHRASES = ["near", "exceed", "exceeded", "exceeding", "at least", "no fewer",
                 "set to", "cannot attain", "minimum attainable", "expected rate",
                 "fewer than", "more than"]
GUARD_BEFORE = re.compile(r"(?:[<>\u2264\u2265]|/\()\s*$")


def is_parameter(window):
    """a threshold, a software setting or a constant in a formula, not a tabulated value"""
    head = str(window).split("[[")[0].rstrip()
    tail = head.lower().split()
    last1 = re.sub(r"[^a-z]", "", tail[-1]) if tail else ""
    near = " ".join(tail[-4:])
    return (any(g in last1 for g in GUARD_NOUNS)
            or any(g in near for g in GUARD_PHRASES)
            or bool(GUARD_BEFORE.search(head[-4:])))


def row_hints(cands, limit=4):
    """identifiers of the rows a column-level check should look at"""
    out = []
    for h in cands[:limit]:
        keys = rowkeys.get((h[0], row_of(h[1])), {})
        label = " ".join(f"{JP_KEY.get(k, k)}={keys[k]}" for k in
                         ("lesion", "cohort", "TF", "compartment", "axis", "version",
                          "pool", "cell_set", "cellset", "donor", "stage")
                         if k in keys)
        if label:
            out.append(f"{h[1]} [{label}]")
    extra = T("m_more", n=len(cands) - limit) if len(cands) > limit else ""
    return "; ".join(out) + extra if out else ""


def decide(cands, window, cited, loose):
    """how precisely this number can be pointed at: a cell, a column, or a sentence"""
    if is_parameter(window):
        return T("how_sent"), "", "", "", "", "", len(cands)
    scored = [(column_score(window, h[0], h[2], cited), h) for h in cands]
    scored.sort(key=lambda t: (-t[0], t[1][4] != "numeric",
                               wb.sheetnames.index(t[1][0]), t[1][1]))
    uniq = [h for _, h in scored]
    best = scored[0][0] if scored else -1
    n_best = sum(1 for sc, _ in scored if sc == best)
    pairs = list(dict.fromkeys((h[0], h[2]) for h in uniq))
    n = len(uniq)
    sh = cell = col = val = ""
    others = ""
    if n == 0:
        how = T("how_sent")
    elif n == 1 or (best >= 1 and n_best <= MAX_CELLS):
        how = T("how_cell")
        sh, cell, col, val = uniq[0][0], uniq[0][1], uniq[0][2], uniq[0][3]
        if n_best > 1:
            others = T("m_ties", list="; ".join(
                f"{a}!{b} ({c})" for a, b, c, *_ in uniq[1:n_best]))
        elif n > 1:
            others = T("m_best_of_n", n=n)
        if loose:
            others = T("m_loose") + others
    else:
        top = [h for h in uniq if (h[0], h[2]) == pairs[0]]
        if len(pairs) <= 4 or (best >= 1 and len(top) <= 6):
            how = T("how_col")
            sh, col = pairs[0]
            hint = row_hints(top)
            others = T("m_col", n=n,
                       rows=T("m_col_rows", hint=hint) if hint else "",
                       cols="; ".join(f"{a} / {c}" for a, c in pairs[:4]))
        else:
            how = T("how_sent")
            others = T("m_pairs",
                       cols="; ".join(f"{a} / {c}" for a, c in pairs[:3]))
    return how, sh, cell, col, val, others, n


RANGE_WORDS = ["range", "mean", "median", "total", "pooled", "across", "each",
               "every", "all ", "respectively", "sum"]


def why_not_a_cell(window, ents, cands):
    """stated reason a sentence-level check is the most that can be specified"""
    low = str(window).lower()
    if is_parameter(window):
        return T("w_param")
    pairs = {(h[0], h[2]) for h in cands}
    if re.search(r"\d\s*[-\u2013]\s*\d", str(window)) or any(
            x in low for x in RANGE_WORDS):
        return T("w_range", n=len(cands), p=len(pairs))
    if not ents:
        return T("w_nokeys", n=len(cands), p=len(pairs))
    return T("w_many", n=len(cands), p=len(pairs))


upgrades = []

HCR, HCR_EV = {}, {}
_h = f"{W}/08_external_review/tables/E12_human_check.csv"
if os.path.exists(_h):
    _hd = pd.read_csv(_h)
    HCR = {int(r.No): str(r.human_check_required) for _, r in _hd.iterrows()}
    HCR_EV = {int(r.No): (r.evidence_ja if LANG == "ja" else r.evidence)
              for _, r in _hd.iterrows() if r.human_check_required == "N"}

rows = []
for i, r in A.iterrows():
    sect, window = ctx[i]
    dp = int(r.decimals)
    hits = index.get((rnd(float(r.value), dp), dp), [])
    loose = False
    if not hits:
        # the manuscript may print a value the table carries at a different precision
        hits = [h for d in range(0, 5)
                for h in index.get((rnd(float(r.value), d), d), [])]
        loose = bool(hits)
    seen, uniq = set(), []
    for h in hits:
        if (h[0], h[1]) not in seen:
            seen.add((h[0], h[1]))
            uniq.append(h)
    # rounding to the printed precision can pull a distant value onto the printed one
    # (0.667 rounds to 1 at zero decimals), so candidates more than 5 % away from the
    # printed value are dropped unless that would leave none.
    close = [h for h in uniq
             if abs(float(h[3]) - abs(float(r.value))) <= 0.05 * abs(float(r.value))]
    if close:
        uniq = close
    if cites[i]:
        narrowed = [h for h in uniq
                    if h[0] == cites[i] or h[0].startswith(cites[i] + "b")]
        if narrowed:
            uniq = narrowed
    ents = context_entities(window)
    narrowed_cands, used_keys = entity_filter(uniq, ents)
    before = decide(uniq, window, cites[i], loose)
    after = decide(narrowed_cands, window, cites[i], loose)
    how, sh, cell, col, val, others, n = after
    if used_keys and after[0] != before[0]:
        others = T("m_narrowed", n=len(uniq),
                   keys=T("sep").join(JP_KEY.get(k, k) for k in used_keys)) + others
    elif used_keys:
        others = T("m_narrowed_done",
                   keys=T("sep").join(JP_KEY.get(k, k) for k in used_keys)) + others
    if n == 0:
        others = T("m_absent", src=(str(r.sources) or T("m_nosrc")))
    if how == T("how_sent") and n > 0:
        others = why_not_a_cell(window, ents, narrowed_cands) + others
    fx = FIXES.get(i + 1)
    if fx is not None and str(fx.action).startswith("cleared"):
        old = f"{sh}!{cell}" if sh and cell else "—"
        _ja = LANG == "ja"
        _why = (fx.get("quantity_ja") if _ja else None) or fx.quantity
        _how = (fx.get("check_instead_ja") if _ja else None) or fx.check_instead
        hw = (T("m_howelse", how=_how) if isinstance(_how, str) and _how else "")
        others = T("m_cleared", old=old, why=_why, how=hw)
        how, sh, cell, col, val = T("how_sent"), "", "", "", ""
        fx = None
    if fx is not None:
        how = T("how_cell")
        sh, cell, col = fx.new_sheet, fx.new_cell, fx.new_column
        try:
            val = float(fx.new_value)
        except (TypeError, ValueError):
            val = fx.new_value
        others = T("m_fixed", old=(fx.old_reference or "—")) + others
    # the agreement table counts a value as matched when a rho read as a percentage
    # happens to equal it; those rows need a closer look than a plain cell lookup.
    srcs = [x for x in str(r.sources or "").split("; ") if x]
    if srcs and all("(as %)" in x for x in srcs):
        others = T("m_pct") + others
    if not bool(r.matched):
        continue
    upgrades.append((i + 1, before[0], how))
    rows.append({
        T("c_no"): i + 1,
        T("c_section"): sect,
        T("c_context"): window,
        T("c_printed"): r.printed,
        T("c_how"): how,
        T("c_hcr"): HCR.get(i + 1, "Y"),
        T("c_file"): os.path.basename(XLSX) if sh else "",
        T("c_sheet"): sh,
        T("c_cell"): cell,
        T("c_col"): col,
        T("c_tabval"): val,
        T("c_ncand"): n,
        T("c_note"): ((T("m_nohuman", ev=HCR_EV[i + 1]) if i + 1 in HCR_EV else "")
                      + others),
        T("c_read"): "",
        T("c_agree"): "",
        T("c_remark"): "",
    })
MAIN = pd.DataFrame(rows)
ORDER = {T("how_sent"): 0, T("how_col"): 1, T("how_cell"): 2}
UPGRADE_PAIRS = {(b, a) for b in ORDER for a in ORDER if ORDER[a] > ORDER[b]}
UPG = pd.DataFrame(upgrades, columns=["No", "entity_filter_off", "entity_filter_on"])
UPG.to_csv(f"{B}/tables/check_worksheet_upgrades"
           f"{'' if LANG == 'ja' else '_' + LANG}.csv", index=False)
print(UPG.groupby(["entity_filter_off", "entity_filter_on"]).size().to_string())

# ------------------------------------------------------------ separate sheet
rows2 = []
for i, r in U.iterrows():
    j = A[(A.printed == r.printed) & (A.context == r.context)].index
    sect, window = ctx[j[0]] if len(j) else (r.section, r.context)
    rows2.append({
        T("c_no"): int(j[0]) + 1 if len(j) else i + 1,
        T("c_section"): sect,
        T("c_context"): window,
        T("c_printed"): r.printed,
        T("c_hcr"): HCR.get(int(j[0]) + 1 if len(j) else i + 1, "Y"),
        T("c_where"): (r.where_the_value_lives if LANG == "ja"
                       else r.where_the_value_lives_en),
        T("c_derivable"): (T("yes") if bool(
            r.get("is_derivable_from_a_supplementary_table", False)) else T("no")),
        T("c_note2"): HCR_EV.get(int(j[0]) + 1 if len(j) else i + 1, ""),
        T("c_read"): "",
        T("c_agree"): "",
        T("c_remark"): "",
    })
OTHER = pd.DataFrame(rows2)

HOWTO = pd.DataFrame({T("sheet_howto"): [x.replace("{u_n}", str(len(U)))
                                          for x in T("steps")]})
HOWTO.insert(0, "", ["" for _ in range(len(HOWTO))])
U_N = len(U)
MC = T("c_how")
META = pd.DataFrame([
    {T("c_item"): T("meta_head"), T("c_value"): ""},
    {T("c_item"): T("meta_ms"), T("c_value"): T("meta_ms_v", v=VER)},
    {T("c_item"): T("meta_n"), T("c_value"): f"{len(A)}{T('unit')}"},
    {T("c_item"): T("meta_matched"), T("c_value"): f"{int(A.matched.sum())}{T('unit')}"},
    {T("c_item"): T("meta_other"), T("c_value"): T("meta_other_v", n=U_N)},
    {T("c_item"): T("meta_cell"),
     T("c_value"): f"{int((MAIN[MC] == T('how_cell')).sum())}{T('unit')}"},
    {T("c_item"): T("meta_col"),
     T("c_value"): f"{int((MAIN[MC] == T('how_col')).sum())}{T('unit')}"},
    {T("c_item"): T("meta_absent"), T("c_value"):
        f"{int(MAIN[T('c_note')].str.startswith(T('m_absent_head')).sum())}{T('unit')}"},
    {T("c_item"): T("meta_upg"), T("c_value"):
        f"{int(sum(1 for _, b, a in upgrades if (b, a) in UPGRADE_PAIRS))}{T('unit')}"},
    {T("c_item"): T("meta_sent"),
     T("c_value"): f"{int((MAIN[MC] == T('how_sent')).sum())}{T('unit')}"},
    {T("c_item"): T("meta_hcr_y"),
     T("c_value"): f"{int((MAIN[T('c_hcr')] == 'Y').sum())}{T('unit')}"},
    {T("c_item"): T("meta_hcr_n"),
     T("c_value"): f"{int((MAIN[T('c_hcr')] == 'N').sum())}{T('unit')}"},
    {T("c_item"): T("meta_src"),
     T("c_value"): "07_submission_bundle/scripts/s5_check_worksheet.py"},
])

MAIN = MAIN.map(safe_text)
OTHER = OTHER.map(safe_text)
HOWTO = HOWTO.map(safe_text)
META = META.map(safe_text)

out = f"{B}/tables/check_worksheet{'' if LANG == 'ja' else '_' + LANG}.xlsx"
with pd.ExcelWriter(out, engine="openpyxl") as xl:
    HOWTO.to_excel(xl, sheet_name=T("sheet_howto"), index=False, startrow=3)
    META.to_excel(xl, sheet_name=T("sheet_howto"), index=False,
                  startrow=len(HOWTO) + 6)
    MAIN.to_excel(xl, sheet_name=T("sheet_main"), index=False)
    OTHER.to_excel(xl, sheet_name=T("sheet_other"), index=False)
    wb2 = xl.book
    ws = wb2[T("sheet_howto")]
    ws["A1"] = T("title")
    ws["A1"].font = Font(bold=True, size=13)
    ws["A2"] = T("subtitle")
    ws.column_dimensions["A"].width = 4
    ws.column_dimensions["B"].width = 110
    for r_ in ws.iter_rows(min_row=1, max_row=ws.max_row):
        for c in r_:
            c.alignment = Alignment(wrap_text=True, vertical="top")
    yellow = PatternFill("solid", fgColor="FFF2CC")
    head = Font(bold=True)
    for name, df in ((T("sheet_main"), MAIN), (T("sheet_other"), OTHER)):
        w = wb2[name]
        widths = {T("c_no"): 5, T("c_section"): 20, T("c_context"): 70,
                  T("c_printed"): 14, T("c_hcr"): 22, T("c_note2"): 60,
                  T("c_file"): 26, T("c_sheet"): 22,
                  T("c_cell"): 9, T("c_col"): 30, T("c_tabval"): 16,
                  T("c_read"): 16, T("c_agree"): 12, T("c_remark"): 28,
                  T("c_how"): 18, T("c_ncand"): 12, T("c_note"): 60,
                  T("c_where"): 70, T("c_derivable"): 16}
        for j, col in enumerate(df.columns, start=1):
            w.cell(row=1, column=j).font = head
            w.column_dimensions[get_column_letter(j)].width = widths.get(col, 16)
            if col in (T("c_read"), T("c_agree"), T("c_remark")):
                for i2 in range(2, len(df) + 2):
                    w.cell(row=i2, column=j).fill = yellow
        w.freeze_panes = "B2"
        w.auto_filter.ref = w.dimensions
print(META.to_string(index=False))
print("written", out)
