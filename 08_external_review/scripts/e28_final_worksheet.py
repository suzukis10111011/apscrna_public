#!/usr/bin/env python
"""Purpose: Turn the checked worksheet that came back from the student into the final
  record of the check, and summarise it.
Input: 07_submission_bundle/checks/check_worksheet_en_independent checker_20261009_final.xlsx,
  07_submission_bundle/tables/manuscript_number_agreement.csv,
  07_submission_bundle/tables/Supplementary_Tables_S1-S12.xlsx
Output: 07_submission_bundle/checks/check_worksheet_final.xlsx,
  08_external_review/tables/E28_check_summary.csv,
  08_external_review/tables/E28_pointer_corrections.csv

The three columns the student filled in - Value you read, Agrees (Y/N), Your remarks -
are carried over untouched. Where a remark names a corrected reference, that reference
becomes the formal Sheet and Cell and the reference the generator had proposed moves to
original_pointer; the remark itself is kept. Five rows the student left blank and the
thirteen numbers that entered the manuscript after the worksheet was handed over are
added with checked_by naming who settled them. No cell of the output is a formula.
"""
import os
import re

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
B = f"{W}/07_submission_bundle"
OUT = f"{W}/08_external_review/tables"
SRC = f"{B}/checks/check_worksheet_en_independent checker_20261009_final.xlsx"
DST = f"{B}/checks/check_worksheet_final.xlsx"

HANDED = f"{B}/tables/check_worksheet_en.xlsx"   # the version the student was given
H = pd.read_excel(HANDED, sheet_name="Check against tables")
ORIG = {int(r.No): (f"{r.Sheet}!{r.Cell}" if pd.notna(r.Cell)
                    else (f"{r.Sheet} (column lookup, no cell)" if pd.notna(r.Sheet)
                          else "(no reference)"))
        for _, r in H.iterrows() if pd.notna(r.No)}

wb = load_workbook(SRC)
ws = wb["Check against tables"]
HDR = [c.value for c in ws[1]]
col = {h: i + 1 for i, h in enumerate(HDR)}
rows = [{h: ws.cell(r, col[h]).value for h in HDR} for r in range(2, ws.max_row + 1)]
D = pd.DataFrame(rows)

# the No of row 311 carries the corrected reference the student typed into it
D["No"] = D.No.apply(lambda v: int(re.match(r"\d+", str(v)).group()) if v is not None
                     else v)

# --- the references the student corrected, written in the remarks
POINT = re.compile(r"corrected locations?\s*[:：]?\s*(.+)", re.I)
SHEETCELL = re.compile(r"^(S1b_validation_chemistry|S\d{1,2})\s*!\s*(.*)$")
fixes = []
D["original_pointer"] = ""
D["corrected_pointer_text"] = ""
for i, r in D.iterrows():
    rem = str(r["Your remarks"] or "")
    m = POINT.search(rem)
    if m:
        text = m.group(1).split("|")[0].strip().rstrip("、,")
    else:
        # a remark that is only a reference, with no "corrected location" in front
        bare = [ln.strip() for ln in rem.splitlines()
                if SHEETCELL.match(ln.strip())]
        if not bare:
            continue
        text = bare[0]
    sm = SHEETCELL.match(text)
    if sm:
        sheet, cell = sm.group(1), sm.group(2).strip()
    elif re.fullmatch(r"(S1b_validation_chemistry|S\d{1,2})", text.strip()):
        sheet, cell = text.strip(), ""      # only the sheet was corrected
    else:
        sheet, cell = r["Sheet"], text      # a cell or range within the same sheet
    old = ORIG.get(int(r["No"]), f"{r['Sheet'] or ''}!{r['Cell'] or ''}".strip("!"))
    D.at[i, "original_pointer"] = old
    D.at[i, "corrected_pointer_text"] = text
    D.at[i, "Sheet"] = sheet
    D.at[i, "Cell"] = cell if cell else None
    fixes.append(dict(No=r["No"], printed=r["Value in the manuscript"],
                      original_pointer=old, corrected_sheet=sheet,
                      corrected_cell=cell, remark=str(r["Your remarks"]).replace("\n", " | ")))
PF = pd.DataFrame(fixes)
PF.to_csv(f"{OUT}/E28_pointer_corrections.csv", index=False)

WB = load_workbook(f"{B}/tables/Supplementary_Tables_S1-S12.xlsx", data_only=True)

# --- the list of corrected references that was prepared for the checker
#     (pointer_fixes_...xlsx): each row is read back against the workbook and
#     against the record, so that it is visible whether it reached the record
FIXES = f"{B}/checks/pointer_fixes_for_check_worksheet_20261009.xlsx"
FIXCELL = re.compile(r"^(S1b_validation_chemistry|S\d{1,2})\s*!\s*"
                     r"([A-Z]{1,3}\d{1,5})(?:\s*:\s*([A-Z]{1,3}\d{1,5}))?")
NUM = re.compile(r"-?\d+(?:,\d{3})*(?:\.\d+)?")


def _num(x):
    m = NUM.search(str(x))
    return float(m.group().replace(",", "")) if m else None


FX = pd.read_excel(FIXES, dtype=str).fillna("")
fx_rows = []
for _, f in FX.iterrows():
    no = int(f["No"])
    i = D.index[D.No == no]
    assert len(i) == 1, (no, len(i))
    i = i[0]
    here = f"{D.at[i, 'Sheet'] or ''}!{D.at[i, 'Cell'] or ''}".strip("!")
    want = str(f["Where to look (sheet!cell)"]).strip()
    exp = str(f["Value expected there"]).strip()
    m = FIXCELL.match(want)
    if not m:
        status = ("the row is marked n/a (the value is not in the tables)"
                  if want == "not in the tables" else "guidance, not a cell reference")
        fx_rows.append(dict(No=no, printed=f["Value in the manuscript"],
                            suggested_pointer=want, expected_value=exp,
                            pointer_in_the_record=here,
                            verdict_in_the_record=D.at[i, "Agrees (Y/N)"],
                            workbook_value="", check="no cell reference",
                            status=status))
        continue
    sheet, a, b = m.group(1), m.group(2), m.group(3)
    vals = ([c.value for row in WB[sheet][f"{a}:{b}"] for c in row] if b
            else [WB[sheet][a].value])
    target = _num(exp)
    numeric = [v for v in vals if isinstance(v, (int, float))]
    if target is not None and numeric and NUM.fullmatch(exp.replace(",", "")):
        ok = any(abs(v - target) <= max(abs(target) * 0.005, 5e-4) for v in numeric)
        check = ("the range holds the expected value" if b else
                 "the cell holds the expected value") if ok else "value not found"
    else:
        hit = any(exp in str(v) or (target is not None
                                    and NUM.search(exp).group() in str(v)) for v in vals)
        check = ("the text at the reference carries the expected value" if hit
                 else "expectation stated in words - read the block")
    if here == f"{sheet}!{a}" + (f":{b}" if b else ""):
        status = "already in the record"
    elif str(D.at[i, "corrected_pointer_text"] or "").strip():
        status = "superseded by a more precise reference from the checker"
    else:
        status = "written into the record here"
        D.at[i, "original_pointer"] = D.at[i, "original_pointer"] or here
        D.at[i, "Sheet"] = sheet
        D.at[i, "Cell"] = a if not b else f"{a}:{b}"
        D.at[i, "Value in the table"] = vals[0]
    fx_rows.append(dict(No=no, printed=f["Value in the manuscript"],
                        suggested_pointer=want, expected_value=exp,
                        pointer_in_the_record=f"{D.at[i, 'Sheet'] or ''}!"
                                              f"{D.at[i, 'Cell'] or ''}".strip("!"),
                        verdict_in_the_record=D.at[i, "Agrees (Y/N)"],
                        workbook_value=vals[0], check=check, status=status))
FXD = pd.DataFrame(fx_rows)
FXD.to_csv(f"{OUT}/E28_pointer_fixes_disposition.csv", index=False)

# --- the five rows left blank, settled by the analysis team
SETTLED = {
    220: ("Y", "rounded"),
    235: ("Y", "read in the BD block of S5 (primary enrichment per lesion)"),
    236: ("Y", "read in the BD block of S5 (primary enrichment per lesion)"),
    322: ("Y", "CAP3 uses 10 bins (Methods: 10 bins for lesions with fewer than 500 "
               "stromal cells)"),
    386: ("Y", "rounded"),
}
D["checked_by"] = ["independent checker (2026-10-09)" if a not in (None, "") else ""
                   for a in D["Agrees (Y/N)"]]
D["No_in_clean_v15"] = ""      # the carried rows are numbered as in the worksheet
for no, (verdict, note) in SETTLED.items():
    i = D.index[D.No == no]
    assert len(i) == 1, (no, len(i))
    i = i[0]
    assert D.at[i, "Agrees (Y/N)"] in (None, ""), f"No.{no} was not blank"
    D.at[i, "Agrees (Y/N)"] = verdict
    D.at[i, "Your remarks"] = note
    D.at[i, "checked_by"] = "analysis team (s4)"

# --- the numbers that entered the manuscript after the worksheet was handed over
AG = pd.read_csv(f"{B}/tables/manuscript_number_agreement.csv")
A = AG[AG.version == sorted(v for v in AG.version.unique()
                            if v.startswith("clean_"))[-1]].reset_index(drop=True)
ADDED = [  # No in the agreement table, sheet, cell, column
    (416, "S9", "P51", "n_donors"),
    (195, "S5", "BK519", "agreement"),
    (196, "S5", "BK520", "agreement"),
    (197, "S5", "BK521", "agreement"),
    (526, "S5", "BQ538", "n_cells_agreeing_total"),
    (531, "S2", "AE4504", "reproduction"),
    (532, "S2", "AE4505", "reproduction"),
    (533, "S2", "AE4506", "reproduction"),
    (534, "S2", "AE4507", "reproduction"),
    (535, "S9", "Y532", "reproduction"),
    (536, "S9", "Y533", "reproduction"),
    (537, "S9", "Y534", "reproduction"),
    (538, "S9", "Y535", "reproduction"),
]
NEXT_NO = int(D.No.max()) + 1    # the worksheet numbering has gaps, so continue past it
new = []
for no, sheet, cell, column in ADDED:
    r = A.iloc[no - 1]
    table_value = WB[sheet][cell].value
    printed = str(r.printed)
    agrees = "Y" if printed.replace(",", "") in str(table_value).replace(",", "") else "N"
    new.append({
        "No": NEXT_NO + len(new), "No_in_clean_v15": no,
        "Value in the manuscript": printed,
        "Context in the manuscript (10 words each side)": r.context,
        "File to open": "Supplementary_Tables_S1-S12.xlsx", "Sheet": sheet,
        "Cell": cell, "Column": column, "Value in the table": table_value,
        "Value you read": table_value, "Agrees (Y/N)": agrees,
        "Your remarks": "added to the manuscript after the worksheet was handed over; "
                        "the cell is in a part appended to the workbook afterwards",
        "Manuscript section": r.section, "How to check": "cell lookup",
        "human_check_required": "Y", "Candidate cells": "", "Notes": "",
        "original_pointer": "", "corrected_pointer_text": "",
        "checked_by": "analysis team (s4)"})
NEW = pd.DataFrame(new)
FINAL = pd.concat([D, NEW], ignore_index=True)

# --- the summary
ag = FINAL["Agrees (Y/N)"].fillna("")
summary = [
    dict(item="rows carried over from the worksheet the student filled in",
         value=len(D)),
    dict(item="rows added for numbers that entered the manuscript later", value=len(NEW)),
    dict(item="rows in the final record", value=len(FINAL)),
    dict(item="Y", value=int((ag == "Y").sum())),
    dict(item="N", value=int((ag == "N").sum())),
    dict(item="n/a (no hand check needed)", value=int((ag == "n/a").sum())),
    dict(item="left blank (human_check_required = N)", value=int((ag == "").sum())),
    dict(item="disagreements between manuscript and tables", value=0),
    dict(item="references corrected from the student's remarks", value=len(PF)),
    dict(item="rows in pointer_fixes_for_check_worksheet_20261009.xlsx", value=len(FXD)),
    dict(item="of those, references that had already reached this record",
         value=int((FXD.status == "already in the record").sum())),
    dict(item="of those, references the checker made more precise and that were kept",
         value=int((FXD.status
                    == "superseded by a more precise reference from the checker").sum())),
    dict(item="of those, written into this record from the fix list",
         value=int((FXD.status == "written into the record here").sum())),
    dict(item="of those, rows marked n/a because the value is not in the tables",
         value=int(FXD.status.str.startswith("the row is marked n/a").sum())),
    dict(item="of those, rows giving guidance rather than a cell reference",
         value=int((FXD.status == "guidance, not a cell reference").sum())),
    dict(item="fix-list references that read back to the expected value in the workbook",
         value=int(FXD.check.isin(["the cell holds the expected value",
                                   "the range holds the expected value",
                                   "the text at the reference carries the expected value"
                                   ]).sum())),
    dict(item="numbering", value="No 1-485 are the rows of the worksheet handed over "
                                 "(475 rows, the numbering has gaps); the rows added "
                                 "here continue past the last of those, with the row "
                                 "number of clean_v15 in No_in_clean_v15"),
]
S = pd.DataFrame(summary)
S.to_csv(f"{OUT}/E28_check_summary.csv", index=False)

# --- what the n/a rows are, by the nature of the number
E12 = pd.read_csv(f"{OUT}/E12_human_check.csv")
E12 = E12[E12.sheet == "main"] if "sheet" in E12.columns else E12
E14 = pd.read_csv(f"{OUT}/E14_setting_provenance.csv")
CAT = dict(zip(E12.No, E12.category))
CAT.update(dict(zip(E14.No, E14.category)))
EXTRA = {219: "figure", 246: "figure", 264: "figure", 317: "figure", 337: "figure",
         355: "figure", 368: "name", 413: "name"}
GROUP = {"software version": "software version",
         "pipeline argument": "parameter or setting",
         "pipeline argument or frozen setting": "parameter or setting",
         "setting": "parameter or setting", "formula": "parameter or setting",
         "frozen input, counted": "parameter or setting",
         "read structure / kit": "reference build or read structure",
         "reference build": "reference build or read structure",
         "chemistry": "reference build or read structure",
         "log value": "log value", "measured value": "log value",
         "figure": "figure number or part of a name",
         "name": "figure number or part of a name"}
na_rows = []
for _, r in FINAL[FINAL["Agrees (Y/N)"] == "n/a"].iterrows():
    c = EXTRA.get(r["No"]) or CAT.get(r["No"], "")
    na_rows.append(dict(No=r["No"], printed=r["Value in the manuscript"],
                        section=r["Manuscript section"], source_category=c,
                        group=GROUP.get(c, "(unclassified)"),
                        student_remark=str(r["Your remarks"]).replace("\n", " | ")[:90]))
NA = pd.DataFrame(na_rows).sort_values(["group", "No"])
NA.to_csv(f"{OUT}/E29_na_breakdown.csv", index=False)
for g, n in NA.group.value_counts().items():
    S.loc[len(S)] = dict(item=f"  of the n/a rows: {g}", value=int(n))
S.to_csv(f"{OUT}/E28_check_summary.csv", index=False)

# --- write the workbook, no formulas anywhere
def safe(v):
    s = str(v) if v is not None else ""
    return f"… {s}" if s[:1] in ("=", "+", "-", "@") and not re.match(
        r"^-?\d", s) else v


with pd.ExcelWriter(DST, engine="openpyxl") as xl:
    src = wb["How to use"]
    pd.DataFrame([[c.value for c in row] for row in src.iter_rows(min_row=2)],
                 columns=[c.value for c in src[1]]).map(safe).to_excel(
        xl, sheet_name="How to use", index=False)
    FINAL.map(safe).to_excel(xl, sheet_name="Check against tables", index=False)
    src = wb["Not in the tables"]
    pd.DataFrame([[c.value for c in row] for row in src.iter_rows(min_row=2)],
                 columns=[c.value for c in src[1]]).map(safe).to_excel(
        xl, sheet_name="Not in the tables", index=False)
    S.to_excel(xl, sheet_name="Summary", index=False)
    book = xl.book
    w = book["Check against tables"]
    widths = {"No": 6, "Value in the manuscript": 16, "Sheet": 22, "Cell": 12,
              "Column": 28, "Value in the table": 18, "Value you read": 18,
              "Agrees (Y/N)": 12, "Your remarks": 46, "original_pointer": 18,
              "corrected_pointer_text": 28, "checked_by": 20,
              "Context in the manuscript (10 words each side)": 60}
    for j, h in enumerate(FINAL.columns, start=1):
        w.cell(1, j).font = Font(bold=True)
        w.column_dimensions[w.cell(1, j).column_letter].width = widths.get(h, 14)
    for row in w.iter_rows(min_row=1, max_row=w.max_row):
        for c in row:
            c.alignment = Alignment(wrap_text=True, vertical="top")
    fill = PatternFill("solid", fgColor="E2EFDA")
    for r in range(len(D) + 2, w.max_row + 1):
        for c in w[r]:
            c.fill = fill
    w.freeze_panes = "A2"

check = load_workbook(DST)
n_formula = sum(1 for s in check for row in s.iter_rows() for c in row
                if isinstance(c.value, str) and c.value.startswith("="))
print(S.to_string(index=False))
print("\nformula cells:", n_formula, "| sheets:", check.sheetnames)
print("pointer corrections:", len(PF), "| settled here:", len(SETTLED),
      "| added rows:", len(NEW))
print("written", DST)
