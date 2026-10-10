#!/usr/bin/env python
"""Purpose: Turn the correct references the external check recorded in its remarks into
  a correction table the worksheet generator applies, and state why the remaining ones
  cannot become a cell address.
Input: 07_submission_bundle/tables/check_worksheet_en_ChatGPT_20261007.xlsx,
  07_submission_bundle/tables/Supplementary_Tables_S1-S12.xlsx,
  07_submission_bundle/tables/manuscript_number_agreement.csv
Output: 08_external_review/tables/E2_worksheet_reference_corrections.csv

Every parsed reference is read back from the workbook and compared with the value the
manuscript prints, at the precision the manuscript prints it, allowing for a fraction
written as a percentage. Only references that read back are handed to the generator.
"""
import os
import re

import pandas as pd
from openpyxl import load_workbook

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
BUN = f"{W}/07_submission_bundle/tables"
OUT = f"{W}/08_external_review/tables"
os.makedirs(OUT, exist_ok=True)

GX = f"{BUN}/check_worksheet_en_ChatGPT_20261007.xlsx"
XLSX = f"{BUN}/Supplementary_Tables_S1-S12.xlsx"
G = pd.ExcelFile(GX).parse("Check against tables")
WB = load_workbook(XLSX, data_only=True)
HEADERS = {s: [c.value for c in WB[s][1]] for s in WB.sheetnames if s != "index"}

CELL = re.compile(r"\b(S\d+[A-Za-z_]*)!([A-Z]{1,3})(\d+)\b")
# a reference that names a range, a minimum or maximum, or an arithmetic step is not a
# single cell; these are listed with the reason rather than written into the worksheet
DERIVED = re.compile(r"\bmin\b|\bmax\b|:[A-Z]{1,3}\d+|gives|derived|abs\(|/ S\d|"
                     r"range|rows \d|has \d+|donor rows")
ABSENT = re.compile(r"not supplied|not found|not tabulated|required|"
                    r"Figure \d heading|cytokine name|threshold", re.I)


def printed_value(s):
    t = (str(s).replace("−", "-").replace("–", "-").replace(",", "")
         .replace(" ", "").rstrip("%").split("/")[0])
    try:
        return abs(float(t)), (len(t.split(".")[1]) if "." in t else 0)
    except ValueError:
        return None, None


def reads_back(printed, value):
    """does the cell hold the printed number, as printed or as a percentage"""
    pv, dp = printed_value(printed)
    if pv is None or value is None:
        return None
    if isinstance(value, (int, float)):
        base = [abs(float(value))]
    else:
        base = [abs(float(x)) for x in re.findall(r"-?\d+(?:\.\d+)?", str(value))]
    for b in base:
        for cand in (b, b * 100, b / 100):
            if round(cand, dp) == round(pv, dp):
                return True
    return False


rows = []
for _, r in G.iterrows():
    note = str(r["Your remarks"] or "")
    m = re.search(r"Correct source:\s*(.*)", note)
    if not m:
        continue
    txt = m.group(1).strip()
    found = CELL.findall(txt)
    old = (f"{r['Sheet']}!{r['Cell']}" if str(r["Sheet"]) not in ("nan", "")
           and str(r["Cell"]) not in ("nan", "") else str(r["Sheet"] or ""))
    new_sheet = new_cell = new_col = ""
    value = None
    ok = None
    if found:
        sh, col_letter, row_no = found[0]
        if sh in HEADERS:
            new_sheet, new_cell = sh, f"{col_letter}{row_no}"
            value = WB[sh][new_cell].value
            idx = 0
            for ch in col_letter:
                idx = idx * 26 + (ord(ch) - 64)
            hdr = HEADERS[sh]
            new_col = hdr[idx - 1] if idx - 1 < len(hdr) else ""
            ok = reads_back(r["Value in the manuscript"], value)
    if ok and not DERIVED.search(txt):
        kind = "applied: single cell, reads back"
        reason = ""
    elif DERIVED.search(txt):
        kind = "not applied: the reference is a range or a derivation"
        reason = ("the correct source is not one cell but a minimum, a maximum, a "
                  "range or an arithmetic step over several cells")
        new_sheet = new_cell = new_col = ""
    elif ABSENT.search(txt):
        kind = "not applied: the quantity is not in S1-S12"
        reason = ("no cell of the supplementary workbook carries this quantity; see "
                  "E3 for the parts added to the tables")
        new_sheet = new_cell = new_col = ""
    elif ok is False:
        kind = "not applied: the manuscript value differs from the named cell"
        reason = "handled as a manuscript-side item in E1/E4"
        new_sheet = new_cell = new_col = ""
    else:
        kind = "not applied: reference could not be resolved"
        reason = "the remark names no cell of the supplementary workbook"
        new_sheet = new_cell = new_col = ""
    rows.append(dict(No=int(r["No"]), printed=str(r["Value in the manuscript"]),
                     how_to_check=r["How to check"], old_reference=old,
                     new_sheet=new_sheet, new_cell=new_cell, new_column=new_col,
                     new_value=value if new_cell else None, reads_back=ok,
                     action=kind, reason=reason, external_remark=txt))

E2 = pd.DataFrame(rows).sort_values("No")
E2.to_csv(f"{OUT}/E2_worksheet_reference_corrections.csv", index=False)
print(E2.action.value_counts().to_string())
print("\napplied rows:", int(E2.action.str.startswith("applied").sum()))
print(E2[E2.action.str.startswith("applied")][["No", "printed", "old_reference",
                                               "new_sheet", "new_cell", "new_column"]]
      .head(12).to_string(index=False))
