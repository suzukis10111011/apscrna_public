#!/usr/bin/env python
"""Purpose: Compare results/tables with results/tables_rerun table by table and column by column.
Input: results/tables/*.csv, results/tables_rerun/*.csv.
Output: results/tables/42_rerun_comparison.csv (per table: shape, numeric-column agreement, worst
  relative difference, and the reason when a table cannot be compared).
"""
import os, sys, glob, math
import numpy as np, pandas as pd

W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A, B = f"{W}/results/tables", f"{W}/results/tables_rerun"
REL_TOL = 1e-9

def numeric_cols(d):
    out = []
    for c in d.columns:
        s = pd.to_numeric(d[c], errors="coerce")
        nn = int(s.notna().sum())
        if nn and nn >= 0.5 * int((d[c].astype(str).str.strip() != "").sum() or 1):
            out.append(c)
    return out

rows = []
for f in sorted(glob.glob(f"{B}/*.csv")):
    name = os.path.basename(f)
    a_path = f"{A}/{name}"
    if not os.path.exists(a_path):
        rows.append(dict(table=name, status="new_in_rerun", note="no counterpart in results/tables"))
        continue
    a = pd.read_csv(a_path, dtype=str, keep_default_na=False)
    b = pd.read_csv(f, dtype=str, keep_default_na=False)
    rec = dict(table=name, rows_orig=len(a), rows_rerun=len(b),
               cols_orig=len(a.columns), cols_rerun=len(b.columns))
    if list(a.columns) != list(b.columns):
        rec.update(status="columns_differ",
                   note=f"only orig: {sorted(set(a.columns)-set(b.columns))[:4]}; only rerun: {sorted(set(b.columns)-set(a.columns))[:4]}")
        rows.append(rec); continue
    if len(a) != len(b):
        rec.update(status="row_count_differs", note="")
        rows.append(rec); continue
    ncols = numeric_cols(a)
    worst, worst_col, bad = 0.0, "", []
    for c in ncols:
        x = pd.to_numeric(a[c], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(b[c], errors="coerce").to_numpy(dtype=float)
        m = ~(np.isnan(x) & np.isnan(y))
        if (np.isnan(x[m]) != np.isnan(y[m])).any():
            bad.append(f"{c}(NaN pattern)"); continue
        k = m & ~np.isnan(x)
        if not k.any():
            continue
        denom = np.maximum(np.abs(x[k]), 1e-300)
        rel = np.abs(x[k] - y[k]) / denom
        mx = float(np.nanmax(rel)) if rel.size else 0.0
        if mx > worst: worst, worst_col = mx, c
        if mx > REL_TOL: bad.append(f"{c}({mx:.3g})")
    txt_diff = [c for c in a.columns if c not in ncols and not a[c].equals(b[c])]
    rec.update(status="match" if not bad else "MISMATCH",
               n_numeric_cols=len(ncols), worst_rel_diff=worst, worst_col=worst_col,
               mismatched_numeric_cols="; ".join(bad), text_cols_differing="; ".join(txt_diff[:6]))
    rows.append(rec)

only_orig = sorted(set(os.path.basename(x) for x in glob.glob(f"{A}/*.csv"))
                   - set(os.path.basename(x) for x in glob.glob(f"{B}/*.csv")))
for name in only_orig:
    rows.append(dict(table=name, status="not_regenerated",
                     note="produced by a stage outside the re-run set (quantification, preprocessing, figure/pptx tooling, or an interactive cell)"))

EXPLAIN = {
    "03c_vs_03b_counts_discovery.csv":
        "state-dependent baseline: 03c records the labels already on the input object as n_03b. The stored "
        "table was written when 03b still had its pre-fix implementation (stromal 9,862 / unassigned 390); "
        "the current 03b gives stromal 6,215 / unassigned 6,962. The result column n_03c is identical in both.",
    "03c_vs_03b_confusion_discovery.csv":
        "same state-dependent baseline as 03c_vs_03b_counts_discovery.csv: the row set is whatever labels the "
        "input object carried, so it cannot be reproduced from the corrected objects.",
    "28_manuscript_numbers.csv":
        "the 7 co-activity rows are appended by an interactive step of the 32 round, not by 28 itself "
        "(recorded in the script docstring); every row 28 does write matches.",
    "28_qc_all_stages.csv":
        "the discovery half comes from qc_cells_by_lesion_discovery.csv, a preprocessing-stage (03) output; "
        "the re-run starts from the processed objects, so that table was seeded rather than regenerated.",
}
for r in rows:
    if r.get("table") in EXPLAIN and r.get("status") in ("MISMATCH", "row_count_differs", "columns_differ"):
        r["explanation"] = EXPLAIN[r["table"]]

D = pd.DataFrame(rows)
D.to_csv(f"{A}/42_rerun_comparison.csv", index=False)
print(D.status.value_counts().to_string())
mis = D[D.status == "MISMATCH"]
print(f"\nMISMATCH: {len(mis)}")
if len(mis):
    print(mis[["table", "mismatched_numeric_cols", "worst_rel_diff"]].to_string(index=False))
cmpd = D[D.status.isin(["match", "MISMATCH"])]
print(f"\ncompared {len(cmpd)} tables | not regenerated {int((D.status=='not_regenerated').sum())} | "
      f"new {int((D.status=='new_in_rerun').sum())}")
print(f"worst relative difference among matching tables: {cmpd[cmpd.status=='match'].worst_rel_diff.max():.3g}"
      if (cmpd.status == "match").any() else "")
