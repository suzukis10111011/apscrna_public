#!/usr/bin/env python
"""Purpose: Redraw Figure S4 with the x axis of panel C named for what the bar shows.
Input: scripts/ms_figS4.py (unchanged on disk) and the result tables it reads
Output: results/figures_ms/FigS4_matrix_vs_requantification.{png,pdf} and the matching
  bare rendering under results/figures_ms/bare/ with its textmeta JSON; the published
  files are moved to results/superseded/figures_ms/ first

The published script is not edited. Its source is read, one string is replaced, and the
patched source is executed with scripts/ on sys.path, so the figure is produced by the
same code path as before apart from the label:

  "Candidates with the same sign of pooled $\\rho$ (%)"
    -> "TFs with the same sign of pooled $\\rho$ (%)"

The bar is same_sign_pct, the percentage of the TFs evaluable in both quantifications
whose pooled rho keeps its sign, so the old label named the wrong set of factors; the
bars, the numbers beside them and every other element are untouched.
"""
import os
import shutil
import sys

W = os.environ.get("APSC_ROOT", "/path/to/Documents/apscrna_v2")
SRC = f"{W}/scripts/ms_figS4.py"
NAME = "FigS4_matrix_vs_requantification"
OLD = 'Candidates with the same sign of pooled $\\\\rho$ (%)'
NEW = 'TFs with the same sign of pooled $\\\\rho$ (%)'

for sub, ext in [("", "png"), ("", "pdf"), ("bare/", "png"), ("bare/", "pdf"),
                 ("bare/textmeta/", "json")]:
    src = f"{W}/results/figures_ms/{sub}{NAME}.{ext}"
    dst = f"{W}/results/superseded/figures_ms/{sub}{NAME}.{ext}"
    if os.path.exists(src) and not os.path.exists(dst):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        print("kept the published file in results/superseded:", dst[len(W) + 1:])

source = open(SRC).read()
assert OLD in source, "the x-axis label of panel C is not where it was"
source = source.replace(OLD, NEW)
sys.path.insert(0, f"{W}/scripts")
g = {"__file__": SRC, "__name__": "__main__"}
exec(compile(source, SRC, "exec"), g)
print("redrawn with the label:", NEW.replace("\\\\", "\\"))
