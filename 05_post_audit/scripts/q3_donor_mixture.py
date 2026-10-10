#!/usr/bin/env python
"""POST HOC (audit item 3).

Every lesion is treated throughout the manuscript as one donor. This counts, per
lesion, cells carrying female-specific and male-specific transcripts in the raw
counts, so that a lesion that is in fact a pool of donors of different sex would be
visible.

  positive        raw count > 0 for that gene in that cell
  XIST_only       XIST positive, none of RPS4Y1 / DDX3Y / UTY / KDM5D positive
  Y_only          at least one Y gene positive, XIST not positive
  both            XIST positive and at least one Y gene positive
  minority_ratio  min(XIST_only, Y_only) / (XIST_only + Y_only)

XIST is absent from the validation matrix (data/processed/all_cells_validation.h5ad,
23,993 symbols), so CAP1-CAP3 are judged on the Y genes alone and the XIST-dependent
columns are left empty for them; judged_on records which rule applied.

Output: 05_post_audit/tables/Q3_donor_mixture.csv
        05_post_audit/tables/Q3_sex_marker_presence.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q0_common import AUD, POSTHOC_NOTE, PROC, log  # noqa: E402

FEM = ["XIST"]
MALE = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D"]
EXTRA = ["RPS4Y2", "EIF1AY"]

rows, pres = [], []
for cohort, h5 in [("discovery_n5", f"{PROC}/all_cells_discovery.h5ad"),
                   ("validation_n3", f"{PROC}/all_cells_validation.h5ad")]:
    A = sc.read_h5ad(h5)
    sym = A.var["symbol"].astype(str).values
    X = A.layers["counts"] if "counts" in A.layers else A.X
    les = A.obs.lesion.astype(str).values

    def pos(genes):
        idx = [i for i, s in enumerate(sym) if s in set(genes)]
        if not idx:
            return None, []
        sub = X[:, idx]
        v = np.asarray(sub.sum(axis=1)).ravel() if sparse.issparse(sub) else sub.sum(axis=1)
        return v > 0, [s for s in genes if s in set(sym)]

    f_pos, f_found = pos(FEM)
    m_pos, m_found = pos(MALE)
    e_pos, e_found = pos(EXTRA)
    for g in FEM + MALE + EXTRA:
        pres.append(dict(cohort=cohort, gene=g, present_in_matrix=bool(g in set(sym)),
                         role=("female" if g in FEM else "male"),
                         used_in_rule=bool(g in set(sym)) and g in FEM + MALE,
                         note=POSTHOC_NOTE))
    log(f"{cohort}: female markers found {f_found}, male markers found {m_found}, "
        f"other Y genes present {e_found}")

    for l in sorted(set(les)):
        m = les == l
        n = int(m.sum())
        r = dict(cohort=cohort, lesion=l, n_cells_qc_passed=n,
                 female_markers_used=",".join(f_found),
                 male_markers_used=",".join(m_found),
                 judged_on=("XIST and Y genes" if f_found else
                            "Y genes only (XIST absent from the matrix)"))
        yp = int(m_pos[m].sum()) if m_pos is not None else np.nan
        r.update(n_any_Y_positive=yp, frac_any_Y_positive=round(yp / n, 4))
        if e_pos is not None:
            r.update(n_RPS4Y2_or_EIF1AY_positive=int(e_pos[m].sum()),
                     frac_RPS4Y2_or_EIF1AY_positive=round(float(e_pos[m].mean()), 4))
        if f_pos is not None:
            xo = int((f_pos[m] & ~m_pos[m]).sum())
            yo = int((~f_pos[m] & m_pos[m]).sum())
            bo = int((f_pos[m] & m_pos[m]).sum())
            ne = int((~f_pos[m] & ~m_pos[m]).sum())
            den = xo + yo
            r.update(n_XIST_positive=int(f_pos[m].sum()), n_XIST_only=xo, n_Y_only=yo,
                     n_both_positive=bo, n_neither=ne,
                     minority_ratio=round(min(xo, yo) / den, 4) if den else np.nan)
        else:
            r.update(n_XIST_positive=np.nan, n_XIST_only=np.nan, n_Y_only=np.nan,
                     n_both_positive=np.nan, n_neither=np.nan, minority_ratio=np.nan)
        r["note"] = POSTHOC_NOTE
        rows.append(r)

T = pd.DataFrame(rows)
T.to_csv(f"{AUD}/tables/Q3_donor_mixture.csv", index=False)
pd.DataFrame(pres).to_csv(f"{AUD}/tables/Q3_sex_marker_presence.csv", index=False)
print(T[["cohort", "lesion", "n_cells_qc_passed", "n_XIST_only", "n_Y_only",
         "n_both_positive", "minority_ratio", "frac_any_Y_positive",
         "judged_on"]].to_string(index=False))
log("Q3 done")
