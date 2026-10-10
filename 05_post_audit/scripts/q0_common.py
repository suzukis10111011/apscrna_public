#!/usr/bin/env python
"""POST HOC (independent audit, added after the pre-registered analysis).

Shared helpers for the 05_post_audit scripts. Read-only with respect to the
repository: nothing here writes into results/tables, data/processed or frozen/.
New intermediates go to data/post_audit, new tables to 05_post_audit/tables.

dctrl_matrix reproduces scripts/04,05 `dctrl` exactly (Spearman of the residuals
of the ranks on [1, logdepth, n_genes]) for a whole activity matrix at once.
"""
import os
import time

import numpy as np
import pandas as pd
from scipy.stats import rankdata

W = "/path/to/Documents/apscrna_v2"
TAB = f"{W}/results/tables"
PROC = f"{W}/data/processed"
OUT = f"{W}/data/post_audit"
AUD = f"{W}/05_post_audit"
AXES = ["lf_A_cytokine_chemokine", "lf_D_tissueforming"]
AXIS_SOURCE = {"lf_A_cytokine_chemokine": "A_cytokine_chemokine",
               "lf_D_tissueforming": "D_tissueforming"}
LESIONS = {"discovery": ["A", "B", "C", "D", "E"], "validation": ["CAP1", "CAP2", "CAP3"]}


def log(msg, fname="run_post_audit.log"):
    s = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(s, flush=True)
    os.makedirs(f"{AUD}/logs", exist_ok=True)
    with open(f"{AUD}/logs/{fname}", "a") as fh:
        fh.write(s + "\n")


def _resid_ranks(V, M, Mpinv):
    """Residuals of the column-wise ranks of V on the design M."""
    R = rankdata(V, axis=0)
    return R - M @ (Mpinv @ R)


def _pearson_cols(A, b):
    """Pearson r between every column of A and the vector b."""
    A = A - A.mean(0)
    b = b - b.mean()
    sa = np.sqrt((A ** 2).sum(0))
    sb = np.sqrt((b ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return (A * b[:, None]).sum(0) / (sa * sb)


def dctrl_matrix(ACT, y, logdepth, n_genes):
    """Vectorised equivalent of dctrl(x, y, dep, ng) over all columns of ACT."""
    M = np.column_stack([np.ones(len(y)), np.asarray(logdepth, float),
                         np.asarray(n_genes, float)])
    Mpinv = np.linalg.pinv(M)
    rx = _resid_ranks(np.asarray(ACT, float), M, Mpinv)
    ry = _resid_ranks(np.asarray(y, float).reshape(-1, 1), M, Mpinv)[:, 0]
    return _pearson_cols(rankdata(rx, axis=0), rankdata(ry))


def rho_by_lesion(ACT, obs, axis_col):
    """TF x lesion matrix of depth-controlled rho, plus a POOLED column."""
    les = sorted(obs["lesion"].astype(str).unique())
    cols = {}
    for l in les:
        m = (obs["lesion"].astype(str) == l).values
        cols[l] = dctrl_matrix(ACT.values[m], obs[axis_col].values[m],
                               obs["logdepth"].values[m], obs["n_genes"].values[m])
    cols["POOLED"] = dctrl_matrix(ACT.values, obs[axis_col].values,
                                  obs["logdepth"].values, obs["n_genes"].values)
    return pd.DataFrame(cols, index=ACT.columns)


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else np.nan


POSTHOC_NOTE = ("post hoc; added in response to the independent audit. "
                "Frozen gene sets, candidate factors and thresholds are unchanged.")
