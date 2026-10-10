#!/usr/bin/env python
"""POST HOC (audit item 3, the methods-reading part).

Whether a lesion is one donor is a claim about the source study, not about the
counts. The primary publications behind the eight lesions were read for any
statement that samples were pooled, multiplexed or combined before capture, and for
the donor sex each one declares; the declared sex is put beside the marker-based
call from Q3_donor_mixture.csv.

Sources read:
  GSE181688 / PRJNA753049 (lesions A, B, C) - Lin et al., Front Cell Dev Biol 2022,
    doi:10.3389/fcell.2021.820274, Materials and Methods 2.1 and the statistics
    paragraph preceding Figure 2
  GSE197680 / PRJNA811571 (lesions D, E) - no publication is linked from the GEO
    record; donor attributes are taken from the GEO sample characteristics
  PRJNA1213680 (lesions CAP1-CAP3) - Front Immunol 2025,
    doi:10.3389/fimmu.2025.1559220, Methods 2.2.1 and 2.2.2

Output: 05_post_audit/tables/Q3b_source_pooling_and_sex.csv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q0_common import AUD, POSTHOC_NOTE, log  # noqa: E402

LIN = ("Lin et al., Front Cell Dev Biol 2022, doi:10.3389/fcell.2021.820274")
ZHANG = ("Front Immunol 2025, doi:10.3389/fimmu.2025.1559220")
NO_POOL_LIN = ("explicit: the statistics paragraph states that each patient's tissue "
               "was processed individually and that each sample's suspension was run "
               "on its own Chromium chip")
NO_POOL_ZHANG = ("no pooling stated: Methods 2.2.1 lists three lesions from three "
                 "named donors and 2.2.2 describes loading the cell suspension on a "
                 "10x Chromium controller; no multiplexing or hashing is mentioned")
NO_PUB = ("no linked publication in the GEO record; GEO sample characteristics give "
          "one donor per sample and no pooling field")

ROWS = [
    dict(lesion="A", cohort="discovery_n5", accession="GSE181688 / GSM5509264",
         publication=LIN, pooling_statement=NO_POOL_LIN,
         declared_donor="1 donor, female, 26 y", declared_sex="Female"),
    dict(lesion="B", cohort="discovery_n5", accession="GSE181688 / GSM5509265",
         publication=LIN, pooling_statement=NO_POOL_LIN,
         declared_donor="1 donor, male, 44 y", declared_sex="Male"),
    dict(lesion="C", cohort="discovery_n5", accession="GSE181688 / GSM5509266",
         publication=LIN, pooling_statement=NO_POOL_LIN,
         declared_donor="1 donor, female, 27 y", declared_sex="Female"),
    dict(lesion="D", cohort="discovery_n5", accession="GSE197680 / GSM5928092",
         publication="(none linked)", pooling_statement=NO_PUB,
         declared_donor="1 donor, female, 26 y (GEO characteristics)",
         declared_sex="Female"),
    dict(lesion="E", cohort="discovery_n5", accession="GSE197680 / GSM5928093",
         publication="(none linked)", pooling_statement=NO_PUB,
         declared_donor="1 donor, female, 26 y (GEO characteristics)",
         declared_sex="Female"),
    dict(lesion="CAP1", cohort="validation_n3", accession="PRJNA1213680 / SRR32119703",
         publication=ZHANG, pooling_statement=NO_POOL_ZHANG,
         declared_donor="3 donors across CAP1-CAP3: male 57 y, female 25 y, "
                        "female 51 y; the paper does not say which run is which",
         declared_sex="(not resolvable per run from the paper)"),
    dict(lesion="CAP2", cohort="validation_n3", accession="PRJNA1213680 / SRR32119704",
         publication=ZHANG, pooling_statement=NO_POOL_ZHANG,
         declared_donor="3 donors across CAP1-CAP3: male 57 y, female 25 y, "
                        "female 51 y; the paper does not say which run is which",
         declared_sex="(not resolvable per run from the paper)"),
    dict(lesion="CAP3", cohort="validation_n3", accession="PRJNA1213680 / SRR32119705",
         publication=ZHANG, pooling_statement=NO_POOL_ZHANG,
         declared_donor="3 donors across CAP1-CAP3: male 57 y, female 25 y, "
                        "female 51 y; the paper does not say which run is which",
         declared_sex="(not resolvable per run from the paper)"),
]

M = pd.read_csv(f"{AUD}/tables/Q3_donor_mixture.csv").set_index("lesion")
out = []
for r in ROWS:
    m = M.loc[r["lesion"]]
    xo, yo = m.n_XIST_only, m.n_Y_only
    if m.judged_on.startswith("XIST"):
        call = ("single donor, female" if (yo == 0 or yo / max(xo, 1) < 0.01)
                else ("single donor, male" if (xo == 0 or xo / max(yo, 1) < 0.01)
                      else "mixed"))
    else:
        call = f"Y-positive fraction {m.frac_any_Y_positive} (XIST not in the matrix)"
    r.update(n_cells_qc_passed=int(m.n_cells_qc_passed),
             n_XIST_only=xo, n_Y_only=yo, n_both_positive=m.n_both_positive,
             minority_ratio=m.minority_ratio,
             frac_any_Y_positive=m.frac_any_Y_positive,
             marker_based_call=call, judged_on=m.judged_on,
             agrees_with_declared=(np.nan if not m.judged_on.startswith("XIST")
                                   else bool(r["declared_sex"].lower() in call)),
             note=POSTHOC_NOTE)
    out.append(r)

T = pd.DataFrame(out)
T.to_csv(f"{AUD}/tables/Q3b_source_pooling_and_sex.csv", index=False)
print(T[["lesion", "accession", "declared_sex", "marker_based_call",
         "agrees_with_declared", "minority_ratio"]].to_string(index=False))
log("Q3b done")
