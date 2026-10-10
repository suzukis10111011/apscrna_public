#!/usr/bin/env python
"""Purpose: Compare GEO matrix and kb requantified data; align key metrics and candidate TFs
Input: none (reads from results/ directories created by upstream analysis)
Output: results/tables/09_geo_vs_kb.csv, 09_tf_agreement_geo_vs_kb.csv, 09_stromal_jaccard.csv

09: Line up the key metrics and candidate TFs of the GEO matrix version and the kb requantified version.
The kb side has two variants:
  kb    = QC identical to 03 applied to all kb barcodes (the QC threshold also serves as cell calling)
  kbgeo = restricted to only the barcodes contained in the GEO matrix (to look only at differences in quantification)
Output: results/tables/09_geo_vs_kb.csv (key metrics)
      results/tables/09_tf_agreement_geo_vs_kb.csv (agreement of candidate TFs)
      results/tables/09_stromal_jaccard.csv (Jaccard of stromal barcodes)
Each block is independently wrapped in try, so that even if one is missing, the rest are produced (for unattended execution).
"""
import os, sys, json, traceback, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Decide the numba cache destination before importing scanpy (site-packages is read-only, and
# if unset, the import itself fails with "cannot cache function ... no locator available").
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".numba_cache"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from _utils import bare_barcode

W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
SRCS = ["", "_kb", "_kbgeo"]
NAME = {"": "geo", "_kb": "kb_allbarcodes", "_kbgeo": "kb_geobarcodes"}


def rd(stem, sfx, **kw):
    p = f"{TAB}/{stem}{sfx}.csv"
    return pd.read_csv(p, **kw) if os.path.exists(p) else None


rows, problems = [], []


def add(quantity, lesion, sfx, value, note=""):
    rows.append(dict(quantity=quantity, lesion=lesion, source=NAME[sfx], value=value, note=note))


# ---- 1. Cell counts and compartment ratios after QC ----
for sfx in SRCS:
    try:
        q = rd("qc_cells_by_lesion_discovery", sfx)
        c = rd("compartment_by_lesion_discovery", sfx, index_col=0)
        if q is not None:
            last = q[q.stage == q.stage.unique()[-1]] if "stage" in q.columns else q
            for r in last.itertuples():
                add("n_cells_after_QC", getattr(r, "lesion", "?"), sfx, int(getattr(r, "n_cells", np.nan)))
        if c is not None:
            for les in c.index:
                for comp in ("stromal", "immune"):
                    col = f"pct_{comp}"
                    if col in c.columns:
                        add(f"pct_{comp}", str(les), sfx, float(c.loc[les, col]))
                if "n_total" in c.columns:
                    add("n_cells_compartment_table", str(les), sfx, int(c.loc[les, "n_total"]))
                if "n_stromal" in c.columns:
                    add("n_stromal", str(les), sfx, int(c.loc[les, "n_stromal"]))
    except Exception:
        problems.append(f"QC/compartment {NAME[sfx]}: {traceback.format_exc(limit=1).strip()}")

# ---- 2. Enrichment, ρ, Definition A axis score ----
for sfx in SRCS:
    try:
        e = rd("04_enrichment_by_lesion_discovery", sfx)
        if e is not None:
            d = e[e.score == "axis_D_tissueforming"]
            for r in d.itertuples():
                add("enrichment_fold_tissueforming", r.lesion, sfx, float(r.enrichment_fold))
                add("stromal_pct_top_decile", r.lesion, sfx, float(r.stromal_pct_top_decile))
        r2 = rd("04_rho_inflam_vs_tissueforming_discovery", sfx)
        if r2 is not None:
            sub = r2[(r2.cellset == "stromal")] if "cellset" in r2.columns else r2
            for defn in sub.definition.unique() if "definition" in sub.columns else [None]:
                s2 = sub[sub.definition == defn] if defn else sub
                for r in s2.itertuples():
                    add(f"depth_controlled_rho[{defn or 'axis'}]", r.lesion, sfx, float(r.depth_controlled_rho))
        a = rd("04_definitionA_discovery_v2", sfx)
        if a is not None:
            for r in a.itertuples():
                add("defA_frac_retained", r.lesion, sfx, float(r.frac_retained))
                add("defA_n_bins", r.lesion, sfx, int(r.n_bins))
    except Exception:
        problems.append(f"enrich/rho/defA {NAME[sfx]}: {traceback.format_exc(limit=1).strip()}")

R = pd.DataFrame(rows)
if len(R):
    P = R.pivot_table(index=["quantity", "lesion"], columns="source", values="value", aggfunc="first").reset_index()
    for c in ("kb_allbarcodes", "kb_geobarcodes"):
        if c in P.columns and "geo" in P.columns:
            P[f"diff_{c}_minus_geo"] = (P[c] - P["geo"]).round(4)
    P.to_csv(f"{TAB}/09_geo_vs_kb.csv", index=False)
    print("=== 09_geo_vs_kb.csv ===")
    print(P.head(40).to_string(index=False))
else:
    print("主要数値が1件も取れなかった")

# ---- 3. Agreement of candidate TFs ----
try:
    tf = {s: rd("04_tf_screen_discovery_candidates", s) for s in SRCS}
    out = []
    base = tf[""]
    if base is not None:
        for sfx in ["_kb", "_kbgeo"]:
            t = tf[sfx]
            if t is None:
                continue
            for ax in sorted(base.axis.unique()):
                b = base[(base.axis == ax) & base.is_candidate.astype(str).str.lower().eq("true")]
                k = t[(t.axis == ax) & t.is_candidate.astype(str).str.lower().eq("true")]
                bs, ks = set(b.TF), set(k.TF)
                m = base[base.axis == ax][["TF", "pooled_rho"]].merge(
                    t[t.axis == ax][["TF", "pooled_rho"]], on="TF", suffixes=("_geo", "_kb"))
                same = int((np.sign(m.pooled_rho_geo) == np.sign(m.pooled_rho_kb)).sum())
                rho = float(m.pooled_rho_geo.corr(m.pooled_rho_kb, method="spearman")) if len(m) > 2 else np.nan
                out.append(dict(source=NAME[sfx], axis=ax, n_candidates_geo=len(bs), n_candidates_kb=len(ks),
                                overlap=len(bs & ks),
                                jaccard=round(len(bs & ks) / max(len(bs | ks), 1), 4),
                                geo_only=";".join(sorted(bs - ks)), kb_only=";".join(sorted(ks - bs)),
                                n_TF_common=len(m), same_sign=same,
                                same_sign_pct=round(100 * same / max(len(m), 1), 2),
                                spearman_rho_of_pooled_rho=round(rho, 4) if rho == rho else np.nan))
    T = pd.DataFrame(out)
    if len(T):
        T.to_csv(f"{TAB}/09_tf_agreement_geo_vs_kb.csv", index=False)
        print("\n=== 候補TFの一致 ===")
        print(T[["source", "axis", "n_candidates_geo", "n_candidates_kb", "overlap", "jaccard",
                 "same_sign_pct", "spearman_rho_of_pooled_rho"]].to_string(index=False))
except Exception:
    problems.append(f"TF agreement: {traceback.format_exc(limit=1).strip()}")

# ---- 4. Jaccard of stromal barcodes ----
try:
    sets = {}
    for sfx in SRCS:
        p = f"{OUT}/stromal_discovery{sfx}.h5ad"
        if os.path.exists(p):
            a = sc.read_h5ad(p, backed="r")
            sets[NAME[sfx]] = pd.DataFrame(dict(bc=[bare_barcode(x) for x in a.obs_names],
                                                lesion=a.obs.lesion.astype(str).values))
            a.file.close()
    jr = []
    if "geo" in sets:
        g = sets["geo"]
        for k in [x for x in sets if x != "geo"]:
            o = sets[k]
            for les in sorted(set(g.lesion) | set(o.lesion)) + ["POOLED"]:
                A = set(g.bc if les == "POOLED" else g.bc[g.lesion == les])
                B = set(o.bc if les == "POOLED" else o.bc[o.lesion == les])
                jr.append(dict(comparison=f"geo vs {k}", lesion=les, n_geo=len(A), n_other=len(B),
                               intersection=len(A & B), union=len(A | B),
                               jaccard=round(len(A & B) / max(len(A | B), 1), 4),
                               pct_of_geo_recovered=round(100 * len(A & B) / max(len(A), 1), 2)))
    J = pd.DataFrame(jr)
    if len(J):
        J.to_csv(f"{TAB}/09_stromal_jaccard.csv", index=False)
        print("\n=== 間質バーコードの Jaccard ===")
        print(J.to_string(index=False))
except Exception:
    problems.append(f"Jaccard: {traceback.format_exc(limit=1).strip()}")

if problems:
    print("\n=== 取れなかったブロック ===")
    for p in problems: print("  -", p)
    open(f"{TAB}/09_problems.txt", "w").write("\n".join(problems))
print("\n09 DONE")
