#!/usr/bin/env python
"""Purpose: Compile manuscript numerical results from GEO discovery cohort; perform LOLO cross-validation,
  enrichment permutation test, and immune decile ranking
Input: all_cells_discovery_obs_scored.parquet, stromal_discovery.h5ad,
  stromal_discovery_tf_activity.parquet, TF screen & validation tables, GSE164241 summary
Output: 28_manuscript_numbers.csv, 28_lolo_grant_factors.csv, 28_enrichment_permutation.csv,
  28_qc_all_stages.csv

28: List of numerical values to fill in the manuscript's brackets + new computations (LOLO, enrichment permutation test, immune top decile).
The main analysis is finalized on the GEO matrix version. The kb requantification is only referenced from the 09 table as a supplementary robustness check.
Output: results/tables/28_manuscript_numbers.csv (item / value / source file)
      results/tables/28_lolo_grant_factors.csv
      results/tables/28_enrichment_permutation.csv
      results/tables/28_qc_all_stages.csv
Note: the 32_* rows of 28_manuscript_numbers.csv are appended by 32_coactivity_ruleset.py,
not by this script; re-running this script alone regenerates only its own rows.
"""
import os, sys, json, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".numba_cache"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "5")
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr, rankdata

W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
OUT = os.environ.get("APSC_OUTDIR", f"{W}/data/processed")
SEED = 0
rng = np.random.default_rng(SEED)
# The factor list in the grant application is private (docs/private/grant_factors.txt). If absent, skip the LOLO block.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _utils import load_grant_factors
_gf = load_grant_factors(W)
INFL5 = _gf["inflammatory"] if _gf else []
TISS5 = _gf["tissueforming"] if _gf else []
AXOF = {**{t: "lf_A_cytokine_chemokine" for t in INFL5}, **{t: "lf_D_tissueforming" for t in TISS5}}
rows = []


def add(item, value, source, note=""):
    rows.append(dict(item=item, value=value, source=source, note=note))


def dctrl(x, y, dep, ng):                                    # Copied verbatim from 04
    M = np.column_stack([np.ones(len(x)), dep, ng])
    def res(v):
        b, *_ = np.linalg.lstsq(M, rankdata(v), rcond=None); return rankdata(v) - M @ b
    return spearmanr(res(x), res(y)).statistic


# ============ 1. leave-one-lesion-out (10 factors from the grant application) ============
St = sc.read_h5ad(f"{OUT}/stromal_discovery.h5ad")
ACT = pd.read_parquet(f"{OUT}/stromal_discovery_tf_activity.parquet")
obs = St.obs
assert len(ACT) == len(obs), (ACT.shape, obs.shape)
les_all = sorted(obs.lesion.astype(str).unique())
print(f"間質 {St.shape} / TF {ACT.shape[1]} / 病変 {les_all}")

# Confirm that the implementation is identical to 04: reproduce one saved pooled_rho
S04 = pd.read_csv(f"{TAB}/04_tf_screen_discovery_candidates.csv")
chk = S04[(S04.axis == "lf_D_tissueforming") & (S04.TF == "SP7")].iloc[0]
got = dctrl(ACT["SP7"].values, obs["lf_D_tissueforming"].values, obs.logdepth.values, obs.n_genes.values)
assert abs(got - chk.pooled_rho) < 1e-6, f"dctrl の再現失敗: {got} vs {chk.pooled_rho}"
print(f"dctrl 再現確認 OK (SP7 lf_D: {got:.6f})")

lolo = []
for held in les_all:
    keep = (obs.lesion.astype(str) != held).values
    o = obs[keep]; A = ACT[keep]
    inc = sorted(o.lesion.astype(str).unique())
    for axis in ("lf_A_cytokine_chemokine", "lf_D_tissueforming"):
        y = o[axis].values
        recs = []
        for tf in A.columns:
            x = A[tf].values
            per = np.array([dctrl(x[(o.lesion.astype(str) == l).values], y[(o.lesion.astype(str) == l).values],
                                  o.logdepth.values[(o.lesion.astype(str) == l).values],
                                  o.n_genes.values[(o.lesion.astype(str) == l).values]) for l in inc])
            pooled = dctrl(x, y, o.logdepth.values, o.n_genes.values)
            recs.append(dict(TF=tf, pooled_rho=pooled, abs_pooled=abs(pooled),
                             same_sign=bool((per > 0).all() or (per < 0).all())))
        R = pd.DataFrame(recs)
        elig = R[R.same_sign].sort_values("abs_pooled", ascending=False)
        cand = set(elig.head(20).TF)
        rank = {t: i + 1 for i, t in enumerate(elig.TF)}
        for tf in (INFL5 if axis.startswith("lf_A") else TISS5):
            r = R[R.TF == tf]
            lolo.append(dict(held_out_lesion=held, n_lesions_used=len(inc), axis=axis, TF=tf,
                             pooled_rho=round(float(r.pooled_rho.iloc[0]), 4) if len(r) else np.nan,
                             sign_consistent=bool(r.same_sign.iloc[0]) if len(r) else False,
                             rank_within_axis=rank.get(tf, np.nan), is_candidate=tf in cand))
    print(f"  LOLO {held} 完了")
L = pd.DataFrame(lolo)
if _gf:
    _priv = os.environ.get("APSC_PRIVATE", f"{W}/docs/private/tables")
    os.makedirs(_priv, exist_ok=True)
    L.to_csv(f"{_priv}/28_lolo_grant_factors.csv", index=False)
per_tf = L.groupby(["axis", "TF"]).is_candidate.agg(["sum", "size"]).reset_index() if len(L) else pd.DataFrame(
    columns=["axis", "TF", "sum", "size"])
per_tf["all5_runs"] = per_tf["sum"] == per_tf["size"] if len(per_tf) else []
n_all = int(per_tf.all5_runs.sum()) if len(per_tf) else 0
if _gf:
    add("leave-one-lesion-out: number of factors among 10 that remained as candidates in all 5 rounds", f"{n_all} / 10",
        "28_lolo_grant_factors.csv",
    "Candidate = at self-axis 'all 4 lesions used have same sign' and top 20 by |pooled ρ|. inflammatory-side 5 factors judged by lf_A, tissue-forming-side 5 factors judged by lf_D")
    for r in per_tf.itertuples():
        add(f"LOLO candidate count {r.TF} ({r.axis})", f"{int(r.sum)}/{int(r.size)}", "28_lolo_grant_factors.csv")
else:
    print("LOLO skipped: docs/private/grant_factors.txt is absent (grant material is not published)")
print(f"\n全5回候補: {n_all}/10\n{per_tf[['axis','TF','sum','size','all5_runs']].to_string(index=False)}")

# ============ 2. Permutation test for enrichment (per lesion, 10,000 times) ============
O = pd.read_parquet(f"{OUT}/all_cells_discovery_obs_scored.parquet")
NPERM = 10000
perm = []
for les, g in list(O.groupby("lesion", observed=True)) + [("POOLED_discovery", O)]:
    sco = g["p_D_tissueforming"].values
    thr = np.quantile(sco, 0.90)
    intop = sco >= thr
    k = int(intop.sum())
    for comp in ("stromal", "immune"):
        lab = (g.compartment == comp).values
        base = lab.mean(); obsf = lab[intop].mean()
        fold = obsf / base if base > 0 else np.nan
        null = np.empty(NPERM)
        for i in range(NPERM):
            null[i] = lab[rng.permutation(len(lab))[:k]].mean()
        p = (1 + int((null >= obsf).sum())) / (1 + NPERM)
        perm.append(dict(lesion=str(les), compartment=comp, n_cells=len(g), n_top_decile=k,
                         pct_overall=round(100 * base, 2), pct_top_decile=round(100 * obsf, 2),
                         enrichment_fold=round(fold, 3), null_mean_pct=round(100 * null.mean(), 2),
                         null_p975_pct=round(100 * np.quantile(null, 0.975), 2),
                         p_permutation=p, n_perm=NPERM, one_sided="greater", seed=SEED))
    print(f"  置換検定 {les} 完了")
PM = pd.DataFrame(perm)
PM.to_csv(f"{TAB}/28_enrichment_permutation.csv", index=False)
for r in PM[PM.compartment == "stromal"].itertuples():
    add(f"Enrichment permutation test p ({r.lesion})", f"{r.p_permutation:.4g} ({r.enrichment_fold}-fold, null mean {r.null_mean_pct}%)",
        "28_enrichment_permutation.csv", f"{NPERM} permutations, one-sided (greater), seed={SEED}")
imm = PM[(PM.compartment == "immune") & (PM.lesion == "POOLED_discovery")].iloc[0]
add("immune share of the top decile (discovery pooled)",
    f"{imm.pct_top_decile}% ({imm.pct_overall}% overall, {imm.enrichment_fold}-fold, permutation p={imm.p_permutation:.4g})",
    "28_enrichment_permutation.csv")
print(f"\n{PM[['lesion','compartment','pct_overall','pct_top_decile','enrichment_fold','p_permutation']].to_string(index=False)}")

# ============ 3. Aggregation from the existing table ============
V25 = pd.read_csv(f"{TAB}/25_validation_states_within_vs_ingest.csv")
w = V25[V25.labelling == "within_cohort"]
for _, r in w.iterrows():
    add(f"Validation state proportions ({r['lesion']}, within-cohort definition)",
        f"inflammatory {r['pct_inflammatory']}% / tissue-forming {r['pct_tissue-forming']}% / "
        f"intermediate {r['pct_intermediate']}% / other {r['pct_other']}%",
        "25_validation_states_within_vs_ingest.csv", f"n={int(r['n_cells'])}")

DF = pd.read_csv(f"{TAB}/17_receptors_detection_fraction.csv")
RECEP = ["IL1R1", "IL1RAP", "OSMR", "LIFR", "TNFRSF1A", "TNFRSF1B", "IL6R", "IL6ST"]
for gene in RECEP:
    d = DF[DF.gene == gene]
    if not len(d): continue
    txt = " / ".join(f"{r.lesion}:{r.state[:4]} {100*r.detection_fraction:.1f}%" for r in d.itertuples())
    add(f"Receptor detection fraction {gene} (by state and lesion)", txt, "17_receptors_detection_fraction.csv",
        f"detection_fraction = fraction of cells with log-normalised expression > 0; {d.lesion.nunique()} lesions x {d.state.nunique()} states")

G13 = pd.read_csv(f"{TAB}/13_gse164241_summary.csv")
for grp, g in G13.groupby("group"):
    pooled = g[g.donor == "POOLED"]
    per = g[g.donor != "POOLED"].dropna(subset=["depth_controlled_rho"])
    add(f"GSE164241 depth-controlled rho ({grp})",
        (f"pooled {float(pooled.depth_controlled_rho.iloc[0]):+.4f}" if len(pooled) else "no pooled value") +
        (f", per donor {per.depth_controlled_rho.min():+.4f} to {per.depth_controlled_rho.max():+.4f}" if len(per) else ""),
        "13_gse164241_summary.csv",
        f"{len(per)} donors, {int(per.n_cells.sum()):,} stromal cells")
    add(f"GSE164241 donor count ({grp})", f"{len(per)} donors ({int(per.n_cells.sum()):,} stromal cells)", "13_gse164241_summary.csv")

RV = pd.read_csv(f"{TAB}/05_rho_inflam_vs_tissueforming_validation.csv")
sv = RV[(RV.cellset == "stromal") & (RV.definition == "axis")]
for r in sv.itertuples():
    add(f"Validation single-axis (lf_) depth-controlled rho ({r.lesion})", f"{r.depth_controlled_rho:+.4f}（raw {r.raw_rho:+.4f}、n={r.n}）",
        "05_rho_inflam_vs_tissueforming_validation.csv")

# ============ 4. Cell counts at all QC stages (both cohorts) ============
qs = []
for f, coh in [("qc_cells_by_lesion_discovery.csv", "discovery_n5"), ("05_qc_cells_by_lesion_validation.csv", "validation_n3")]:
    p = f"{TAB}/{f}"
    if os.path.exists(p):
        d = pd.read_csv(p); d["cohort"] = coh; d["source_file"] = f; qs.append(d)
if qs:
    Q = pd.concat(qs, ignore_index=True)
    Q.to_csv(f"{TAB}/28_qc_all_stages.csv", index=False)
    for coh, g in Q.groupby("cohort"):
        stages = list(dict.fromkeys(g.stage))
        add(f"Cell count by QC stage ({coh})",
            " → ".join(f"{s}: {int(g[g.stage == s].n_cells.sum()):,}" for s in stages),
            "28_qc_all_stages.csv", f"per sample, see the table ({g.lesion.nunique()} samples x {len(stages)} stages)")
    print(f"\nQC 表: {Q.shape}, 段階 {sorted(Q.stage.unique())}")

M = pd.DataFrame(rows)
M.to_csv(f"{TAB}/28_manuscript_numbers.csv", index=False)
print(f"\n28_manuscript_numbers.csv: {len(M)} 項目")
print(M[["item", "value"]].head(12).to_string(index=False))
print("\n28 DONE")
