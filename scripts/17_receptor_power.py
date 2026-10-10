import os
#!/usr/bin/env python
"""Purpose: Statistical power analysis for receptor expression across lesion states using linear mixed models
Input: none (loads preprocessed receptor expression data from prior analysis)
Output: {TAB}/17_receptors_detection_fraction.csv, {TAB}/17_receptors_mixedmodel.csv

17: 受容体比較の検定力。細胞単位の log 正規化発現 ~ state + (1|lesion) の線形混合モデル。
両状態が10細胞以上ある病変のみを対象にし、探索と検証を合わせて当てはめる。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("17")
# 25: On the validation side, use the within-cohort definition (state_within) as default, rather than ingest projection.
# because agreement with the projected labels is only 40.9%, and projection is sensitive to gene space and unstable.
_VAL_STATE = os.environ.get("APSC_VAL_STATE", "state_within")
import scanpy as sc, pandas as pd, numpy as np, scipy.sparse as sp
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
REC = ["TNFRSF1A","TNFRSF1B","IL1R1","IL1RAP","IL6R","IL6ST","OSMR","LIFR","CCR2","CXCR4"]
PAIR = ("inflammatory", "tissue-forming")
frames, det = [], []
for path, coh, col in [(f"{W}/data/processed/stromal_discovery_labelled.h5ad", "discovery_n5", "state"),
                       (f"{W}/data/processed/stromal_validation.h5ad", "validation_n3", _VAL_STATE)]:
    if not os.path.exists(path): continue
    A = sc.read_h5ad(path)
    if col not in A.obs: continue
    s2i = {}
    for gid, s in zip(A.var_names, A.var["symbol"].astype(str)): s2i.setdefault(s, gid)
    genes = [g for g in REC if g in s2i]
    X = A[:, [s2i[g] for g in genes]].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    d = pd.DataFrame(X, columns=genes)
    d["state"] = A.obs[col].astype(str).values
    d["lesion"] = coh[:4] + "_" + A.obs["lesion"].astype(str).values
    d["cohort"] = coh
    d = d[d.state.isin(PAIR)]
    keep = [l for l, g in d.groupby("lesion")
            if (g.state == PAIR[0]).sum() >= 10 and (g.state == PAIR[1]).sum() >= 10]
    print(f"{coh}: 両状態10細胞以上の病変 {len(keep)} / {d.lesion.nunique()} -> {keep}")
    for l, g in d.groupby("lesion"):
        for s in PAIR:
            gg = g[g.state == s]
            if len(gg) == 0: continue
            for gene in genes:
                det.append(dict(cohort=coh, lesion=l, state=s, gene=gene, n_cells=len(gg),
                                detection_fraction=round(float((gg[gene] > 0).mean()), 4),
                                mean_lognorm=round(float(gg[gene].mean()), 4),
                                included_in_model=l in keep))
    frames.append(d[d.lesion.isin(keep)])
DET = pd.DataFrame(det); DET.to_csv(f"{TAB}/17_receptors_detection_fraction.csv", index=False)
Dall = pd.concat(frames, ignore_index=True)
print(f"\nモデル対象: {len(Dall)} 細胞 / 病変 {Dall.lesion.nunique()} "
      f"({Dall.groupby('cohort').lesion.nunique().to_dict()})")
rows = []
for gene in REC:
    if gene not in Dall.columns: 
        rows.append(dict(gene=gene, status="absent from the matrix")); continue
    d = Dall[[gene, "state", "lesion"]].rename(columns={gene: "y"}).copy()
    d["is_infl"] = (d.state == PAIR[0]).astype(float)
    try:
        m = smf.mixedlm("y ~ is_infl", d, groups=d["lesion"]).fit(reml=True, method="lbfgs")
        ci = m.conf_int().loc["is_infl"]
        rows.append(dict(gene=gene, status="OK", n_cells=len(d), n_lesions=d.lesion.nunique(),
                         coef_inflammatory_minus_tissueforming=round(float(m.params["is_infl"]), 4),
                         ci_low=round(float(ci[0]), 4), ci_high=round(float(ci[1]), 4),
                         pvalue=float(m.pvalues["is_infl"]),
                         group_var=round(float(m.cov_re.iloc[0, 0]), 5), converged=bool(m.converged)))
    except Exception as e:
        rows.append(dict(gene=gene, status=f"失敗 {type(e).__name__}"))
R = pd.DataFrame(rows)
ok = R.pvalue.notna() if "pvalue" in R.columns else pd.Series(False, index=R.index)
R["padj_BH"] = np.nan
if ok.sum(): R.loc[ok, "padj_BH"] = multipletests(R.loc[ok, "pvalue"], method="fdr_bh")[1]
R["model"] = "MixedLM: log-normalised expression ~ state, random intercept = lesion (discovery + validation pooled)"
R.to_csv(f"{TAB}/17_receptors_mixedmodel.csv", index=False)
print("\n=== 線形混合モデル ===")
print(R[["gene","n_cells","n_lesions","coef_inflammatory_minus_tissueforming","ci_low","ci_high","pvalue","padj_BH"]].to_string(index=False))
require_heartbeat("17 完了時")
print("\n17 DONE")
