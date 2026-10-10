import os
#!/usr/bin/env python
"""Purpose: Bilateral TF co-activity within single cells; SMAD3 sensitivity, 11×11 correlation matrix, 25-pair
  lesion sign consistency in discovery and validation cohorts
Input: discovery_stromal.h5ad, validation_stromal.h5ad, collectri_regulons.csv
Output: 11_coactivity_smad3_sensitivity.csv, 11_tf_activity_spearman_matrix.csv,
  11_tf_crossside_pairs_sign_consistency.csv

11 (図5C/D): 同一細胞内の両側TF共活性。SMAD3 を含む/除くの感度、11x11 相関行列、
群間25ペアの病変間符号一致を、探索と検証の両方で出す。
使い方: python 11_coactivity.py [discovery_stromal.h5ad] [validation_stromal.h5ad]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("11")
# 25: On the validation side, use the within-cohort definition (state_within) as default, rather than ingest projection.
# because agreement with the projected labels is only 40.9%, and projection is sensitive to gene space and unstable.
_VAL_STATE = os.environ.get("APSC_VAL_STATE", "state_within")
import scanpy as sc, pandas as pd, numpy as np
from scipy.stats import spearmanr
from _utils import ulm_on_symbols
try:
    import diptest; HAVE_DIP = True
except Exception:
    HAVE_DIP = False

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
DISC = sys.argv[1] if len(sys.argv) > 1 else f"{W}/data/processed/stromal_discovery_labelled.h5ad"
VAL  = sys.argv[2] if len(sys.argv) > 2 else f"{W}/data/processed/stromal_validation.h5ad"
NET  = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")[["source","target","weight"]]
INFL = ["NFKB2","RELB","REL","ATF1","MTA1"]
TISS5 = ["SP7","DLX5","RUNX2","SMAD3","MYRF"]
TISS4 = [g for g in TISS5 if g != "SMAD3"]
ALL11 = INFL + TISS5 + ["MYC"]

def shape2d(x, y):
    Z = np.column_stack([(x - x.mean())/x.std(ddof=1), (y - y.mean())/y.std(ddof=1)])
    ev, evec = np.linalg.eigh(np.cov(Z.T)); o = np.argsort(ev)[::-1]; ev, evec = ev[o], evec[:, o]
    r = spearmanr(x, y)
    d = dict(n=len(x), spearman_rho=round(r.statistic, 4), spearman_p=r.pvalue,
             pc1_var_frac=round(float(ev[0]/ev.sum()), 4),
             pc1_angle_deg=round(float(np.degrees(np.arctan2(evec[1,0], evec[0,0])) % 180), 1))
    if HAVE_DIP:
        d["dip_p_infl"] = round(float(diptest.diptest(x)[1]), 4)
        d["dip_p_tiss"] = round(float(diptest.diptest(y)[1]), 4)
    return d

sens, mats, pairs = [], [], []
for path, coh, statecol in [(DISC, "discovery", "state"), (VAL, "validation", _VAL_STATE)]:
    if not os.path.exists(path): print(f"{coh}: {path} 無し — スキップ"); continue
    A = sc.read_h5ad(path)
    print(f"{coh}: {A.shape} X max {float(A.X.max()):.2f}")
    ACT = ulm_on_symbols(A, NET, tmin=5)
    have = [t for t in ALL11 if t in ACT.columns]
    print(f"  活性が得られたTF: {len(have)}/11", [t for t in ALL11 if t not in have] or "")
    D = pd.DataFrame(index=A.obs_names)
    D["lesion"] = A.obs["lesion"].astype(str).values
    D["state"] = A.obs[statecol].astype(str).values if statecol in A.obs else "NA"
    for t in have: D[t] = ACT[t].values
    D["infl5"] = D[[t for t in INFL if t in have]].mean(1)
    D["tiss5"] = D[[t for t in TISS5 if t in have]].mean(1)
    D["tiss4"] = D[[t for t in TISS4 if t in have]].mean(1)
    for variant, ycol in [("tissueforming_5TF_incl_SMAD3", "tiss5"), ("tissueforming_4TF_excl_SMAD3", "tiss4")]:
        for grp, lab in [("lesion", None), ("state", None)]:
            for k, g in D.groupby(grp, observed=True):
                if len(g) < 30: continue
                sens.append(dict(cohort=coh, variant=variant, group=grp, label=str(k),
                                 **shape2d(g.infl5.values, g[ycol].values)))
            sens.append(dict(cohort=coh, variant=variant, group=grp, label="POOLED",
                             **shape2d(D.infl5.values, D[ycol].values)))
    M = pd.DataFrame(spearmanr(D[have].values)[0], index=have, columns=have)
    M.index.name = "TF"; M = M.reset_index(); M.insert(0, "cohort", coh)
    mats.append(M)
    for a in [t for t in INFL if t in have]:
        for b in [t for t in TISS5 if t in have]:
            rs = [spearmanr(g[a], g[b]).statistic for _, g in D.groupby("lesion", observed=True) if len(g) >= 30]
            pairs.append(dict(cohort=coh, infl_TF=a, tiss_TF=b,
                              pooled_rho=round(float(spearmanr(D[a], D[b]).statistic), 4),
                              n_lesions=len(rs), n_positive=int(sum(r > 0 for r in rs)),
                              rho_min=round(min(rs), 4) if rs else np.nan,
                              rho_max=round(max(rs), 4) if rs else np.nan,
                              all_same_sign=bool(rs and (all(r > 0 for r in rs) or all(r < 0 for r in rs)))))
S = pd.DataFrame(sens).drop_duplicates(); S.to_csv(f"{TAB}/11_coactivity_smad3_sensitivity.csv", index=False)
pd.concat(mats).to_csv(f"{TAB}/11_tf_activity_spearman_matrix.csv", index=False)
PR = pd.DataFrame(pairs); PR.to_csv(f"{TAB}/11_tf_crossside_pairs_sign_consistency.csv", index=False)
print("\n=== SMAD3 感度（病変別ρ）===")
print(S[(S.group == "lesion")][["cohort","variant","label","n","spearman_rho","pc1_var_frac"]].round(4).to_string(index=False))
print("\n=== 群間ペア: 全病変で符号一致 ===")
print(PR.groupby("cohort").all_same_sign.agg(["sum","size"]).to_string())
require_heartbeat("11 完了時")
print("\n11 DONE")
