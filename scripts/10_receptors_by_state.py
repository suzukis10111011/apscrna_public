import os
#!/usr/bin/env python
"""Purpose: Compare receptor expression across disease states using pseudobulk analysis with paired Wilcoxon
  test (5 discovery / 3 validation lesions)
Input: discovery_stromal.h5ad, validation_stromal.h5ad (or kb-quantified equivalents)
Output: 10_receptors_by_state.csv, 10_receptors_dotplot_input.csv, 10_receptors_pseudobulk_log2cpm.csv

10 (図5A/B): 受容体発現を state 間で比較する。
pseudobulk = 病変×状態で counts layer を合計 -> CPM -> log2(CPM+1)。
病変を単位とした対応のある Wilcoxon（探索5病変 / 検証3病変）、BH 補正。
使い方: python 10_receptors_by_state.py [discovery_stromal.h5ad] [validation_stromal.h5ad]
再定量後は kb 版のオブジェクトを引数で渡せば同じ表が出る。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("10")
# 25: On the validation side, use the within-cohort definition (state_within) as default, rather than ingest projection.
# because agreement with the projected labels is only 40.9%, and projection is sensitive to gene space and unstable.
_VAL_STATE = os.environ.get("APSC_VAL_STATE", "state_within")
import scanpy as sc, pandas as pd, numpy as np, scipy.sparse as sp
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
DISC = sys.argv[1] if len(sys.argv) > 1 else f"{W}/data/processed/stromal_discovery_labelled.h5ad"
VAL  = sys.argv[2] if len(sys.argv) > 2 else f"{W}/data/processed/stromal_validation.h5ad"
REC = ["TNFRSF1A","TNFRSF1B","IL1R1","IL1RAP","IL6R","IL6ST","OSMR","LIFR","CCR2","CXCR4"]

def run(path, cohort, statecol, pair=("inflammatory","tissue-forming")):
    A = sc.read_h5ad(path)
    if statecol not in A.obs:
        print(f"{cohort}: 状態列 {statecol} が無い — スキップ"); return None, None
    lay = "counts" if "counts" in A.layers else None
    if lay is None:
        print(f"{cohort}: counts layer が無い — pseudobulk は正規化値から作る（注記）")
    s2i = {}
    for gid, s in zip(A.var_names, A.var["symbol"].astype(str)): s2i.setdefault(s, gid)
    genes = [g for g in REC if g in s2i]
    miss = [g for g in REC if g not in s2i]
    ids = [s2i[g] for g in genes]
    X = A[:, ids].layers[lay] if lay else A[:, ids].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    Xall = A.layers[lay] if lay else A.X
    tot = np.asarray(Xall.sum(1)).ravel()
    st = A.obs[statecol].astype(str).values; les = A.obs["lesion"].astype(str).values
    pb, cell = [], []
    for l in sorted(set(les)):
        for s in pair:
            m = (les == l) & (st == s)
            if m.sum() < 10:
                pb.append(dict(cohort=cohort, lesion=l, state=s, n_cells=int(m.sum()),
                               **{g: np.nan for g in genes})); continue
            summed = X[m].sum(0); lib = tot[m].sum()
            cpm = summed / lib * 1e6
            pb.append(dict(cohort=cohort, lesion=l, state=s, n_cells=int(m.sum()),
                           **{g: float(np.log2(c + 1)) for g, c in zip(genes, cpm)}))
            for j, g in enumerate(genes):
                v = X[m][:, j]
                cell.append(dict(cohort=cohort, lesion=l, state=s, gene=g, n_cells=int(m.sum()),
                                 pct_expressing=round(100 * float((v > 0).mean()), 2),
                                 mean_expr_log2cpm=round(float(np.log2(v.sum()/lib*1e6 + 1)), 4)))
    P = pd.DataFrame(pb)
    rows = []
    for g in genes:
        a = P[P.state == pair[0]].set_index("lesion")[g]
        b = P[P.state == pair[1]].set_index("lesion")[g]
        idx = a.index.intersection(b.index); a, b = a[idx].dropna(), b[idx].dropna()
        idx = a.index.intersection(b.index); a, b = a[idx], b[idx]
        if len(a) >= 3:
            try: st_, p = wilcoxon(a.values, b.values)
            except Exception: st_, p = np.nan, np.nan
        else: st_, p = np.nan, np.nan
        rows.append(dict(cohort=cohort, gene=g, n_lesions=len(a),
                         mean_log2cpm_inflammatory=round(float(a.mean()), 4) if len(a) else np.nan,
                         mean_log2cpm_tissueforming=round(float(b.mean()), 4) if len(b) else np.nan,
                         delta_infl_minus_tiss=round(float((a - b).mean()), 4) if len(a) else np.nan,
                         wilcoxon_stat=st_, pvalue=p))
    R = pd.DataFrame(rows)
    ok = R.pvalue.notna()
    R["padj_BH"] = np.nan
    if ok.sum(): R.loc[ok, "padj_BH"] = multipletests(R.pvalue[ok], method="fdr_bh")[1]
    R["pseudobulk"] = "counts summed per lesion x state -> CPM -> log2(CPM+1)" + ("" if lay else "(no counts layer; normalised values used)")
    R["genes_absent_from_matrix"] = ",".join(miss) if miss else ""
    return R, (P, pd.DataFrame(cell))

out_stats, out_pb, out_cell = [], [], []
for path, coh, col in [(DISC, "discovery_n5", "state"), (VAL, "validation_n3", _VAL_STATE)]:
    if not os.path.exists(path): print(f"{coh}: {path} が無い — スキップ"); continue
    R, extra = run(path, coh, col)
    if R is None: continue
    out_stats.append(R); out_pb.append(extra[0]); out_cell.append(extra[1])
    print(f"\n=== {coh} ===")
    print(R[["gene","n_lesions","mean_log2cpm_inflammatory","mean_log2cpm_tissueforming",
             "delta_infl_minus_tiss","pvalue","padj_BH"]].round(4).to_string(index=False))
pd.concat(out_stats).to_csv(f"{TAB}/10_receptors_by_state.csv", index=False)
pd.concat(out_pb).to_csv(f"{TAB}/10_receptors_pseudobulk_log2cpm.csv", index=False)
pd.concat(out_cell).to_csv(f"{TAB}/10_receptors_dotplot_input.csv", index=False)
require_heartbeat("10 完了時")
print("\n10 DONE")
