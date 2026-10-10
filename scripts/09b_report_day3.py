#!/usr/bin/env python
"""Purpose: Write REPORT_day3.md from 09 results; prepend section to dev_log_ja.md; append kb-row to methods
  CSV
Input: 16_methods_numbers.csv
Output: REPORT_day3.md, docs/dev_log_ja.md, 16_methods_numbers.csv

09b: Writes REPORT_day3.md from the results of 09, adds a section at the top of docs/dev_log_ja.md,
appends a kb-version row to 16_methods_numbers.csv. Does not create figures.
"""
import os, sys, json, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".numba_cache"))
import pandas as pd, numpy as np, pathlib

W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
now = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%MZ")


def rd(n):
    p = f"{TAB}/{n}"
    return pd.read_csv(p) if os.path.exists(p) else None


P, T, J = rd("09_geo_vs_kb.csv"), rd("09_tf_agreement_geo_vs_kb.csv"), rd("09_stromal_jaccard.csv")
if P is None:
    print("09_geo_vs_kb.csv が無い — 09 が走っていない。中止"); sys.exit(1)
have_kb = "kb_allbarcodes" in P.columns
have_kbgeo = "kb_geobarcodes" in P.columns
cols = ["quantity", "lesion", "geo"] + [c for c in ("kb_allbarcodes", "kb_geobarcodes") if c in P.columns]


def md(df, keep=None, r=4):
    d = df if keep is None else df[df.quantity.isin(keep)]
    d = d[cols].copy()
    for c in d.columns[2:]:
        d[c] = pd.to_numeric(d[c], errors="coerce").round(r)
    head = "| " + " | ".join(d.columns) + " |\n|" + "---|" * len(d.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
                                   for row in d.itertuples(index=False))


parts = [f"""# REPORT_day3 — kb 再定量版と GEO 行列版の対照

{now}。08（探索5検体の kb 再定量）完了後に `99_watch_09.sh` が自動起動した工程の結果。
図は作っていない（指示どおり）。

## 比較の枠組みと自己判断

kb の出力は `counts_unfiltered`（全バーコード）で、GEO 行列は作成者側で cell calling 済み。
そのため kb をそのまま使うと「定量の差」と「cell calling の差」が混ざる。2通り並べて切り分けた。

| 列 | 意味 |
|---|---|
| `geo` | GEO の filtered 行列（本解析の主系列） |
| `kb_allbarcodes` | kb 全バーコードに 03 と同一の QC（>=500 UMI, >=200 遺伝子, MT<=20%）。cell calling は QC 閾値が兼ねる |
| `kb_geobarcodes` | kb のうち GEO 行列に含まれるバーコードだけ。**定量の差だけを見る系列** |

**自己判断⑯**: 指示の監視対象は `data/counts/Sample*/verify.json` だったが、08 の出力先は
`data/counts_discovery/`（`data/counts/` は検証コホート CAP1–3 用）。後者を監視した。

**自己判断⑰**: GEO のバーコードは CellRanger の `-1` 付き、kb は素の16mer。接尾辞を外して突合する
（付けたままでは一致 0 件。SampleA で 8,107/8,107 = 100% 一致を事前確認した）。

**自己判断⑱**: kb 側は連結前に細胞側 QC 閾値で絞った。QC は細胞単位なので 03 本来の順序
（連結 → QC → filter_genes）と結果は同一で、491k バーコード×5 をメモリに抱えずに済む。

---

## 1. 細胞数と区画比率
"""]
parts.append(md(P, ["n_cells_after_QC", "n_stromal", "pct_stromal", "pct_immune"], 2))
parts.append("\n## 2. 濃縮・ρ・定義A\n")
parts.append(md(P, ["enrichment_fold_tissueforming", "stromal_pct_top_decile",
                    "depth_controlled_rho[axis]", "depth_controlled_rho[composite]",
                    "defA_frac_retained"], 4))
if T is not None:
    parts.append("\n## 3. 候補TFの一致\n\n| 系列 | 軸 | 候補 geo | 候補 kb | 重複 | Jaccard | 符号一致% | pooled ρ の Spearman | geo のみ | kb のみ |\n|---|---|---|---|---|---|---|---|---|---|\n")
    parts.append("\n".join(
        f"| {r.source} | {r.axis} | {r.n_candidates_geo} | {r.n_candidates_kb} | {r.overlap} | {r.jaccard} | "
        f"{r.same_sign_pct} | {r.spearman_rho_of_pooled_rho} | {r.geo_only if isinstance(r.geo_only,str) else ''} | "
        f"{r.kb_only if isinstance(r.kb_only,str) else ''} |" for r in T.itertuples()))
if J is not None:
    parts.append("\n\n## 4. 間質バーコードの一致（Jaccard）\n\n| 比較 | 病変 | n geo | n 他 | 共通 | Jaccard | geo の回収率% |\n|---|---|---|---|---|---|---|\n")
    parts.append("\n".join(f"| {r.comparison} | {r.lesion} | {r.n_geo} | {r.n_other} | {r.intersection} | "
                           f"{r.jaccard} | {r.pct_of_geo_recovered} |" for r in J.itertuples()))

# How to read (values are mechanically pulled from the table)
def val(q, les, col):
    m = P[(P.quantity == q) & (P.lesion == les)]
    return float(m[col].iloc[0]) if len(m) and col in m.columns and pd.notna(m[col].iloc[0]) else np.nan

notes = ["\n\n## 5. 読み方\n"]
if have_kbgeo:
    ge, kg = val("enrichment_fold_tissueforming", "POOLED_discovery", "geo"), val("enrichment_fold_tissueforming", "POOLED_discovery", "kb_geobarcodes")
    if kg == kg:
        notes.append(f"- 同一バーコードでの濃縮倍率は geo {ge:.3f} に対し kb {kg:.3f}（差 {kg-ge:+.3f}）。\n")
if J is not None and len(J):
    p = J[(J.lesion == "POOLED")]
    for r in p.itertuples():
        notes.append(f"- 間質集合の一致（{r.comparison}, pooled）: Jaccard {r.jaccard}、geo 側 {r.n_geo} のうち {r.intersection} を回収。\n")
if T is not None and len(T):
    for r in T[T.axis.isin(["lf_A_cytokine_chemokine", "lf_D_tissueforming"])].itertuples():
        notes.append(f"- 候補TF（{r.source}, {r.axis}）: 重複 {r.overlap}/{r.n_candidates_geo}、符号一致 {r.same_sign_pct}%。\n")
notes.append("\n数値はすべて `results/tables/09_*.csv` から機械的に引いており、転記していない。\n")
parts += notes
if os.path.exists(f"{TAB}/09_problems.txt"):
    parts.append("\n## 取れなかったブロック\n\n```\n" + open(f"{TAB}/09_problems.txt").read() + "\n```\n")
pathlib.Path(f"{W}/REPORT_day3.md").write_text("".join(parts))

# ---- Append to 16_methods_numbers.csv ----
MN = pd.read_csv(f"{TAB}/16_methods_numbers.csv")
add = []
for col, lab in [("kb_allbarcodes", "kb全BC"), ("kb_geobarcodes", "kbGEO-BC")]:
    if col not in P.columns: continue
    for q, ph in [("n_cells_after_QC", "[kb 細胞数]"), ("enrichment_fold_tissueforming", "[kb 濃縮倍率]"),
                  ("depth_controlled_rho[axis]", "[kb 単一軸ρ]"), ("defA_frac_retained", "[kb 定義A]")]:
        d = P[(P.quantity == q) & P[col].notna()]
        if not len(d): continue
        add.append(dict(placeholder=f"{ph[:-1]} {lab}]", source="09_geo_vs_kb.csv",
                        value=" / ".join(f"{r.lesion} {float(getattr(r, col)):.4g}（geo {float(r.geo):.4g}）"
                                         for r in d.itertuples() if pd.notna(r.geo))))
if J is not None and len(J):
    add.append(dict(placeholder="[kb 間質バーコード Jaccard]", source="09_stromal_jaccard.csv",
                    value=" / ".join(f"{r.comparison} {r.lesion} {r.jaccard}" for r in J[J.lesion == "POOLED"].itertuples())))
if T is not None and len(T):
    add.append(dict(placeholder="[kb 候補TFの一致]", source="09_tf_agreement_geo_vs_kb.csv",
                    value=" / ".join(f"{r.source} {r.axis} 重複{r.overlap} 符号一致{r.same_sign_pct}%"
                                     for r in T.itertuples())))
if add:
    MN = MN[~MN.placeholder.isin([a["placeholder"] for a in add])]
    MN = pd.concat([MN, pd.DataFrame(add)], ignore_index=True)
    MN.to_csv(f"{TAB}/16_methods_numbers.csv", index=False)

sec = f"""## 追記（{now}）: kb 再定量版と GEO 版の対照（09、自動実行）

08 が5検体を終えたのを `99_watch_09.sh` が検知し、03→04 を `APSC_SOURCE=kb` と `APSC_SOURCE=kbgeo`
の2系列で回して 09 で対照した（GEO 系列の出力は接尾辞なしのまま、上書きしていない）。
結果は `REPORT_day3.md`、表は `results/tables/09_*.csv`。

- 比較は2系列。`kb_allbarcodes`（kb 全バーコードに同一QC、cell calling は QC 閾値が兼ねる）と
  `kb_geobarcodes`（GEO のバーコードに限定、**定量の差だけ**を見る系列）。
- 自己判断⑯ 監視対象は `data/counts_discovery/`（指示の `data/counts/` は検証コホート用）。
- 自己判断⑰ GEO バーコードの `-1` 接尾辞を外して突合（外さないと一致 0 件）。
- 自己判断⑱ kb 側は連結前に細胞側 QC で絞った（順序を変えても結果は同一、メモリ対策）。

"""
q = pathlib.Path(f"{W}/docs/dev_log_ja.md"); q.write_text(sec + q.read_text())
print("REPORT_day3.md と docs/dev_log_ja.md を更新 |", len(pathlib.Path(f'{W}/REPORT_day3.md').read_text().encode()), "bytes")
print("09b DONE")
