import os
W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#!/usr/bin/env python
"""Purpose: Extract and classify text elements from paper figures; preserve kept text, record deleted text
  verbatim
Input: none (reads from results/tables/)
Output: results/tables/32_figure_text_moved.csv

35: 論文図の文字要素の棚卸しと分類記録。
KEEP は現行スクリプトから機械抽出し、MOVE は削除した原文を逐語で持つ。
数値は results/tables/ から読み、丸めない。
"""
import re, json, pathlib, pandas as pd

SCRIPTS = pathlib.Path("scripts")
NOISE = {"center", "right", "left", "italic", "offset points", "_", "\u2014", "-|>", "bold",
         "both", "round,pad=0.06", "upper right", "lower center", "upper center", "data"}

PATTERNS = [(r'\.set_title\(', "panel_title"), (r'\.set_xlabel\(', "axis_label"),
            (r'\.set_ylabel\(', "axis_label"), (r'\.set_xticklabels\(', "tick_labels"),
            (r'\.set_yticklabels\(', "tick_labels"), (r'\.annotate\(', "annotation"),
            (r'\.text\(', "annotation"), (r'\bpanel\(', "panel_letter"),
            (r'\.set_label\(', "colorbar_label"), (r'suptitle\(', "suptitle"),
            (r'label=', "legend_key")]
STR = re.compile(r'''f?"((?:[^"\\]|\\.)*)"|f?'((?:[^'\\]|\\.)*)'”''', re.X)
STR = re.compile(r'f?"((?:[^"\\]|\\.)*)"' + r"|f?'((?:[^'\\]|\\.)*)'")


def logical_lines(src):
    out, buf, start, depth = [], "", None, 0
    for i, ln in enumerate(src.split("\n"), 1):
        if start is None:
            start = i
        buf += ln + " "
        depth += ln.count("(") + ln.count("[") + ln.count("{")
        depth -= ln.count(")") + ln.count("]") + ln.count("}")
        if depth <= 0:
            out.append((start, buf.strip())); buf, start, depth = "", None, 0
    if buf.strip():
        out.append((start, buf.strip()))
    return out


def extract(fn, fig):
    rows, src = [], (SCRIPTS / fn).read_text()
    for lineno, ll in logical_lines(src):
        if ll.lstrip().startswith("#"):
            continue
        for pat, etype in PATTERNS:
            for m in re.finditer(pat, ll):
                sm = STR.search(ll[m.end():])
                if not sm:
                    continue
                txt = sm.group(1) if sm.group(1) is not None else sm.group(2)
                if not txt.strip() or txt.strip() in NOISE or len(txt) > 130:
                    continue
                rows.append(dict(figure=fig, panel="", element_type=etype, original_text=txt,
                                 action="KEEP", borderline="FALSE", script=fn, line=lineno))
                break
    return rows


FIG = {}
for f in sorted(SCRIPTS.glob("ms_fig*.py")):
    m = re.search(r'save\(f2?,\s*"([A-Za-z0-9_]+)"\)', f.read_text())
    if m:
        FIG[f.name] = m.group(1).split("_")[0]

keep = []
for fn, fig in FIG.items():
    keep += extract(fn, fig)

# Manually classified MOVE (the deleted original text is kept verbatim. line is the line number before deletion)
MOVE = [
 ("Fig1","a","annotation","one-way: definitions fixed in discovery and never revised from validation","FALSE","ms_fig1.py",37),
 ("Fig1","b,c","panel_title","{lab}  n={len(U):,}  （n 部分のみ）","TRUE","ms_fig1.py",51),
 ("Fig1","b,c","panel_title","Discovery, 5 lesions / Validation, 3 lesions  （病変数のみ）","TRUE","ms_fig1.py",44),
 ("Fig1","b,c","legend_key","{comp(c)} ({len(g):,})  （細胞数のみ）","TRUE","ms_fig1.py",47),
 ("Fig1","d","panel_title","Compartment markers (discovery)\\ndot size = % cells detected  （2行目と (discovery)）","FALSE","ms_fig1.py",76),
 ("Fig2","b","panel_title","Stromal share of the top decile\\n(permutation p = 1.0 $\\\\times$ 10$^{-4}$ (floor) in every lesion)","TRUE","ms_fig2.py",64),
 ("Fig2","c","legend_key","Discovery stromal ({len(sd):,}) / Validation, {s} ({len(g):,})  （細胞数のみ）","TRUE","ms_fig2.py",71),
 ("Fig2","c","panel_title","Validation stromal cells projected onto\\nthe discovery embedding","TRUE","ms_fig2.py",79),
 ("Fig3","a,b","panel_title","{prog(axis)}\\n{len(d)} factors; candidates same sign {same}/{len(cd)}  （2行目）","FALSE","ms_fig3.py",45),
 ("Fig3","e","annotation","{r.pooled}/{r.n} と {r.all3}/{r.n}（棒の値ラベル）","TRUE","ms_fig3.py",90),
 ("Fig3","f","panel_title","{n_all}/{len(cnt.dropna())} factors retained in all {nr} runs","FALSE","ms_fig3.py",111),
 ("Fig4","b","annotation","{r.frac_retained:.3f}（棒の値ラベル）","TRUE","ms_fig4.py",56),
 ("Fig4","b","annotation","ns（帰無帯と重なる病変の印）","TRUE","ms_fig4.py",59),
 ("Fig4","b","panel_title","Dispersion left after matching on the\\nupstream regulator's activity","FALSE","ms_fig4.py",69),
 ("Fig4","c","tick_labels","{r.lesion}\\n({int(r.n_cells):,})  （細胞数のみ）","TRUE","ms_fig4.py",85),
 ("Fig4","c","panel_title","Cell states, defined within each cohort  （説明句のみ）","FALSE","ms_fig4.py",89),
 ("Fig4","d","panel_title","{l}  n={len(g):,}  （n 部分のみ）","TRUE","ms_fig4.py",105),
 ("Fig4","e","annotation","p = 0.05（閾値線のラベル）","TRUE","ms_fig4.py",123),
 ("Fig4","e","panel_title","Both score distributions are unimodal in every lesion ({reject}/{n} reject)","FALSE","ms_fig4.py",128),
 ("Fig5","a","panel_title","Discovery, pseudobulk (5 lesions)  （病変数のみ）","TRUE","ms_fig5.py",19),
 ("Fig5","b","axis_label","dot size = % cells expressing","FALSE","ms_fig5.py",33),
 ("Fig5","d","legend_key","NF-κB block | osteo-mesenchymal block\\nNFKB1 and RELA share 83% of the smaller regulon","FALSE","ms_fig5.py",64),
 ("Fig6","a","panel_title","dots = lesions, bar = pooled","FALSE","ms_fig6.py",21),
 ("Fig2","a","legend_key","Discovery (5 lesions) / Validation (3 lesions)  （病変数のみ）","TRUE","ms_fig2.py",36),
 ("Fig4","a","tick_labels","Discovery (n=5) / Validation (n=3)  （n のみ）","TRUE","ms_fig4.py",38),
 ("Fig5","c","tick_labels","Discovery (n=5) / Validation (n=3)  （n のみ）","TRUE","ms_fig5.py",48),
 ("Fig7","a","annotation","pooled {y:+.4f}（群ごとの pooled 値）","TRUE","ms_fig7.py",40),
 ("Fig7","a","tick_labels","{GRPN.get(g, g)}\\n({n} donor{'s'})  （ドナー数のみ）","TRUE","ms_fig7.py",46),
 ("Fig7","a","panel_title","GSE164241 gingiva, {n:,} stromal cells\\nhorizontal bar = pooled estimate per group  （細胞数と2行目）","FALSE","ms_fig7.py",49),
 ("Fig7","b","annotation","{r.rep}/{r.n} と {r.cons}/{r.n}（棒の値ラベル）","TRUE","ms_fig7.py",65),
 ("FigS1","a,b","panel_title","Discovery (5 lesions) / Validation (3 lesions)  （病変数のみ）","TRUE","ms_figS1.py",14),
 ("FigS2","a","annotation","{v:,}（棒の値ラベル）","TRUE","ms_figS2.py",14),
 ("FigS2","a","panel_title","Jaccard {I.stromal_jaccard:.3f}","FALSE","ms_figS2.py",18),
 ("FigS2","b","panel_title","Agreement over {n} factors: Spearman {rho:.3f};\\ntop-20 overlap {k}/20","FALSE","ms_figS2.py",30),
 ("FigS2","b","annotation","{n_tf} factors here vs {n_screen} in the main screen: activities for this panel were computed on all cells, the screen on stromal cells only, so more regulons pass the minimum-target filter ({extra} extra).","FALSE","ms_figS2.py",32),
 ("FigS3","a","panel_title","Negative in {neg}/{n} gene-set × scoring combinations","FALSE","ms_figS3.py",19),
 ("FigS3","a","tick_labels","(this study)（集合名に付く注記の行）","TRUE","ms_figS3.py",13),
 ("FigS3","b","panel_title","{min:.3f}–{max:.3f}","FALSE","ms_figS3.py",26),
 ("FigS4","c","annotation","{overlap}/{n_candidates_geo} cand.（棒の値ラベル）","TRUE","ms_figS4.py",46),
 ("FigS6","a","panel_title","Rank on its own program axis\\n(top 5 inflammatory / bottom 5 tissue-forming)","FALSE","ms_figS6.py",26),
 ("FigS6","b","panel_title","{n_all}/{n} retained in all {runs} runs","FALSE","ms_figS6.py",37),
 ("FigS7","a","annotation","{agreement_pct:.1f}% agree","TRUE","ms_figS7.py",26),
 ("FigS7","a","panel_title","State labels are sensitive to the labelling route","FALSE","ms_figS7.py",31),
 ("FigS8","a","annotation","{r.rep}/{r.n}（棒の値ラベル）","TRUE","ms_figS8.py",25),
 ("FigS8","a","panel_title","Validation cohort (3 lesions)  （病変数のみ）","TRUE","ms_figS8.py",28),
 ("FigS8","b","annotation","{r.rep}/{r.n}（棒の値ラベル）","TRUE","ms_figS8.py",41),
]
mv = [dict(figure=f, panel=p, element_type=e, original_text=o, action="MOVE",
           borderline=b, script=s, line=l) for f, p, e, o, b, s, l in MOVE]

ORD = ["Fig1","Fig2","Fig3","Fig4","Fig5","Fig6","Fig7",
       "FigS1","FigS2","FigS3","FigS4","FigS5","FigS6","FigS7","FigS8"]
D = pd.DataFrame(keep + mv)
D["_o"] = D.figure.map({k: i for i, k in enumerate(ORD)})
D = D.sort_values(["_o", "action", "element_type", "line"]).drop(columns="_o")
D.to_csv(f'{os.environ.get("APSC_TABDIR", f"{W}/results/tables")}/32_figure_text_moved.csv', index=False)
print(f"32_figure_text_moved.csv: {len(D)} 行 "
      f"(KEEP {int((D.action=='KEEP').sum())} / MOVE {int((D.action=='MOVE').sum())}, "
      f"borderline {int((D.borderline=='TRUE').sum())})")
print(D.groupby(["figure", "action"]).size().unstack(fill_value=0)
       .reindex(ORD).to_string())
