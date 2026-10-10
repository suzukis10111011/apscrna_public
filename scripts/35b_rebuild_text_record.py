import os
W = os.environ.get("APSC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
"""Purpose: Rebuild text record by comparing old script (tag v0.9-draft-figtext) and current extraction; mark
  existing elements as KEEP, missing-only elements as MOVE; inherit borderline status from v32
Input: results/tables/32_figure_text_moved.csv
Output: results/tables/32_figure_text_moved.csv

旧スクリプト（tag v0.9-draft-figtext）と現行を同じ抽出器にかけ、
現行に在る要素を KEEP、旧にしか無い要素を MOVE として記録を作り直す。
borderline は 32 の前版の判定を引き継ぐ。
"""
import re, pathlib, pandas as pd, sys
sys.path.insert(0, "scripts")
NOISE = {"center","right","left","italic","offset points","_","\u2014","-|>","bold","both",
         "round,pad=0.06","upper right","lower center","upper center","data"}
PAT = [(r'\.set_title\(',"panel_title"),(r'\.set_xlabel\(',"axis_label"),(r'\.set_ylabel\(',"axis_label"),
       (r'\.set_xticklabels\(',"tick_labels"),(r'\.set_yticklabels\(',"tick_labels"),
       (r'\.annotate\(',"annotation"),(r'\.text\(',"annotation"),(r'\bpanel\(',"panel_letter"),
       (r'\.set_label\(',"colorbar_label"),(r'suptitle\(',"suptitle"),(r'label=',"legend_key")]
STR = re.compile(r'f?"((?:[^"\\]|\\.)*)"' + r"|f?'((?:[^'\\]|\\.)*)'")

def logical_lines(src):
    out, buf, st, d = [], "", None, 0
    for i, l in enumerate(src.split("\n"), 1):
        if st is None: st = i
        buf += l + " "; d += l.count("(")+l.count("[")+l.count("{"); d -= l.count(")")+l.count("]")+l.count("}")
        if d <= 0: out.append((st, buf.strip())); buf, st, d = "", None, 0
    if buf.strip(): out.append((st, buf.strip()))
    return out

def extract(path, fig, script):
    """1論理行に同じ呼び出しが複数あっても全て拾う。
    以前は最初の一致で break していたため、1行に2つ label= を書いた凡例
    （図2a の Discovery / Validation）の2つ目が記録から落ちていた。"""
    rows, seen = [], set()
    for lineno, ll in logical_lines(pathlib.Path(path).read_text()):
        if ll.lstrip().startswith("#"): continue
        for pat, et in PAT:
            for m in re.finditer(pat, ll):
                sm = STR.search(ll[m.end():])
                if not sm: continue
                t = sm.group(1) if sm.group(1) is not None else sm.group(2)
                if not t.strip() or t.strip() in NOISE or len(t) > 130: continue
                k = (et, t, lineno)
                if k in seen: continue
                seen.add(k)
                rows.append(dict(figure=fig, element_type=et, original_text=t, script=script, line=lineno))
    return rows

FIG = {}
for f in sorted(pathlib.Path("scripts").glob("ms_fig*.py")):
    m = re.search(r'save\(f2?,\s*"([A-Za-z0-9_]+)"\)', f.read_text())
    if m: FIG[f.name] = m.group(1).split("_")[0]

pre, now = [], []
for fn, fig in FIG.items():
    pre += extract(f"/tmp/orig/{fn}", fig, fn)
    now += extract(f"scripts/{fn}", fig, fn)
key = lambda r: (r["figure"], r["element_type"], r["original_text"])
now_keys = {key(r) for r in now}
moved = [r for r in pre if key(r) not in now_keys]

# Carry over the borderline / panel determinations from the previous version
prev = pd.read_csv("results/tables/32_figure_text_moved.csv", dtype=str)
prev["borderline"] = prev.borderline.str.upper()
BL = {(r.figure, r.element_type) for r in prev[prev.borderline == "TRUE"].itertuples()}
PANEL = {(r.figure, r.element_type, str(r.original_text)[:40]): r.panel for r in prev.itertuples()}

def panel_of(r):
    for (f, e, txt), p in PANEL.items():
        if f == r["figure"] and e == r["element_type"] and (txt.startswith(r["original_text"][:20]) or
                                                            r["original_text"][:20] in txt):
            return p if isinstance(p, str) else ""
    return ""

rows = []
for r in now:
    rows.append(dict(**r, panel=panel_of(r), action="KEEP",
                     borderline="TRUE" if (r["figure"], r["element_type"]) in BL else "FALSE",
                     restored_2026_09_29="TRUE" if (r["figure"], r["element_type"]) in BL else "FALSE"))
for r in moved:
    rows.append(dict(**r, panel=panel_of(r), action="MOVE",
                     borderline="TRUE" if (r["figure"], r["element_type"]) in BL else "FALSE",
                     restored_2026_09_29="FALSE"))

ORD = ["Fig1","Fig2","Fig3","Fig4","Fig5","Fig6","Fig7","FigS1","FigS2","FigS3","FigS4","FigS5","FigS6","FigS7","FigS8"]
D = pd.DataFrame(rows)[["figure","panel","element_type","original_text","action","borderline",
                        "restored_2026_09_29","script","line"]]
D["_o"] = D.figure.map({k: i for i, k in enumerate(ORD)})
D = D.sort_values(["_o","action","element_type","line"]).drop(columns="_o")
D.to_csv(f'{os.environ.get("APSC_TABDIR", f"{W}/results/tables")}/32_figure_text_moved.csv', index=False)
print(f"{len(D)} 行 | KEEP {int((D.action=='KEEP').sum())} / MOVE {int((D.action=='MOVE').sum())}")
print("MOVE のまま:")
print(D[D.action=="MOVE"][["figure","element_type","original_text"]].to_string(index=False, max_colwidth=58))
