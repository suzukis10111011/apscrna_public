"""Purpose: Define manuscript figure specifications and styling for results/figures_ms/; reads data from
  results/tables/*.csv and applies filtering
Input: none (loads styling parameters and CSV source paths)
Output: {OUT}/{name}.pdf, {OUT}/{name}.png, {TAB}/_msn_{name}.csv

Common specification for manuscript figures (results/figures_ms/).
Figures must always read numerical values from results/tables/*.csv. "grant" / "previous set" /
"previous stromal set" / "v2" terms and series must not appear in the main-text figures.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
BARE = os.environ.get("BARE") == "1"          # Raw figure version (PowerPoint base): hide descriptive text
OUT = f"{W}/results/" + os.environ.get("APSC_MSFIGDIR", "figures_ms") + ("/bare" if BARE else "")
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
os.makedirs(OUT, exist_ok=True)
MM = 1 / 25.4

DISC, VAL = "#2166ac", "#e08214"                     # Discovery = blue tones / Validation = orange tones
# "other" is changed to tan because gray would be indistinguishable from the discovery background (gray) in Fig2c.
# To avoid the same label having different colors across figures, replace STATE itself
# (the "other" in Fig4c and FigS7 also becomes the same tan).
STATE = {"inflammatory": "#440154", "intermediate": "#31688e",
         "tissue-forming": "#35b779", "other": "#C8A96E"}   # viridis-based palette
COMPC = {"stromal": "#e08214", "immune": "#2166ac", "endothelial": "#1b7837",
         "pericyte_smc": "#c994c7", "epithelial": "#80cdc1", "schwann_neural": "#8c6bb1",
         "immune_myeloid": "#4393c3", "immune_T": "#92c5de", "immune_plasma": "#053061",
         "unassigned": "#bdbdbd", "erythrocyte": "#f4a582"}

# Notation for readers (axes, panels, gene sets)
PROG = {"lf_A_cytokine_chemokine": "Inflammatory programme",
        "lf_D_tissueforming": "Tissue-forming programme",
        "A_cytokine_chemokine": "Inflammatory programme",
        "D_tissueforming": "Tissue-forming programme",
        "p_A_cytokine_chemokine": "Inflammatory programme",
        "p_D_tissueforming": "Tissue-forming programme",
        "lf_B_survival_antiapoptotic": "Survival / anti-apoptotic programme",
        "lf_C1_proliferation": "Proliferation programme",
        "lf_C2_migration": "Migration programme",
        "lf_C_prolif_migration": "Proliferation + migration programme",
        "B_survival_antiapoptotic": "Survival / anti-apoptotic programme",
        "C1_proliferation": "Proliferation programme",
        "C2_migration": "Migration programme",
        "C_prolif_migration": "Proliferation + migration programme"}
COMP = {"stromal": "Stromal cells", "immune": "Immune cells", "endothelial": "Endothelium",
        "pericyte_smc": "Pericytes/SMC", "epithelial": "Epithelium",
        "schwann_neural": "Schwann / neural", "immune_myeloid": "Myeloid cells",
        "immune_T": "T cells", "immune_plasma": "Plasma cells",
        "unassigned": "Unassigned", "erythrocyte": "Erythrocytes"}
GSET = {"current_lf_A_cytokine_chemokine": "Inflammatory programme\n(this study)",
        "current_lf_D_tissueforming": "Tissue-forming programme\n(this study)",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "HALLMARK_TNFA_\nSIGNALING_VIA_NFKB",
        "GOBP_INFLAMMATORY_RESPONSE": "GOBP_INFLAMMATORY_\nRESPONSE",
        "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION": "REACTOME_EXTRACELLULAR_\nMATRIX_ORGANIZATION"}
OKABE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#999999"]
LESION_ORDER = ["A", "B", "C", "D", "E", "CAP1", "CAP2", "CAP3"]


def prog(x):
    return PROG.get(str(x), str(x).replace("NF-kB", "NF-\u03baB"))


def comp(x):
    return COMP.get(str(x), str(x).replace("_", " "))


def gset(x):
    return GSET.get(str(x), str(x))


def setup():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.dpi": 110,
        "pdf.fonttype": 42, "ps.fonttype": 42})


setup()


def fig(width_mm=180, height_mm=100):
    f, a = plt.subplots(figsize=(width_mm * MM, height_mm * MM))
    return f, a


def gridfig(width_mm=180, height_mm=120, **kw):
    f = plt.figure(figsize=(width_mm * MM, height_mm * MM))
    return f, f.add_gridspec(**kw)


def panel(ax, letter, dx=-0.13, dy=1.06):
    # Panel labels are unified in lowercase (matching the notation in the main text and legend). Position, font, and weight remain as before.
    # In BARE mode, create the Text first and then make it invisible (take the pptx text box coordinates from here).
    return ax.text(dx, dy, str(letter).lower(), transform=ax.transAxes, fontsize=9,
                   fontweight="bold", va="top", ha="left", _kind="panel_letter")



# ===== Plain figure version (BARE=1): "hide without deleting" the explanatory text =====
# Policy: rather than deleting the text, create Text/Legend and then set_visible(False).
#   - Non-BARE output does not change at all (the patch only removes bare_keep / _kind).
#   - Even when invisible, window_extent can still be obtained, so the pptx text box coordinates can be output precisely.
# For text to keep ([1]), the caller passes bare_keep=True: ticks, axis labels, and colorbars are
# Not subject to patching in the first place; numeric values within heatmap cells and bar value labels, ns, p=0.05 are bare_keep.
import re as _re
import matplotlib.axes as _maxes
import matplotlib.figure as _mfigure

_TEXTMETA = []          # Record for 1 figure's worth. Flush to JSON with save() and clear.
_PARENS = _re.compile(r"\s*[\n ]?\((?:n=)?[^()]*\)$")     # "Discovery (n=5)" / "A\n(2,314)"


def move_ticklabels(ax, which="y"):
    """その軸の目盛ラベルを素図版では消し、pptx のテキストボックスに移す。
    位置は _dump_textmeta で各ラベルの外形から取る（先に消すと外形がゼロになる）。"""
    _TEXTMETA.append(dict(artist=None, kind="tick_moved", ax=ax, keep=False,
                          extra=dict(which=which)))
    return ax


def _rec(artist, kind, ax=None, keep=False, extra=None):
    """記録だけする。非表示にするのは save() の直前（_dump_textmeta の最後）。
    不可視の Text は get_window_extent がゼロ矩形を返すため、先に隠すと
    pptx のテキストボックス座標が全て (0,0) になってしまう。"""
    _TEXTMETA.append(dict(artist=artist, kind=kind, ax=ax, keep=keep, extra=extra or {}))
    return artist


def _patch():
    ax_text, ax_ann = _maxes.Axes.text, _maxes.Axes.annotate
    ax_title, ax_leg = _maxes.Axes.set_title, _maxes.Axes.legend
    ax_xtl, ax_ytl = _maxes.Axes.set_xticklabels, _maxes.Axes.set_yticklabels
    fig_text, fig_sup = _mfigure.Figure.text, _mfigure.Figure.suptitle

    def text(self, *a, **kw):
        keep = kw.pop("bare_keep", False); kind = kw.pop("_kind", "annotation")
        if kw.pop("numbox", False):
            # Numeric values (bar value labels, heatmap cells) are removed in the plain figure, and on the pptx side
            # placed in a centered, zero-margin box. The reason kind is separated is for that formatting.
            keep, kind = False, "number"
        return _rec(ax_text(self, *a, **kw), kind, self, keep)

    def annotate(self, *a, **kw):
        keep = kw.pop("bare_keep", False); kind = kw.pop("_kind", "annotation")
        art = ax_ann(self, *a, **kw)
        # The arrow itself is retained even in the plain figure ([2] removes only the "arrow label").
        # Emptying the label is done after measuring the coordinates.
        return _rec(art, kind, self, keep, dict(keep_arrow=getattr(art, "arrow_patch", None) is not None))

    def set_title(self, label, *a, **kw):
        # set_title returns Text. loc="left"/"right" are artists separate from ax.title
        # (_left_title / _right_title), so failing to use the return value will fail to remove them.
        kw.pop("bare_keep", None)
        art = ax_title(self, label, *a, **kw)
        return _rec(art if art is not None else self.title, "panel_title", self, False)

    def legend(self, *a, **kw):
        lg = ax_leg(self, *a, **kw)
        ent = []
        try:
            for h, lab in zip(lg.legend_handles, [x.get_text() for x in lg.get_texts()]):
                ent.append(dict(label=lab, colour=_hex_of(h), marker=_marker_of(h)))
        except Exception:
            pass
        return _rec(lg, "legend", self, False, dict(entries=ent))

    def _strip(self, labels, setter, axis):
        """BARE では群名に付く n / 病変数の括弧だけを落とす（群名は目盛として残す）。"""
        out, dropped = [], []
        for s in labels:
            s0 = s.get_text() if hasattr(s, "get_text") else str(s)
            s1 = _PARENS.sub("", s0)
            if s1 != s0:
                dropped.append(dict(full=s0, kept=s1))
            out.append(s1 if BARE else s0)
        if dropped:
            _TEXTMETA.append(dict(artist=None, kind=f"tick_paren_{axis}", ax=self,
                                  keep=False, extra=dict(dropped=dropped)))
        return setter(self, out, *(), **{})

    def _ticks(self, labels, setter, which, a, kw):
        """BARE では群名に付く n / 病変数の括弧だけを落とす（群名は目盛として残す）。
        落とした括弧は目盛の index つきで記録し、pptx はその目盛位置に置く。
        fontsize 等の kwargs 付きで呼ばれても記録する（以前は素通りしていた）。"""
        txt = [(s.get_text() if hasattr(s, "get_text") else str(s)) for s in labels]
        if not BARE:
            return setter(self, labels, *a, **kw)
        strip = [_PARENS.sub("", s) for s in txt]
        drop = [dict(i=i, full=f, kept=k) for i, (f, k) in enumerate(zip(txt, strip)) if f != k]
        if drop:
            _TEXTMETA.append(dict(artist=None, kind=f"tick_paren_{which}", ax=self,
                                  keep=False, extra=dict(dropped=drop, which=which)))
        return setter(self, strip, *a, **kw)

    def set_xticklabels(self, labels, *a, **kw):
        return _ticks(self, labels, ax_xtl, "x", a, kw)

    def set_yticklabels(self, labels, *a, **kw):
        return _ticks(self, labels, ax_ytl, "y", a, kw)

    def f_text(self, *a, **kw):
        keep = kw.pop("bare_keep", False); kind = kw.pop("_kind", "figure_text")
        return _rec(fig_text(self, *a, **kw), kind, None, keep)

    def f_sup(self, *a, **kw):
        kw.pop("bare_keep", None)
        return _rec(fig_sup(self, *a, **kw), "suptitle", None, False)

    _maxes.Axes.text, _maxes.Axes.annotate = text, annotate
    _maxes.Axes.set_title, _maxes.Axes.legend = set_title, legend
    _maxes.Axes.set_xticklabels, _maxes.Axes.set_yticklabels = set_xticklabels, set_yticklabels
    _mfigure.Figure.text, _mfigure.Figure.suptitle = f_text, f_sup


def _hex_of(h):
    from matplotlib.colors import to_hex
    for g in ("get_facecolor", "get_color", "get_markerfacecolor"):
        try:
            c = getattr(h, g)()
            c = c[0] if hasattr(c, "__len__") and len(c) and hasattr(c[0], "__len__") else c
            return to_hex(c)
        except Exception:
            continue
    return ""


def _marker_of(h):
    try:
        m = h.get_marker()
        return "" if m in (None, "None", "") else str(m)
    except Exception:
        return "patch"


_patch()


def cell_text_colour(cmap, v, vmin, vmax):
    """ヒートマップのセル内文字を、その塗り色の相対輝度で黒/白に切り替える。
    淡い単色系（Blues 等）では大半が黒になり、濃い側だけ白になる。"""
    import numpy as _np
    from matplotlib.colors import Normalize
    if v is None or not _np.isfinite(v):
        return "0.1"
    r, g, b, _ = plt.get_cmap(cmap)(float(Normalize(vmin, vmax)(v)))
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "white" if lum < 0.55 else "0.1"


MIN_PT = 7.0
_RAISED = []


def _enforce_min_pt(fig, floor=MIN_PT):
    """全テキストを下限 floor pt に引き上げる。
    Tick のラベルは描画時に Axis から再設定されるため、Text を直接書き換えるだけでは
    効かない（実測で 5.2-6.5 pt が PDF に残った）。軸は tick_params、凡例は Text、
    カラーバーもその軸として同じ経路で処理する。"""
    import matplotlib.text as mtext
    raised = []
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            labs = axis.get_ticklabels()
            cur = min([l.get_fontsize() for l in labs], default=floor)
            if cur < floor:
                raised.append((f"ticks:{axis.axis_name}", round(cur, 2)))
                ax.tick_params(axis=axis.axis_name, labelsize=floor)
        lg = ax.get_legend()
        if lg is not None:
            for txt in lg.get_texts():
                if txt.get_fontsize() < floor:
                    raised.append((txt.get_text()[:30], round(txt.get_fontsize(), 2)))
                    txt.set_fontsize(floor)
    for o in fig.findobj(mtext.Text):
        s = o.get_text()
        if not s or not str(s).strip():
            continue
        fs = o.get_fontsize()
        if fs < floor:
            raised.append((str(s)[:30], round(fs, 2)))
            o.set_fontsize(floor)
    return raised


def save(f, name):
    _RAISED[:] = _enforce_min_pt(f)
    if _RAISED:
        print(f"  [{name}] 7pt 未満を引き上げ: {len(_RAISED)} 要素")
    if BARE:
        # Pass the tight region measured before removing the text to bbox_inches. In the same region as the version with text
        # Since it is cropped out, the pptx text box coordinates correspond 1:1 with the image,
        # Even the off-axis legend fits within the image (with automatic tight layout, the region would shrink and shift by the amount that was removed).
        bb = _dump_textmeta(f, name)
        f.savefig(f"{OUT}/{name}.pdf", bbox_inches=bb)
        f.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches=bb)
    else:
        f.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
        f.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    plt.close(f)
    print(f"saved {name}.pdf / .png")


def _dump_textmeta(f, name):
    """不可視にした文字の位置・字体を図座標（0-1）で JSON に出す。pptx がこれを読む。"""
    import json
    from matplotlib.transforms import Bbox
    f.canvas.draw()
    r = f.canvas.get_renderer()
    # Coordinate reference is not the canvas, but "the same tight saved area as the version with text."
    # If placed on canvas coordinates, the legend positioned outside (below) the axes extends beyond the figure image (y0 > 1 for frac),
    # The box was clamped to the slide edge on the pptx side, causing overlap (Fig4/Fig7/FigS4).
    PAD = 0.1                                  # Default pad_inches for savefig(bbox_inches="tight")
    _tb = f.get_tightbbox(r)                   # in inches. Not yet hidden, so same area as the version with text
    # Union of canvas and tight bounding boxes (both width and height). With canvas alone, legends placed outside the axes, or
    # long tick labels that do not fit in the left margin (FigS8a) get clipped, causing the corresponding coordinates to go negative,
    # which get clamped on the pptx side, making right-aligned text extend beyond the axis and appear on top of the bars.
    # With tight alone, the width falls short of canvas in figures with little margin. Taking the union preserves at least canvas,
    # while still including the overflowing portion in the image. The layout of the figure itself is not changed.
    _wi, _hi = f.get_size_inches()
    BB = Bbox.from_extents(min(0, _tb.x0 - PAD), min(0, _tb.y0 - PAD),
                           max(_wi, _tb.x1 + PAD), max(_hi, _tb.y1 + PAD))
    _dpi, _W, _H = f.dpi, BB.width, BB.height
    rows = []

    def frac(bb):
        return dict(x0=(bb.x0 / _dpi - BB.x0) / _W, y0=1 - (bb.y1 / _dpi - BB.y0) / _H,
                    x1=(bb.x1 / _dpi - BB.x0) / _W, y1=1 - (bb.y0 / _dpi - BB.y0) / _H)

    kept = []
    for m in _TEXTMETA:
        if m["keep"]:
            a = m["artist"]
            try:
                kept.append(dict(text=str(a.get_text())[:24], **frac(a.get_window_extent(r))))
            except Exception:
                pass
            continue
        a, kind = m["artist"], m["kind"]
        d = dict(kind=kind, panel=_panel_of(m["ax"], f))
        if kind == "legend" and a is not None:
            try:
                d.update(frac(a.get_window_extent(r)))
            except Exception:
                continue
            # Record the actual position of each item. Previously, only the overall outline of the entire legend was passed, and the pptx side
            # performed a repositioning that "stacked them into a single column," so multi-column legends (ncol>1) or legends placed outside the axes
            # flowed into the adjacent panel, causing the text boxes to overlap with each other.
            ents = [dict(e) for e in m["extra"].get("entries", [])]
            try:
                for e, tx in zip(ents, a.get_texts()):
                    e["pos"] = frac(tx.get_window_extent(r))
            except Exception:
                pass
            d["entries"] = ents
            d["fontsize"] = 7.0
            d["text"] = " | ".join(e["label"] for e in d["entries"])
        elif kind == "tick_moved":
            ax, which = m["ax"], m["extra"]["which"]
            labs = ax.get_xticklabels() if which == "x" else ax.get_yticklabels()
            for lab in labs:
                s = str(lab.get_text())
                if not s.strip():
                    continue
                dd = dict(kind="tick_moved", axis=which, panel=_panel_of(ax, f),
                          fontsize=7.0, text=s)
                try:
                    dd.update(frac(lab.get_window_extent(r)))
                except Exception:
                    continue
                rows.append(dd)
                lab.set_visible(False)          # Hide after measuring the outline
            continue
        elif kind.startswith("tick_paren"):
            ax, which = m["ax"], m["extra"]["which"]
            labs = (ax.get_xticklabels() if which == "x" else ax.get_yticklabels())
            for x in m["extra"]["dropped"]:
                if x["i"] >= len(labs):
                    continue
                dd = dict(kind="tick_paren", axis=which, panel=_panel_of(ax, f), fontsize=7.0,
                          text=x["full"][len(x["kept"]):].strip(), full=x["full"], kept=x["kept"])
                try:
                    dd.update(frac(labs[x["i"]].get_window_extent(r)))
                except Exception:
                    continue
                rows.append(dd)
            continue
        else:
            if a is None or not str(getattr(a, "get_text", lambda: "")()).strip():
                continue
            try:
                d.update(frac(a.get_window_extent(r)))
            except Exception:
                continue
            d["text"] = a.get_text()
            d["fontsize"] = float(a.get_fontsize())
            d["bold"] = a.get_fontweight() in ("bold", 700)
            d["colour"] = _hex_of(a)
        rows.append(d)
    os.makedirs(f"{OUT}/textmeta", exist_ok=True)
    # Coordinates have been measured, so only now set to hidden
    n_hidden = 0
    for m in _TEXTMETA:
        if m["keep"] or m["artist"] is None:
            continue
        a = m["artist"]
        try:
            if m["extra"].get("keep_arrow"):
                # For labeled arrows (leader lines in Fig3a/b), remove the leader line along with its label.
                # Unlabeled arrows (flow diagram in Fig1a) do not go into _rec, so they remain.
                a.set_text("")
                ap = getattr(a, "arrow_patch", None)
                if ap is not None:
                    ap.set_visible(False)
            else:
                a.set_visible(False)
                # Even if the Legend is hidden, its child Text remains visible, so remove them individually
                # (If not removed, it will contaminate the "remaining text" tally and overlap check of the plain figure plate)
                if m["kind"] == "legend":
                    for _t in a.get_texts():
                        _t.set_visible(False)
            n_hidden += 1
        except Exception:
            pass
    # Also output rectangles for axes, ticks, and colorbar (text drawn even in BARE): used for the pptx overlap check
    import matplotlib.text as _mt
    drawn = []
    for o in f.findobj(_mt.Text):
        s = str(o.get_text()).strip()
        if not s or not o.get_visible():
            continue
        try:
            bb = o.get_window_extent(r)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        drawn.append(dict(text=s[:24], **frac(bb)))
    json.dump(dict(figure=name, width_mm=_W * 25.4, height_mm=_H * 25.4, items=rows,
                   kept=kept, drawn=drawn),
              open(f"{OUT}/textmeta/{name}.json", "w"), indent=1)
    _TEXTMETA.clear()
    print(f"  [{name}] 非表示 {n_hidden} 件 / 記録 {len(rows)} 件 / 素図版に残る文字 {len(drawn)} 件")
    return BB


def _panel_of(ax, f):
    if ax is None:
        return ""
    try:
        return str(f.axes.index(ax))
    except Exception:
        return ""


def paneln(rows, name):
    import pandas as pd
    pd.DataFrame(rows).to_csv(f"{TAB}/_msn_{name}.csv", index=False)
