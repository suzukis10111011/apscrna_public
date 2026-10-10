#!/usr/bin/env python
"""Purpose: Create editable PowerPoint from base PNG with text elements as movable textboxes
Input: none (generates from BARE=1 rendering metadata in JSON)
Output: {BARE}/colour_keys.csv, {W}/results/tables/36_pptx_layout_checks.csv

36: Using the bare figure PNG as a base, create a pptx with the removed text placed as editable text boxes.

Coordinates are not estimated. At BARE=1 rendering time, `_msstyle._dump_textmeta()` outputs
the window_extent of the hidden Text / Legend objects in figure coordinates (0-1, origin at top-left) to JSON,
so we simply linearly transform that into slide coordinates. Since the bare figure is saved without bbox_inches,
figure coordinates, PNG pixels, and slide dimensions correspond 1:1.
"""
import json, glob, os, pathlib
import pandas as pd
from pptx import Presentation
from pptx.util import Mm, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.enum.shapes import MSO_SHAPE

import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
BARE, OUTP = f"{W}/results/figures_ms/bare", f"{W}/results/figures_ms"
ORDER = ["Fig1_design_composition", "Fig2_enrichment", "Fig3_tf_screen",
         "Fig4_coupling_coexistence", "Fig5_receptors_coactivity", "Fig6_compartments",
         "Fig7_independent_cohort", "FigS1_qc", "FigS2_integration", "FigS3_genesets_scoring",
         "FigS4_matrix_vs_requantification", "FigS5_raw_vs_depth", "FigS6_leave_one_out",
         "FigS7_state_labelling", "FigS8_other_programs"]
KIND_JP = {"panel_letter": "パネル記号", "panel_title": "タイトル", "legend": "凡例",
           "annotation": "注記", "tick_paren": "群サイズ", "figure_text": "図中の文字",
           "suptitle": "図タイトル"}


MATH = [("$\\rho$", "\u03c1"), ("$\\times$", "\u00d7"), ("$\\Delta$", "\u0394"),
        ("$-$", "\u2212"), ("$\\pm$", "\u00b1"), ("$\\leq$", "\u2264"), ("$\\geq$", "\u2265"),
        ("log$_{10}$", "log\u2081\u2080"), ("$^{-4}$", "\u207b\u2074"), ("$_{10}$", "\u2081\u2080")]
TEXT_COLOUR = "#000000"        # All pptx text boxes are black ([1])


def demath(s):
    """図では mathtext で組んでいる記号を、pptx のテキストボックス用に unicode に置く。
    残った $ は素の文字として見えてしまうので取り除き、呼び出し側で残存を検査する。"""
    s = str(s)
    for a, b in MATH:
        s = s.replace(a, b)
    if "$" in s:
        import re as _re
        s = _re.sub(r"\$\\?([A-Za-z]+)\$", r"\1", s).replace("$", "")
    return s


def rgb(h):
    h = (h or "#333333").lstrip("#")
    return RGBColor(int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))


def box(slide, x_mm, y_mm, w_mm, h_mm, text, pt=7.0, bold=False, colour="#1a1a1a",
        align=PP_ALIGN.LEFT, slide_mm=None, wrap=False):
    # Positions derived from figure coordinates may fall outside the slide (since panel labels are placed outside the axes,
    # in plain figures not using bbox_inches they can fall outside the canvas). Clamp them to the inside.
    w_mm, h_mm = max(w_mm, 3.0), max(h_mm, 2.6)
    if slide_mm:
        x_mm = min(max(x_mm, 0.4), slide_mm[0] - w_mm - 0.4)
        y_mm = min(max(y_mm, 0.4), slide_mm[1] - h_mm - 0.4)
    tb = slide.shapes.add_textbox(Mm(x_mm), Mm(y_mm), Mm(w_mm), Mm(h_mm))
    tf = tb.text_frame
    tf.word_wrap = wrap
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = MSO_ANCHOR.TOP
    lines = demath(text).split("\n")
    for i, ln in enumerate(lines):
        par = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        par.alignment = align
        r = par.add_run(); r.text = ln
        r.font.size = Pt(max(pt, 7.0)); r.font.bold = bold
        r.font.name = "Arial"; r.font.color.rgb = rgb(TEXT_COLOUR)   # Colors are unified to black
    return tb


def build(names, out_path, per_figure_size):
    """per_figure_size=True: 図ごとに1スライドだけの pptx を、その図の実寸で作る。
    False: 全図を 1 ファイルにまとめ、スライド寸法は最大の図に合わせる（PowerPoint は
    ファイル単位でしかスライド寸法を持てない）。"""
    prs = Presentation()
    metas = {n: json.load(open(f"{BARE}/textmeta/{n}.json")) for n in names}
    # PowerPoint can only hold slide dimensions on a per-file basis (self-determined).
    # The width is 180 mm for all figures, and the height is 78-330 mm per figure, so the slide is unified to 180 x 330 mm,
    # and each figure is placed at actual size (1:1) in the top-left. Below shorter figures there is margin.
    # Dimensions are determined from the PNG pixels (rounding the mm in textmeta causes a 1 mm offset,
    # resulted in a mismatch where the image does not fit on the slide / leaves a remainder)
    from PIL import Image as _Im

    def _png_mm(n):
        im = _Im.open(f"{BARE}/{n}.png")
        return im.width / 600 * 25.4, im.height / 600 * 25.4

    _sz = {n: _png_mm(n) for n in names}
    SW = max(v[0] for v in _sz.values())
    SH = max(v[1] for v in _sz.values())
    prs.slide_width, prs.slide_height = Mm(SW), Mm(SH)
    keys, placed, checks = [], 0, []
    put = []          # The rectangle of the placed text box. Used for the judgment to avoid the scale bracket
    for name in names:
        mp = f"{BARE}/textmeta/{name}.json"
        png = f"{BARE}/{name}.png"
        assert os.path.exists(mp) and os.path.exists(png), name
        M = metas[name]
        Wm, Hm = _sz[name]          # Dimensions derived from PNG pixels. Keep image, slide, and coordinates aligned to the same value
        s = prs.slides.add_slide(prs.slide_layouts[6])
        s.shapes.add_picture(png, 0, 0, Mm(Wm), Mm(Hm))       # Raw figure at 1:1, full-page
        # The clamp upper limit of the box is "the dimensions of that figure's image". Since the combined version has the slide
        # matched to the largest figure, cutting with SH causes the legend of a shorter figure to escape outside the image.
        CLAMP = (Wm, Hm)
        panels, n_box = [], 0
        # Place in 2 passes: the scale bracket is placed "directly below the group name", but in a figure where the legend is placed there (Fig4a)
        # it overlaps with the legend. Place other boxes such as the legend first, and have the bracket check for collisions and avoid them last.
        for it in sorted(M["items"], key=lambda q: q.get("kind") == "tick_paren"):
            x, y = it["x0"] * Wm, it["y0"] * Hm
            w, h = (it["x1"] - it["x0"]) * Wm, (it["y1"] - it["y0"]) * Hm
            k = it["kind"]
            if k == "legend":
                ent = it.get("entries", [])
                step = 3.4 if not ent else min(3.4, max(2.7, h / max(len(ent), 1)))
                for j, e in enumerate(ent):
                    # If an actual position exists, place it there (so multi-column legends/off-axis legends don't shift).
                    # Only when absent, fall back to the conventional arrangement of stacking into a single column based on the overall outline of the entire legend.
                    pos = e.get("pos")
                    if pos:
                        lx, ly = pos["x0"] * SW, pos["y0"] * SH
                        lh = max((pos["y1"] - pos["y0"]) * SH, 2.6)
                        lw = max((pos["x1"] - pos["x0"]) * SW + 1.0, 8.0)
                    else:
                        lx, ly = x + 3.0, y + j * step
                        lh, lw = max(step - 0.2, 2.6), max(len(e["label"]) * 1.5, 16)
                    sx = min(max(lx - 3.0, 0.4), SW - 2.6)
                    sy = min(max(ly + 0.4, 0.4), SH - 2.6)
                    sq = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Mm(sx), Mm(sy),
                                            Mm(2.2), Mm(2.2))
                    sq.fill.solid(); sq.fill.fore_color.rgb = rgb(e.get("colour"))
                    sq.line.fill.background(); sq.shadow.inherit = False
                    box(s, lx, ly, lw, lh, e["label"], 7.0, slide_mm=CLAMP)
                    put.append((lx, ly, lx + lw, ly + lh))
                    n_box += 1
                    keys.append(dict(figure=name.split("_")[0], panel=it.get("panel", ""),
                                     colour_hex=e.get("colour", ""), marker=e.get("marker", ""),
                                     label=e["label"]))
            elif k == "tick_paren":
                # Since the group names remain as tick marks, shift the parentheses directly below for the x axis, or to the left for the y axis
                pt = 7.0
                if it.get("axis") == "x":
                    y = y + h + 0.3
                else:
                    x = max(x - (len(it["text"]) * 1.5 + 0.6), 0.4)
                bw, bh = len(it["text"]) * 1.6 + 1.0, h + 0.4
                # The placement of the parentheses is "directly below the group name", but in figures (Fig4a) where the legend is placed there,
                # it overlaps with the legend box. Only when it collides with an already-placed box, escape downward.
                for _ in range(6):
                    hit = [q for q in put if not (x + bw <= q[0] or q[2] <= x
                                                  or y + bh <= q[1] or q[3] <= y)]
                    if not hit or y + bh + 0.4 > SH:
                        break
                    y = max(h[3] for h in [(0, 0, 0, max(q[3] for q in hit))]) + 0.3
                box(s, x, y, bw, bh, it["text"], pt, False, "#1a1a1a", slide_mm=CLAMP)
                put.append((x, y, x + bw, y + bh))
                n_box += 1
            elif k in ("number", "tick_moved"):
                # The recorded outline is based on the original font size (5.2-6 pt). Since the box text is 7 pt,
                # it does not fit as is. Estimate the width needed at 7 pt, and for the numbers fix the center,
                # for the y-axis tick labels fix the right edge (axis side) and expand. The format is shared.
                lines = str(it["text"]).split("\n")
                need = max(len(ln) for ln in lines) * 1.55 + 0.6
                bw = max(w, need, 3.0)
                bh = max(h, 2.6 * len(lines))
                wrap = False
                if k == "number":
                    al = PP_ALIGN.CENTER
                    bx = x + w / 2 - bw / 2
                elif it.get("axis") == "y":
                    al = PP_ALIGN.RIGHT
                    avail = (x + w) - 0.4
                    if bw > avail:
                        # Labels too long to fit within the original margin (for 5.2-5.8 pt) at 7 pt.
                        # If it overflows to the right it sits on top of the bar, so wrap it within the margin width.
                        bw, bx, wrap = avail, 0.4, True
                        bh = max(bh, 2.6 * (1 + int(need // max(avail, 1))))
                    else:
                        bx = x + w - bw
                else:
                    al = PP_ALIGN.CENTER
                    bx = x + w / 2 - bw / 2
                by = y + h / 2 - bh / 2
                box(s, bx, by, bw, bh, it["text"], 7.0, False, "#1a1a1a",
                    align=al, slide_mm=CLAMP, wrap=wrap)
                put.append((bx, by, bx + bw, by + bh))
                n_box += 1
            else:
                pt = it.get("fontsize", 7.0)
                bold = bool(it.get("bold")) or k == "panel_letter"
                if k == "panel_letter":
                    pt = 8.0
                    panels.append(it["text"])
                pad = 0.0 if k == "panel_letter" else 2.0
                bw2, bh2 = w + pad, h + (0.0 if k == "panel_letter" else 0.6)
                box(s, x, y, bw2, bh2,
                    it["text"], pt, bold, it.get("colour") or "#1a1a1a", slide_mm=CLAMP)
                put.append((x, y, x + bw2, y + bh2))
                n_box += 1
            placed += 1
        # Notes field: colour correspondence and panel list
        kk = [x for x in keys if x["figure"] == name.split("_")[0]]
        note = [f"{name}  figure {Wm:.0f} x {Hm:.0f} mm / slide {SW:.0f} x {SH:.0f} mm",
                "パネル: " + (", ".join(panels) if panels else "（記号なし）"),
                f"テキストボックス: {n_box} 個", "", "colour keys:"]
        note += [f"  {x['colour_hex']}  {x['marker'] or 'patch'}  {x['label']}" for x in kk] or ["  （凡例なし）"]
        s.notes_slide.notes_text_frame.text = "\n".join(note)
        # Validation: overlap with outside-slide area / figure text
        drawn = M.get("drawn", [])
        EMM = Mm(1).emu
        for sh in s.shapes:
            if not sh.has_text_frame and sh.shape_type != MSO_SHAPE.RECTANGLE:
                continue                                   # Background image is excluded
            x0, y0 = sh.left / EMM, sh.top / EMM            # mm
            x1, y1 = (sh.left + sh.width) / EMM, (sh.top + sh.height) / EMM
            out = x0 < -0.2 or y0 < -0.2 or x1 > SW + 0.2 or y1 > SH + 0.2
            TOL = 0.3            # mm: encroachment below this is considered adjacent
            ov, ovmm = [], 0.0
            for d in drawn:
                dx = min(x1, d["x1"] * Wm) - max(x0, d["x0"] * Wm)
                dy = min(y1, d["y1"] * Hm) - max(y0, d["y0"] * Hm)
                if dx > TOL and dy > TOL:
                    ov.append(d["text"]); ovmm = max(ovmm, min(dx, dy))
            if out or ov:
                checks.append(dict(figure=name,
                                   shape=(sh.text_frame.text[:28] if sh.has_text_frame else "swatch"),
                                   x_mm=round(x0, 1), y_mm=round(y0, 1),
                                   outside=bool(out), overlaps=len(ov),
                                   overlap_min_mm=round(ovmm, 2),
                                   overlap_with="; ".join(ov[:3])))
        print(f"  {name}: テキストボックス {n_box} 個, 凡例キー {len(kk)} 件")
    prs.save(out_path)
    return keys, checks, placed, prs, SW, SH


def main():
    keys, checks, placed, prs, SW, SH = build(ORDER, f"{OUTP}/figures_editable.pptx", False)
    os.makedirs(f"{OUTP}/pptx", exist_ok=True)
    for n in ORDER:
        _, _, _, _, w, h = build([n], f"{OUTP}/pptx/{n.split('_')[0]}.pptx", True)
        print(f"  pptx/{n.split('_')[0]}.pptx: スライド {w:.0f} x {h:.0f} mm")
    # Append abbreviated symbol (FigS3) to colour_keys
    sc = f"{W}/results/tables/_figS3_setcodes.csv"
    if os.path.exists(sc):
        for r in pd.read_csv(sc).itertuples():
            keys.append(dict(figure="FigS3", panel="a,b", colour_hex="", marker="code",
                             label=f"{r.code} = {r.gene_set} ({r.side})"))
    pd.DataFrame(keys).to_csv(f"{BARE}/colour_keys.csv", index=False)
    pd.DataFrame(checks).to_csv(f'{os.environ.get("APSC_TABDIR", f"{W}/results/tables")}/36_pptx_layout_checks.csv', index=False)
    print(f"\nfigures_editable.pptx: {len(prs.slides._sldIdLst)} スライド "
          f"({SW:.0f} x {SH:.0f} mm) / 図形計 {placed} 個")
    print(f"colour_keys.csv: {len(keys)} 行")
    print(f"検査で引っかかった図形: {len(checks)} 件"
          + (f" -> results/tables/36_pptx_layout_checks.csv" if checks else ""))


main()
