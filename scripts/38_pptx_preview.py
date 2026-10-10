#!/usr/bin/env python
"""Purpose: Reconstruct saved PPTX slides as matplotlib plots to audit colors, dollar signs, and overflow
Input: none (reads saved PPTX files)
Output: results/figures_ms/pptx/preview/*.png, results/tables/38_pptx_text_audit.csv

38: Rasterize each slide for verification.

Since this environment has no LibreOffice / PowerPoint, **read back the saved pptx**
and redraw the background images and text boxes in matplotlib at the same coordinates, same text, and same colors,
outputting results/figures_ms/pptx/preview/*.png. What is checked is
(1) non-black text such as blue text, (2) residual "$", (3) overflow beyond the slide.
"""
import os, glob, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
from pptx import Presentation
from pptx.util import Mm

import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
PD = f"{W}/results/figures_ms/pptx"
OUT = f"{PD}/preview"
os.makedirs(OUT, exist_ok=True)
EMM = Mm(1).emu
rows = []

for fp in sorted(glob.glob(f"{PD}/*.pptx")):
    name = pathlib.Path(fp).stem
    prs = Presentation(fp)
    SW, SH = prs.slide_width / EMM, prs.slide_height / EMM
    s = list(prs.slides)[0]
    fig = plt.figure(figsize=(SW / 25.4, SH / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, SW); ax.set_ylim(SH, 0); ax.axis("off")
    for sh in s.shapes:
        x, y = sh.left / EMM, sh.top / EMM
        w, h = sh.width / EMM, sh.height / EMM
        if sh.shape_type == 13:                      # Background image
            img = mpimg.imread(sh.image.blob and __import__("io").BytesIO(sh.image.blob))
            ax.imshow(img, extent=[x, x + w, y + h, y], zorder=0, aspect="auto")
            continue
        if not sh.has_text_frame or not sh.text_frame.text.strip():
            ax.add_patch(plt.Rectangle((x, y), w, h, facecolor="#888888", lw=0, zorder=3))
            continue
        # PowerPoint stacks paragraphs vertically. Previously all paragraphs were drawn at the same y, so
        # the two-line tick labels (FigS4c) appeared overlapping. Add line spacing.
        line = 0
        for par in sh.text_frame.paragraphs:
            for r in par.runs:
                hexc = "#000000"
                try:
                    if r.font.color and r.font.color.type is not None:
                        hexc = "#" + str(r.font.color.rgb)
                except Exception:
                    pass
                pt = r.font.size.pt if r.font.size else 7.0
                # Reflect paragraph alignment (previously always drawn from the left edge, so right-aligned
                # tick labels appeared to overflow to the right of the axis)
                al = str(par.alignment)
                if "CENTER" in al:
                    tx, ha = x + w / 2, "center"
                elif "RIGHT" in al:
                    tx, ha = x + w, "right"
                else:
                    tx, ha = x, "left"
                txt = r.text
                if sh.text_frame.word_wrap:
                    # For boxes with wrap specified, PowerPoint wraps at the box width. Show it wrapped in the
                    # preview as well (without wrapping it looks like it doesn't fit). From the 7 pt measurement,
                    # approximately 1.36 mm per character.
                    import textwrap
                    txt = "\n".join(textwrap.wrap(txt, max(int(w / 1.36), 8)))
                nlines = txt.count(chr(10)) + 1
                ax.text(tx, y + 1.0 + line * pt * 0.36, txt, fontsize=pt, color=hexc,
                        va="top", ha=ha,
                        fontweight="bold" if r.font.bold else "normal", zorder=4,
                        linespacing=1.15)
                line += nlines
                rows.append(dict(figure=name, text=r.text[:40], colour=hexc.lower(), pt=round(pt, 1),
                                 has_dollar="$" in r.text,
                                 outside=bool(x < -0.2 or y < -0.2 or x + w > SW + 0.2 or y + h > SH + 0.2)))
                break                                 # Paragraph 1, line 1
    fig.savefig(f"{OUT}/{name}.png", dpi=200)
    plt.close(fig)

R = pd.DataFrame(rows)
R.to_csv(f'{os.environ.get("APSC_TABDIR", f"{W}/results/tables")}/38_pptx_text_audit.csv', index=False)
nonblack = R[R.colour != "#000000"]
print(f"テキスト run: {len(R)} 件 / {R.figure.nunique()} 図")
print(f"  黒以外の文字: {len(nonblack)} 件", nonblack.colour.unique().tolist() if len(nonblack) else "")
print(f"  '$' を含む: {int(R.has_dollar.sum())} 件", R[R.has_dollar].text.tolist()[:4])
print(f"  スライド外: {int(R.outside.sum())} 件")
print(f"  7 pt 未満: {int((R.pt < 7).sum())} 件")
print(f"preview: {len(glob.glob(f'{OUT}/*.png'))} 枚")
