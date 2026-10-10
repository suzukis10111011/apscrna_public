"""Purpose: Enforce shared figure style (size, font, color, export formats)
Input: none detected
Output: {OUT}/{name}.pdf, {OUT}/{name}.png

Common figure specification: width 180mm (2-column)/85mm (1-column), Arial 7-8pt, colorblind-friendly, PDF + 600dpi PNG.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os as _os
W = _os.environ.get("APSC_ROOT") or _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
# The output destination is switched via APSC_FIGDIR (default figures_draft, final version is figures_final)
OUT = f"{W}/results/" + os.environ.get("APSC_FIGDIR", "figures_draft")
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
os.makedirs(OUT, exist_ok=True)
MM = 1 / 25.4
DISC, VAL = "#2166ac", "#e08214"          # Discovery = blue tones / Validation = orange tones
STATE = {"inflammatory": "#440154", "intermediate": "#31688e",
         "tissue-forming": "#35b779", "other": "#bdbdbd"}
OKABE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#999999"]

def setup():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "savefig.bbox": "tight", "pdf.fonttype": 42})

def fig(width_mm=180, height_mm=100):
    setup()
    return plt.subplots(figsize=(width_mm * MM, height_mm * MM))

def gridfig(width_mm=180, height_mm=120, **kw):
    setup()
    f = plt.figure(figsize=(width_mm * MM, height_mm * MM))
    return f, f.add_gridspec(**kw)

def panel(ax, letter, dx=-0.13, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")

def save(f, name):
    os.makedirs(OUT, exist_ok=True)
    f.savefig(f"{OUT}/{name}.pdf")
    f.savefig(f"{OUT}/{name}.png", dpi=600)
    print(f"saved {name}.pdf / .png")
