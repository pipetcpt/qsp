"""Shared matplotlib style for the 2026-tox talk figures.

Palette and typography are chosen to match (a) the uploaded 2026-tox.pptx
template (blue title master, white content master, QR + PIPET logos) and
(b) the author's docs/ reveal decks (navy / teal / gold, Inter typeface).

All figures render on a transparent/white canvas at 300 dpi so they drop
cleanly onto the white content slides.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from pathlib import Path

# ---- register Inter -------------------------------------------------------
for _f in fm.findSystemFonts():
    if "/inter/" in _f.lower() and "display" not in _f.lower():
        try:
            fm.fontManager.addfont(_f)
        except Exception:
            pass
_FAMILY = "Inter" if any("Inter" == f.name for f in fm.fontManager.ttflist) else "DejaVu Sans"

# ---- palette (content-informed: benefit vs risk) --------------------------
NAVY   = "#163A5F"   # primary dark / titles / structure
BLUE   = "#2B6CB0"   # secondary
TEAL   = "#1F8A8C"   # accent / benefit, efficacy
GOLD   = "#C8922A"   # sharp accent / caution
RED    = "#C0392B"   # harm / risk / toxicity
GREEN  = "#2E8B6F"   # benefit / safe
PURPLE = "#6B4E8C"   # fourth series
GREY   = "#8A98A6"   # muted / reference
LGREY  = "#D7DEE5"   # gridlines / fills
INK    = "#1A222C"   # body text
BG     = "#FFFFFF"

# ordered categorical cycle (max visual separation, colour-blind aware)
CYCLE = [NAVY, TEAL, GOLD, RED, BLUE, GREEN, PURPLE, GREY]

# sequential (navy ramp) for heatmaps
def navy_ramp():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("navy", ["#EAF2F8", "#9FC2DE", BLUE, NAVY])

def risk_ramp():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("risk", ["#EAF4F0", "#BFE0D2", GOLD, RED])


def apply():
    plt.rcParams.update({
        "font.family": _FAMILY,
        "font.size": 13,
        "axes.titlesize": 15,
        "axes.titleweight": "bold",
        "axes.titlecolor": NAVY,
        "axes.labelsize": 13,
        "axes.labelcolor": INK,
        "axes.edgecolor": "#AEB9C4",
        "axes.linewidth": 1.0,
        "axes.facecolor": "none",
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": LGREY,
        "grid.linewidth": 0.8,
        "grid.alpha": 0.9,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": 11.5,
        "ytick.labelsize": 11.5,
        "legend.fontsize": 11.5,
        "legend.frameon": False,
        "figure.facecolor": "none",
        "savefig.facecolor": "none",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.06,
        "lines.linewidth": 2.6,
        "lines.solid_capstyle": "round",
    })
    for s in ("top", "right"):
        plt.rcParams[f"axes.spines.{s}"] = False


def despine(ax, left=True, bottom=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)


FIGDIR = Path(__file__).resolve().parent / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    p = FIGDIR / name
    fig.savefig(p, dpi=300, transparent=True)
    plt.close(fig)
    print("wrote", p)
    return p
