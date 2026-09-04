"""Figure S1: the pipeline, as implemented in src/nucleoside_analogues.

    uv run --extra figures python figures/make_workflow_figure.py
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

matplotlib.use("Agg")

OUT = Path(__file__).parent / "workflow"
INK = "#1a1a1a"
FILL = {"in": "#eef1f5", "mid": "#ffffff", "out": "#e4ebe3"}

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
})


def box(ax, x, y, w, h, text, kind="mid", bold_first=True):
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        linewidth=0.8, edgecolor=INK, facecolor=FILL[kind], zorder=2))
    head, _, rest = text.partition("\n")
    ax.text(x, y + (0.018 if rest else 0), head, ha="center", va="center",
            fontsize=8, fontweight="bold" if bold_first else "normal",
            color=INK, zorder=3)
    if rest:
        ax.text(x, y - 0.026, rest, ha="center", va="center",
                fontsize=6.8, color="#4a4a4a", zorder=3, linespacing=1.35)
    return x, y - h / 2, y + h / 2


def arrow(ax, start, end, label=None):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=8,
        linewidth=0.8, color=INK, shrinkA=0, shrinkB=0, zorder=1))
    if label:
        ax.text((start[0] + end[0]) / 2 + 0.012, (start[1] + end[1]) / 2, label,
                ha="left", va="center", fontsize=6.3, color="#4a4a4a", zorder=3)


fig, ax = plt.subplots(figsize=(7.0, 4.3))
ax.set_xlim(0, 1); ax.set_ylim(-0.015, 1.015); ax.axis("off")

W, H, C = 0.40, 0.10, 0.70
xl, xr = 0.245, 0.755

a1 = box(ax, xl, 0.945, W, H, "Nucleoside analogue library\nMOLGEN enumeration, CHO and CHNO", "in")
a2 = box(ax, xr, 0.945, W, H, "Chemical reaction network\nMØD graph grammar, seeds expanded G0–G6", "in")
b1 = box(ax, xl, 0.775, W, H, "Substitute and flatten\nCl → OH / NH$_2$; stereoisomers collapsed")
b2 = box(ax, xr, 0.775, W, H, "Products and reaction relations\nspecies per generation; reagents → products")
c = box(ax, 0.5, 0.585, C, H, "Match on InChIKey first block\nconstitution and charge; stereochemistry not encoded")
d = box(ax, 0.5, 0.405, C, H, "Reaction free energies\neQuilibrator component contribution: ΔrG′° ± σ")
e = box(ax, 0.5, 0.225, C, H, "Three-way spontaneity call\nspontaneous / non-spontaneous / undetermined")
f = box(ax, 0.5, 0.060, C, H, "Minimum-weight hyperpath (Knuth)\nsteps to each target, and minimum-cost route count", "out")

arrow(ax, (xl, a1[1]), (xl, b1[2]))
arrow(ax, (xr, a2[1]), (xr, b2[2]))
for x in (xl, xr):
    ax.plot([x, x], [b1[1], 0.680], lw=0.8, color=INK, zorder=1)
ax.plot([xl, xr], [0.680, 0.680], lw=0.8, color=INK, zorder=1)
arrow(ax, (0.5, 0.680), (0.5, c[2]))
arrow(ax, (0.5, c[1]), (0.5, d[2]), "matched analogues per network")
arrow(ax, (0.5, d[1]), (0.5, e[2]), "uncertainty retained")
arrow(ax, (0.5, e[1]), (0.5, f[2]), "spontaneous reactions only")

fig.tight_layout(pad=0.3)
OUT.mkdir(exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure_S1_workflow.{ext}", dpi=600, bbox_inches="tight")
print("wrote", *(OUT / f"Figure_S1_workflow.{e}" for e in ("png", "pdf")))
