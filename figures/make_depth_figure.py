"""Why chain depth and reaction count differ, on one real route.

    uv run --extra figures python figures/make_depth_figure.py

The formose route to deoxyribose is convergent: two branches meet at the last
step. Its longest chain of consecutive reactions is three, while the whole
derivation uses four. Both are exact minima of their own objective, which is
why the SI tables report them separately.
"""

from pathlib import Path

from acs_style import DOUBLE, figure, save
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "figures" / "output"

CHAIN, OFF, RXN = "#1f3d63", "#9bb0c9", "#b8860b"
HW, HH = 0.31, 0.10          # node half width and half height, data units

#: (x, y, label, on the longest chain). From figures/routes/Formose_Deoxyribose.yaml.
NODES = {
    "w":    (0.0,  1.25, "water", True),
    "gly1": (0.0,  0.98, "glycolaldehyde", True),
    "f1":   (0.0,  0.71, "formaldehyde", True),
    "eg":   (1.0,  0.98, "ethylene glycol", True),
    "ac":   (2.0,  0.98, "acetaldehyde", True),
    "gly2": (0.0,  0.25, "glycolaldehyde", False),
    "f2":   (0.0, -0.02, "formaldehyde", False),
    "ga":   (1.0,  0.115, "glyceraldehyde", False),
    "dr":   (3.0,  0.55, "2-deoxyribose", True),
}
#: (x, y, id, reagents, product, on the longest chain)
STEPS = [
    (0.5, 0.98,  "r1", ["w", "gly1", "f1"], "eg", True),
    (1.5, 0.98,  "r2", ["eg"], "ac", True),
    (0.5, 0.115, "r3", ["gly2", "f2"], "ga", False),
    (2.5, 0.55,  "r4", ["ac", "ga"], "dr", True),
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    xlo, xhi, ylo, yhi = -0.45, 3.55, -0.45, 1.62
    # match the figure box to the data box, so aspect "equal" costs no space
    fig, ax = figure(DOUBLE, DOUBLE * (yhi - ylo) / (xhi - xlo))

    for key, (x, y, label, on) in NODES.items():
        colour = CHAIN if on else OFF
        ax.add_patch(FancyBboxPatch((x - HW, y - HH), 2 * HW, 2 * HH,
                                    boxstyle="round,pad=0,rounding_size=0.035",
                                    linewidth=1.1 if key == "dr" else 0.7,
                                    edgecolor=colour, facecolor="white", zorder=3))
        ax.text(x, y, label, ha="center", va="center", fontsize=6.4, color="#111",
                zorder=4, fontweight="bold" if key == "dr" else "normal")

    for sx, sy, rid, reagents, product, on in STEPS:
        colour = CHAIN if on else OFF
        ax.add_patch(FancyBboxPatch((sx - 0.03, sy - 0.03), 0.06, 0.06,
                                    boxstyle="square,pad=0", linewidth=0,
                                    facecolor=RXN if on else OFF, zorder=3))
        ax.text(sx, sy + 0.055, rid, ha="center", va="bottom", fontsize=6, color=colour)
        for r in reagents:
            rx, ry, *_ = NODES[r]
            ax.add_patch(FancyArrowPatch((rx + HW + 0.01, ry), (sx - 0.035, sy),
                                         arrowstyle="-", linewidth=0.7, color=colour,
                                         shrinkA=0, shrinkB=0, zorder=2))
        px, py, *_ = NODES[product]
        ax.add_patch(FancyArrowPatch((sx + 0.035, sy), (px - HW - 0.02, py),
                                     arrowstyle="-|>", mutation_scale=5, linewidth=0.9,
                                     color=colour, shrinkA=0, shrinkB=0, zorder=2))

    ax.plot([0.36, 2.64], [1.44, 1.44], linewidth=0.8, color=CHAIN,
            solid_capstyle="butt", zorder=2)
    for x in (0.36, 2.64):
        ax.plot([x, x], [1.41, 1.47], linewidth=0.8, color=CHAIN, zorder=2)
    ax.text(1.5, 1.50, "longest chain of consecutive reactions: 3", ha="center",
            va="bottom", fontsize=6.6, color=CHAIN)
    ax.text(1.5, -0.34, "whole derivation: 4 reactions, counting r1–r4 once each",
            ha="center", va="center", fontsize=6.6, color="#111")

    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_aspect("equal")
    ax.axis("off")
    save(fig, str(OUT / "Figure_S4_depth_vs_count"))
    print("wrote figures/output/Figure_S4_depth_vs_count.png")


if __name__ == "__main__":
    main()
