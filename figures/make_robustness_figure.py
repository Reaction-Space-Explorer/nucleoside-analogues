"""Figure S5: which routes survive removing a disputed mechanism.

    uv run --extra figures python figures/make_robustness_figure.py

One cell per target and CRNR, giving the shortest spontaneous chain with the
whole rule set and the same figure once a family of rules is removed. Colour
carries the outcome, so what survives can be read at a glance.

Numbers from ProcessedData/SI/rule_dependence.csv.
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from acs_style import DOUBLE, panel, save, use
from matplotlib.patches import Patch, Rectangle

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "figures" / "output"

NETS = ["Formose", "FormoseAmm", "Glucose", "GlucoseAmm", "PyruvicAcid"]
HEAD = ["F", "FA", "G", "GA", "PA"]
TARGETS = ["Threose", "Ribose", "Glycerol", "Deoxyribose", "Glyoxylate"]
PANELS = [
    ("carbonyl migration", "carbonyl migration removed"),
    ("Cannizzaro", "Cannizzaro reaction removed"),
]
COLOUR = {"unchanged": "#cfe0d2", "longer": "#f6e3bf", "lost": "#f2cdcd", "none": "#f4f4f4"}
EDGE = "#ffffff"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = {}
    with (REPO / "ProcessedData" / "SI" / "rule_dependence.csv").open() as handle:
        for r in csv.DictReader(handle):
            data[(r["removed"], r["network"], r["target"])] = (
                r["depth_all_rules"],
                r["depth_after_removal"],
            )

    use()
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.5))
    for ax, (family, title) in zip(axes, PANELS, strict=True):
        for col, net in enumerate(NETS):
            for row, target in enumerate(TARGETS):
                before, after = data[(family, net, target)]
                if not before:
                    kind, text = "none", "—"
                elif not after:
                    kind, text = "lost", f"{before} → lost"
                elif int(after) > int(before):
                    kind, text = "longer", f"{before} → {after}"
                else:
                    kind, text = "unchanged", f"{before} → {after}"
                ax.add_patch(
                    Rectangle(
                        (col, row), 0.94, 0.9, facecolor=COLOUR[kind], edgecolor=EDGE, linewidth=1.2
                    )
                )
                ax.text(
                    col + 0.47,
                    row + 0.45,
                    text,
                    ha="center",
                    va="center",
                    fontsize=6.2,
                    color="#222222",
                )
        ax.set_xlim(0, len(NETS))
        ax.set_ylim(0, len(TARGETS))
        ax.set_xticks([i + 0.47 for i in range(len(NETS))])
        ax.set_xticklabels(HEAD, fontsize=7)
        ax.set_yticks([i + 0.45 for i in range(len(TARGETS))])
        ax.set_yticklabels(TARGETS, fontsize=7)
        ax.invert_yaxis()
        ax.set_title(title, fontsize=7.5, pad=4)
        for side in ("top", "right", "bottom", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(length=0)
    for ax, tag in zip(axes, "ab", strict=True):
        panel(ax, tag, x=-0.22, y=1.04)
    fig.legend(
        handles=[
            Patch(facecolor=COLOUR[k], label=v)
            for k, v in (
                ("unchanged", "reached in the same number of steps"),
                ("longer", "still reached, but by a longer route"),
                ("lost", "no longer reachable"),
                ("none", "not reachable with the whole rule set"),
            )
        ],
        loc="lower center",
        ncol=4,
        fontsize=6.5,
        frameon=False,
        handlelength=1.3,
        columnspacing=1.3,
        bbox_to_anchor=(0.5, -0.04),
    )
    fig.tight_layout(w_pad=3.0, rect=(0, 0.10, 1, 1))
    save(fig, str(OUT / "Figure_S5_robustness"))
    print("wrote figures/output/Figure_S5_robustness.png")


if __name__ == "__main__":
    main()
