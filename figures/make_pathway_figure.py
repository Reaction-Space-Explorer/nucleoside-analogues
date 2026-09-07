"""Figure 7: how many spontaneous routes there are, and how favourable they are.

    uv run --extra figures python figures/make_pathway_figure.py

Two panels against route length: how many matched analogues are reached at
that length, and the free energy of the route divided by the reactions it
uses. The route total is not plotted: it falls close to linearly with length
in every CRNR, which is what summing more favourable steps must do, and the
per-reaction figure is what carries information.

Numbers from ProcessedData/SI/figure_pathway_energetics.csv, in kJ/mol.
"""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from acs_style import DOUBLE, save, use

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "figures" / "output"

STYLE = {
    "Formose": ("Formose (F)", "#1f3d63", "o", "-"),
    "FormoseAmm": ("Formose ammonia (FA)", "#4f7cac", "s", "-"),
    "Glucose": ("Glucose (G)", "#b8860b", "^", "-"),
    "GlucoseAmm": ("Glucose ammonia (GA)", "#d9a441", "v", "-"),
    "PyruvicAcid": ("Pyruvic acid (PA)", "#a03623", "D", ":"),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = defaultdict(list)
    with (REPO / "ProcessedData" / "SI" / "figure_pathway_energetics.csv").open() as handle:
        for r in csv.DictReader(handle):
            data[r["network"]].append(
                (int(r["length"]), int(r["pathways"]), float(r["mean_dg_per_step_kJ_mol"]))
            )

    use()
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 2.75))
    for ax, col, ylabel, logy in (
        (axes[0], 1, "routes to a matched analogue", True),
        (axes[1], 2, "$\\Delta_\\mathrm{r}G'^{\\circ}$ per reaction  (kJ mol$^{-1}$)", False),
    ):
        for net, (label, colour, marker, ls) in STYLE.items():
            pts = sorted(data.get(net, []))
            if not pts:
                continue
            ax.plot(
                [p[0] for p in pts],
                [p[col] for p in pts],
                ls,
                color=colour,
                marker=marker,
                markersize=3.0,
                linewidth=1.0,
                label=label,
                clip_on=False,
                markerfacecolor="white" if net == "PyruvicAcid" else colour,
                markeredgewidth=0.8,
            )
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel("route length (reactions in the longest chain)")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0, 21)
        ax.set_xticks([1, 5, 10, 15, 20])

    axes[1].axhline(0, color="#999", linewidth=0.5, zorder=0)
    axes[1].set_ylim(-45, 2)
    # PA contributes a single spontaneous route; say so rather than let one marker pass
    axes[0].annotate(
        "PA: 1 route",
        xy=(3, 1),
        xytext=(7.5, 1.8),
        fontsize=6,
        color=STYLE["PyruvicAcid"][1],
        arrowprops=dict(
            arrowstyle="-", linewidth=0.5, color=STYLE["PyruvicAcid"][1], shrinkA=1, shrinkB=2
        ),
    )
    for ax, tag in zip(axes, "ab", strict=True):
        ax.text(
            -0.16,
            1.03,
            f"({tag})",
            transform=ax.transAxes,
            fontsize=8,
            fontweight="bold",
            va="bottom",
        )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=5,
        fontsize=6.5,
        handlelength=1.8,
        columnspacing=1.4,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.tight_layout(w_pad=2.0, rect=(0, 0.07, 1, 1))
    save(fig, str(OUT / "Figure_7_pathway_energetics"))
    print("wrote figures/output/Figure_7_pathway_energetics.png")


if __name__ == "__main__":
    main()
