"""Figure 7: how many spontaneous routes there are, and how favourable they are.

    uv run --extra figures python figures/make_pathway_figure.py

Three panels against route length: how many matched analogues are reached by a
route of that length, the free energy of the whole route, and the free energy
per reaction. The third is the one that carries information: the second falls
with length largely because a longer route sums more negative steps.

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
    with (REPO / "ProcessedData" / "SI" / "figure_pathway_energetics.csv").open() as h:
        for r in csv.DictReader(h):
            data[r["network"]].append(
                (int(r["length"]), int(r["pathways"]),
                 float(r["mean_dg_kJ_mol"]), float(r["mean_dg_per_step_kJ_mol"]))
            )

    use()
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 2.35))
    panels = [
        (0, 1, "routes to a matched analogue", True),
        (1, 2, "route $\\Delta_\\mathrm{r}G'^{\\circ}$  (kJ mol$^{-1}$)", False),
        (2, 3, "per reaction  (kJ mol$^{-1}$)", False),  # divided by distinct reactions
    ]
    for ax, col, ylabel, logy in panels:
        a = axes[ax]
        for net, (label, colour, marker, ls) in STYLE.items():
            pts = sorted(data.get(net, []))
            if not pts:
                continue
            xs = [p[0] for p in pts]
            ys = [p[col] for p in pts]
            a.plot(xs, ys, ls, color=colour, marker=marker, markersize=2.6,
                   linewidth=0.9, label=label, clip_on=False,
                   markerfacecolor="white" if net == "PyruvicAcid" else colour,
                   markeredgewidth=0.7)
        if logy:
            a.set_yscale("log")
        a.set_xlabel("route length (reactions in the longest chain)", fontsize=6.6)
        a.set_ylabel(ylabel, fontsize=6.8)
        a.set_xlim(0, 21)
        a.set_xticks([1, 5, 10, 15, 20])
        a.tick_params(labelsize=6.4)

    axes[2].axhline(0, color="#999", linewidth=0.5, zorder=0)
    # PA contributes a single spontaneous route; say so rather than let one marker pass
    axes[0].annotate("PA: 1 route", xy=(3, 1), xytext=(6.4, 1.9),
                     fontsize=5.8, color=STYLE["PyruvicAcid"][1],
                     arrowprops=dict(arrowstyle="-", linewidth=0.5,
                                     color=STYLE["PyruvicAcid"][1], shrinkA=1, shrinkB=2))
    axes[0].legend(loc="upper right", fontsize=5.8, handlelength=1.6, borderpad=0.3)
    for a, tag in zip(axes, "abc", strict=True):
        a.text(-0.20, 1.04, f"({tag})", transform=a.transAxes, fontsize=7.6,
               fontweight="bold", va="bottom")
    fig.tight_layout(w_pad=1.6)
    save(fig, str(OUT / "Figure_7_pathway_energetics"))
    print("wrote figures/output/Figure_7_pathway_energetics.png")


if __name__ == "__main__":
    main()
