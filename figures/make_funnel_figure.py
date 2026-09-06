"""Reachability funnel: how many analogues survive each filter.

    uv run --extra figures python figures/make_funnel_figure.py

Products -> matched analogues -> spontaneously reachable -> reachable using
only reactions whose free energy can be estimated. Numbers from
ProcessedData/SI/figure_funnel.csv.
"""

import csv
from pathlib import Path

import numpy as np
from acs_style import DOUBLE, save, use

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "figures" / "output"
LABEL = {"Formose": "Formose\n(F)", "FormoseAmm": "Formose\nAmmonia (FA)",
         "Glucose": "Glucose\n(G)", "GlucoseAmm": "Glucose\nAmmonia (GA)",
         "PyruvicAcid": "Pyruvic\nAcid (PA)"}
STAGES = [("products", "network products", "#c9d1d9"),
          ("matched", "matched analogues", "#8fa8c8"),
          ("reachable_with_unestimable", "spontaneously reachable", "#4f7cac"),
          ("reachable_estimable_only", "…using estimable reactions only", "#1f3d63")]

rows = list(csv.DictReader((REPO / "ProcessedData" / "SI" / "figure_funnel.csv").open()))
use()
import matplotlib.pyplot as plt  # noqa: E402

fig, (ax, axr) = plt.subplots(1, 2, figsize=(DOUBLE, 3.1), gridspec_kw={"width_ratios": [1.35, 1]})

x = np.arange(len(rows))
width = 0.2
for k, (key, label, colour) in enumerate(STAGES):
    vals = [int(r[key]) for r in rows]
    ax.bar(x + (k - 1.5) * width, vals, width, color=colour, label=label, edgecolor="none")
    for xi, v in zip(x + (k - 1.5) * width, vals, strict=True):
        ax.text(xi, v * 1.15, f"{v:,}", ha="center", va="bottom", fontsize=4.6,
                rotation=90, color="#333333")
ax.set_yscale("log")
ax.set_ylim(0.6, 5e5)
ax.set_xticks(x)
ax.set_xticklabels([LABEL[r["network"]] for r in rows], fontsize=6)
ax.set_ylabel("structures")
ax.legend(fontsize=5.6, frameon=False, loc="upper center", ncol=2, columnspacing=0.9,
          handlelength=1.1, handletextpad=0.4, borderpad=0.2,
          bbox_to_anchor=(0.5, 1.30))
ax.set_title("A", loc="left", fontweight="bold", fontsize=8, y=1.16)

# right panel: the same as survival fractions, which is what the funnel is about
for r in rows:
    m, s, e = int(r["matched"]), int(r["reachable_with_unestimable"]), int(r["reachable_estimable_only"])
    axr.plot([0, 1, 2], [100, 100 * s / m, 100 * e / m], marker="o", markersize=3,
             linewidth=1, label=LABEL[r["network"]].replace("\n", " "))
axr.set_xticks([0, 1, 2])
axr.set_xticklabels(["matched", "spontaneously\nreachable", "estimable\nreactions only"], fontsize=6)
axr.set_ylabel("% of matched analogues")
axr.set_ylim(-4, 108)
axr.legend(fontsize=5.6, frameon=False, loc="lower left", handlelength=1.2, handletextpad=0.4)
axr.set_title("B", loc="left", fontweight="bold", fontsize=8, y=1.16)

fig.tight_layout(pad=0.4, w_pad=1.4, rect=(0, 0, 1, 0.90))
OUT.mkdir(parents=True, exist_ok=True)
save(fig, str(OUT / "Figure_reachability_funnel"))
print("wrote", OUT / "Figure_reachability_funnel.png")
