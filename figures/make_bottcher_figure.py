"""Bottcher complexity of matched analogues, by generation.

    uv run --extra figures python figures/make_bottcher_figure.py

Reads the complexity values deposited in ProcessedData/ComplexityData, which
were computed by notebooks/pipeline/BottcherComplexity.ipynb.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from acs_style import DOUBLE, save, use

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "ProcessedData" / "ComplexityData"
OUT = REPO / "figures" / "workflow"
NETWORKS = [
    ("Formose", "Formose (F)"),
    ("FormoseAmm", "Formose Ammonia (FA)"),
    ("Glucose", "Glucose (G)"),
    ("GlucoseAmm", "Glucose Ammonia (GA)"),
    ("PyruvicAcid", "Pyruvic Acid (PA)"),
]
#: Bottcher Cm of the four target nucleosides, from the same data.
TARGETS = {"Glycerol": 36.68, "Threose": 92.04, "Deoxyribose": 98.04, "Ribose": 111.55}
INK, BAND = "#1a1a1a", "#cfd8e3"

use()
fig, axes = plt.subplots(1, len(NETWORKS), figsize=(DOUBLE, 2.5), sharey=True)
lo, hi = min(TARGETS.values()), max(TARGETS.values())

for ax, (key, label) in zip(axes, NETWORKS, strict=True):
    frame = pd.read_csv(DATA / f"{key}ComplexityData.tsv", sep="\t", index_col=0)
    gens = sorted(frame["Generation"].unique(), key=lambda g: int(g[1:]))
    groups = [frame.loc[frame["Generation"] == g, "Complexity"].to_numpy() for g in gens]

    ax.axhspan(lo, hi, color=BAND, zorder=0, linewidth=0)
    parts = ax.boxplot(
        groups, positions=range(len(gens)), widths=0.62, showfliers=False,
        patch_artist=True, medianprops={"color": INK, "linewidth": 1.1},
        boxprops={"facecolor": "white", "edgecolor": INK, "linewidth": 0.6},
        whiskerprops={"color": INK, "linewidth": 0.6},
        capprops={"color": INK, "linewidth": 0.6},
    )
    del parts
    for i, values in enumerate(groups):
        ax.text(i, 205, f"{len(values):,}", ha="center", va="center", fontsize=5.2,
                color="#4a4a4a")
    ax.set_xticks(range(len(gens)))
    ax.set_xticklabels([g[1:] for g in gens])
    ax.set_title(label, fontsize=7, pad=3)
    ax.set_xlabel("generation")
    ax.set_ylim(20, 215)

axes[0].set_ylabel("Böttcher complexity $C_m$")
axes[0].set_yticks([40, 60, 80, 100, 120, 140, 160, 180])
handle = plt.Rectangle((0, 0), 1, 1, facecolor=BAND, edgecolor="none")
fig.legend([handle], ["range spanned by the four target nucleosides"], loc="lower center",
           bbox_to_anchor=(0.5, -0.04), fontsize=6, frameon=False)
fig.tight_layout(pad=0.4, w_pad=0.5, rect=(0, 0.04, 1, 1))
OUT.mkdir(parents=True, exist_ok=True)
save(fig, str(OUT / "Figure_Bottcher_complexity"))
print("wrote", OUT / "Figure_Bottcher_complexity.png")
