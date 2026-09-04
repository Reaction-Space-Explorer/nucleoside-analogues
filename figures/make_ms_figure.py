"""Mass-spectrometry validation figure, all five networks.

    uv run --extra figures python figures/make_ms_figure.py

Experimental in blue, recovered by the network in magenta, following the
Chem. Sci. ESI for this workflow. One representative spectrum per network.
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from acs_style import DOUBLE, save, use
from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from ms_validation import (  # noqa: E402
    MS,
    PRODUCTS,
    SAMPLES,
    network_formulas,
    neutral_formula,
    read_midas,
)

from nucleoside_analogues.rels import read_products  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = REPO / "figures" / "workflow"
CUTOFF = 200.0        # the networks build nothing heavier
FLOOR = 160.0         # no sample has an assigned CHNO peak much below this
BLUE, MAGENTA, GREY, BAND = "#1f5fa8", "#c2258a", "#b9bfc7", "#f2e6ef"
#: One representative spectrum per network, best-matching where there is a choice.
PANELS = [("50", "Formose (F)"), ("40", "Formose Ammonia (FA)"), ("38", "Glucose (G)"),
          ("37", "Glucose Ammonia (GA)"), ("46", "Pyruvic Acid (PA)")]


def counts(formula: str) -> dict[str, int]:
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula) if el}


use()
fig, axes = plt.subplots(2, len(PANELS), figsize=(DOUBLE, 3.9))

for column, (number, title) in enumerate(PANELS):
    network, products_file, _ = SAMPLES[number]
    peaks = [p for p in read_midas(next(MS.glob(f"*_{number}_*")))
             if p["organic"] and p["mass"] < CUTOFF]
    net = network_formulas(products_file)
    frame = read_products(PRODUCTS / products_file)
    model = [ExactMolWt(m) for s in frame["Smiles"]
             if (m := Chem.MolFromSmiles(str(s))) and ExactMolWt(m) < CUTOFF]

    # ---- top row: mirror of experiment against the network's product masses
    ax = axes[0, column]
    ax.axvspan(FLOOR, CUTOFF, color=BAND, zorder=0, linewidth=0)
    if peaks:
        top = max(p["abundance"] for p in peaks)
        for p in peaks:
            hit = neutral_formula(p["ion"]) in net
            ax.vlines(p["mass"], 0, 100 * p["abundance"] / top,
                      color=MAGENTA if hit else BLUE, linewidth=0.5,
                      zorder=3 if hit else 2)
    else:
        ax.text(125, 55, "no assigned peak\nbelow 200 Da", ha="center", va="center",
                fontsize=5.4, color="#8a5878")
    hist, edges = np.histogram(model, bins=np.arange(0, CUTOFF + 5, 5))
    ax.bar(edges[:-1], -100 * hist / hist.max(), width=4.3, align="edge",
           color=GREY, edgecolor="none")
    ax.axhline(0, color="#1a1a1a", linewidth=0.6)
    ax.set_xlim(50, CUTOFF)
    ax.set_ylim(-118, 118)
    ax.set_yticks([-100, 0, 100])
    ax.set_yticklabels(["100", "0", "100"])
    ax.set_title(title, fontsize=7, pad=3)
    if column == 0:
        ax.set_ylabel("relative abundance")
        ax.text(0.04, 0.96, "experimental", transform=ax.transAxes, fontsize=5.2,
                color=BLUE, va="top")
        ax.text(0.04, 0.04, "network", transform=ax.transAxes, fontsize=5.2,
                color="#6b7280", va="bottom")

    # ---- bottom row: van Krevelen of the assigned formulas
    ax = axes[1, column]
    hits = 0
    for p in peaks:
        formula = neutral_formula(p["ion"])
        c = counts(formula)
        if not c.get("C"):
            continue
        hit = formula in net
        hits += hit
        ax.scatter(c.get("O", 0) / c["C"], c.get("H", 0) / c["C"], s=4,
                   color=MAGENTA if hit else GREY, linewidths=0,
                   zorder=3 if hit else 2, alpha=0.9 if hit else 0.5)
    ax.set_xlim(0, 1.5)
    ax.set_ylim(0.4, 2.7)
    ax.set_xlabel("O/C")
    if column == 0:
        ax.set_ylabel("H/C")
    pct = f"{100 * hits / len(peaks):.0f}%" if peaks else "n/a"
    ax.text(0.96, 0.95, f"{hits}/{len(peaks)}\n{pct}", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.8, color=MAGENTA if peaks else "#8a8f98")

for ax in axes[0]:
    ax.set_xlabel("neutral mass (Da)", fontsize=6.5)
axes[0, 0].text(-0.34, 1.12, "A", transform=axes[0, 0].transAxes, fontsize=8,
                fontweight="bold", va="top")
axes[1, 0].text(-0.34, 1.10, "B", transform=axes[1, 0].transAxes, fontsize=8,
                fontweight="bold", va="top")
fig.tight_layout(pad=0.4, w_pad=0.7, h_pad=1.0)
OUT.mkdir(parents=True, exist_ok=True)
save(fig, str(OUT / "Figure_MS_validation"))
print("wrote", OUT / "Figure_MS_validation.png")
