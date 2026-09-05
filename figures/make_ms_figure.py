"""Mass-spectrometry validation figure, all five networks.

    uv run --extra figures python figures/make_ms_figure.py

Conventions follow plots/ in Reaction-Space-Explorer/reac-space-exp, which
produced the Chem. Sci. ESI for this workflow:

  mirror plot    experimental in cornflower blue above, network in deep pink
                 below, both on a log abundance axis, over 155-205 Da
  van Krevelen   experimental as open black circles, network filled and
                 coloured by generation

The network's y axis is the number of structures sharing each exact mass,
normalised to 100, as in mirror_plot_spectrum.ipynb -- not a binned histogram.
"""

import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from acs_style import DOUBLE, save, use
from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from ms_validation import MS, PRODUCTS, SAMPLES, neutral_formula, read_midas  # noqa: E402

from nucleoside_analogues.rels import read_products  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = REPO / "figures" / "output"
#: Upper bound per network. Four were generated under a 200 amu cutoff on
#: product mass; the ammonia-seeded formose network was not, and reaches 312 Da,
#: so comparing it at 200 would manufacture a gap that is not there.
CEILING = {"Formose": 200.0, "FormoseAmm": 315.0, "Glucose": 200.0,
           "GlucoseAmm": 200.0, "PyruvicAcid": 201.0}
#: Peak lists were exported from m/z 150 upward, so nothing below this is seen.
FLOOR = 155.0

EXP, MODEL = "cornflowerblue", "deeppink"
PANELS = [("50", "Formose (F)"), ("40", "Formose Ammonia (FA)"), ("38", "Glucose (G)"),
          ("37", "Glucose Ammonia (GA)"), ("46", "Pyruvic Acid (PA)")]
GEN_COLOURS = plt.get_cmap("Spectral")


def counts(formula: str) -> dict[str, int]:
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula) if el}


use()
fig = plt.figure(figsize=(DOUBLE, 4.8))
outer = fig.add_gridspec(2, 1, height_ratios=[2.0, 1.7], hspace=0.42)
mirror = outer[0].subgridspec(2, len(PANELS), hspace=0, wspace=0.34)
lower = outer[1].subgridspec(1, len(PANELS), wspace=0.34)
axes = [[fig.add_subplot(mirror[0, c]) for c in range(len(PANELS))],
        [fig.add_subplot(mirror[1, c]) for c in range(len(PANELS))],
        [fig.add_subplot(lower[0, c]) for c in range(len(PANELS))]]

for column, (number, title) in enumerate(PANELS):
    network, products_file, _ = SAMPLES[number]
    LOW, HIGH = FLOOR, CEILING[network]
    peaks = [p for p in read_midas(next(MS.glob(f"*_{number}_*")))
             if p["organic"] and FLOOR <= p["mass"] <= HIGH]
    frame = read_products(PRODUCTS / products_file)

    masses, generations = [], []
    for smiles, generation in zip(frame["Smiles"], frame["Generation"], strict=True):
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is not None:
            masses.append(round(ExactMolWt(mol), 4))
            generations.append(int(generation))
    frequency = Counter(masses)
    peak_of = max(frequency.values())

    # ---- experimental, above
    ax = axes[0][column]
    window = peaks
    if window:
        top = max(p["abundance"] for p in window)
        ax.vlines([p["mass"] for p in window], 0.875,
                  [100 * p["abundance"] / top for p in window], color=EXP, linewidth=0.6)
    else:
        ax.text((LOW + HIGH) / 2, 11, "no assigned peak in range",
                ha="center", va="center", fontsize=5.2, color="#7a8290")
    ax.set_yscale("log"); ax.set_ylim(0.875, 125); ax.set_xlim(LOW, HIGH)
    ax.set_xticklabels([]); ax.set_title(title, fontsize=7, pad=3)
    ax.tick_params(labelsize=5.6)

    # ---- network, below, inverted
    ax = axes[1][column]
    inside = [(m, f) for m, f in frequency.items() if LOW <= m <= HIGH]
    if inside:
        ax.vlines([m for m, _ in inside], 0.875,
                  [100 * f / peak_of for _, f in inside], color=MODEL, linewidth=0.6)
    ax.set_yscale("log"); ax.set_ylim(0.875, 125); ax.set_xlim(LOW, HIGH)
    ax.invert_yaxis(); ax.set_xlabel("exact mass (Da)", fontsize=6.2)
    ax.tick_params(labelsize=5.6)

    # ---- van Krevelen
    ax = axes[2][column]
    for generation in sorted(set(generations), reverse=True):
        xs, ys = [], []
        for mass, gen in zip(masses, generations, strict=True):
            del mass
            if gen != generation:
                continue
        sel = [i for i, g in enumerate(generations) if g == generation]
        for i in sel:
            mol = Chem.MolFromSmiles(str(frame["Smiles"].iloc[i]))
            c = counts(Chem.rdMolDescriptors.CalcMolFormula(mol).replace("+", "").replace("-", ""))
            if c.get("C"):
                xs.append(c.get("O", 0) / c["C"]); ys.append(c.get("H", 0) / c["C"])
        ax.scatter(xs, ys, s=6, color=GEN_COLOURS(generation / max(generations)),
                   alpha=0.55, linewidths=0, zorder=2, label=f"G{generation}")
    xs, ys = [], []
    for p in peaks:
        c = counts(neutral_formula(p["ion"]))
        if c.get("C"):
            xs.append(c.get("O", 0) / c["C"]); ys.append(c.get("H", 0) / c["C"])
    ax.scatter(xs, ys, s=5, facecolors="none", edgecolors="black",
               linewidths=0.18, alpha=0.12, zorder=3)
    ax.set_xlim(0, 1.6); ax.set_ylim(0.3, 2.8)
    ax.set_xlabel("O/C"); ax.tick_params(labelsize=5.6)
    if column == 0:
        ax.set_ylabel("H/C")

handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=3,
                      color=GEN_COLOURS(g / 6), label=f"G{g}") for g in range(1, 7)]
handles.append(plt.Line2D([], [], marker="o", linestyle="", markersize=3,
                          markerfacecolor="none", markeredgecolor="black",
                          markeredgewidth=0.4, label="FT-ICR-MS"))
fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=5.8,
           frameon=False, bbox_to_anchor=(0.5, -0.015), handletextpad=0.3,
           columnspacing=1.1)
axes[0][0].set_ylabel("experimental", fontsize=6, color=EXP)
axes[1][0].set_ylabel("network", fontsize=6, color=MODEL)
for row, letter in ((0, "A"), (2, "B")):
    axes[row][0].text(-0.40, 1.12, letter, transform=axes[row][0].transAxes,
                      fontsize=8, fontweight="bold", va="top")
OUT.mkdir(parents=True, exist_ok=True)
save(fig, str(OUT / "Figure_MS_validation"))
print("wrote", OUT / "Figure_MS_validation.png")
