"""Figure 2: how far the five CRNRs make the same products.

    uv run --extra figures python figures/make_overlap_figure.py

Each cell is the share of the two networks' combined product sets that both
reach, so the matrix is symmetric and the diagonal carries no information.
Species are taken from the reaction listings rather than the deposited product
files, which for Formose stop a generation short.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from acs_style import SINGLE, save, use

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from make_si_tables import PRODUCTS, species_generations

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "figures" / "output"
SHORT = {"Formose": "F", "FormoseAmm": "FA", "Glucose": "G",
         "GlucoseAmm": "GA", "PyruvicAcid": "PA"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    nets = list(PRODUCTS)
    species = {n: set(species_generations(n)) for n in nets}
    n = len(nets)
    m = np.full((n, n), np.nan)
    for i, a in enumerate(nets):
        for j, b in enumerate(nets):
            if i != j:
                m[i, j] = 100 * len(species[a] & species[b]) / len(species[a] | species[b])

    use()
    fig, ax = plt.subplots(figsize=(SINGLE, 2.9))
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#f4f4f4")
    im = ax.imshow(np.ma.masked_invalid(m), cmap=cmap, vmin=0, vmax=35)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            ax.text(j, i, f"{m[i, j]:.1f}", ha="center", va="center", fontsize=6.4,
                    color="white" if m[i, j] > 20 else "#222222")
    ax.set_xticks(range(n), [SHORT[x] for x in nets], fontsize=7)
    ax.set_yticks(range(n), [SHORT[x] for x in nets], fontsize=7)
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.4)
    ax.tick_params(which="both", length=0)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    bar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    bar.set_label("products reached by both, % of their union", fontsize=6.5)
    bar.ax.tick_params(labelsize=6.5)
    bar.outline.set_visible(False)
    fig.tight_layout(pad=0.4)
    save(fig, str(OUT / "Figure_2_network_overlap"))
    print("wrote figures/output/Figure_2_network_overlap.png")


if __name__ == "__main__":
    main()
