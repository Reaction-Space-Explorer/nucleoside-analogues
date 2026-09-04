"""Mass-spectrometry validation figure.

    uv run --extra figures python figures/make_ms_figure.py

Convention follows the Chem. Sci. ESI for this workflow: experimental in blue,
model in magenta.
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
CUTOFF = 200.0          # the networks' own mass ceiling
BLUE, MAGENTA, GREY = "#1f5fa8", "#c2258a", "#b9bfc7"
SHOWCASE = "50"         # formose 85 C


def counts(formula: str) -> dict[str, int]:
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula) if el}


use()
fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 2.4))

# ---------------------------------------------------------------- A: mirror
network, products_file, label = SAMPLES[SHOWCASE]
peaks = [p for p in read_midas(next(MS.glob(f"*_{SHOWCASE}_*"))) if p["organic"]]
peaks = [p for p in peaks if p["mass"] < CUTOFF]
top = max(p["abundance"] for p in peaks)
net = network_formulas(products_file)
frame = read_products(PRODUCTS / products_file)
model = sorted({round(ExactMolWt(m), 4) for s in frame["Smiles"]
                if (m := Chem.MolFromSmiles(str(s))) and ExactMolWt(m) < CUTOFF})

ax = axes[0]
for p in peaks:
    matched = neutral_formula(p["ion"]) in net
    ax.vlines(p["mass"], 0, 100 * p["abundance"] / top,
              color=MAGENTA if matched else BLUE, linewidth=0.5,
              zorder=3 if matched else 2)
hist, edges = np.histogram(model, bins=np.arange(0, CUTOFF + 4, 4))
ax.bar(edges[:-1], -100 * hist / hist.max(), width=3.4, align="edge",
       color=GREY, edgecolor="none")
ax.axvspan(160, CUTOFF, color="#f0e4ee", zorder=0, linewidth=0)
ax.axhline(0, color="#1a1a1a", linewidth=0.6)
ax.text(180, 108, "overlap", ha="center", fontsize=5.4, color="#8a5878", va="top")
ax.set_xlim(50, CUTOFF)
ax.set_ylim(-115, 115)
ax.set_yticks([-100, -50, 0, 50, 100])
ax.set_yticklabels(["100", "50", "0", "50", "100"])
ax.set_xlabel("neutral mass (Da)")
ax.set_ylabel("relative abundance")
ax.text(0.03, 0.95, "experimental", transform=ax.transAxes, fontsize=5.8,
        color=BLUE, va="top")
ax.text(0.03, 0.05, "network products", transform=ax.transAxes, fontsize=5.8,
        color="#6b7280", va="bottom")
ax.set_title(label, fontsize=7, pad=3)

# ---------------------------------------------------------------- B: van Krevelen
ax = axes[1]
for p in peaks:
    c = counts(neutral_formula(p["ion"]))
    if not c.get("C"):
        continue
    matched = neutral_formula(p["ion"]) in net
    ax.scatter(c.get("O", 0) / c["C"], c.get("H", 0) / c["C"], s=3.5,
               color=MAGENTA if matched else GREY, linewidths=0,
               zorder=3 if matched else 2, alpha=0.9 if matched else 0.55)
ax.set_xlabel("O/C"); ax.set_ylabel("H/C")
ax.set_xlim(0, 1.6); ax.set_ylim(0, 2.6)
ax.set_title("van Krevelen, < 200 Da", fontsize=7, pad=3)

# ---------------------------------------------------------------- C: recovery
ax = axes[2]
bars = []
for number, (net_name, pfile, lab) in SAMPLES.items():
    ps = [p for p in read_midas(next(MS.glob(f"*_{number}_*")))
          if p["organic"] and p["mass"] < CUTOFF]
    formulas = {neutral_formula(p["ion"]) for p in ps}
    if len(formulas) < 10:
        bars.append((lab, None, len(formulas)))
        continue
    hit = formulas & set(network_formulas(pfile))
    bars.append((lab, 100 * len(hit) / len(formulas), len(formulas)))
bars.sort(key=lambda b: (b[1] is None, -(b[1] or 0)))
ypos = range(len(bars))
ax.barh(list(ypos), [b[1] or 0 for b in bars],
        color=[GREY if b[1] is None else MAGENTA for b in bars], height=0.62)
for y, (lab, pct, n) in zip(ypos, bars, strict=True):
    ax.text(2, y, f"{lab}  (n={n})" if pct is not None else f"{lab}  (n={n}, too few)",
            va="center", fontsize=5.3, color="white" if pct and pct > 25 else "#333333")
ax.set_yticks([]); ax.invert_yaxis()
ax.set_xlabel("% of experimental formulas\nrecovered by the network")
ax.set_xlim(0, 80)
ax.set_title("recovery below 200 Da", fontsize=7, pad=3)

for ax, letter in zip(axes, "ABC", strict=True):
    ax.text(-0.20, 1.06, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="top")

fig.tight_layout(pad=0.4, w_pad=1.0)
OUT.mkdir(parents=True, exist_ok=True)
save(fig, str(OUT / "Figure_MS_validation"))
print("wrote", OUT / "Figure_MS_validation.png")
