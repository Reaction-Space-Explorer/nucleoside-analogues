"""Figure S3: what each reaction rule contributes energetically, and how far pH moves it.

    uv run --extra figures python figures/make_rule_figure.py

The rule, not the network, is the unit: a rule is a general transformation, so
its reactions are pooled across the five CRNRs. Only the twenty rules firing
most often are drawn; the rest are in the accompanying file.
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from acs_style import DOUBLE, panel, save, use

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "figures" / "output"
SHOW = 20


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (REPO / "ProcessedData" / "SI" / "rule_energetics.csv").open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["q1_dG_pH7.4"]]
    rows = sorted(rows, key=lambda r: -int(r["reactions"]))[:SHOW]
    rows = sorted(rows, key=lambda r: float(r["median_dG_pH7.4"]))

    labels = [r["rule"].replace(" for ", " ").replace(" membered rings", "-ring") for r in rows]
    med = np.array([float(r["median_dG_pH7.4"]) for r in rows])
    q1 = np.array([float(r["q1_dG_pH7.4"]) for r in rows])
    q3 = np.array([float(r["q3_dG_pH7.4"]) for r in rows])
    shift = np.array([float(r["median_shift_pH7_to_11"] or 0) for r in rows])
    y = np.arange(len(rows))

    use()
    fig, axes = plt.subplots(
        1, 2, figsize=(DOUBLE, 3.6), sharey=True, gridspec_kw={"width_ratios": [1.55, 1]}
    )
    ax = axes[0]
    ax.axvline(0, color="#999", linewidth=0.6, zorder=0)
    ax.hlines(y, q1, q3, color="#9bb0c9", linewidth=3.4, zorder=2)
    ax.plot(med, y, "o", markersize=3.4, color="#1f3d63", zorder=3)
    ax.set_yticks(y, labels, fontsize=6.2)
    ax.set_xlabel("$\\Delta_\\mathrm{r}G'^{\\circ}$ at pH 7.4  (kJ mol$^{-1}$)")
    ax.set_ylim(-0.7, len(rows) - 0.3)

    ax2 = axes[1]
    ax2.axvline(0, color="#999", linewidth=0.6, zorder=0)
    ax2.barh(
        y,
        shift,
        height=0.62,
        zorder=2,
        color=["#a03623" if abs(s) > 5 else "#cfd8e3" for s in shift],
    )
    ax2.set_xlabel("median shift, pH 7 to 11  (kJ mol$^{-1}$)")
    for ax_, tag in zip(axes, "ab", strict=True):
        panel(ax_, tag, x=-0.02 if tag == "b" else -0.62, y=1.03)
    fig.tight_layout(w_pad=1.0)
    save(fig, str(OUT / "Figure_S3_rule_energetics"))
    print("wrote figures/output/Figure_S3_rule_energetics.png")


if __name__ == "__main__":
    main()
