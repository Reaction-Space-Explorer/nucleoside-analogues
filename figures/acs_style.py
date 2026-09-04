"""Matplotlib defaults sized for ACS figures.

    from acs_style import figure, save
    fig, ax = figure("single", 2.4)
    save(fig, "figures/out/Figure_5")

Arial leads the font stack because Helvetica on macOS lacks the arrow and
subscript glyphs these figures use.
"""

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")

#: ACS one- and two-column widths, inches.
SINGLE, DOUBLE = 3.3, 7.0

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.frameon": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
}


def use() -> None:
    plt.rcParams.update(RC)


def figure(width: str | float = "single", height: float = 2.4, **kwargs):
    """A figure whose width is an ACS column."""
    use()
    inches = {"single": SINGLE, "double": DOUBLE}.get(width, width)
    return plt.subplots(figsize=(float(inches), height), **kwargs)


def save(fig, stem: str) -> None:
    """Write PNG and PDF beside each other."""
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
