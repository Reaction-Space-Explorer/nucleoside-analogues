"""Figure 7: the hexose branch point, against the experiment of Yi et al.

    uv run --extra figures python figures/make_ketohexose_figure.py

Two routes over the spontaneous Formose network, from glycolaldehyde alone,
which is the feedstock Yi et al. used. They differ only in where the carbonyl
migration falls, and that is the whole result: migrate then add, and the route
is the one their labelling supports; add then migrate, and it runs through the
aldohexose they did not detect. Our search returns the second, being a step
shorter, so minimum step count picks the route experiment argues against.

Panels come from the same specs as figures/routes/case_*.yaml, rendered without
their own legends, so the figure cannot drift from the traced data. Arrow width
scales with the magnitude of the free energy, which is why the migrations, at
-2.0 kJ/mol, are nearly invisible beside the aldols.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from acs_style import DOUBLE, save, use
from make_figure6 import trimmed

REPO = Path(__file__).resolve().parent.parent
PANELS_DIR = REPO / "figures" / "routes" / "panels"
OUT = REPO / "figures" / "output"

PANELS = [
    ("a", "case_3ketohexose", "aldol, migration, aldol: the route Yi et al. assign"),
    ("b", "case_2ketohexose", "aldol, aldol, migration: the route minimum step count returns"),
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    use()
    images = [trimmed(PANELS_DIR / f"{stem}.png") for _, stem, _ in PANELS]
    ratios = [im.shape[0] / im.shape[1] for im in images]
    fig, axes = plt.subplots(
        2, 1, figsize=(DOUBLE, DOUBLE * sum(ratios) + 0.5), gridspec_kw={"height_ratios": ratios}
    )
    for ax, image, (letter, _, title) in zip(axes, images, PANELS, strict=True):
        ax.imshow(image)
        ax.set_axis_off()
        ax.set_title(title, fontsize=7.5, pad=3)
        ax.text(
            -0.015,
            1.0,
            f"({letter})",
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            va="top",
            ha="right",
        )
    fig.tight_layout(h_pad=1.0)
    save(fig, str(OUT / "Figure_7_ketohexose_branch"))
    print("wrote figures/output/Figure_7_ketohexose_branch.png")


if __name__ == "__main__":
    main()
