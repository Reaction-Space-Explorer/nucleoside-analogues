"""Figure S4: the two routes to glyoxylate, carbon against nitrogen.

    uv run --extra figures python figures/make_glyoxylate_figure.py

Glyoxylate is the one target the ammonia-seeded formose network reaches more
directly than its parent, and it does so by different chemistry rather than by a
shortcut: F builds it in six steps of carbon chemistry, ending in a benzilic acid
rearrangement and a retro-aldol, while FA reaches it in three through glycine,
formed by a Cannizzaro reduction and carried to glyoxylate by transamination.

Panels are the autocycle renderings of figures/routes/*_Glyoxylate.yaml drawn
without their own legends, so the figure cannot drift from the traced data. Both
are placed at the same scale, not stretched to a common height, so a bond in one
is a bond in the other; the panels differ in shape because the routes do.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from acs_style import DOUBLE, save, use
from make_figure6 import trimmed

REPO = Path(__file__).resolve().parent.parent
PANELS_DIR = REPO / "figures" / "routes" / "panels"
OUT = REPO / "figures" / "output"

PANELS = [
    ("a", "Formose_Glyoxylate", "F: six steps, carbon chemistry throughout"),
    ("b", "FormoseAmm_Glyoxylate", "FA: three steps, through glycine"),
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    use()
    images = [trimmed(PANELS_DIR / f"{stem}.png") for _, stem, _ in PANELS]
    widths = [im.shape[1] for im in images]
    # equal scale: the width each panel gets is proportional to its own width,
    # so both are reduced by the same factor and structures stay comparable
    scale = DOUBLE / sum(widths)
    height = max(im.shape[0] for im in images) * scale
    fig, axes = plt.subplots(
        1, 2, figsize=(DOUBLE, height + 0.45), gridspec_kw={"width_ratios": widths}
    )
    for ax, image, (letter, _, title) in zip(axes, images, PANELS, strict=True):
        ax.imshow(image)
        ax.set_axis_off()
        ax.set_anchor("N")  # top-align; the shorter panel keeps its own height
        ax.set_title(title, fontsize=7.5, pad=3)
        ax.text(
            -0.02,
            1.02,
            f"({letter})",
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            va="top",
            ha="right",
        )
    fig.tight_layout(w_pad=1.2)
    save(fig, str(OUT / "Figure_S4_glyoxylate_routes"))
    print("wrote figures/output/Figure_S4_glyoxylate_routes.png")


if __name__ == "__main__":
    main()
