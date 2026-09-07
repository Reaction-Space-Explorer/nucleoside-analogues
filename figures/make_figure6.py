"""Figure 6: the four formose routes, tiled from the autocycle renderings.

    uv run --extra figures python figures/make_figure6.py

The panels are the spontaneous routes to threose, ribose, glycerol and
deoxyribose in the Formose CRNR. FormoseAmm reaches all four by the same
chemistry, so only one network is shown. Panels are rendered from the same
specs as figures/routes/ but without their own legends and titles, which the
composite supplies once; the figure therefore cannot drift from the data.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from acs_style import DOUBLE, save, use
from PIL import Image

REPO = Path(__file__).resolve().parent.parent
PANELS_DIR = REPO / "figures" / "routes" / "panels"
OUT = REPO / "figures" / "output"

PANELS = [
    ("a", "Formose_Threose", "threose, one step"),
    ("b", "Formose_Ribose", "ribose, two steps"),
    ("c", "Formose_Glycerol", "glycerol, two steps"),
    ("d", "Formose_Deoxyribose", "2-deoxyribose, three steps"),
]


def trimmed(path: Path, pad: int = 12) -> np.ndarray:
    """The panel with its white margin removed, so a sparse route and a dense
    one are drawn at comparable scale rather than one being padded out."""
    image = Image.open(path).convert("RGB")
    grey = np.asarray(image.convert("L"))
    ink = np.argwhere(grey < 250)
    top, left = ink.min(axis=0)
    bottom, right = ink.max(axis=0)
    box = (max(0, left - pad), max(0, top - pad),
           min(image.width, right + pad), min(image.height, bottom + pad))
    return np.asarray(image.crop(box))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    use()
    images = {stem: trimmed(PANELS_DIR / f"{stem}.png") for _, stem, _ in PANELS}
    # the top pair are wide and the bottom pair tall, so equal rows would squash
    # the tall ones; give each row the height its own panels need
    cell = DOUBLE / 2
    rows = [
        max(images[PANELS[i][1]].shape[0] / images[PANELS[i][1]].shape[1] for i in pair)
        for pair in ((0, 1), (2, 3))
    ]
    fig, axes = plt.subplots(
        2, 2, figsize=(DOUBLE, cell * sum(rows) + 0.45),
        gridspec_kw={"height_ratios": rows},
    )
    for ax, (letter, stem, title) in zip(axes.flat, PANELS, strict=True):
        ax.imshow(images[stem])
        ax.set_axis_off()
        ax.set_title(title, fontsize=7.5, pad=2)
        ax.text(-0.02, 1.0, f"({letter})", transform=ax.transAxes, fontsize=9,
                fontweight="bold", va="top", ha="right")
    fig.tight_layout(w_pad=0.4, h_pad=0.8)
    save(fig, str(OUT / "Figure_6_formose_routes"))
    print("wrote figures/output/Figure_6_formose_routes.png")


if __name__ == "__main__":
    main()
