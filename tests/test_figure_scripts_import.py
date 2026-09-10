"""Every figure script must at least import.

The gates lint and type-check figures/ but never execute it, so a rename in
scripts/ can leave a figure script importing a name that no longer exists.
That happened: replacing ms_validation.MS with ms_dir() broke
make_ms_figure.py, and nothing caught it until the figure was rebuilt by hand.

Importing is enough to catch it. Several of these scripts render at import
time rather than under a __main__ guard, so importing also exercises them;
make_ms_figure needs MS_DATA and is skipped without it.
"""

import importlib
import os
import sys
from pathlib import Path

import pytest

FIGURES = Path(__file__).resolve().parent.parent / "figures"
SCRIPTS = sorted(p.stem for p in FIGURES.glob("make_*.py"))


@pytest.mark.parametrize("module", SCRIPTS)
def test_figure_script_imports(module: str) -> None:
    pytest.importorskip("matplotlib")
    if module == "make_ms_figure" and not os.environ.get("MS_DATA"):
        pytest.skip("needs MS_DATA; the peak lists are not in the repository")
    for path in (str(FIGURES), str(FIGURES.parent / "scripts")):
        if path not in sys.path:
            sys.path.insert(0, path)
    importlib.import_module(module)
