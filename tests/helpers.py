"""Paths and guards for tests that read the deposited data."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ORIGINAL = REPO / "OriginalData" / "OriginalNetworkData"
PROCESSED = REPO / "ProcessedData"

#: Networks reported in the manuscript, with their product listings.
NETWORKS: dict[str, str] = {
    "Formose": "formose_output.tsv",
    "FormoseAmm": "formose_amm_output.tsv",
    "Glucose": "glucose_degradation_output.tsv",
    "GlucoseAmm": "glucose_amm_output.tsv",
    "PyruvicAcid": "pyruvic_output.tsv",
}


def products_path(network: str) -> Path:
    return ORIGINAL / "Products" / NETWORKS[network]


def processed_rels_path(network: str, generation: int) -> Path:
    return PROCESSED / "RelsFiles" / network / f"{network}G{generation}ProcessedRels.tsv"


def spontaneous_rels_path(network: str) -> Path:
    return (
        PROCESSED / "SpontaneousRelsWithThermoFiles" / f"Spontaneous{network}G3RelsWithThermo.tsv"
    )


def requires(path: Path) -> None:
    """Skip rather than fail when the deposited data is not checked out."""
    if not path.exists():
        pytest.skip(f"deposited data not present: {path.relative_to(REPO)}")
