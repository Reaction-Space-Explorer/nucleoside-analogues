"""The MS comparison numbers the manuscript quotes, pinned to the deposited table.

The raw FT-ICR peak lists are not redistributable and live outside the repo, so
this checks the deposited summary rather than recomputing it. It exists because
the caption once quoted a ceiling of 315 Da, a range of 45 to 70% and a count of
16 of 511, none of which any script produced; all three were wrong.
"""

import csv

from helpers import REPO

CUTOFF_NETWORKS = {"Formose", "Glucose", "GlucoseAmm", "PyruvicAcid"}


def rows() -> list[dict]:
    with (REPO / "ProcessedData" / "SI" / "ms_validation.csv").open() as handle:
        return list(csv.DictReader(handle))


def test_manuscript_ms_numbers() -> None:
    table = rows()
    ceilings = {r["network"]: float(r["ceiling_Da"]) for r in table}
    # four CRNRs were built under a 200 amu cutoff; formose-ammonia was not
    assert ceilings["FormoseAmm"] == 312.15
    for network in CUTOFF_NETWORKS:
        assert 198.0 <= ceilings[network] <= 200.0, network

    within = [float(r["percent_lt_ceiling"]) for r in table if r["network"] in CUTOFF_NETWORKS]
    assert round(min(within)) == 50, min(within)
    assert round(max(within)) == 71, max(within)

    amm = next(r for r in table if r["network"] == "FormoseAmm")
    assert (int(amm["matched_lt_ceiling"]), int(amm["formulas_lt_ceiling"])) == (16, 492)
    assert float(amm["percent_lt_ceiling"]) == 3.3
