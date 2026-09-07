"""The control-comparison statistics quoted in the Results.

These were computed ad hoc for an earlier draft and deposited nowhere, which is
how the permutation p came to be quoted as 0.010 when the exact value over all
3,268,760 splits is 0.0074.
"""

import csv

import pytest
from helpers import REPO

SI = REPO / "ProcessedData" / "SI"


def stats() -> dict[str, str]:
    with (SI / "control_statistics.csv").open() as handle:
        return {r["statistic"]: r["value"] for r in csv.DictReader(handle)}


def test_per_pair_p_is_the_phipson_smyth_estimator() -> None:
    with (SI / "matched_controls.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 25
    for row in rows:
        faster, controls = int(row["controls_strictly_faster"]), int(row["controls"])
        assert float(row["p_value"]) == pytest.approx((faster + 1) / (controls + 1), abs=5e-5)
        # the estimator cannot return zero, which is the reason for citing it
        assert float(row["p_value"]) > 0


def test_fisher_matches_an_independent_implementation() -> None:
    scipy_stats = pytest.importorskip("scipy.stats")
    summary = stats()
    expected = scipy_stats.chi2.sf(float(summary["fisher_chi2"]), int(summary["fisher_df"]))
    assert float(summary["fisher_combined_p"]) == pytest.approx(expected, rel=1e-3)


def test_permutation_is_exact_and_separates_the_formose_networks() -> None:
    summary = stats()
    assert int(summary["permutation_splits"]) == 3_268_760
    assert float(summary["permutation_p_exact"]) == pytest.approx(0.00742, abs=5e-5)
    assert float(summary["mean_p_formose"]) < float(summary["mean_p_other"])
