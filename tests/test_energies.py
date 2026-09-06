def test_null_estimates_are_not_spontaneous():
    """Component contribution returns 0 +/- 0 where reagents and products share
    a decomposition. A point value of -1e-05 must not pass the 95% test on
    rounding, so such reactions count as unestimable."""
    from make_si_tables import is_null

    assert is_null({"dG_prime_kJ_mol": "-1.018336507740969e-05", "sigma_kJ_mol": "0.0"})
    assert is_null({"dG_prime_kJ_mol": "0.0", "sigma_kJ_mol": "0.0"})
    # a real estimate that happens to be small keeps its uncertainty
    assert not is_null({"dG_prime_kJ_mol": "-0.0001", "sigma_kJ_mol": "1.4"})
    # a real estimate that is large and certain is not null
    assert not is_null({"dG_prime_kJ_mol": "-13.4", "sigma_kJ_mol": "0.0"})
