"""Phase 2: Dryden gust model statistics (checked against closed-form theory)."""

import numpy as np
import pytest
from scipy.linalg import expm

from params import aerosonde as P
from sim import wind

VA0 = 25.0
CHANNELS = [("u", P.dryden_Lu, P.dryden_sigma_u), ("v", P.dryden_Lv, P.dryden_sigma_v),
            ("w", P.dryden_Lw, P.dryden_sigma_w)]


@pytest.mark.parametrize("ch, L, sigma", CHANNELS)
def test_filter_normalisation_gives_variance_sigma_squared(ch, L, sigma):
    A, B, C = wind.dryden_state_space(VA0, L, sigma, ch)
    _, var = wind.dryden_output_covariance(A, B, C)
    assert var == pytest.approx(sigma**2, rel=1e-10)


@pytest.mark.parametrize("ch, L, sigma", CHANNELS)
def test_autocorrelation_matches_dryden_closed_form(ch, L, sigma):
    """R_u = s^2 exp(-xi/L);  R_v,w = s^2 (1 - xi/(2L)) exp(-xi/L),  xi = Va tau.  Independent of the filter algebra."""
    A, B, C = wind.dryden_state_space(VA0, L, sigma, ch)
    Pcov, _ = wind.dryden_output_covariance(A, B, C)
    for tau in (0.0, 0.5, 2.0, 8.0, 20.0):
        xi = VA0 * tau
        R_model = (C @ expm(A * tau) @ Pcov @ C.T).item()
        R_theory = sigma**2 * np.exp(-xi / L) * (1.0 if ch == "u" else 1 - xi / (2 * L))
        assert R_model == pytest.approx(R_theory, abs=1e-9)


def test_discretisation_is_exact_and_independent_of_Ts():
    A, B, C = wind.dryden_state_space(VA0, 50.0, 0.7, "w")
    Pcov, _ = wind.dryden_output_covariance(A, B, C)
    for Ts in (0.002, 0.02, 0.2):
        Phi, Qd = wind.discretize_noise_driven(A, B, Ts)
        np.testing.assert_allclose(Phi @ Pcov @ Phi.T + Qd, Pcov, rtol=1e-8, atol=1e-12)  # stationary


def test_simulated_gust_statistics():
    rng = np.random.default_rng(42)
    g = wind.DrydenGust(Ts=0.1, Va0=VA0, rng=rng)
    N = 200_000
    X = np.array([g.step() for _ in range(N)])
    assert abs(X.mean(axis=0)).max() < 0.15
    for i, (ch, L, sigma) in enumerate(CHANNELS):
        assert X[:, i].var() == pytest.approx(sigma**2, rel=0.15), ch
    # u-channel correlation at lag L/Va is exp(-1)
    lag = int(round(P.dryden_Lu / VA0 / 0.1))
    x = X[:, 0]
    rho = np.mean(x[:-lag] * x[lag:]) / x.var()
    assert rho == pytest.approx(np.exp(-1), abs=0.08)
    # the three channels are independent
    assert abs(np.corrcoef(X.T)[0, 1]) < 0.05 and abs(np.corrcoef(X.T)[1, 2]) < 0.05


def test_gust_reproducible_with_seed():
    a = [wind.DrydenGust(0.1, VA0, np.random.default_rng(7)).step() for _ in range(3)]
    b = [wind.DrydenGust(0.1, VA0, np.random.default_rng(7)).step() for _ in range(3)]
    np.testing.assert_array_equal(a, b)
