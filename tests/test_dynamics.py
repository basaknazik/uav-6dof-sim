"""Angular-dynamics algebra tests. None depend on any book value or on specific parameters."""

import numpy as np

from sim import dynamics as dyn

RNG = np.random.default_rng(2)


def random_body():
    """Random physically valid inertia: Jxz small enough that J is positive definite."""
    Jx, Jy, Jz = RNG.uniform(0.5, 3.0, 3)
    Jxz = RNG.uniform(-0.9, 0.9) * np.sqrt(Jx * Jz)
    return dyn.RigidBody(mass=RNG.uniform(1, 20), Jx=Jx, Jy=Jy, Jz=Jz, Jxz=Jxz)


def test_gamma_equations_match_numerical_inverse_random_inertia():
    """Gamma1..8 form == J^-1 (m - w x Jw) to machine precision, for random J, w, m."""
    for _ in range(200):
        body = random_body()
        omega, moment = RNG.normal(scale=3.0, size=3), RNG.normal(scale=5.0, size=3)
        a = dyn.omega_dot_gamma(omega, moment, body)
        b = dyn.omega_dot_numeric(omega, moment, body)
        np.testing.assert_allclose(a, b, rtol=1e-11, atol=1e-12 * max(1.0, np.abs(b).max()))


def test_gamma_equations_match_numerical_inverse_aerosonde():
    for _ in range(200):
        omega, moment = RNG.normal(scale=3.0, size=3), RNG.normal(scale=5.0, size=3)
        np.testing.assert_allclose(dyn.omega_dot_gamma(omega, moment, dyn.AEROSONDE),
                                   dyn.omega_dot_numeric(omega, moment, dyn.AEROSONDE),
                                   rtol=1e-11, atol=1e-12)


def test_gamma_equations_detect_wrong_coefficient():
    """Negative control: corrupting one coefficient must break the match."""
    body = dyn.AEROSONDE
    omega, moment = np.array([1.0, 2.0, -1.0]), np.array([0.3, -0.2, 0.1])
    good = dyn.omega_dot_numeric(omega, moment, body)
    p, q, r = omega
    g = body.gammas
    bad_pdot = g.g1 * p * q + g.g2 * q * r + g.g3 * moment[0] + g.g4 * moment[2]  # sign flip on Gamma2 term
    assert not np.isclose(bad_pdot, good[0], rtol=1e-6)


def test_inertia_matrix_symmetric_positive_definite():
    J = dyn.AEROSONDE.J
    np.testing.assert_allclose(J, J.T)
    assert np.all(np.linalg.eigvalsh(J) > 0)
