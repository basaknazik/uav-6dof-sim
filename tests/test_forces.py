"""Phase 2: air data, aerodynamic forces/moments (symmetry, signs, lift model)."""

import numpy as np
import pytest

from params import aerosonde as P
from sim import aero, wind
from sim import rotations as rot

RNG = np.random.default_rng(3)
Z3 = np.zeros(3)


# ------------------------------------------------------------------ air data
def test_airdata_no_wind():
    v = np.array([22.0, 3.0, 4.0])
    Va, alpha, beta, v_air = wind.airdata(v, np.eye(3))
    assert Va == pytest.approx(np.linalg.norm(v))
    assert alpha == pytest.approx(np.arctan2(4.0, 22.0))
    assert beta == pytest.approx(np.arcsin(3.0 / np.linalg.norm(v)))
    np.testing.assert_allclose(v_air, v)


@pytest.mark.parametrize("wind_ned, Va_expected", [([5.0, 0, 0], 20.0),     # tailwind
                                                   ([-5.0, 0, 0], 30.0)])   # headwind
def test_head_and_tail_wind_change_airspeed_not_groundspeed(wind_ned, Va_expected):
    Va, alpha, beta, _ = wind.airdata([25.0, 0, 0], np.eye(3), wind_ned)
    assert Va == pytest.approx(Va_expected)
    assert alpha == pytest.approx(0.0) and beta == pytest.approx(0.0)


def test_crosswind_sign_of_beta_and_crab_angle():
    # Wind blowing toward East (+5). Heading North, ground speed 25: relative wind comes from the LEFT -> beta<0.
    Va, _, beta, v_air = wind.airdata([25.0, 0, 0], np.eye(3), [0, 5.0, 0])
    assert beta < 0 and v_air[1] == pytest.approx(-5.0)
    # Crab: to track due North with Va=25 and 5 m/s wind toward East, point into the wind (heading west of north).
    psi = -np.arcsin(5.0 / 25.0)
    Vg = np.sqrt(25.0**2 - 5.0**2)
    R = rot.euler_to_R(0, 0, psi)
    v_b = R @ np.array([Vg, 0, 0])  # ground velocity (due North) in body axes
    Va, alpha, beta, _ = wind.airdata(v_b, R, [0, 5.0, 0])
    assert Va == pytest.approx(25.0) and beta == pytest.approx(0.0, abs=1e-12) and alpha == pytest.approx(0.0, abs=1e-12)


def test_airspeed_equals_inertial_computation_for_random_states():
    for _ in range(50):
        ang = (RNG.uniform(-1, 1), RNG.uniform(-1, 1), RNG.uniform(-3, 3))
        R = rot.euler_to_R(*ang)
        v_b, w = RNG.normal(size=3) * 10 + [20, 0, 0], RNG.normal(size=3) * 4
        Va, *_ = wind.airdata(v_b, R, w)
        assert Va == pytest.approx(np.linalg.norm(R.T @ v_b - w))  # compare in the inertial frame


def test_gust_is_added_in_body_frame():
    g = np.array([1.0, -2.0, 3.0])
    R = rot.euler_to_R(0.3, 0.2, 1.0)
    _, _, _, v_air = wind.airdata([20.0, 0, 0], R, Z3, g)
    np.testing.assert_allclose(v_air, [20.0, 0, 0] - g)


def test_alpha_positive_when_sinking_relative_to_air():
    _, alpha, _, _ = wind.airdata([20.0, 0, 3.0], np.eye(3))   # w>0 is down
    assert alpha > 0
    _, alpha_up, _, _ = wind.airdata([20.0, 0, 0.0], np.eye(3), [0, 0, -3.0])  # air moving UP (wind d<0)
    assert alpha_up > 0  # same relative flow: aircraft sinks relative to the air


# ------------------------------------------------------------------ lift model
def test_stall_blend_properties():
    assert aero.stall_blend(P.alpha0) == pytest.approx(0.5, abs=1e-12)
    assert aero.stall_blend(-P.alpha0) == pytest.approx(0.5, abs=1e-12)
    a = np.linspace(-1.5, 1.5, 301)
    np.testing.assert_allclose(aero.stall_blend(a), aero.stall_blend(-a), atol=1e-12)
    assert aero.stall_blend(0.0) < 1e-6 and aero.stall_blend(1.2) > 1 - 1e-6


def test_lift_coefficient_regimes_and_continuity():
    for a in (-0.1, 0.0, 0.1):
        assert aero.lift_coefficient(a) == pytest.approx(P.C_L_0 + P.C_L_alpha * a, rel=1e-5)
    # sigma(0.25) ~ e^{-M (alpha0-0.25)} ~ 2e-5: still linear to 1e-3
    assert aero.lift_coefficient(0.25) == pytest.approx(P.C_L_0 + P.C_L_alpha * 0.25, rel=1e-3)
    for a in (0.9, 1.2, 1.5):
        assert aero.lift_coefficient(a) == pytest.approx(2 * np.sin(a) ** 2 * np.cos(a), rel=1e-5)
    grid = np.arange(-1.5, 1.5, 1e-4)
    CL = np.array([aero.lift_coefficient(a) for a in grid])
    assert np.abs(np.diff(CL)).max() < 0.02          # no jumps
    # stall: lift peaks near alpha0 then falls; post-stall max is far below the linear extrapolation
    assert CL.max() < P.C_L_0 + P.C_L_alpha * P.alpha0 + 1e-9
    assert aero.lift_coefficient(0.8) < 0.7 * (P.C_L_0 + P.C_L_alpha * 0.8)


# ------------------------------------------------------------------ forces / moments
def test_stability_axis_rotation_preserves_lift_drag_magnitude_and_alpha_zero_case():
    Va = 25.0
    qbarS = 0.5 * P.rho * Va**2 * P.S_wing
    f, _ = aero.aero_forces_moments(Va, 0.0, 0.0, Z3, Z3)
    assert f[0] == pytest.approx(-qbarS * aero.drag_coefficient(0.0))   # alpha=0: fx=-D
    assert f[2] == pytest.approx(-qbarS * aero.lift_coefficient(0.0))   # fz=-L
    for a in np.linspace(-0.4, 0.6, 11):
        f, _ = aero.aero_forces_moments(Va, a, 0.0, Z3, Z3)
        L, D = qbarS * aero.lift_coefficient(a), qbarS * aero.drag_coefficient(a)
        assert np.hypot(f[0], f[2]) == pytest.approx(np.hypot(L, D))
        # lift is perpendicular to the relative wind: wind direction in body is [cos a, 0, sin a]
        wind_dir = np.array([np.cos(a), 0, np.sin(a)])
        assert f @ wind_dir == pytest.approx(-D)         # force component along relative wind = -drag
        assert f @ np.array([-np.sin(a), 0, np.cos(a)]) == pytest.approx(-L)


def test_zero_lateral_inputs_give_zero_lateral_aero():
    f, m = aero.aero_forces_moments(25.0, 0.05, 0.0, [0.0, 0.3, 0.0], [0.0, 0.1, 0.0])
    assert f[1] == 0 and m[0] == 0 and m[2] == 0


def test_left_right_mirror_symmetry():
    for _ in range(100):
        Va, a, be = RNG.uniform(10, 40), RNG.uniform(-0.3, 0.5), RNG.uniform(-0.3, 0.3)
        p, q, r = RNG.normal(size=3)
        da, de, dr = RNG.uniform(-0.5, 0.5, 3)
        f1, m1 = aero.aero_forces_moments(Va, a, be, [p, q, r], [da, de, dr])
        f2, m2 = aero.aero_forces_moments(Va, a, -be, [-p, q, -r], [-da, de, -dr])
        np.testing.assert_allclose(f2 * [1, -1, 1], f1, atol=1e-10)   # fx, fz even; fy odd
        np.testing.assert_allclose(m2 * [-1, 1, -1], m1, atol=1e-10)  # m even; l, n odd


def _delta(fm_fn, base, key, h):
    f0, m0 = fm_fn(**base)
    f1, m1 = fm_fn(**{**base, key: h})
    return f1 - f0, m1 - m0


def _fm(Va=25.0, alpha=0.05, beta=0.0, pqr=(0, 0, 0), delta=(0, 0, 0)):
    return aero.aero_forces_moments(Va, alpha, beta, np.array(pqr, float), np.array(delta, float))


def test_control_surface_signs():
    df, dm = _delta(_fm, {}, "delta", (0, 0.1, 0))      # positive elevator
    assert dm[1] < 0 and df[2] < 0                        # nose-down moment; lift up (fz more negative)
    df, dm = _delta(_fm, {}, "delta", (0.1, 0, 0))      # positive aileron
    assert dm[0] > 0                                      # positive roll (right wing down)
    df, dm = _delta(_fm, {}, "delta", (0, 0, 0.1))      # positive rudder
    assert dm[2] < 0 and df[1] > 0                        # nose-left yaw moment per C_n_delta_r<0; C_Y_delta_r>0


def test_static_stability_and_damping_signs():
    _, dm = _delta(_fm, {}, "alpha", 0.15)               # nose-up disturbance -> nose-down moment
    assert dm[1] < 0
    df, dm = _delta(_fm, {}, "beta", 0.1)                # wind from the right
    assert dm[2] > 0 and df[1] < 0 and dm[0] < 0         # weathercock; side force opposes; dihedral roll away
    _, dm = _delta(_fm, {}, "pqr", (1.0, 0, 0))
    assert dm[0] < 0                                      # roll damping
    _, dm = _delta(_fm, {}, "pqr", (0, 1.0, 0))
    assert dm[1] < 0                                      # pitch damping
    _, dm = _delta(_fm, {}, "pqr", (0, 0, 1.0))
    assert dm[2] < 0                                      # yaw damping


def test_forces_scale_with_dynamic_pressure():
    f1, m1 = _fm(Va=20.0, beta=0.1, delta=(0.1, 0.1, 0.1))
    f2, m2 = _fm(Va=40.0, beta=0.1, delta=(0.1, 0.1, 0.1))
    np.testing.assert_allclose(f2, 4 * f1, rtol=1e-10)
    np.testing.assert_allclose(m2, 4 * m1, rtol=1e-10)


def test_cruise_lift_is_consistent_with_weight():
    """At 25 m/s and CL needed for L = W, alpha is a small positive angle (sanity of C_L0, C_La, S, rho)."""
    qbarS = 0.5 * P.rho * 25.0**2 * P.S_wing
    CL_needed = P.mass * P.gravity / qbarS
    assert 0.3 < CL_needed < 0.8
    alpha = (CL_needed - P.C_L_0) / P.C_L_alpha
    assert 0.0 < alpha < np.deg2rad(8)
