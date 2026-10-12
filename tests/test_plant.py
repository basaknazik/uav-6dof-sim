"""Phase 2: full plant (rigid body + aero + propulsion + actuators) consistency tests."""

import numpy as np
import pytest

from sim import aircraft as ac
from sim import integrators as integ
from sim import rotations as rot
from sim import wind

CMD = np.array([0.05, -0.04, 0.02, 0.6])


def _x0(v_air_b, wind_ned, angles=(0.1, 0.05, 0.7), pqr=(0.02, 0.0, 0.01)):
    R = rot.euler_to_R(*angles)
    v_ground_b = np.asarray(v_air_b) + R @ np.asarray(wind_ned)
    return np.concatenate([[0, 0, -150], v_ground_b, rot.euler_to_quat(*angles), pqr]), R


def test_galilean_invariance_steady_wind():
    """A uniform steady wind only translates the air mass: air-relative motion must be IDENTICAL.

    Same commands, same air-relative initial state; with wind the ground position drifts by w*t and
    ground velocity = air velocity + w, but attitude, rates, V_a, alpha, beta must not change."""
    w = np.array([6.0, -4.0, 1.5])
    v_air = np.array([24.0, 0.5, 1.5])
    x_a, _ = _x0(v_air, np.zeros(3))
    x_b, _ = _x0(v_air, w)
    A = ac.Aircraft(x_a, [0, 0, 0, 0.6])
    B = ac.Aircraft(x_b, [0, 0, 0, 0.6])
    dt, T = 0.01, 8.0
    for k in range(int(T / dt)):
        A.step(dt, CMD)
        B.step(dt, CMD, wind_ned=w)
    t = A.t
    np.testing.assert_allclose(B.x[0:3] - w * t, A.x[0:3], atol=1e-7)          # position relative to air mass
    np.testing.assert_allclose(B.x[6:], A.x[6:], atol=1e-9)                    # attitude, rates, actuators
    assert A.air_data() == pytest.approx(B.air_data(wind_ned=w), abs=1e-8)
    # ...and it is a non-trivial manoeuvre (attitude actually changed)
    assert np.abs(A.x[6:10] - x_a[6:10]).max() > 1e-3


def test_without_wind_correction_the_result_differs_negative_control():
    """If we forget to give the wind to the plant, the dynamics must change (so the test above has teeth)."""
    w = np.array([6.0, -4.0, 1.5])
    v_air = np.array([24.0, 0.5, 1.5])
    x_a, _ = _x0(v_air, np.zeros(3))
    x_b, _ = _x0(v_air, w)
    A, B = ac.Aircraft(x_a, [0, 0, 0, 0.6]), ac.Aircraft(x_b, [0, 0, 0, 0.6])
    for _ in range(400):
        A.step(0.01, CMD)
        B.step(0.01, CMD)  # wind ignored
    assert np.abs(B.x[6:10] - A.x[6:10]).max() > 1e-3 or abs(B.x[3] - A.x[3]) > 1.0


def test_gust_in_body_frame_equals_equivalent_ground_velocity_shift():
    """Gust g (body frame) is equivalent to body-frame air velocity v - g: same forces."""
    x13, _ = _x0([24.0, 0.5, 1.5], np.zeros(3))
    g = np.array([1.0, -2.0, 0.5])
    d = np.array([0.0, 0.0, 0.0, 0.5])
    f1, m1, ad1 = ac.forces_moments(x13, d, gust_b=g)
    x13b = x13.copy()
    x13b[3:6] -= g
    f0, m0, ad0 = ac.forces_moments(x13b, d)
    # aero + prop identical (gravity depends only on attitude, unchanged)
    np.testing.assert_allclose(f1, f0, atol=1e-12)
    np.testing.assert_allclose(m1, m0, atol=1e-12)
    assert ad1 == pytest.approx(ad0)


def test_aircraft_rk4_matches_scipy_reference():
    x_a, _ = _x0([24.0, 0.5, 1.5], np.zeros(3))
    a = ac.Aircraft(x_a, [0, 0, 0, 0.6])
    dt, T = 0.005, 4.0
    for _ in range(int(T / dt)):
        a.step(dt, CMD)
    f = lambda t, x: a.deriv(x, CMD, np.zeros(3), np.zeros(3))  # noqa: E731
    ref = integ.simulate_scipy(f, np.concatenate([x_a, [0, 0, 0, 0.6]]), np.array([0.0, T]), rtol=1e-11, atol=1e-11)
    assert np.abs(a.x - ref[-1]).max() < 5e-7   # (a.x is renormalised each step; reference is not)


def test_plant_runs_with_dryden_gusts_without_blowing_up():
    x_a, _ = _x0([25.0, 0.0, 1.0], np.zeros(3), angles=(0, 0.04, 0.0), pqr=(0, 0, 0))
    a = ac.Aircraft(x_a, [0, 0, 0, 0.6])
    g = wind.DrydenGust(0.01, 25.0, np.random.default_rng(0))
    for _ in range(2000):
        a.step(0.01, [0, 0, 0, 0.6], wind_ned=[3.0, 0, 0], gust_b=g.step())
    assert np.all(np.isfinite(a.x)) and np.linalg.norm(a.x[6:10]) == pytest.approx(1.0)
    assert 10 < a.air_data()[0] < 40
