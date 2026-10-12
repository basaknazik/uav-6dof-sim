"""Phase 2: motor/propeller model."""

import numpy as np
import pytest
from scipy.optimize import brentq

from params import aerosonde as P
from sim import aero, propulsion as pr


def _prop_torque_direct(w, Va):
    """Q_p from the definition (J = Va/(nD)), independent of the polynomial expansion in propulsion.py."""
    n = w / (2 * np.pi)
    J = Va / (n * P.D_prop)
    return P.rho * n**2 * P.D_prop**5 * (P.C_Q2 * J**2 + P.C_Q1 * J + P.C_Q0)


def _motor_torque(w, delta_t):
    return P.KQ * ((P.V_max * delta_t - P.KV * w) / P.R_motor - P.i0)


@pytest.mark.parametrize("Va", [5.0, 15.0, 25.0, 35.0])
@pytest.mark.parametrize("delta_t", [0.5, 0.7, 1.0])
def test_torque_balance_matches_independent_root_find(Va, delta_t):
    w = pr.prop_speed(Va, delta_t)
    w_ref = brentq(lambda x: _motor_torque(x, delta_t) - _prop_torque_direct(x, Va), 1.0, 2000.0, xtol=1e-12)
    assert w == pytest.approx(w_ref, rel=1e-9)
    assert _motor_torque(w, delta_t) == pytest.approx(_prop_torque_direct(w, Va), rel=1e-8)


def test_thrust_matches_definition():
    Va, dt = 20.0, 0.8
    w = pr.prop_speed(Va, dt)
    n = w / (2 * np.pi)
    J = Va / (n * P.D_prop)
    T_def = P.rho * n**2 * P.D_prop**4 * (P.C_T2 * J**2 + P.C_T1 * J + P.C_T0)
    T, Q = pr.thrust_torque(Va, dt)
    assert T == pytest.approx(T_def, rel=1e-10)
    assert Q == pytest.approx(_prop_torque_direct(w, Va), rel=1e-10)


def test_thrust_monotone_in_throttle_and_decreasing_in_airspeed():
    for Va in (15.0, 25.0):
        T = [pr.thrust_torque(Va, d)[0] for d in np.linspace(0.4, 1.0, 13)]
        assert all(b > a for a, b in zip(T, T[1:]))
    T = [pr.thrust_torque(Va, 0.8)[0] for Va in (0.0, 10.0, 20.0, 30.0)]
    assert all(b < a for a, b in zip(T, T[1:]))


def test_static_thrust_and_limits():
    T_full, Q_full = pr.thrust_torque(0.0, 1.0)
    assert 30 < T_full < 200 and Q_full > 0        # tens of N for an 11 kg UAV (thrust-to-weight order 1)
    assert pr.thrust_torque(0.0, 0.0) == (0.0, 0.0)  # no voltage, no airflow: nothing
    assert pr.thrust_torque(25.0, 1.5) == pr.thrust_torque(25.0, 1.0)   # throttle is clipped to [0, 1]
    assert pr.thrust_torque(25.0, -0.3) == pr.thrust_torque(25.0, 0.0)


def test_cruise_throttle_is_plausible_and_efficiency_below_one():
    """Level flight at 25 m/s needs T = D. The required throttle must exist in (0.3, 0.9) and the
    propulsive efficiency T*Va / (V_in * i) must lie in (0, 1)."""
    Va = 25.0
    qbarS = 0.5 * P.rho * Va**2 * P.S_wing
    alpha = brentq(lambda a: qbarS * aero.lift_coefficient(a) * np.cos(a) - P.mass * P.gravity, 0.0, 0.2)
    drag = qbarS * aero.drag_coefficient(alpha) / np.cos(alpha) * 1.0
    # exact force balance along the flight path in level, unaccelerated flight: T cos(a) = D_body_x term
    f, _ = aero.aero_forces_moments(Va, alpha, 0.0, np.zeros(3), np.zeros(3))
    T_needed = -f[0]
    dt = brentq(lambda d: pr.thrust_torque(Va, d)[0] - T_needed, 0.2, 1.0)
    assert 0.3 < dt < 0.9, dt
    w = pr.prop_speed(Va, dt)
    i = (P.V_max * dt - P.KV * w) / P.R_motor
    eta = T_needed * Va / (P.V_max * dt * i)
    assert 0.2 < eta < 1.0, eta
    assert 0 < i < 120.0     # battery current plausible
    assert drag > 0
