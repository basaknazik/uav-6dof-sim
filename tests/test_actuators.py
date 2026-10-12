"""Phase 2: actuator saturation and first-order dynamics."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from params import aerosonde as P
from sim import actuators as act
from sim.aircraft import Aircraft
from sim import rotations as rot


def _simulate(cmd, T, delta0=np.zeros(4)):
    sol = solve_ivp(lambda t, x: act.actuator_rate(x, cmd), (0, T), delta0,
                    t_eval=np.linspace(0, T, 400), rtol=1e-11, atol=1e-12)
    return sol.t, sol.y.T


def test_first_order_step_response_matches_analytic():
    cmd = np.array([0.2, -0.3, 0.1, 0.6])
    t, X = _simulate(cmd, 1.0)
    expected = cmd * (1 - np.exp(-t[:, None] / act.TAU))
    np.testing.assert_allclose(X, expected, atol=1e-8)
    # 63.2 % at t = tau (property of first-order systems)
    t2, X2 = _simulate(cmd, act.TAU[0], )
    assert X2[-1, 0] == pytest.approx(0.632120558 * cmd[0], rel=1e-6)


def test_saturation_and_state_never_exceeds_limits():
    cmd = np.array([5.0, -5.0, 10.0, 3.0])               # far beyond limits
    t, X = _simulate(cmd, 3.0)
    assert np.all(X <= act.LIMITS_HI + 1e-9) and np.all(X >= act.LIMITS_LO - 1e-9)
    np.testing.assert_allclose(X[-1], act.saturate(cmd), atol=1e-6)   # settles at the clipped command
    assert act.saturate(np.array([0, 0, 0, -1.0]))[3] == 0.0


def test_limits_values():
    assert act.LIMITS_HI[0] == pytest.approx(np.deg2rad(45))
    assert act.LIMITS_HI[3] == 1.0 and act.LIMITS_LO[3] == 0.0


def test_aircraft_actuator_states_follow_commands():
    """Through the full plant: elevator state follows the first-order response (plant forces do not feed back into it)."""
    x13 = np.concatenate([[0, 0, -100, 25, 0, 0], rot.euler_to_quat(0, 0.05, 0), [0, 0, 0]])
    ac = Aircraft(x13, np.array([0, 0, 0, 0.6]))
    cmd = np.array([0, 0.1, 0, 0.6])
    for _ in range(100):
        ac.step(0.002, cmd)          # 0.2 s = 4 tau
    expected = 0.1 * (1 - np.exp(-0.2 / P.tau_servo))
    assert ac.x[14] == pytest.approx(expected, rel=1e-6)
    assert np.linalg.norm(ac.x[6:10]) == pytest.approx(1.0, abs=1e-12)
