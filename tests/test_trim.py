"""Phase 3: trim solutions and simulation from trim."""

import numpy as np
import pytest

from params import aerosonde as P
from sim import actuators, aircraft, dynamics as dyn, rotations as rot, trim

CASES = {
    "level": dict(Va=25.0),
    "slow": dict(Va=18.0),
    "fast": dict(Va=30.0),
    "climb5": dict(Va=25.0, gamma=np.deg2rad(5)),
    "descent5": dict(Va=25.0, gamma=np.deg2rad(-5)),
    "right_turn": dict(Va=25.0, R=150.0),
    "left_turn": dict(Va=25.0, R=-150.0),
    "climbing_turn": dict(Va=22.0, gamma=np.deg2rad(3), R=200.0),
}


@pytest.fixture(scope="module")
def trims():
    return {k: trim.compute_trim(**kw) for k, kw in CASES.items()}


@pytest.mark.parametrize("name", CASES)
def test_residual_is_zero_and_inputs_not_saturated(trims, name):
    t = trims[name]
    assert t.cost < 1e-10
    assert np.all(t.delta[:3] > actuators.LIMITS_LO[:3] + 1e-3) and np.all(t.delta[:3] < actuators.LIMITS_HI[:3] - 1e-3)
    assert 0.0 < t.delta[3] < 1.0
    assert np.linalg.norm(t.x13[6:10]) == pytest.approx(1.0)


@pytest.mark.parametrize("name", CASES)
def test_trim_verified_with_independent_euler_model(trims, name):
    """Re-check with the 12-state EULER dynamics (different attitude parameterisation than the solver):
    all body accelerations and Euler rates zero, psi_dot = Va cos(gamma)/R, h_dot = Va sin(gamma),
    ground track is a circle of radius R at speed Va cos(gamma)."""
    t = trims[name]
    f_b, m_b, (Va, alpha, beta) = aircraft.forces_moments(t.x13, t.delta)
    assert (Va, beta) == pytest.approx((t.Va, t.beta), abs=1e-12)
    xd = dyn.f_euler(t.x12, f_b, m_b)
    np.testing.assert_allclose(xd[3:6], 0, atol=1e-9)         # u_dot, v_dot, w_dot
    np.testing.assert_allclose(xd[9:12], 0, atol=1e-9)        # p_dot, q_dot, r_dot
    np.testing.assert_allclose(xd[6:8], 0, atol=1e-9)         # phi_dot, theta_dot
    psi_dot = 0.0 if not np.isfinite(t.R) else t.Va * np.cos(t.gamma) / t.R
    assert xd[8] == pytest.approx(psi_dot, abs=1e-9)
    assert -xd[2] == pytest.approx(t.Va * np.sin(t.gamma), abs=1e-9)
    assert np.hypot(xd[0], xd[1]) == pytest.approx(t.Va * np.cos(t.gamma), abs=1e-9)


@pytest.mark.parametrize("name", CASES)
def test_simulation_from_trim_stays_at_trim(trims, name):
    """Full nonlinear plant (RK4 dt=0.01 s, actuators at trim): body velocity, rates, attitude, V_a, alpha, beta constant."""
    t = trims[name]
    a = aircraft.Aircraft(t.x13, t.delta)
    for _ in range(1000):  # 10 s
        a.step(0.01, t.delta)
    np.testing.assert_allclose(a.x[3:6], t.x13[3:6], atol=1e-8)
    np.testing.assert_allclose(a.x[10:13], t.x13[10:13], atol=1e-8)
    np.testing.assert_allclose(a.x[14:17], t.delta[1:], atol=1e-12)
    Va, alpha, beta = a.air_data()
    assert (Va, alpha, beta) == pytest.approx((t.Va, t.alpha, t.beta), abs=1e-8)
    # position evolves as predicted: climb rate and (for turns) a circle of radius |R|
    T = a.t
    assert -(a.x[2] - t.x13[2]) == pytest.approx(t.Va * np.sin(t.gamma) * T, abs=1e-6)
    if np.isfinite(t.R):
        # circle centre lies |R| to the right (left for R<0) of the initial GROUND COURSE (not the heading:
        # with bank, w != 0 makes course differ from heading by a fraction of a degree)
        v_ned = rot.quat_to_R(t.x13[6:10]).T @ t.x13[3:6]
        chi = np.arctan2(v_ned[1], v_ned[0])
        centre = t.x13[0:2] + t.R * np.array([-np.sin(chi), np.cos(chi)])
        assert np.linalg.norm(a.x[0:2] - centre) == pytest.approx(abs(t.R), rel=1e-7)


def test_level_flight_is_nearly_wings_level_and_coordinated(trims):
    t = trims["level"]
    assert abs(np.rad2deg(t.phi)) < 0.2            # tiny bank offsets prop-torque/rudder side force
    assert abs(t.delta[0]) < 0.05 and abs(t.delta[2]) < 0.02
    assert t.theta == pytest.approx(t.alpha, abs=1e-4)   # gamma = theta - alpha when phi ~ 0, beta = 0


@pytest.mark.parametrize("name", ["climb5", "descent5"])
def test_flight_path_angle_geometry(trims, name):
    t = trims[name]
    assert t.theta - t.alpha == pytest.approx(t.gamma, abs=2e-4)


def test_coordinated_turn_bank_angle_formula(trims):
    """tan(phi) ~ V_a psi_dot / g for a level coordinated turn (side-force terms are small)."""
    t = trims["right_turn"]
    psi_dot = t.Va / t.R
    assert np.tan(t.phi) == pytest.approx(t.Va * psi_dot / P.gravity, rel=0.04)
    assert trims["left_turn"].phi < 0 < t.phi


def test_trends_with_airspeed_climb_and_turn(trims):
    lv, slow, fast = trims["level"], trims["slow"], trims["fast"]
    assert slow.alpha > lv.alpha > fast.alpha                  # slower -> larger alpha
    assert slow.delta[1] < lv.delta[1] < fast.delta[1]         # slower -> more nose-up (negative) elevator
    assert trims["climb5"].delta[3] > lv.delta[3] > trims["descent5"].delta[3]   # climb needs more throttle
    assert trims["right_turn"].delta[3] > lv.delta[3]          # extra induced drag in the turn
    assert trims["right_turn"].alpha > lv.alpha                # load factor 1/cos(phi) -> more alpha


def test_left_right_turns_have_mirrored_longitudinal_trim(trims):
    r, l = trims["right_turn"], trims["left_turn"]
    assert r.phi == pytest.approx(-l.phi, abs=np.deg2rad(0.3))
    assert r.alpha == pytest.approx(l.alpha, abs=2e-3)
    assert r.delta[1] == pytest.approx(l.delta[1], abs=3e-3)
    assert r.delta[3] == pytest.approx(l.delta[3], abs=3e-3)


def test_infeasible_trim_raises():
    with pytest.raises(RuntimeError):
        trim.compute_trim(Va=60.0)          # more speed than the 11 kg airframe can hold with the available thrust
    with pytest.raises(RuntimeError):
        trim.compute_trim(Va=25.0, R=10.0)  # 10 m radius: load factor far beyond the lift available
