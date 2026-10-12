"""Phase 5: autopilot design (from linear models), linear closed-loop analysis, nonlinear mission tests."""

import dataclasses

import numpy as np
import pytest
from scipy import signal

from params import aerosonde as P
from sim import actuators, mission, propulsion, trim
from sim import autopilot as ap
from sim import dynamics as dyn
from sim import linearize as lin
from sim import rotations as rot

V = 25.0
G = dyn.AEROSONDE.gammas
D = np.deg2rad
DT = 0.02          # controller/plant step used in the mission tests (the design is rate independent)


@pytest.fixture(scope="module")
def setup():
    t = trim.compute_trim(V)
    m = lin.linearize(t)
    g, dm = ap.design_gains(m)
    return t, m, g, dm


# ------------------------------------------------------------------ design models and gains
def test_design_model_coefficients_match_analytic_expressions(setup):
    t, m, g, dm = setup
    qS = 0.5 * P.rho * V**2 * P.S_wing
    assert dm.a_phi2 == pytest.approx(qS * P.b * (G.g3 * P.C_ell_delta_a + G.g4 * P.C_n_delta_a), rel=1e-6)
    assert dm.a_theta3 == pytest.approx(qS * P.c * P.C_m_delta_e / P.Jy, rel=1e-6)
    assert dm.a_theta2 == pytest.approx(-qS * P.c * P.C_m_alpha / P.Jy, rel=0.02)   # alpha ~ w/V: small CL/CD cross terms
    assert dm.a_V3 == pytest.approx(P.gravity * np.cos(t.gamma), rel=1e-6)           # Va_dot = -g sin(theta - alpha) ...
    h = 1e-6
    dT = (propulsion.thrust_torque(V, t.delta[3] + h)[0] - propulsion.thrust_torque(V, t.delta[3] - h)[0]) / (2 * h)
    assert dm.a_V2 == pytest.approx(np.cos(t.alpha) * dT / P.mass, rel=1e-5)
    assert dm.a_V1 > 0       # drag grows with speed at fixed throttle: speed-stable


def test_gain_signs_follow_plant_gain_signs(setup):
    _, _, g, dm = setup
    assert dm.a_phi2 > 0 and g.kp_phi > 0 and g.kd_phi >= 0
    assert dm.a_theta3 < 0 and g.kp_theta < 0                 # positive elevator is nose-DOWN -> negative gain
    assert g.kp_chi > 0 and g.ki_chi > 0 and g.K_theta_DC > 0
    assert g.kp_h > 0 and g.ki_h > 0 and g.kp_V > 0 and g.ki_V > 0
    assert g.kp_V2 < 0 and g.ki_V2 < 0                         # pitch down to gain speed
    assert g.kr_yaw > 0                                        # positive rudder yaws nose left (C_n_delta_r < 0)


def test_design_equations_reproduce_the_targets_on_the_design_models(setup):
    """Algebra check: on the 2nd-order design models, the gains give exactly the specified wn and zeta."""
    _, _, g, dm = setup
    s = g.spec
    # roll: s^2 + (a1 + a2 kd) s + a2 kp
    assert np.sqrt(dm.a_phi2 * g.kp_phi) == pytest.approx(s.wn_phi)
    # pitch: s^2 + (a1 + a3 kd) s + (a2 + a3 kp)
    wn = np.sqrt(dm.a_theta2 + dm.a_theta3 * g.kp_theta)
    zeta = (dm.a_theta1 + dm.a_theta3 * g.kd_theta) / (2 * wn)
    assert (wn, zeta) == pytest.approx((s.wn_theta, s.zeta_theta))
    # throttle speed loop: s^2 + (a1 + a2 kp) s + a2 ki
    assert np.sqrt(dm.a_V2 * g.ki_V) == pytest.approx(s.wn_V)
    assert (dm.a_V1 + dm.a_V2 * g.kp_V) / (2 * s.wn_V) == pytest.approx(s.zeta_V)
    # bandwidth separation between nested loops
    assert s.wn_phi / np.sqrt(g.ki_chi * P.gravity / V) == pytest.approx(s.W_chi)
    assert s.wn_theta / np.sqrt(g.ki_h * g.K_theta_DC * V) == pytest.approx(s.W_h)


# ------------------------------------------------------------------ linear closed-loop analysis
def test_closed_loop_is_stable_and_spiral_is_stabilised(setup):
    t, m, g, _ = setup
    Acl, _, _ = ap.closed_loop_hold_zone(m, g)
    ev = np.linalg.eigvals(Acl)
    assert lin.modal_summary(m)["spiral"].eigenvalue.real > 0       # open loop: unstable spiral
    assert ev.real.max() < -0.1                                      # closed loop: everything stable
    cplx = ev[np.abs(ev.imag) > 1e-6]
    assert (-cplx.real / np.abs(cplx)).min() > 0.25


def test_yaw_damper_improves_dutch_roll_damping(setup):
    t, m, g, _ = setup
    def min_zeta(gains):
        ev = np.linalg.eigvals(ap.closed_loop_hold_zone(m, gains)[0])
        c = ev[np.abs(ev.imag) > 1e-6]
        return (-c.real / np.abs(c)).min()
    assert min_zeta(g) > min_zeta(dataclasses.replace(g, kr_yaw=0.0)) + 0.15
    assert g.spec.zeta_dr - 0.1 < 0.7                                # spec itself is sane


def test_fixed_gains_remain_stable_over_the_airspeed_range(setup):
    _, _, g, _ = setup
    for Va in (18.0, 22.0, 28.0, 30.0):
        mm = lin.linearize(trim.compute_trim(Va))
        assert np.linalg.eigvals(ap.closed_loop_hold_zone(mm, g)[0]).real.max() < -0.05, Va


def test_linear_closed_loop_step_specifications(setup):
    _, m, g, _ = setup
    Acl, Bcl, C = ap.closed_loop_hold_zone(m, g)
    t = np.linspace(0, 60, 3001)
    for ref, amp, row, os_max, t90_max in ((0, D(5), 0, 0.30, 2.0), (1, 5.0, 1, 0.35, 3.0), (2, 1.0, 2, 0.15, 3.0)):
        sys = signal.StateSpace(Acl, Bcl[:, [ref]] * amp, C, np.zeros((10, 1)))
        _, y = signal.step(sys, T=t)
        r = y[:, row]
        assert r[-1] == pytest.approx(amp, rel=1e-3)                  # integral action: zero steady-state error
        assert r.max() / amp - 1 < os_max
        assert t[np.argmax(r > 0.9 * amp)] < t90_max


@pytest.mark.parametrize("ref, amp, name", [(0, D(2), "chi"), (1, 2.0, "h"), (2, 0.5, "Va")])
def test_nonlinear_closed_loop_matches_linear_closed_loop_for_small_steps(setup, ref, amp, name):
    """Validates the design analysis AND the controller implementation: same small step, linear model vs. full plant."""
    t, m, g, _ = setup
    Acl, Bcl, C = ap.closed_loop_hold_zone(m, g)
    cmd = {"chi": lambda a: (0, V, 100.0, a), "h": lambda a: (0, V, 100.0 + a, 0.0), "Va": lambda a: (0, V + a, 100.0, 0.0)}[name]
    o = mission.simulate_mission(t, g, [cmd(amp)], 20.0, dt=0.01)
    sys = signal.StateSpace(Acl, Bcl[:, [ref]] * amp, C, np.zeros((10, 1)))
    _, y, _ = signal.lsim(sys, np.ones_like(o["t"]), o["t"])
    nl = {"chi": o["chi"] - o["chi"][0], "h": o["h"] - 100.0, "Va": o["Va"] - V}[name]
    assert np.abs(nl - y[:, {"chi": 0, "h": 1, "Va": 2}[name]]).max() / amp < 0.03


# ------------------------------------------------------------------ building blocks
def test_pi_anti_windup_limits_integrator_and_recovers_quickly():
    pi, naive = ap.PI(1.0, 2.0), 0.0
    for _ in range(1000):                                  # 10 s of a large error with the output stuck at its limit
        pi.update(5.0, 0.01, -1.0, 1.0)
        naive += 2.0 * 5.0 * 0.01
    assert abs(pi.I) < 1.0 + 5.0 + 1e-9                    # bounded by (limit + kp*error)
    assert naive > 90
    steps = 0
    while pi.update(-1.0, 0.01, -1.0, 1.0) > -0.99 and steps < 1000:   # error reverses: how long to leave the limit
        steps += 1
    assert steps < 150                                      # naive integrator would need > 40 s (4000 steps)


def test_pi_bumpless_reset_sets_requested_output():
    pi = ap.PI(2.0, 1.0)
    pi.reset(output=0.3, error=0.1, ff=0.05)
    assert pi.update(0.1, 0.0, -10, 10, ff=0.05) == pytest.approx(0.3)


def test_wrap_pi():
    assert ap.wrap_pi(D(190)) == pytest.approx(D(-170))
    assert ap.wrap_pi(D(-190)) == pytest.approx(D(170))
    assert ap.wrap_pi(0.3) == pytest.approx(0.3)
    assert ap.wrap_pi(D(3 * 180)) == pytest.approx(np.pi) or ap.wrap_pi(D(3 * 180)) == pytest.approx(-np.pi)


# ------------------------------------------------------------------ nonlinear mission tests
CMDS = [(0, 25.0, 100.0, 0.0), (5, 25.0, 150.0, D(90)), (60, 28.0, 150.0, D(90)),
        (90, 28.0, 120.0, D(45)), (130, 25.0, 120.0, D(45))]
SEGMENT_ENDS = [(55, 25.0, 150.0, 90.0), (88, 28.0, 150.0, 90.0), (125, 28.0, 120.0, 45.0), (169.9, 25.0, 120.0, 45.0)]


@pytest.fixture(scope="module")
def nominal(setup):
    t, _, g, _ = setup
    return mission.simulate_mission(t, g, CMDS, 170.0, dt=DT)


@pytest.fixture(scope="module")
def windy(setup):
    t, _, g, _ = setup
    return mission.simulate_mission(t, g, CMDS, 170.0, dt=DT, wind_ned=(-4.0, 6.0, 0.0), gusts=True, seed=1)


def _at(o, tq):
    return np.searchsorted(o["t"], tq)


def test_mission_tracks_commands_without_wind(nominal):
    o = nominal
    for tq, Va_c, h_c, chi_c in SEGMENT_ENDS:
        i = _at(o, tq)
        assert abs(o["h"][i] - h_c) < 0.3, (tq, o["h"][i])
        assert abs(np.rad2deg(o["chi"][i]) - chi_c) < 0.3
        assert abs(o["Va"][i] - Va_c) < 0.2
        assert abs(np.rad2deg(o["beta"][i])) < 0.5               # coordinated
    assert np.all(np.isfinite(o["x"]))


def test_mission_zones_and_actuator_limits(nominal):
    o = nominal
    zones = list(dict.fromkeys(o["zone"]))
    assert zones[0] == "hold" and "climb" in zones and "descend" in zones
    n_switch = int(np.sum(o["zone"][1:] != o["zone"][:-1]))
    assert n_switch <= 6                                          # no zone chattering
    x = o["x"]
    assert np.all(x[:, 13:] <= actuators.LIMITS_HI + 1e-9) and np.all(x[:, 13:] >= actuators.LIMITS_LO - 1e-9)
    assert np.rad2deg(np.abs(o["phi"]).max()) < 40
    assert np.rad2deg(np.abs(o["theta"]).max()) < 32
    in_climb = o["zone"] == "climb"
    assert np.all(o["delta_cmd"][in_climb, 3] == 1.0)
    assert np.all(o["delta_cmd"][o["zone"] == "descend", 3] == 0.0)


def test_climb_performance_and_speed_hold_during_climb(nominal):
    o = nominal
    i = np.where(o["zone"] == "climb")[0]
    roc = (o["h"][i[-1]] - o["h"][i[0]]) / (o["t"][i[-1]] - o["t"][i[0]])
    assert roc > 4.0                                              # m/s at full throttle
    assert np.abs(o["Va"][i] - 25.0).max() < 3.0                  # speed held with the elevator while climbing


def test_mission_with_wind_and_gusts_tracks_within_tolerance(windy):
    o = windy
    for tq, Va_c, h_c, chi_c in SEGMENT_ENDS:
        i = _at(o, tq)
        j = slice(max(i - 150, 0), i + 1)                          # average the last 3 s (gusts)
        assert abs(o["h"][j].mean() - h_c) < 1.5, tq
        assert abs(np.rad2deg(o["chi"][j].mean()) - chi_c) < 1.5, tq
        assert abs(o["Va"][j].mean() - Va_c) < 1.0, tq
    assert np.all(np.isfinite(o["x"]))
    assert np.rad2deg(np.abs(o["beta"]).max()) < 20
    x = o["x"]
    assert np.all(x[:, 13:] <= actuators.LIMITS_HI + 1e-9) and np.all(x[:, 13:] >= actuators.LIMITS_LO - 1e-9)


def test_course_is_tracked_not_heading_crab_angle_matches_wind_triangle(setup):
    """Steady 8 m/s crosswind, no gusts: the aircraft must hold GROUND course; heading differs by the crab angle."""
    t, _, g, _ = setup
    w = np.array([0.0, 8.0, 0.0])                                  # wind toward East
    o = mission.simulate_mission(t, g, [(0, 25.0, 100.0, 0.0)], 80.0, dt=DT, wind_ned=w)
    assert abs(np.rad2deg(o["chi"][-1])) < 0.2                     # tracking due North
    # wind triangle: Vg e_chi = Va e_psi + w  with chi = 0:  0 = Va sin(psi) + w_E  ->  psi = -asin(w_E / Va)
    psi_expected = -np.arcsin(w[1] / 25.0)
    assert o["psi"][-1] == pytest.approx(psi_expected, abs=D(0.3))
    assert abs(np.rad2deg(o["psi"][-1])) > 15                      # a large crab angle: heading hold would have failed
    assert abs(o["Va"][-1] - 25.0) < 0.2 and abs(o["h"][-1] - 100.0) < 0.3


def test_course_command_across_plus_minus_180_takes_the_short_way(setup):
    t, _, g, _ = setup
    x0 = t.x13.copy()
    phi, theta, _ = rot.quat_to_euler(x0[6:10])
    x0[6:10] = rot.euler_to_quat(phi, theta, D(170.0))
    o = mission.simulate_mission(t, g, [(0, 25.0, 100.0, D(-170.0))], 25.0, dt=DT, x0=x0)
    chi_deg = np.rad2deg(o["chi"])
    assert np.abs(ap.wrap_pi(o["chi"][-1] - D(-170))) < D(1.0)
    assert np.all(np.abs(chi_deg) > 150)                           # never swung through the long way round (through 0)
    assert o["phi"][_at(o, 1.0):_at(o, 6.0)].max() > 0             # the short way from +170 to -170 is a RIGHT turn


def test_unreachable_altitude_command_saturates_gracefully_and_recovers(setup):
    """Hold a 2 km altitude command for 20 s (climb zone, throttle saturated, pitch/throttle commands railed),
    then command 160 m. The aircraft is then ~60 m too high; it must descend and settle without a wound-up
    integrator causing a large undershoot, and without exceeding the pitch limit."""
    t, _, g, _ = setup
    cmds = [(0, 25.0, 100.0 + 2000.0, 0.0), (20, 25.0, 160.0, 0.0)]
    o = mission.simulate_mission(t, g, cmds, 90.0, dt=DT)
    assert np.all(np.isfinite(o["x"])) and np.rad2deg(np.abs(o["theta"]).max()) < 32
    k = _at(o, 20.0)
    assert o["h"][k] > 200.0                                      # it really was climbing hard
    assert o["h"][k:].min() > 160.0 - 15.0                        # modest undershoot after the descent
    i = _at(o, 89.9)
    assert abs(o["h"][i] - 160.0) < 1.0 and abs(o["Va"][i] - 25.0) < 0.3
