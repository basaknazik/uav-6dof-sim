"""Phase 4: linearisation, modal analysis, linear-vs-nonlinear comparison."""

import numpy as np
import pytest

from params import aerosonde as P
from sim import aero, propulsion, trim
from sim import dynamics as dyn
from sim import linearize as lin
from sim import response as rs

V = 25.0
G = dyn.AEROSONDE.gammas


@pytest.fixture(scope="module")
def level():
    t = trim.compute_trim(V)
    return t, lin.linearize(t)


# ------------------------------------------------------------------ the Jacobian machinery itself
def test_numerical_jacobian_exact_for_linear_and_accurate_for_nonlinear():
    rng = np.random.default_rng(0)
    M = rng.normal(size=(5, 4))
    J = lin.numerical_jacobian(lambda x: M @ x, rng.normal(size=4), 1e-6)
    np.testing.assert_allclose(J, M, atol=1e-8)
    x0 = np.array([0.3, -1.2])
    J = lin.numerical_jacobian(lambda x: np.array([np.sin(x[0]) * x[1], x[0] ** 3]), x0, 1e-5)
    np.testing.assert_allclose(J, [[np.cos(0.3) * -1.2, np.sin(0.3)], [3 * 0.09, 0.0]], atol=1e-9)


def test_jacobian_independent_of_step_size(level):
    t, m = level
    m2 = lin.linearize(t, h_x=1e-5, h_u=1e-5)
    np.testing.assert_allclose(m2.A, m.A, atol=2e-5)
    np.testing.assert_allclose(m2.B, m.B, atol=2e-5)


# ------------------------------------------------------------------ entries checked against hand derivatives
def test_kinematic_and_gravity_entries(level):
    t, m = level
    A = m.A
    phi, theta = t.phi, t.theta
    assert A[lin.IDX["phi"], lin.IDX["p"]] == pytest.approx(1.0, abs=1e-7)
    assert A[lin.IDX["theta"], lin.IDX["q"]] == pytest.approx(np.cos(phi), abs=1e-7)
    assert A[lin.IDX["psi"], lin.IDX["r"]] == pytest.approx(np.cos(phi) / np.cos(theta), abs=1e-7)
    assert A[lin.IDX["u"], lin.IDX["theta"]] == pytest.approx(-P.gravity * np.cos(theta), abs=1e-6)   # aero indep. of theta
    assert A[lin.IDX["w"], lin.IDX["theta"]] == pytest.approx(-P.gravity * np.sin(theta) * np.cos(phi), abs=1e-6)
    # altitude does not enter the model (constant density)
    assert np.abs(A[3:, lin.IDX["pd"]]).max() == 0.0


def test_control_derivatives_against_analytic_aero(level):
    t, m = level
    B = m.B
    Va, a = t.Va, t.alpha
    qS = 0.5 * P.rho * Va**2 * P.S_wing
    # elevator: pitch moment and normal force
    assert B[lin.IDX["q"], lin.IN["de"]] == pytest.approx(qS * P.c * P.C_m_delta_e / P.Jy, rel=1e-6)
    fz_de = -np.cos(a) * qS * P.C_L_delta_e - np.sin(a) * qS * P.C_D_delta_e
    fx_de = -np.cos(a) * qS * P.C_D_delta_e + np.sin(a) * qS * P.C_L_delta_e
    assert B[lin.IDX["w"], lin.IN["de"]] == pytest.approx(fz_de / P.mass, rel=1e-6)
    assert B[lin.IDX["u"], lin.IN["de"]] == pytest.approx(fx_de / P.mass, rel=1e-6)
    # aileron / rudder through the Gamma coefficients
    for ch, Cl, Cn, CY in (("da", P.C_ell_delta_a, P.C_n_delta_a, P.C_Y_delta_a),
                           ("dr", P.C_ell_delta_r, P.C_n_delta_r, P.C_Y_delta_r)):
        l, n = qS * P.b * Cl, qS * P.b * Cn
        assert B[lin.IDX["p"], lin.IN[ch]] == pytest.approx(G.g3 * l + G.g4 * n, rel=1e-6)
        assert B[lin.IDX["r"], lin.IN[ch]] == pytest.approx(G.g4 * l + G.g8 * n, rel=1e-6)
        assert B[lin.IDX["v"], lin.IN[ch]] == pytest.approx(qS * CY / P.mass, rel=1e-6)
    # throttle: thrust slope from the propulsion model, independent central difference
    h = 1e-6
    dT = (propulsion.thrust_torque(Va, t.delta[3] + h)[0] - propulsion.thrust_torque(Va, t.delta[3] - h)[0]) / (2 * h)
    assert B[lin.IDX["u"], lin.IN["dt"]] == pytest.approx(dT / P.mass, rel=1e-5)


def test_elevator_sign_in_state_space(level):
    _, m = level
    assert m.B[lin.IDX["q"], lin.IN["de"]] < 0       # positive elevator -> nose-down pitch acceleration
    assert m.B[lin.IDX["w"], lin.IN["de"]] < 0       # and small lift increase (fz more negative)
    assert m.B[lin.IDX["p"], lin.IN["da"]] > 0
    assert m.B[lin.IDX["u"], lin.IN["dt"]] > 0


# ------------------------------------------------------------------ structure
def test_longitudinal_lateral_decoupling_at_level_trim(level):
    """Longitudinal states depend on lateral states only through the (tiny) trim bank angle phi0, e.g.
    d(w_dot)/d(phi) = -g sin(phi0), d(theta_dot)/d(r) = -sin(phi0); so the 8-state eigenvalues equal the union of the
    decoupled 4+4 eigenvalues up to O(phi0^2)."""
    _, m = level
    lon = [lin.IDX[s] for s in ("u", "w", "q", "theta")]
    lat = [lin.IDX[s] for s in ("v", "p", "r", "phi")]
    coupling = np.abs(m.A[np.ix_(lon, lat)]).max()
    assert coupling == pytest.approx(P.gravity * abs(np.sin(level[0].phi)), rel=0.05)
    assert coupling < 0.01
    A8 = m.A[np.ix_(lon + lat, lon + lat)]
    full = np.sort_complex(np.linalg.eigvals(A8))
    parts = np.sort_complex(np.concatenate([np.linalg.eigvals(m.longitudinal(False)[0]),
                                            np.linalg.eigvals(m.lateral(False)[0])]))
    np.testing.assert_allclose(full, parts, atol=1e-4)


def test_reduced_models_follow_book_ordering(level):
    _, m = level
    Alon, Blon, names = m.longitudinal()
    assert names == ["u", "w", "q", "theta", "h"] and Alon.shape == (5, 5) and Blon.shape == (5, 2)
    # h_dot = u sin(theta) - w cos(theta) (wings level)  ->  dh_dot/dtheta = u cos(theta) + w sin(theta) = V cos(gamma)
    assert Alon[4, 3] == pytest.approx(V * np.cos(level[0].theta - level[0].alpha), rel=1e-6)
    Alat, Blat, names = m.lateral()
    assert names == ["v", "p", "r", "phi", "psi"] and Alat.shape == (5, 5) and Blat.shape == (5, 2)


# ------------------------------------------------------------------ modes
def test_modes_are_the_expected_types(level):
    _, m = level
    S = lin.modal_summary(m)
    for k in ("short_period", "phugoid", "dutch_roll"):
        assert S[k].eigenvalue.imag != 0 and S[k].eigenvalue.real < 0, k          # stable oscillatory
    assert S["roll"].eigenvalue.imag == 0 and S["roll"].eigenvalue.real < 0
    assert S["spiral"].eigenvalue.imag == 0
    # separation of time scales: roll fastest lateral, short period >> phugoid, spiral slowest
    assert S["short_period"].wn > 10 * S["phugoid"].wn
    assert S["roll"].wn > S["dutch_roll"].wn > S["spiral"].wn
    assert S["spiral"].wn < 0.5 and 0.02 < S["roll"].tau < 0.2
    assert 0.1 < S["dutch_roll"].zeta < 0.6 and 0.1 < S["phugoid"].zeta < 0.6 and 0.3 < S["short_period"].zeta < 0.8
    print({k: (round(v.wn, 3), round(v.zeta, 3), round(v.tau, 3)) for k, v in S.items()})


def test_short_period_matches_two_dof_approximation(level):
    """Classic 2-DOF (alpha, q) approximation with constant speed: char. poly s^2 - (Za+Mq)s + (Za Mq - Ma)."""
    t, m = level
    qS = 0.5 * P.rho * V**2 * P.S_wing
    Za = -qS * (P.C_L_alpha + aero.drag_coefficient(t.alpha)) / (P.mass * V)
    Ma = qS * P.c * P.C_m_alpha / P.Jy
    Mq = qS * P.c**2 * P.C_m_q / (2 * V * P.Jy)
    approx = np.linalg.eigvals(np.array([[Za, 1.0], [Ma, Mq]]))
    sp = lin.modal_summary(m)["short_period"]
    assert sp.wn == pytest.approx(abs(approx[0]), rel=0.05)
    assert sp.zeta == pytest.approx(-approx[0].real / abs(approx[0]), abs=0.03)


def test_roll_subsidence_matches_single_dof_roll(level):
    t, m = level
    qS = 0.5 * P.rho * V**2 * P.S_wing
    Lp = qS * P.b * (P.b / (2 * V)) * (G.g3 * P.C_ell_p + G.g4 * P.C_n_p)
    assert lin.modal_summary(m)["roll"].eigenvalue.real == pytest.approx(Lp, rel=0.05)


def test_dutch_roll_matches_two_dof_approximation(level):
    """(beta, r) approximation: [[Yb, -1], [Nb V, Nr]] with Gamma-weighted yaw derivatives."""
    t, m = level
    qS = 0.5 * P.rho * V**2 * P.S_wing
    Yb = qS * P.C_Y_beta / (P.mass * V)
    Nr = qS * P.b * (P.b / (2 * V)) * (G.g4 * P.C_ell_r + G.g8 * P.C_n_r)
    Nb = qS * P.b * (G.g4 * P.C_ell_beta + G.g8 * P.C_n_beta) / V
    approx = np.linalg.eigvals(np.array([[Yb, -1.0], [Nb * V, Nr]]))
    dr = lin.modal_summary(m)["dutch_roll"]
    assert dr.wn == pytest.approx(abs(approx[0]), rel=0.12)
    assert dr.zeta == pytest.approx(-approx[0].real / abs(approx[0]), abs=0.05)


def test_phugoid_frequency_lanchester_estimate(level):
    """Lanchester: wn = sqrt(2) g / V  (neglects thrust variation and pitch damping coupling) -- loose bound."""
    sp = lin.modal_summary(level[1])["phugoid"]
    assert sp.wn == pytest.approx(np.sqrt(2) * P.gravity / V, rel=0.2)


def test_eigenvector_content_identifies_modes(level):
    _, m = level
    S = lin.modal_summary(m)

    def share(mode, scale):
        v = np.abs(mode.vector * scale)
        return v / v.max()

    # [u, w, q, theta] with u, w scaled by 1/V -> [u/V, alpha, q, theta]; angles/rates in rad, rad/s
    sp, ph = share(S["short_period"], [1 / V, 1 / V, 1, 1]), share(S["phugoid"], [1 / V, 1 / V, 1, 1])
    assert sp[2] == 1.0 and sp[0] < 0.05                     # short period: pitch-rate dominated, speed constant
    assert ph[3] == 1.0 and ph[0] > 0.5 and ph[1] < 0.05    # phugoid: theta and speed, alpha ~ const
    ro, sprl, dr = (share(S["roll"], [1 / V, 1, 1, 1]), share(S["spiral"], [1 / V, 1, 1, 1]),
                    share(S["dutch_roll"], [1 / V, 1, 1, 1]))
    assert ro[1] == 1.0 and ro[0] < 0.01 and ro[2] < 0.05    # roll subsidence: almost pure roll rate
    assert sprl[3] == 1.0                                    # spiral: bank angle dominated
    assert dr[0] > 0.1 and dr[2] > 0.5                       # Dutch roll: sideslip and yaw rate involved


def test_mode_trends_with_airspeed():
    modes = {}
    for Va in (18.0, 25.0, 30.0):
        modes[Va] = lin.modal_summary(lin.linearize(trim.compute_trim(Va)))
    # short-period and Dutch-roll frequencies and roll damping scale ~ V; phugoid frequency falls with V
    assert modes[30.0]["short_period"].wn / modes[18.0]["short_period"].wn == pytest.approx(30 / 18, rel=0.08)
    assert modes[30.0]["roll"].wn / modes[18.0]["roll"].wn == pytest.approx(30 / 18, rel=0.08)
    assert modes[18.0]["phugoid"].wn > modes[25.0]["phugoid"].wn > modes[30.0]["phugoid"].wn


# ------------------------------------------------------------------ linear vs nonlinear
CASES = [("de", [3, 5, 10, 7]), ("dt", [3, 5, 10, 7]), ("da", [4, 9, 11, 6]), ("dr", [4, 9, 11, 6])]


def _rel_error(t, m, channel, amp, states, T=3.0):
    f = rs.step_input(m, channel, amp)
    tt, Xn = rs.simulate_nonlinear(t, T, 0.005, f)
    Xl = rs.simulate_linear(m, tt, f)
    dn, dl = Xn - t.x12, Xl - t.x12
    return max(np.abs(dn[:, i] - dl[:, i]).max() / np.abs(dl[:, i]).max() for i in states)


@pytest.mark.parametrize("channel, states", CASES)
def test_linear_matches_nonlinear_for_small_inputs(level, channel, states):
    t, m = level
    amp = 0.005 if channel != "dt" else 0.01
    assert _rel_error(t, m, channel, amp, states) < 0.06


@pytest.mark.parametrize("channel, states", CASES)
def test_linearisation_error_scales_with_amplitude(level, channel, states):
    """First-order accuracy: halving the input amplitude halves the *relative* error (second-order remainder)."""
    t, m = level
    a = 0.008 if channel != "dt" else 0.02
    e1, e2 = _rel_error(t, m, channel, a, states), _rel_error(t, m, channel, a / 2, states)
    assert e1 / e2 > 1.7, (e1, e2)     # error vanishes at least linearly with amplitude


def test_linear_model_invalid_for_large_inputs(level):
    """Negative control: a big elevator step makes the linear prediction clearly wrong."""
    t, m = level
    assert _rel_error(t, m, "de", 0.3, [3, 5, 10, 7]) > 0.3


def test_wrong_linear_model_has_error_floor_negative_control(level):
    """A wrong A (5 % too stiff) leaves a relative-error floor that does NOT shrink with amplitude."""
    t, m = level
    bad = lin.LinearModel(m.A * 1.05, m.B, m.x_trim, m.u_trim, m.trim)
    e1, e2 = _rel_error(t, bad, "de", 0.008, [3, 5, 10, 7]), _rel_error(t, bad, "de", 0.004, [3, 5, 10, 7])
    assert e1 / e2 < 1.3 and e2 > 0.02
