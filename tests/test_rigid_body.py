"""Physics verification of the 6-DOF rigid-body model against analytic/conservation results."""

import dataclasses

import numpy as np
import pytest

from params import aerosonde
from sim import dynamics as dyn
from sim import integrators as integ
from sim import rotations as rot

BODY = dyn.AEROSONDE
ZERO3 = np.zeros(3)


def quat_state(pos=(0, 0, 0), vel=(0, 0, 0), angles=(0, 0, 0), omega=(0, 0, 0)):
    return np.concatenate([pos, vel, rot.euler_to_quat(*angles), omega]).astype(float)


def closure_quat(f_b=ZERO3, m_b=ZERO3, body=BODY):
    return lambda t, x: dyn.f_quat(x, f_b, m_b, body)


# ---------------------------------------------------------------- tumbling
def test_torque_free_tumbling_conserves_energy_and_angular_momentum():
    omega0 = np.array([1.3, -2.1, 0.8])
    x0 = quat_state(angles=(0.2, -0.4, 0.9), omega=omega0)
    t, X = integ.simulate_rk4(closure_quat(), x0, 20.0, 1e-3, post_step=dyn.normalize_quat_state)
    J = BODY.J
    w = X[:, 10:13]
    Jw = w @ J.T
    T = 0.5 * np.einsum("ij,ij->i", w, Jw)
    h_b_norm = np.linalg.norm(Jw, axis=1)
    h_i = np.array([rot.quat_to_R(e).T @ jw for e, jw in zip(X[:, 6:10], Jw)])  # R_b_to_v J w

    assert np.abs(T / T[0] - 1).max() < 1e-9
    assert np.abs(h_b_norm / h_b_norm[0] - 1).max() < 1e-9
    # all three inertial components conserved (not just the norm)
    np.testing.assert_allclose(h_i, np.tile(h_i[0], (len(t), 1)), rtol=0, atol=1e-9 * np.linalg.norm(h_i[0]))
    # sanity: the BODY-frame momentum vector does change (so the test is non-trivial)
    assert np.abs(Jw - Jw[0]).max() > 0.1


def test_inertial_angular_momentum_check_detects_wrong_kinematics():
    """Negative control: using R instead of R^T for h_i must NOT be conserved."""
    x0 = quat_state(angles=(0.2, -0.4, 0.9), omega=(1.3, -2.1, 0.8))
    _, X = integ.simulate_rk4(closure_quat(), x0, 10.0, 1e-3, post_step=dyn.normalize_quat_state)
    Jw = X[:, 10:13] @ BODY.J.T
    h_wrong = np.array([rot.quat_to_R(e) @ jw for e, jw in zip(X[:, 6:10], Jw)])
    assert np.abs(h_wrong - h_wrong[0]).max() > 0.1


# ---------------------------------------------------------------- intermediate axis
def principal_axes():
    vals, vecs = np.linalg.eigh(BODY.J)  # ascending
    return vals, vecs


def test_intermediate_axis_theorem_using_principal_axes():
    vals, vecs = principal_axes()
    assert vals[0] < vals[1] < vals[2]
    Omega, eps = 5.0, 1e-3
    results = {}
    for k in range(3):
        axis = vecs[:, k]
        perp = vecs[:, (k + 1) % 3] + 0.7 * vecs[:, (k + 2) % 3]
        omega0 = Omega * axis + eps * Omega * perp
        x0 = quat_state(omega=omega0)
        _, X = integ.simulate_rk4(closure_quat(), x0, 30.0, 2e-3, post_step=dyn.normalize_quat_state)
        w = X[:, 10:13]
        results[k] = (w @ axis) / np.linalg.norm(w, axis=1)  # cosine of angle to the principal axis
    assert results[0].min() > 0.99      # min-inertia axis: stable
    assert results[2].min() > 0.99      # max-inertia axis: stable
    assert results[1].min() < -0.9      # intermediate axis: flips (instability)


def test_body_y_is_principal_but_x_z_are_not():
    """With Jxz != 0 only body-y is principal; spinning about body x wobbles."""
    _, vecs = principal_axes()
    assert np.abs(vecs[1, :]).max() > 1 - 1e-12  # one eigenvector is exactly body-y
    ex = np.array([1.0, 0, 0])
    assert np.linalg.norm(np.cross(ex, BODY.J @ ex)) > 0.1


# ---------------------------------------------------------------- free fall
def test_free_fall_matches_analytic_with_nonzero_attitude():
    g = aerosonde.gravity
    angles = (0.3, -0.2, 1.0)
    R = rot.euler_to_R(*angles)
    v_i0 = np.array([10.0, -3.0, 2.0])
    x0 = quat_state(pos=(1, 2, -50), vel=R @ v_i0, angles=angles)

    def deriv(t, x):  # gravity depends on attitude -> evaluated inside the derivative
        return dyn.f_quat(x, dyn.gravity_force_body(rot.quat_to_R(x[6:10]), BODY.mass), ZERO3, BODY)

    T = 5.0
    t, X = integ.simulate_rk4(deriv, x0, T, 1e-2, post_step=dyn.normalize_quat_state)
    g_vec = np.array([0, 0, g])
    p_exact = np.array([1, 2, -50]) + v_i0 * T + 0.5 * g_vec * T**2
    np.testing.assert_allclose(X[-1, 0:3], p_exact, atol=1e-9)
    v_i_end = rot.quat_to_R(X[-1, 6:10]).T @ X[-1, 3:6]
    np.testing.assert_allclose(v_i_end, v_i0 + g_vec * T, atol=1e-9)
    np.testing.assert_allclose(X[-1, 10:13], 0, atol=1e-14)  # no torque, no spin


# ---------------------------------------------------------------- transport theorem
@pytest.mark.parametrize("jxz_zero", [False, True])
def test_steady_turn_is_a_circle_of_radius_u_over_r(jxz_zero):
    """Constant yaw rate r, speed u, f = m (w x v) = [0, m r u, 0]  ->  circle of radius u/r."""
    body = dataclasses.replace(BODY, Jxz=0.0) if jxz_zero else BODY
    u, r = 25.0, 0.2
    omega = np.array([0.0, 0.0, r])
    f_b = body.mass * np.cross(omega, [u, 0, 0])
    m_b = np.cross(omega, body.J @ omega)  # keeps omega constant (== 0 if Jxz == 0)
    if jxz_zero:
        np.testing.assert_allclose(m_b, 0, atol=1e-15)
    else:
        assert m_b[1] == pytest.approx(-body.Jxz * r**2)
    np.testing.assert_allclose(f_b, [0, body.mass * r * u, 0], atol=1e-12)

    x0 = quat_state(vel=(u, 0, 0), omega=omega)
    T = 2 * np.pi / r  # one full lap
    t, X = integ.simulate_rk4(closure_quat(f_b, m_b, body), x0, T, 1e-2, post_step=dyn.normalize_quat_state)

    R_circle = u / r
    centre = np.array([0.0, R_circle])  # started at origin heading North, turning right (East)
    radius = np.linalg.norm(X[:, 0:2] - centre, axis=1)
    assert np.abs(radius - R_circle).max() < 1e-5 * R_circle
    np.testing.assert_allclose(X[:, 3], u, atol=1e-8)          # speed constant
    np.testing.assert_allclose(X[:, 4:6], 0, atol=1e-8)        # v = w = 0
    # full analytic trajectory (heading North at t=0, turning right): not just the radius
    p_exact = np.stack([R_circle * np.sin(r * t), R_circle * (1 - np.cos(r * t))], axis=1)
    np.testing.assert_allclose(X[:, 0:2], p_exact, atol=1e-5)
    phi, theta, psi = zip(*[rot.quat_to_euler(e) for e in X[:, 6:10]])
    np.testing.assert_allclose(phi, 0, atol=1e-8)
    np.testing.assert_allclose(theta, 0, atol=1e-8)
    psi_unwrapped = np.unwrap(psi)
    np.testing.assert_allclose(psi_unwrapped, r * t, atol=1e-8)


def test_steady_turn_without_centripetal_force_is_not_a_circle():
    """Negative control: with f_y = 0 the omega x v term makes v drift and the path is not that circle."""
    u, r = 25.0, 0.2
    omega = np.array([0.0, 0.0, r])
    m_b = np.cross(omega, BODY.J @ omega)
    x0 = quat_state(vel=(u, 0, 0), omega=omega)
    _, X = integ.simulate_rk4(closure_quat(ZERO3, m_b), x0, 10.0, 1e-2, post_step=dyn.normalize_quat_state)
    assert np.abs(X[:, 4]).max() > 1.0  # v grows large; in the circle it stays 0
    centre = np.array([0.0, u / r])
    assert np.abs(np.linalg.norm(X[:, 0:2] - centre, axis=1) - u / r).max() > 1.0


# ---------------------------------------------------------------- quaternion norm
def test_quaternion_norm_drift_with_and_without_renormalisation():
    x0 = quat_state(angles=(0.1, 0.2, 0.3), omega=(3.0, -4.0, 5.0))
    f = closure_quat()
    _, X_free = integ.simulate_rk4(f, x0, 120.0, 0.02)
    _, X_norm = integ.simulate_rk4(f, x0, 120.0, 0.02, post_step=dyn.normalize_quat_state)
    drift_free = np.abs(np.linalg.norm(X_free[:, 6:10], axis=1) - 1).max()
    drift_norm = np.abs(np.linalg.norm(X_norm[:, 6:10], axis=1) - 1).max()
    print(f"quaternion norm drift: free={drift_free:.3e}, renormalised={drift_norm:.3e}")
    assert drift_norm < 1e-14
    assert drift_free < 1e-3          # RK4 keeps it small...
    assert drift_free > 10 * drift_norm  # ...but it is nonzero, hence renormalisation


# ---------------------------------------------------------------- Euler vs quaternion
def test_euler_and_quaternion_models_agree_away_from_gimbal_lock():
    f_b, m_b = np.array([3.0, 1.0, -2.0]), np.array([0.5, -0.3, 0.2])
    xe0 = np.array([0, 0, 0, 20, 1, -1, 0.1, 0.2, 0.3, 0.4, -0.3, 0.5])
    xq0 = dyn.euler_state_to_quat_state(xe0)
    _, XE = integ.simulate_rk4(lambda t, x: dyn.f_euler(x, f_b, m_b, BODY), xe0, 5.0, 1e-3)
    _, XQ = integ.simulate_rk4(closure_quat(f_b, m_b), xq0, 5.0, 1e-3, post_step=dyn.normalize_quat_state)
    assert np.abs(XE[:, 7]).max() < 1.4  # never approached theta = 90 deg
    XQ_as_euler = np.array([dyn.quat_state_to_euler_state(x) for x in XQ])
    np.testing.assert_allclose(XE[:, [0, 1, 2, 3, 4, 5, 9, 10, 11]], XQ_as_euler[:, [0, 1, 2, 3, 4, 5, 9, 10, 11]], atol=1e-8)
    np.testing.assert_allclose(np.unwrap(XE[:, 6:9], axis=0), np.unwrap(XQ_as_euler[:, 6:9], axis=0), atol=1e-8)


# ---------------------------------------------------------------- integrators
def _manoeuvre():
    f_b, m_b = np.array([4.0, 1.0, -3.0]), np.array([0.4, -0.2, 0.3])
    x0 = quat_state(vel=(20, 1, 0), angles=(0.2, 0.1, 0.0), omega=(0.5, 0.3, -0.4))
    return closure_quat(f_b, m_b), x0


def test_rk4_matches_scipy_dop853():
    deriv, x0 = _manoeuvre()
    t, X = integ.simulate_rk4(deriv, x0, 10.0, 1e-3, post_step=dyn.normalize_quat_state)
    Xref = integ.simulate_scipy(deriv, x0, t)
    assert np.abs(X - Xref).max() < 1e-8


def test_rk4_convergence_order_is_four():
    deriv, x0 = _manoeuvre()
    T = 4.0
    errs = []
    for dt in (0.2, 0.1, 0.05):
        t, X = integ.simulate_rk4(deriv, x0, T, dt)  # no renormalisation: pure integrator error
        Xref = integ.simulate_scipy(deriv, x0, t, rtol=1e-13, atol=1e-13)
        errs.append(np.abs(X[-1] - Xref[-1]).max())
    ratios = [errs[0] / errs[1], errs[1] / errs[2]]
    print("RK4 error ratios per halving:", ratios)
    assert all(10 < r < 24 for r in ratios)  # ~16 expected for 4th order
