"""Rotation / attitude tests. None depend on any book value."""

import numpy as np
import pytest

from sim import rotations as rot

RNG = np.random.default_rng(1)


def random_angles(n):
    """Euler angles away from gimbal lock: phi, psi in (-pi, pi), theta in (-80, 80) deg."""
    return [(RNG.uniform(-np.pi, np.pi), RNG.uniform(-1.39, 1.39), RNG.uniform(-np.pi, np.pi))
            for _ in range(n)]


def test_rotation_matrix_is_proper_orthogonal():
    for a in random_angles(50):
        R = rot.euler_to_R(*a)
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
        assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-14)


def test_pure_yaw_90deg_points_east():
    """After psi=90 deg the nose points East, so inertial North is on the body's left (-y)."""
    R = rot.euler_to_R(0.0, 0.0, np.pi / 2)
    np.testing.assert_allclose(R @ [1, 0, 0], [0, -1, 0], atol=1e-15)
    np.testing.assert_allclose(R @ [0, 1, 0], [1, 0, 0], atol=1e-15)  # East appears on the nose


def test_euler_quaternion_round_trip_and_same_matrix():
    for a in random_angles(50):
        e = rot.euler_to_quat(*a)
        assert np.linalg.norm(e) == pytest.approx(1.0, abs=1e-14)
        np.testing.assert_allclose(rot.quat_to_R(e), rot.euler_to_R(*a), atol=1e-14)
        np.testing.assert_allclose(rot.quat_to_euler(e), a, atol=1e-12)


@pytest.mark.parametrize("angles", random_angles(20))
def test_euler_rate_matrix_via_finite_difference(angles):
    """S(phi,theta) verified independently of the book.

    Moving the angles along S(phi,theta)@omega must rotate R_b_to_v with
    dR_b_to_v/dt = R_b_to_v @ [omega]x   (Poisson equation, omega in body frame).
    Central finite difference in the Euler angles; also checks the equivalent form
    dR_v_to_b/dt = -[omega]x @ R_v_to_b.
    """
    omega = RNG.normal(size=3)
    angle_dot = rot.euler_rate_matrix(angles[0], angles[1]) @ omega
    h = 1e-6
    a_plus = np.array(angles) + h * angle_dot
    a_minus = np.array(angles) - h * angle_dot
    R_plus, R_minus = rot.euler_to_R(*a_plus), rot.euler_to_R(*a_minus)

    Rdot_v_to_b = (R_plus - R_minus) / (2 * h)
    R = rot.euler_to_R(*angles)
    np.testing.assert_allclose(Rdot_v_to_b.T, R.T @ rot.skew(omega), atol=1e-8)
    np.testing.assert_allclose(Rdot_v_to_b, -rot.skew(omega) @ R, atol=1e-8)


def test_euler_rate_matrix_detects_wrong_matrix():
    """Negative control: a deliberately wrong S (swap sin/cos) must fail the same check."""
    angles = (0.4, 0.3, 1.0)
    omega = np.array([0.7, -0.2, 0.5])
    S_wrong = rot.euler_rate_matrix(angles[0], angles[1]).copy()
    S_wrong[1, 1], S_wrong[1, 2] = -S_wrong[1, 2], S_wrong[1, 1]
    h = 1e-6
    ad = S_wrong @ omega
    Rdot = (rot.euler_to_R(*(np.array(angles) + h * ad))
            - rot.euler_to_R(*(np.array(angles) - h * ad))) / (2 * h)
    R = rot.euler_to_R(*angles)
    assert not np.allclose(Rdot.T, R.T @ rot.skew(omega), atol=1e-6)


def test_gimbal_lock_euler_rates_blow_up():
    norms = [np.linalg.norm(rot.euler_rate_matrix(0.3, np.deg2rad(t)), 2)
             for t in (80.0, 89.0, 89.9, 89.99)]
    assert all(b > a for a, b in zip(norms, norms[1:]))
    assert norms[2] > 500  # 1/cos(89.9 deg) ~ 573
    # The quaternion kinematics matrix stays bounded for the same body rate.
    omega = np.array([1.0, 1.0, 1.0])
    assert np.linalg.norm(rot.quat_omega_matrix(omega), 2) < 2.0
