"""Attitude representations: Euler angles (3-2-1) and unit quaternions.

Conventions (Beard & McLain):
    * Rotation sequence: yaw psi (z), pitch theta (y), roll phi (x).
    * R_v_to_b maps components of a vector from the inertial/vehicle NED frame
      into the body frame:  x_b = R_v_to_b @ x_v.   R_b_to_v is its transpose.
    * Quaternion e = [e0, ex, ey, ez], scalar first, unit norm, rotating the
      vehicle frame into the body frame (same R_v_to_b as the Euler matrix).
All angles in rad, rates in rad/s.
"""

import numpy as np


def skew(w):
    """Cross-product matrix [w]x such that skew(w) @ a == np.cross(w, a)."""
    return np.array([[0.0, -w[2], w[1]],
                     [w[2], 0.0, -w[0]],
                     [-w[1], w[0], 0.0]])


def euler_to_R(phi, theta, psi):
    """Rotation matrix R_v_to_b (vehicle NED -> body) from 3-2-1 Euler angles [rad]."""
    cphi, sphi = np.cos(phi), np.sin(phi)
    cth, sth = np.cos(theta), np.sin(theta)
    cpsi, spsi = np.cos(psi), np.sin(psi)
    return np.array([
        [cth * cpsi, cth * spsi, -sth],
        [sphi * sth * cpsi - cphi * spsi, sphi * sth * spsi + cphi * cpsi, sphi * cth],
        [cphi * sth * cpsi + sphi * spsi, cphi * sth * spsi - sphi * cpsi, cphi * cth],
    ])


def euler_rate_matrix(phi, theta):
    """S(phi, theta) with [phi_dot, theta_dot, psi_dot] = S @ [p, q, r].

    Input angles [rad]; body rates [rad/s] in the body frame. Singular at
    theta = +-90 deg (gimbal lock): entries contain tan(theta), 1/cos(theta).
    """
    cphi, sphi = np.cos(phi), np.sin(phi)
    cth = np.cos(theta)
    tth = np.tan(theta)
    return np.array([
        [1.0, sphi * tth, cphi * tth],
        [0.0, cphi, -sphi],
        [0.0, sphi / cth, cphi / cth],
    ])


def euler_to_quat(phi, theta, psi):
    """Unit quaternion [e0, ex, ey, ez] (scalar first) from Euler angles [rad]."""
    cphi, sphi = np.cos(phi / 2), np.sin(phi / 2)
    cth, sth = np.cos(theta / 2), np.sin(theta / 2)
    cpsi, spsi = np.cos(psi / 2), np.sin(psi / 2)
    return np.array([
        cpsi * cth * cphi + spsi * sth * sphi,
        cpsi * cth * sphi - spsi * sth * cphi,
        cpsi * sth * cphi + spsi * cth * sphi,
        spsi * cth * cphi - cpsi * sth * sphi,
    ])


def quat_to_euler(e):
    """Euler angles (phi, theta, psi) [rad] from a quaternion [e0, ex, ey, ez]."""
    e0, ex, ey, ez = e
    phi = np.arctan2(2 * (e0 * ex + ey * ez), e0**2 + ez**2 - ex**2 - ey**2)
    theta = np.arcsin(np.clip(2 * (e0 * ey - ex * ez), -1.0, 1.0))
    psi = np.arctan2(2 * (e0 * ez + ex * ey), e0**2 + ex**2 - ey**2 - ez**2)
    return phi, theta, psi


def quat_to_R(e):
    """Rotation matrix R_v_to_b (vehicle NED -> body) from quaternion [e0, ex, ey, ez]."""
    e0, ex, ey, ez = e
    return np.array([
        [e0**2 + ex**2 - ey**2 - ez**2, 2 * (ex * ey + e0 * ez), 2 * (ex * ez - e0 * ey)],
        [2 * (ex * ey - e0 * ez), e0**2 - ex**2 + ey**2 - ez**2, 2 * (ey * ez + e0 * ex)],
        [2 * (ex * ez + e0 * ey), 2 * (ey * ez - e0 * ex), e0**2 - ex**2 - ey**2 + ez**2],
    ])


def quat_normalize(e):
    """Return e / ||e||."""
    return e / np.linalg.norm(e)


def quat_omega_matrix(omega):
    """Omega(omega) in e_dot = 0.5 * Omega @ e, omega = [p, q, r] [rad/s] in body frame."""
    p, q, r = omega
    return np.array([
        [0.0, -p, -q, -r],
        [p, 0.0, r, -q],
        [q, -r, 0.0, p],
        [r, q, -p, 0.0],
    ])
