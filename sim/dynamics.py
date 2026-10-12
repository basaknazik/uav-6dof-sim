"""Rigid-body 6-DOF equations of motion (flat, non-rotating Earth, constant mass).

State layouts (all SI; position in the inertial NED frame, velocities and rates
expressed in the body frame):
    Euler state, 12 : [pn, pe, pd, u, v, w, phi, theta, psi, p, q, r]
    Quaternion state, 13 : [pn, pe, pd, u, v, w, e0, ex, ey, ez, p, q, r]

Inputs are a body-frame force f_b = [fx, fy, fz] [N] (including gravity, which
the caller supplies) and a body-frame moment m_b = [l, m, n] [N m] about the CG.
No aerodynamics here -- those arrive in Phase 2.
"""

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

from params import aerosonde
from sim import rotations as rot


class Gammas(NamedTuple):
    """Inertia coefficients Gamma1..Gamma8 [units: 1/(kg m^2) or dimensionless]."""
    g1: float
    g2: float
    g3: float
    g4: float
    g5: float
    g6: float
    g7: float
    g8: float


@dataclass(frozen=True)
class RigidBody:
    """Mass and inertia. mass [kg]; Jx, Jy, Jz, Jxz [kg m^2] about the CG, body axes."""
    mass: float
    Jx: float
    Jy: float
    Jz: float
    Jxz: float

    @property
    def J(self):
        """Inertia matrix [kg m^2] in body axes (note the -Jxz off-diagonals)."""
        return np.array([[self.Jx, 0.0, -self.Jxz],
                         [0.0, self.Jy, 0.0],
                         [-self.Jxz, 0.0, self.Jz]])

    @property
    def gammas(self):
        """Gamma1..Gamma8 from the inversion of J (docs/01_rigid_body.md, section 5)."""
        Jx, Jy, Jz, Jxz = self.Jx, self.Jy, self.Jz, self.Jxz
        G = Jx * Jz - Jxz**2
        return Gammas(
            Jxz * (Jx - Jy + Jz) / G,
            (Jz * (Jz - Jy) + Jxz**2) / G,
            Jz / G,
            Jxz / G,
            (Jz - Jx) / Jy,
            Jxz / Jy,
            ((Jx - Jy) * Jx + Jxz**2) / G,
            Jx / G,
        )


AEROSONDE = RigidBody(aerosonde.mass, aerosonde.Jx, aerosonde.Jy,
                      aerosonde.Jz, aerosonde.Jxz)


def omega_dot_gamma(omega, moment, body):
    """Angular acceleration [rad/s^2] via the Gamma form. omega, moment in body frame."""
    p, q, r = omega
    l, m, n = moment
    g = body.gammas
    return np.array([
        g.g1 * p * q - g.g2 * q * r + g.g3 * l + g.g4 * n,
        g.g5 * p * r - g.g6 * (p**2 - r**2) + m / body.Jy,
        g.g7 * p * q - g.g1 * q * r + g.g4 * l + g.g8 * n,
    ])


def omega_dot_numeric(omega, moment, body):
    """Angular acceleration [rad/s^2] as solve(J, m - omega x J omega); cross-check of the Gamma form."""
    omega = np.asarray(omega, dtype=float)
    return np.linalg.solve(body.J, np.asarray(moment) - np.cross(omega, body.J @ omega))


def _translational_dot(v_b, omega, f_b, mass):
    """[u_dot, v_dot, w_dot] [m/s^2] in body frame: transport theorem, -omega x V + f/m."""
    return -np.cross(omega, v_b) + np.asarray(f_b) / mass


def gravity_force_body(R_v_to_b, mass, g=aerosonde.gravity):
    """Gravity force [N] in body frame; R_v_to_b maps NED -> body, NED gravity = [0, 0, g]."""
    return mass * g * R_v_to_b[:, 2]


def f_euler(x, f_b, m_b, body=AEROSONDE):
    """Time derivative of the 12-state Euler model. See module docstring for layout."""
    v_b, angles, omega = x[3:6], x[6:9], x[9:12]
    R_v_to_b = rot.euler_to_R(*angles)
    return np.concatenate([
        R_v_to_b.T @ v_b,
        _translational_dot(v_b, omega, f_b, body.mass),
        rot.euler_rate_matrix(angles[0], angles[1]) @ omega,
        omega_dot_gamma(omega, m_b, body),
    ])


def f_quat(x, f_b, m_b, body=AEROSONDE):
    """Time derivative of the 13-state quaternion model. See module docstring for layout."""
    v_b, e, omega = x[3:6], x[6:10], x[10:13]
    R_v_to_b = rot.quat_to_R(e)
    return np.concatenate([
        R_v_to_b.T @ v_b,
        _translational_dot(v_b, omega, f_b, body.mass),
        0.5 * rot.quat_omega_matrix(omega) @ e,
        omega_dot_gamma(omega, m_b, body),
    ])


def normalize_quat_state(x):
    """Return a copy of a 13-state vector with its quaternion renormalised to unit norm."""
    x = np.array(x, dtype=float)
    x[6:10] = rot.quat_normalize(x[6:10])
    return x


def euler_state_to_quat_state(x):
    """Convert a 12-state Euler vector to the 13-state quaternion vector."""
    return np.concatenate([x[0:6], rot.euler_to_quat(*x[6:9]), x[9:12]])


def quat_state_to_euler_state(x):
    """Convert a 13-state quaternion vector to the 12-state Euler vector."""
    return np.concatenate([x[0:6], rot.quat_to_euler(x[6:10]), x[10:13]])
