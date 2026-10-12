"""Trim: steady straight / climbing / turning flight found by numerical optimisation.

A trim condition is specified by (V_a, gamma, R, beta):
    V_a   air-relative speed [m/s]            (no wind during trim)
    gamma flight-path angle [rad]             (positive = climbing)
    R     signed turn radius [m]              (+ right turn, np.inf = straight)
    beta  sideslip [rad]                      (default 0 = coordinated)
beta is a specification, not an unknown: without it steady sideslip with bank gives a one-parameter
family of solutions (docs/03_trim.md).

Unknowns z = [alpha, theta, phi, delta_a, delta_e, delta_r, delta_t]. Steady helical flight means
phi_dot = theta_dot = 0 and psi_dot = V_a cos(gamma)/R, which fixes the body rates
    p = -psi_dot sin(theta),  q = psi_dot sin(phi) cos(theta),  r = psi_dot cos(phi) cos(theta).
Residuals (7): u_dot, v_dot, w_dot, p_dot, q_dot, r_dot = 0 and h_dot - V_a sin(gamma) = 0.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from params import aerosonde as P
from sim import actuators, aircraft, rotations as rot

N_UNKNOWN = 7


@dataclass(frozen=True)
class Trim:
    Va: float
    gamma: float
    R: float
    beta: float
    alpha: float
    theta: float
    phi: float
    delta: np.ndarray       # [da, de, dr, dt]
    x13: np.ndarray         # quaternion-model trim state (position [0, 0, -h])
    cost: float             # max |residual|

    @property
    def x12(self):
        """12-state Euler vector [pn pe pd u v w phi theta psi p q r]."""
        e = self.x13
        return np.concatenate([e[0:6], rot.quat_to_euler(e[6:10]), e[10:13]])


def _turn_rate(Va, gamma, R):
    return 0.0 if not np.isfinite(R) else Va * np.cos(gamma) / R


def build_state(Va, alpha, beta, theta, phi, psi_dot, h=100.0, psi=0.0):
    """Quaternion-model state for steady flight (uvw from V_a, alpha, beta; rates from psi_dot)."""
    u = Va * np.cos(alpha) * np.cos(beta)
    v = Va * np.sin(beta)
    w = Va * np.sin(alpha) * np.cos(beta)
    p = -psi_dot * np.sin(theta)
    q = psi_dot * np.sin(phi) * np.cos(theta)
    r = psi_dot * np.cos(phi) * np.cos(theta)
    return np.concatenate([[0.0, 0.0, -h], [u, v, w], rot.euler_to_quat(phi, theta, psi), [p, q, r]])


def residuals(z, Va, gamma, R, beta):
    """Trim residual vector (7,) for unknowns z; zero at a trim."""
    alpha, theta, phi = z[0:3]
    delta = z[3:7][[0, 1, 2, 3]]
    x13 = build_state(Va, alpha, beta, theta, phi, _turn_rate(Va, gamma, R))
    xdot = aircraft.plant_deriv(x13, delta)
    h_dot = -xdot[2]
    return np.concatenate([xdot[3:6], xdot[10:13], [h_dot - Va * np.sin(gamma)]])


def compute_trim(Va, gamma=0.0, R=np.inf, beta=0.0, h=100.0, tol=1e-10):
    """Solve for trim. Raises RuntimeError if the residual does not reach ``tol``."""
    psi_dot = _turn_rate(Va, gamma, R)
    phi0 = np.arctan(Va * psi_dot / P.gravity)
    z0 = np.array([0.05, 0.05 + gamma, phi0, 0.0, 0.0, 0.0, 0.5])
    lo = np.concatenate([[-0.4, -1.2, -1.2], actuators.LIMITS_LO])
    hi = np.concatenate([[0.4, 1.2, 1.2], actuators.LIMITS_HI])
    z0 = np.clip(z0, lo + 1e-9, hi - 1e-9)
    sol = least_squares(residuals, z0, args=(Va, gamma, R, beta), bounds=(lo, hi),
                        xtol=1e-15, ftol=1e-15, gtol=1e-15, x_scale=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5])
    cost = np.abs(sol.fun).max()
    if cost > tol:
        raise RuntimeError(f"trim did not converge: max residual {cost:.2e} (Va={Va}, gamma={gamma}, R={R})")
    alpha, theta, phi = sol.x[:3]
    x13 = build_state(Va, alpha, beta, theta, phi, psi_dot, h=h)
    return Trim(Va, gamma, R, beta, alpha, theta, phi, sol.x[3:7].copy(), x13, cost)
