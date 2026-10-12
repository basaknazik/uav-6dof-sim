"""Aerodynamic forces and moments (body frame) from stability/control derivatives.

All forces in N, moments in N m, about the CG, expressed in the BODY frame.
Excludes gravity (sim.dynamics) and propulsion (sim.propulsion).
Model: Beard & McLain Ch. 4 (nonlinear lift with stall blending; linear lateral model);
coefficients from params/aerosonde.py.
"""

import numpy as np

from params import aerosonde as P
from sim.wind import V_A_MIN


def stall_blend(alpha):
    """sigma(alpha) in [0,1]: 0 in the linear regime, 1 post-stall. sigma(+-alpha0) = 0.5."""
    e1 = np.exp(-P.M_blend * (alpha - P.alpha0))
    e2 = np.exp(P.M_blend * (alpha + P.alpha0))
    return (1 + e1 + e2) / ((1 + e1) * (1 + e2))


def lift_coefficient(alpha):
    """C_L(alpha) [-]: blend of linear (C_L0 + C_La alpha) and flat-plate 2 sign(a) sin^2(a) cos(a)."""
    s = stall_blend(alpha)
    linear = P.C_L_0 + P.C_L_alpha * alpha
    flat_plate = 2.0 * np.sign(alpha) * np.sin(alpha) ** 2 * np.cos(alpha)
    return (1 - s) * linear + s * flat_plate


def drag_coefficient(alpha):
    """C_D(alpha) [-]: parasitic + induced, C_Dp + (C_L0 + C_La alpha)^2 / (pi e AR)."""
    return P.C_D_p + (P.C_L_0 + P.C_L_alpha * alpha) ** 2 / (np.pi * P.e_oswald * P.AR)


def aero_forces_moments(Va, alpha, beta, pqr, delta_aer):
    """Aerodynamic force [fx, fy, fz] and moment [l, m, n] in body frame.

    Va [m/s], alpha, beta [rad], pqr = body rates [rad/s], delta_aer = [delta_a, delta_e, delta_r] [rad].
    """
    p, q, r = pqr
    da, de, dr = delta_aer
    Va_s = max(Va, V_A_MIN)
    qbar = 0.5 * P.rho * Va**2 * P.S_wing                  # dynamic pressure * S  [N]
    cq = P.c / (2 * Va_s)                                  # q-hat scaling [s]
    bq = P.b / (2 * Va_s)                                  # p-hat, r-hat scaling [s]

    CL = lift_coefficient(alpha)
    CD = drag_coefficient(alpha)
    F_lift = qbar * (CL + P.C_L_q * cq * q + P.C_L_delta_e * de)
    F_drag = qbar * (CD + P.C_D_q * cq * q + P.C_D_delta_e * de)

    # stability axes -> body axes: rotate (-drag, -lift) by alpha about body y
    ca, sa = np.cos(alpha), np.sin(alpha)
    fx = -ca * F_drag + sa * F_lift
    fz = -sa * F_drag - ca * F_lift
    fy = qbar * (P.C_Y_0 + P.C_Y_beta * beta + P.C_Y_p * bq * p + P.C_Y_r * bq * r
                 + P.C_Y_delta_a * da + P.C_Y_delta_r * dr)

    l = qbar * P.b * (P.C_ell_0 + P.C_ell_beta * beta + P.C_ell_p * bq * p + P.C_ell_r * bq * r
                      + P.C_ell_delta_a * da + P.C_ell_delta_r * dr)
    m = qbar * P.c * (P.C_m_0 + P.C_m_alpha * alpha + P.C_m_q * cq * q + P.C_m_delta_e * de)
    n = qbar * P.b * (P.C_n_0 + P.C_n_beta * beta + P.C_n_p * bq * p + P.C_n_r * bq * r
                      + P.C_n_delta_a * da + P.C_n_delta_r * dr)
    return np.array([fx, fy, fz]), np.array([l, m, n])
