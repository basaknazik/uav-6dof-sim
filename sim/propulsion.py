"""Electric motor + propeller model (McLain addendum form used by the authors' parameter file).

Model (derivation in docs/02_forces_moments.md section 5):
    V_in = V_max * delta_t                             [V]
    motor : Q_m = K_Q ((V_in - K_V w)/R - i0)          [N m],  w = prop speed [rad/s]
    prop  : T = rho n^2 D^4 C_T(J),  Q = rho n^2 D^5 C_Q(J),  n = w/(2 pi),  J = Va/(n D)
            C_T, C_Q quadratic in J (coefficients in params)
    steady state Q_m = Q  ->  quadratic a w^2 + b w + c = 0, take the positive root.
"""

import numpy as np

from params import aerosonde as P


def prop_coefficients(Va, delta_t):
    """(a, b, c) of the motor/prop torque balance a w^2 + b w + c = 0."""
    V_in = P.V_max * delta_t
    D = P.D_prop
    a = P.rho * D**5 * P.C_Q0 / (4 * np.pi**2)
    b = P.rho * D**4 * Va * P.C_Q1 / (2 * np.pi) + P.KQ * P.KV / P.R_motor
    c = P.rho * D**3 * Va**2 * P.C_Q2 - P.KQ * V_in / P.R_motor + P.KQ * P.i0
    return a, b, c


def prop_speed(Va, delta_t):
    """Steady-state propeller speed w [rad/s] (>= 0). Va = airspeed along the prop axis [m/s]."""
    a, b, c = prop_coefficients(Va, delta_t)
    if c >= 0:  # motor cannot overcome no-load current / aero load: prop stalled
        return 0.0
    return (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)


def thrust_torque(Va, delta_t):
    """Propeller thrust T [N, along +body x] and aerodynamic torque Q [N m, opposing rotation].

    Written in polynomial form in w so it stays finite at w = 0.
    """
    delta_t = float(np.clip(delta_t, P.delta_t_min, P.delta_t_max))
    w = prop_speed(Va, delta_t)
    D, rho = P.D_prop, P.rho
    T = (rho * D**2 * Va**2 * P.C_T2 + rho * D**3 * Va * P.C_T1 * w / (2 * np.pi)
         + rho * D**4 * P.C_T0 * w**2 / (4 * np.pi**2))
    Q = (rho * D**3 * Va**2 * P.C_Q2 + rho * D**4 * Va * P.C_Q1 * w / (2 * np.pi)
         + rho * D**5 * P.C_Q0 * w**2 / (4 * np.pi**2))
    return T, Q
