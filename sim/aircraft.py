"""Full plant: rigid body (sim.dynamics) + gravity + aero + propulsion + actuators.

State (17): quaternion model x13 = [pn pe pd u v w e0 ex ey ez p q r] followed by the
actuator positions [delta_a delta_e delta_r delta_t] (rad, rad, rad, -).
"""

import numpy as np

from params import aerosonde as P
from sim import actuators, aero, dynamics as dyn, propulsion, rotations as rot, wind
from sim.integrators import rk4_step

N_X = 13


def forces_moments(x13, delta, wind_ned=(0, 0, 0), gust_b=(0, 0, 0)):
    """Total body-frame force [N] and moment [N m] (gravity + aero + propulsion) and air data.

    x13   : quaternion-model state; delta: ACTUAL [da, de, dr, dt]
    Returns f_b, m_b, (Va, alpha, beta).
    """
    v_b, e, pqr = x13[3:6], x13[6:10], x13[10:13]
    R = rot.quat_to_R(e)
    Va, alpha, beta, _ = wind.airdata(v_b, R, wind_ned, gust_b)
    f_aero, m_aero = aero.aero_forces_moments(Va, alpha, beta, pqr, delta[:3])
    T, Q = propulsion.thrust_torque(Va, delta[3])
    f_b = dyn.gravity_force_body(R, P.mass) + f_aero + np.array([T, 0.0, 0.0])
    m_b = m_aero + np.array([-Q, 0.0, 0.0])  # prop reaction torque on the airframe about body x
    return f_b, m_b, (Va, alpha, beta)


def plant_deriv(x13, delta, wind_ned=(0, 0, 0), gust_b=(0, 0, 0)):
    """dx13/dt for given ACTUAL surface/throttle positions (no actuator dynamics)."""
    f_b, m_b, _ = forces_moments(x13, delta, wind_ned, gust_b)
    return dyn.f_quat(x13, f_b, m_b, dyn.AEROSONDE)


class Aircraft:
    """Stateful simulator: x17 with fixed-step RK4; wind and commands held constant over a step."""

    def __init__(self, x13, delta0):
        self.x = np.concatenate([x13, delta0]).astype(float)
        self.t = 0.0

    def deriv(self, x, delta_cmd, wind_ned, gust_b):
        return np.concatenate([plant_deriv(x[:N_X], x[N_X:], wind_ned, gust_b),
                               actuators.actuator_rate(x[N_X:], delta_cmd)])

    def step(self, dt, delta_cmd, wind_ned=(0, 0, 0), gust_b=(0, 0, 0)):
        d = lambda t, x: self.deriv(x, delta_cmd, wind_ned, gust_b)  # noqa: E731
        self.x = rk4_step(d, self.t, self.x, dt)
        self.x[6:10] = rot.quat_normalize(self.x[6:10])
        self.t += dt
        return self.x

    def air_data(self, wind_ned=(0, 0, 0), gust_b=(0, 0, 0)):
        return forces_moments(self.x[:N_X], self.x[N_X:], wind_ned, gust_b)[2]
