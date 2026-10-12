"""Actuators: command saturation followed by a first-order lag (rate = (cmd_sat - x)/tau).

Order of ops: saturate the COMMAND, then lag, so the state can never leave the limits.
Command order everywhere: [delta_a, delta_e, delta_r, delta_t].
"""

import numpy as np

from params import aerosonde as P

LIMITS_LO = np.array([-P.delta_surface_max] * 3 + [P.delta_t_min])
LIMITS_HI = np.array([P.delta_surface_max] * 3 + [P.delta_t_max])
TAU = np.array([P.tau_servo] * 3 + [P.tau_throttle])


def saturate(delta_cmd):
    """Clip [da, de, dr, dt] to actuator limits (rad for surfaces, 0..1 throttle)."""
    return np.clip(delta_cmd, LIMITS_LO, LIMITS_HI)


def actuator_rate(delta, delta_cmd):
    """d(delta)/dt [rad/s or 1/s] of the first-order actuators."""
    return (saturate(delta_cmd) - delta) / TAU
