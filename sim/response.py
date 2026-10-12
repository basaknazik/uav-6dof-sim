"""Open-loop time responses: nonlinear plant vs. linearised model (inputs applied directly, no actuator lag)."""

import numpy as np
from scipy import signal

from sim import linearize as lin
from sim.integrators import rk4_step


def simulate_nonlinear(t_trim, T, dt, delta_fn):
    """RK4 simulation of the 12-state Euler plant from the trim state.

    delta_fn(t) -> actual [da, de, dr, dt] (held constant over each step). Returns (t, X) with X (N+1, 12).
    """
    n = int(round(T / dt))
    t = dt * np.arange(n + 1)
    X = np.empty((n + 1, 12))
    X[0] = t_trim.x12
    for k in range(n):
        d = delta_fn(t[k])
        X[k + 1] = rk4_step(lambda tt, x: lin.plant_deriv_euler(x, d), t[k], X[k], dt)
    return t, X


def simulate_linear(model, t, delta_fn):
    """Linear response x_trim + dx(t) of the 12-state model (exact matrix-exponential solver, linear input interpolation)."""
    U = np.array([delta_fn(tk) - model.u_trim for tk in t])
    sys = signal.StateSpace(model.A, model.B, np.eye(12), np.zeros((12, 4)))
    _, y, _ = signal.lsim(sys, U, t, X0=np.zeros(12))
    return model.x_trim + y


def step_input(model, channel, amplitude, t_on=0.0):
    """delta_fn applying a step of `amplitude` [rad or throttle units] on channel in {da, de, dr, dt}."""
    k = lin.IN[channel]

    def f(t):
        d = model.u_trim.copy()
        if t >= t_on:
            d[k] += amplitude
        return d
    return f
