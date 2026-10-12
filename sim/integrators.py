"""Fixed-step RK4 and a scipy reference integrator.

``deriv(t, x) -> dx/dt`` is any callable; inputs (forces, moments) are
captured by the closure. Fixed step is deliberate: SIL/HIL loops run at a
fixed rate with deterministic timing.
"""

import numpy as np
from scipy.integrate import solve_ivp


def rk4_step(deriv, t, x, dt):
    """One classical RK4 step of size dt [s]; local error O(dt^5)."""
    k1 = deriv(t, x)
    k2 = deriv(t + dt / 2, x + dt / 2 * k1)
    k3 = deriv(t + dt / 2, x + dt / 2 * k2)
    k4 = deriv(t + dt, x + dt * k3)
    return x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate_rk4(deriv, x0, t_final, dt, post_step=None, t0=0.0):
    """Integrate with fixed-step RK4 from t0 to t_final [s].

    post_step(x) -> x is applied after every step (e.g. quaternion renormalisation).
    Returns (t, X): t shape (N+1,), X shape (N+1, n_states).
    """
    n_steps = int(round((t_final - t0) / dt))
    t = t0 + dt * np.arange(n_steps + 1)
    X = np.empty((n_steps + 1, len(x0)))
    X[0] = x0
    for k in range(n_steps):
        x_next = rk4_step(deriv, t[k], X[k], dt)
        X[k + 1] = post_step(x_next) if post_step is not None else x_next
    return t, X


def simulate_scipy(deriv, x0, t_eval, method="DOP853", rtol=1e-12, atol=1e-12):
    """Reference solution with scipy.solve_ivp at the times t_eval [s]. Returns X (len(t_eval), n)."""
    sol = solve_ivp(deriv, (t_eval[0], t_eval[-1]), x0, t_eval=t_eval,
                    method=method, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T
