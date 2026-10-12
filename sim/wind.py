"""Air-relative velocity (V_a, alpha, beta) and the Dryden gust model.

Frames: steady wind is given in the inertial NED frame [m/s]; gusts are
generated and applied in the BODY frame (Beard & McLain convention).
"""

import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov

from params import aerosonde as P

V_A_MIN = 1e-3  # m/s; floor to avoid dividing by zero in alpha/beta/nondimensional rates


def airdata(v_b, R_v_to_b, wind_ned=(0.0, 0.0, 0.0), gust_b=(0.0, 0.0, 0.0)):
    """Air-relative velocity and flow angles.

    v_b      : ground-relative CG velocity [u, v, w] in body frame [m/s]
    R_v_to_b : rotation NED -> body
    wind_ned : steady wind velocity of the air mass in NED [m/s]
    gust_b   : gust velocity in body frame [m/s]
    Returns (Va [m/s], alpha [rad], beta [rad], v_air_b [m/s, body frame]) with
        v_air_b = v_b - (R_v_to_b @ wind_ned + gust_b)        = [u_r, v_r, w_r]
        Va = |v_air_b|, alpha = atan2(w_r, u_r), beta = asin(v_r / Va)
    """
    v_air_b = np.asarray(v_b, float) - (R_v_to_b @ np.asarray(wind_ned, float) + np.asarray(gust_b, float))
    ur, vr, wr = v_air_b
    Va = np.linalg.norm(v_air_b)
    alpha = np.arctan2(wr, ur)
    beta = np.arcsin(np.clip(vr / max(Va, V_A_MIN), -1.0, 1.0))
    return Va, alpha, beta, v_air_b


# ----------------------------------------------------------------- Dryden
def dryden_state_space(Va0, L, sigma, channel):
    """Continuous state-space (A, B, C) of one Dryden shaping filter driven by unit white noise.

    channel 'u':  H = sigma*sqrt(2 Va/L) / (s + Va/L)
    channel 'v'/'w': H = sigma*sqrt(3 Va/L) (s + Va/(sqrt(3) L)) / (s + Va/L)^2
    Va0 [m/s] nominal airspeed, L [m] scale length, sigma [m/s] intensity.
    """
    a = Va0 / L
    if channel == "u":
        return (np.array([[-a]]), np.array([[1.0]]), np.array([[sigma * np.sqrt(2 * a)]]))
    K = sigma * np.sqrt(3 * a)
    A = np.array([[0.0, 1.0], [-a * a, -2 * a]])
    B = np.array([[0.0], [1.0]])
    C = K * np.array([[a / np.sqrt(3), 1.0]])
    return A, B, C


def discretize_noise_driven(A, B, Ts):
    """Exact discretisation x+ = Phi x + w, w ~ N(0, Qd), of dx = A x dt + B dW (unit white noise).

    Uses the Van Loan matrix-exponential method. Returns (Phi, Qd).
    """
    n = A.shape[0]
    M = np.zeros((2 * n, 2 * n))
    M[:n, :n] = -A
    M[:n, n:] = B @ B.T
    M[n:, n:] = A.T
    E = expm(M * Ts)
    Phi = E[n:, n:].T
    Qd = Phi @ E[:n, n:]
    return Phi, 0.5 * (Qd + Qd.T)


def dryden_output_covariance(A, B, C):
    """Stationary state covariance P (A P + P A^T + B B^T = 0) and output variance C P C^T."""
    Pcov = solve_continuous_lyapunov(A, -B @ B.T)
    return Pcov, (C @ Pcov @ C.T).item()


class DrydenGust:
    """Dryden turbulence [u_g, v_g, w_g] in the BODY frame, sampled every Ts seconds.

    Filters use a fixed nominal airspeed Va0 (as in the book) and are discretised exactly,
    so stationary statistics do not depend on Ts. Starts from the stationary distribution.
    """

    def __init__(self, Ts, Va0, rng=None, L=None, sigma=None):
        L = L or (P.dryden_Lu, P.dryden_Lv, P.dryden_Lw)
        sigma = sigma or (P.dryden_sigma_u, P.dryden_sigma_v, P.dryden_sigma_w)
        self.rng = rng or np.random.default_rng()
        self.filters = []
        for ch, Li, si in zip("uvw", L, sigma):
            A, B, C = dryden_state_space(Va0, Li, si, ch)
            Phi, Qd = discretize_noise_driven(A, B, Ts)
            Pcov, _ = dryden_output_covariance(A, B, C)
            x0 = np.linalg.cholesky(Pcov) @ self.rng.standard_normal(len(A))
            self.filters.append([Phi, np.linalg.cholesky(Qd), C, x0])

    def step(self):
        """Advance one sample; returns gust [u_g, v_g, w_g] in body frame [m/s]."""
        out = np.empty(3)
        for i, f in enumerate(self.filters):
            Phi, Lq, C, x = f
            f[3] = Phi @ x + Lq @ self.rng.standard_normal(len(x))
            out[i] = (C @ f[3]).item()
        return out
