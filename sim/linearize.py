"""Numerical linearisation about a trim and modal analysis.

Full model: 12-state Euler plant  x = [pn pe pd u v w phi theta psi p q r],  inputs d = [da de dr dt]
(actual surface/throttle positions; actuator lag is NOT part of the linear design model, as in the book).
Linear model: d(dx)/dt = A dx + B dd with dx = x - x_trim, dd = d - d_trim, from central differences.

Reduced models (Beard & McLain Ch. 5 ordering):
    longitudinal x = [u w q theta h],     inputs [de dt]      (h = -pd)
    lateral      x = [v p r phi psi],     inputs [da dr]
Modal analysis drops the pure integrators h and psi (nothing depends on them).
"""

from dataclasses import dataclass

import numpy as np

from sim import aircraft, dynamics as dyn, rotations as rot, trim as trim_mod

IDX = dict(pn=0, pe=1, pd=2, u=3, v=4, w=5, phi=6, theta=7, psi=8, p=9, q=10, r=11)
IN = dict(da=0, de=1, dr=2, dt=3)

# state steps for central differences (SI units; velocities m/s, angles rad, rates rad/s)
_H_X = 1e-6
_H_U = 1e-6


def plant_deriv_euler(x12, delta, wind_ned=(0, 0, 0), gust_b=(0, 0, 0)):
    """dx12/dt of the full nonlinear plant (gravity + aero + prop) with the Euler-angle model.

    x12 in the layout above; delta = ACTUAL [da, de, dr, dt]; returns dx12/dt.
    """
    x13 = np.concatenate([x12[0:6], rot.euler_to_quat(*x12[6:9]), x12[9:12]])
    f_b, m_b, _ = aircraft.forces_moments(x13, delta, wind_ned, gust_b)
    return dyn.f_euler(x12, f_b, m_b)


def numerical_jacobian(fun, x0, h):
    """Central-difference Jacobian d fun / d x at x0 (fun: R^n -> R^m); h scalar or per-variable step."""
    x0 = np.asarray(x0, float)
    h = np.broadcast_to(h, x0.shape).astype(float)
    cols = []
    for i in range(len(x0)):
        e = np.zeros_like(x0)
        e[i] = h[i]
        cols.append((fun(x0 + e) - fun(x0 - e)) / (2 * h[i]))
    return np.array(cols).T


@dataclass(frozen=True)
class LinearModel:
    A: np.ndarray      # (12, 12)
    B: np.ndarray      # (12, 4), columns [da de dr dt]
    x_trim: np.ndarray  # 12-state Euler trim
    u_trim: np.ndarray  # [da de dr dt]
    trim: trim_mod.Trim

    def sub(self, states, inputs, flip=()):
        """Sub-model with the given state names / input names; states named in `flip` change sign (h = -pd)."""
        si = [IDX[s] for s in states]
        ii = [IN[i] for i in inputs]
        A, B = self.A[np.ix_(si, si)].copy(), self.B[np.ix_(si, ii)].copy()
        for k, s in enumerate(states):
            if s in flip:
                A[k, :] *= -1
                A[:, k] *= -1
                B[k, :] *= -1
        return A, B

    def longitudinal(self, with_h=True):
        """(A, B, state_names) for [u w q theta (h)] and inputs [de dt]."""
        names = ["u", "w", "q", "theta"] + (["pd"] if with_h else [])
        A, B = self.sub(names, ["de", "dt"], flip=("pd",))
        return A, B, names[:-1] + ["h"] if with_h else names

    def lateral(self, with_psi=True):
        """(A, B, state_names) for [v p r phi (psi)] and inputs [da dr]."""
        names = ["v", "p", "r", "phi"] + (["psi"] if with_psi else [])
        A, B = self.sub(names, ["da", "dr"])
        return A, B, names


def linearize(t, h_x=_H_X, h_u=_H_U):
    """Linearise the nonlinear plant about the trim ``t`` (a sim.trim.Trim)."""
    x0, u0 = t.x12, t.delta
    A = numerical_jacobian(lambda x: plant_deriv_euler(x, u0), x0, h_x)
    B = numerical_jacobian(lambda u: plant_deriv_euler(x0, u), u0, h_u)
    return LinearModel(A, B, x0, u0, t)


# ------------------------------------------------------------------ modal analysis
@dataclass(frozen=True)
class Mode:
    name: str
    eigenvalue: complex
    vector: np.ndarray          # eigenvector (state order of the model it came from)

    @property
    def wn(self):
        """Natural frequency [rad/s] = |lambda| (undamped frequency for a complex pair)."""
        return abs(self.eigenvalue)

    @property
    def zeta(self):
        """Damping ratio -Re(lambda)/|lambda| (for real eigenvalues: +-1 by this definition)."""
        return -self.eigenvalue.real / abs(self.eigenvalue)

    @property
    def tau(self):
        """Time constant -1/Re(lambda) [s] (negative = unstable growth time constant magnitude)."""
        return -1.0 / self.eigenvalue.real

    @property
    def wd(self):
        """Damped frequency Im(lambda) [rad/s] (0 for real modes)."""
        return abs(self.eigenvalue.imag)


def _pick(eigs, vecs, k):
    return Mode("", eigs[k], vecs[:, k])


def longitudinal_modes(A4):
    """Short-period and phugoid from the 4-state [u w q theta] matrix: complex pairs, higher |lambda| = short period."""
    eigs, vecs = np.linalg.eig(A4)
    pairs = [k for k in range(4) if eigs[k].imag > 1e-9]
    if len(pairs) != 2:
        raise RuntimeError(f"expected two oscillatory longitudinal modes, got eigenvalues {eigs}")
    pairs.sort(key=lambda k: abs(eigs[k]))
    ph, sp = pairs
    return {"phugoid": Mode("phugoid", eigs[ph], vecs[:, ph]), "short_period": Mode("short_period", eigs[sp], vecs[:, sp])}


def lateral_modes(A4):
    """Roll subsidence, spiral and Dutch roll from the 4-state [v p r phi] matrix."""
    eigs, vecs = np.linalg.eig(A4)
    cplx = [k for k in range(4) if eigs[k].imag > 1e-9]
    real = [k for k in range(4) if abs(eigs[k].imag) <= 1e-9]
    if len(cplx) != 1 or len(real) != 2:
        raise RuntimeError(f"expected one oscillatory + two real lateral modes, got {eigs}")
    real.sort(key=lambda k: abs(eigs[k].real))
    spiral, roll = real
    return {"dutch_roll": Mode("dutch_roll", eigs[cplx[0]], vecs[:, cplx[0]]),
            "roll": Mode("roll", eigs[roll], vecs[:, roll]),
            "spiral": Mode("spiral", eigs[spiral], vecs[:, spiral])}


def modal_summary(model):
    """Dict mode name -> Mode for the decoupled longitudinal and lateral 4-state models."""
    A_lon, _, _ = model.longitudinal(with_h=False)
    A_lat, _, _ = model.lateral(with_psi=False)
    return {**longitudinal_modes(A_lon), **lateral_modes(A_lat)}
