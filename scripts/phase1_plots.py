"""Generate Phase 1 result figures into results/. Run: python scripts/phase1_plots.py"""

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from sim import dynamics as dyn  # noqa: E402
from sim import integrators as integ  # noqa: E402
from sim import rotations as rot  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parents[1] / "results"
BODY = dyn.AEROSONDE


def qstate(vel=(0, 0, 0), angles=(0, 0, 0), omega=(0, 0, 0)):
    return np.concatenate([[0, 0, 0], vel, rot.euler_to_quat(*angles), omega]).astype(float)


def fig_tumbling():
    x0 = qstate(angles=(0.2, -0.4, 0.9), omega=(1.3, -2.1, 0.8))
    t, X = integ.simulate_rk4(lambda t, x: dyn.f_quat(x, np.zeros(3), np.zeros(3), BODY), x0, 20, 1e-3,
                              post_step=dyn.normalize_quat_state)
    Jw = X[:, 10:13] @ BODY.J.T
    T = 0.5 * np.einsum("ij,ij->i", X[:, 10:13], Jw)
    h_i = np.array([rot.quat_to_R(e).T @ j for e, j in zip(X[:, 6:10], Jw)])
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    ax[0].plot(t, X[:, 10:13]); ax[0].set(title="Body rates (change)", xlabel="t [s]", ylabel="rad/s")
    ax[0].legend(["p", "q", "r"])
    ax[1].plot(t, Jw); ax[1].set(title="Body-frame J·ω (changes)", xlabel="t [s]", ylabel="kg m²/s")
    ax[2].plot(t, h_i - h_i[0]); ax[2].set(title="Inertial h = R_b→v J ω minus h(0) (conserved)", xlabel="t [s]",
                                           ylabel="kg m²/s")
    ax[2].legend(["N", "E", "D"])
    fig.tight_layout(); fig.savefig(OUT / "phase1_tumbling.png", dpi=130); plt.close(fig)


def fig_turn():
    u, r = 25.0, 0.2
    w = np.array([0, 0, r])
    f_b, m_b = BODY.mass * np.cross(w, [u, 0, 0]), np.cross(w, BODY.J @ w)
    x0 = qstate(vel=(u, 0, 0), omega=w)
    t, X = integ.simulate_rk4(lambda t, x: dyn.f_quat(x, f_b, m_b, BODY), x0, 2 * np.pi / r, 1e-2,
                              post_step=dyn.normalize_quat_state)
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    th = np.linspace(0, 2 * np.pi, 400)
    # Heading North at t=0, turning right: North = R sin(rt), East = R (1 - cos(rt)); plot x=East, y=North.
    R = u / r
    east_exact, north_exact = R * (1 - np.cos(r * t)), R * np.sin(r * t)
    err = np.hypot(X[:, 1] - east_exact, X[:, 0] - north_exact).max()
    if err > 1e-5:  # plot-script guard: analytic circle and RK4 path must overlap
        raise AssertionError(f"turn plot: RK4 path and analytic circle differ by {err:.3e} m")
    ax.plot(R * (1 - np.cos(th)), R * np.sin(th), "k--", label=f"analytic circle R=u/r={R:.0f} m, centre (E={R:.0f}, N=0)")
    ax.plot(X[:, 1], X[:, 0], label="RK4 trajectory", lw=1)
    ax.set(xlabel="East [m]", ylabel="North [m]", title="Steady turn: f = m(ω×v)"); ax.axis("equal")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(OUT / "phase1_turn.png", dpi=130); plt.close(fig)


def fig_convergence():
    f_b, m_b = np.array([4.0, 1.0, -3.0]), np.array([0.4, -0.2, 0.3])
    deriv = lambda t, x: dyn.f_quat(x, f_b, m_b, BODY)  # noqa: E731
    x0 = qstate(vel=(20, 1, 0), angles=(0.2, 0.1, 0), omega=(0.5, 0.3, -0.4))
    dts = np.array([0.2, 0.1, 0.05, 0.025, 0.0125])
    errs = []
    for dt in dts:
        t, X = integ.simulate_rk4(deriv, x0, 4.0, dt)
        ref = integ.simulate_scipy(deriv, x0, t, rtol=1e-13, atol=1e-13)
        errs.append(np.abs(X[-1] - ref[-1]).max())
    fig, ax = plt.subplots(figsize=(5, 3.8))
    ax.loglog(dts, errs, "o-", label="RK4 vs DOP853 (tol 1e-13)")
    ax.loglog(dts, errs[0] * (dts / dts[0]) ** 4, "k--", label="slope 4")
    ax.set(xlabel="dt [s]", ylabel="max state error at t=4 s", title="RK4 convergence"); ax.legend()
    fig.tight_layout(); fig.savefig(OUT / "phase1_convergence.png", dpi=130); plt.close(fig)


def fig_gimbal():
    th = np.linspace(0, 89.99, 400)
    n = [np.linalg.norm(rot.euler_rate_matrix(0.3, np.deg2rad(a)), 2) for a in th]
    fig, ax = plt.subplots(figsize=(5, 3.6))
    ax.semilogy(th, n, label="‖S(φ,θ)‖ (Euler rates)")
    ax.axhline(np.linalg.norm(rot.quat_omega_matrix([1, 1, 1]), 2) / 2, color="C1", label="‖½Ω‖ (quaternion), |ω|=√3")
    ax.set(xlabel="pitch θ [deg]", title="Gimbal lock"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "phase1_gimbal.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for f in (fig_tumbling, fig_turn, fig_convergence, fig_gimbal):
        f()
    print("wrote figures to", OUT)
