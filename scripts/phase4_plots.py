"""Phase 4 figures and mode table -> results/. Run: python scripts/phase4_plots.py"""
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sim import linearize as lin, response as rs, trim  # noqa: E402

OUT = ROOT / "results"
V = 25.0


def mode_table():
    lines = ["| Va [m/s] | mode | eigenvalue [1/s] | wn [rad/s] | zeta | time const / period [s] |", "|---|---|---|---|---|---|"]
    for Va in (18.0, 25.0, 30.0):
        for name, m in lin.modal_summary(lin.linearize(trim.compute_trim(Va))).items():
            lam = m.eigenvalue
            if lam.imag:
                extra = f"period {2 * np.pi / m.wd:.2f}"
                ev = f"{lam.real:.3f} ± {m.wd:.3f}j"
            else:
                extra = f"tau {m.tau:.3f}" + (" (UNSTABLE)" if lam.real > 0 else "")
                ev = f"{lam.real:.3f}"
            lines.append(f"| {Va:.0f} | {name} | {ev} | {m.wn:.3f} | {m.zeta:.3f} | {extra} |")
    (OUT / "phase4_modes.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


def fig_eigs():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for Va, mk in ((18.0, "^"), (25.0, "o"), (30.0, "s")):
        m = lin.linearize(trim.compute_trim(Va))
        for a, (A, _, _) in zip(ax, (m.longitudinal(False), m.lateral(False))):
            ev = np.linalg.eigvals(A)
            a.plot(ev.real, ev.imag, mk, label=f"Va={Va:.0f}", fillstyle="none", ms=8)
    for a, t in zip(ax, ("Longitudinal poles", "Lateral poles")):
        a.axvline(0, color="k", lw=0.5); a.set(xlabel="Re [1/s]", ylabel="Im [rad/s]", title=t, xscale="symlog",
                                               xlim=(-60, 1)); a.grid(alpha=0.3); a.legend()
    fig.tight_layout(); fig.savefig(OUT / "phase4_poles.png", dpi=130)


def fig_step_compare():
    t0 = trim.compute_trim(V); m = lin.linearize(t0)
    cases = [("de", 0.01, [3, 5, 10, 7], ["u", "w", "q", "theta"], "elevator +0.57°"),
             ("dt", 0.02, [3, 5, 10, 7], ["u", "w", "q", "theta"], "throttle +0.02"),
             ("da", 0.01, [4, 9, 11, 6], ["v", "p", "r", "phi"], "aileron +0.57°"),
             ("dr", 0.01, [4, 9, 11, 6], ["v", "p", "r", "phi"], "rudder +0.57°")]
    fig, axs = plt.subplots(2, 4, figsize=(16, 6))
    for col, (ch, amp, idx, names, title) in enumerate(cases):
        f = rs.step_input(m, ch, amp)
        t, Xn = rs.simulate_nonlinear(t0, 10.0, 0.005, f)
        Xl = rs.simulate_linear(m, t, f)
        for row, sel in enumerate((idx[:2], idx[2:])):
            for i in sel:
                line, = axs[row, col].plot(t, Xn[:, i] - t0.x12[i], label=f"{names[idx.index(i)]} nonlin")
                axs[row, col].plot(t, Xl[:, i] - t0.x12[i], "--", color=line.get_color(), label="linear")
            axs[row, col].legend(fontsize=7); axs[row, col].grid(alpha=0.3)
            axs[row, col].set(xlabel="t [s]")
        axs[0, col].set_title(title)
    fig.suptitle("Step responses about trim (Va=25 m/s): nonlinear plant (solid) vs linear model (dashed)")
    fig.tight_layout(); fig.savefig(OUT / "phase4_step_compare.png", dpi=120)


def fig_impulse():
    """2 deg doublets (+2 deg for 0.25 s, then -2 deg for 0.25 s, starting at t = 1 s) on the NONLINEAR plant,
    linear model dashed. Angles in deg, rates in deg/s."""
    t0 = trim.compute_trim(V); m = lin.linearize(t0)
    amp, t_on, T = np.deg2rad(2.0), 1.0, 15.0

    def doublet(channel):
        k = lin.IN[channel]
        def f(t):
            d = m.u_trim.copy()
            if t_on <= t < t_on + 0.25:
                d[k] += amp
            elif t_on + 0.25 <= t < t_on + 0.5:
                d[k] -= amp
            return d
        return f

    def angles(X):
        u, v, w = X[:, 3], X[:, 4], X[:, 5]
        return {"alpha": np.arctan2(w, u), "beta": np.arcsin(v / np.sqrt(u**2 + v**2 + w**2))}

    fig, axs = plt.subplots(3, 2, figsize=(12, 8), sharex=True, gridspec_kw={"height_ratios": [2, 2, 1]})
    for col, (ch, title, sig_a, sig_r) in enumerate(
            (("de", "Elevator doublet 2° (short period; phugoid barely excited)", [("theta", 7), ("alpha", None)], [("q", 10)]),
             ("dr", "Rudder doublet 2° (Dutch roll; spiral grows slowly)", [("phi", 6), ("beta", None)], [("r", 11), ("p", 9)]))):
        f = doublet(ch)
        t, Xn = rs.simulate_nonlinear(t0, T, 0.005, f)
        Xl = rs.simulate_linear(m, t, f)
        for row, sigs in enumerate((sig_a, sig_r)):
            for name, idx in sigs:
                for X, ls in ((Xn, "-"), (Xl, "--")):
                    y = X[:, idx] - t0.x12[idx] if idx is not None else angles(X)[name] - angles(t0.x12[None, :])[name][0]
                    axs[row, col].plot(t, np.rad2deg(y), ls, color=f"C{[n for n, _ in sigs].index(name)}",
                                       label=f"{name}" if ls == "-" else None)
            axs[row, col].grid(alpha=0.3); axs[row, col].legend(fontsize=8, loc="upper right")
        axs[0, col].set(title=title, ylabel="perturbation [deg]")
        axs[1, col].set(ylabel="rate perturbation [deg/s]")
        d = np.array([f(tk)[lin.IN[ch]] - m.u_trim[lin.IN[ch]] for tk in t])
        axs[2, col].plot(t, np.rad2deg(d), "k"); axs[2, col].set(ylabel="input [deg]", xlabel="t [s]"); axs[2, col].grid(alpha=0.3)
    fig.suptitle("Small doublet responses about trim, Va = 25 m/s: nonlinear plant (solid) vs linear model (dashed)")
    fig.tight_layout(); fig.savefig(OUT / "phase4_impulse.png", dpi=130)


if __name__ == "__main__":
    mode_table(); fig_eigs(); fig_step_compare(); fig_impulse()
