"""Phase 2 figures -> results/. Run: python scripts/phase2_plots.py"""

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from params import aerosonde as P  # noqa: E402
from sim import aero, propulsion, wind  # noqa: E402

OUT = ROOT / "results"


def main():
    OUT.mkdir(exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
    a = np.linspace(-0.6, 1.6, 800)
    ax[0].plot(np.rad2deg(a), [aero.lift_coefficient(x) for x in a], label="blended $C_L$")
    ax[0].plot(np.rad2deg(a), P.C_L_0 + P.C_L_alpha * a, "--", label="linear")
    ax[0].plot(np.rad2deg(a), 2 * np.sign(a) * np.sin(a) ** 2 * np.cos(a), ":", label="flat plate")
    ax[0].plot(np.rad2deg(a), aero.stall_blend(a), label=r"$\sigma(\alpha)$")
    ax[0].set(xlabel="alpha [deg]", title="Lift with stall blending", ylim=(-1.5, 3.5)); ax[0].legend(fontsize=8)
    for Va in (0, 10, 20, 30):
        d = np.linspace(0, 1, 60)
        ax[1].plot(d, [propulsion.thrust_torque(Va, x)[0] for x in d], label=f"Va={Va} m/s")
    ax[1].set(xlabel="throttle", ylabel="thrust [N]", title="Propeller thrust (note extrapolation at idle)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    g = wind.DrydenGust(0.05, 25.0, np.random.default_rng(1))
    G = np.array([g.step() for _ in range(2400)])
    ax[2].plot(0.05 * np.arange(len(G)), G, lw=0.8)
    ax[2].set(xlabel="t [s]", ylabel="gust [m/s]", title="Dryden gusts (body frame)"); ax[2].legend(["u", "v", "w"])
    fig.tight_layout(); fig.savefig(OUT / "phase2_models.png", dpi=130)
    print("wrote", OUT / "phase2_models.png")


if __name__ == "__main__":
    main()
