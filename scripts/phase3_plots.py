"""Phase 3 figure -> results/phase3_trim_map.png"""
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sim import trim  # noqa: E402


def main():
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
    Vas = np.linspace(14, 32, 25)
    for g in (-5, 0, 5):
        rows = [trim.compute_trim(V, gamma=np.deg2rad(g)) for V in Vas if _ok(V, g)]
        V = [t.Va for t in rows]
        ax[0].plot(V, [np.rad2deg(t.alpha) for t in rows], label=f"gamma={g} deg")
        ax[1].plot(V, [t.delta[3] for t in rows], label=f"gamma={g} deg")
    ax[0].set(xlabel="Va [m/s]", ylabel="alpha [deg]", title="Trim angle of attack"); ax[0].legend()
    ax[1].set(xlabel="Va [m/s]", ylabel="throttle", title="Trim throttle (saturates near 33 m/s)"); ax[1].legend()
    Rs = np.linspace(80, 400, 20)
    ts = [trim.compute_trim(25.0, R=R) for R in Rs]
    ax[2].plot(Rs, [np.rad2deg(t.phi) for t in ts], label="bank [deg]")
    ax[2].plot(Rs, [100 * t.delta[3] for t in ts], label="throttle [%]")
    ax[2].set(xlabel="turn radius [m]", title="Level coordinated turn, Va=25 m/s"); ax[2].legend()
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(ROOT / "results/phase3_trim_map.png", dpi=130)


def _ok(V, g):
    try:
        trim.compute_trim(V, gamma=np.deg2rad(g))
        return True
    except RuntimeError:
        return False


if __name__ == "__main__":
    main()
