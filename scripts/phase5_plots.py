"""Phase 5 figures -> results/. Run: python scripts/phase5_plots.py"""
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sim import autopilot as ap, linearize as lin, mission, trim  # noqa: E402

D = np.deg2rad
CMDS = [(0, 25.0, 100.0, 0.0), (5, 25.0, 150.0, D(90)), (60, 28.0, 150.0, D(90)),
        (90, 28.0, 120.0, D(45)), (130, 25.0, 120.0, D(45))]


def main():
    t0 = trim.compute_trim(25.0)
    g, _ = ap.design_gains(lin.linearize(t0))
    runs = {"no wind": mission.simulate_mission(t0, g, CMDS, 170.0, dt=0.01),
            "wind (-4, 6, 0) m/s + Dryden gusts": mission.simulate_mission(
                t0, g, CMDS, 170.0, dt=0.01, wind_ned=(-4.0, 6.0, 0.0), gusts=True, seed=1)}
    fig, ax = plt.subplots(2, 3, figsize=(16, 7))
    for label, o in runs.items():
        ax[0, 0].plot(o["t"], o["h"], label=label)
        ax[0, 1].plot(o["t"], np.rad2deg(np.unwrap(o["chi"])), label=label)
        ax[0, 2].plot(o["t"], o["Va"], label=label)
        ax[1, 0].plot(o["t"], np.rad2deg(o["phi"]), label=label)
        ax[1, 1].plot(o["t"], np.rad2deg(o["theta"]), label=label)
        ax[1, 2].plot(o["t"], o["x"][:, 16], label=label)
    c = runs["no wind"]
    ax[0, 0].plot(c["t"], c["cmd"][:, 1], "k--", label="command")
    ax[0, 1].plot(c["t"], np.rad2deg(c["cmd"][:, 2]), "k--")
    ax[0, 2].plot(c["t"], c["cmd"][:, 0], "k--")
    for a, (title, yl) in zip(ax.ravel(), [("Altitude", "m"), ("Ground course", "deg"), ("Airspeed", "m/s"),
                                           ("Bank angle", "deg"), ("Pitch angle", "deg"), ("Throttle (actual)", "-")]):
        a.set(title=title, ylabel=yl, xlabel="t [s]"); a.grid(alpha=0.3)
    ax[0, 0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(ROOT / "results/phase5_mission.png", dpi=120)
    # ground track
    fig, a = plt.subplots(figsize=(5.5, 5.5))
    for label, o in runs.items():
        a.plot(o["x"][:, 1], o["x"][:, 0], label=label)
    a.set(xlabel="East [m]", ylabel="North [m]", title="Ground track"); a.axis("equal"); a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(ROOT / "results/phase5_track.png", dpi=120)
    np.savez(ROOT / "results/phase5_mission_nowind.npz", t=c["t"], x=c["x"], cmd=c["cmd"])


if __name__ == "__main__":
    main()
