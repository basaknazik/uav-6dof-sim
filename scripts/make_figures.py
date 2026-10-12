"""Regenerate ALL result figures and the animation into results/.  Run: python scripts/make_figures.py

Takes a few minutes (mission simulations at dt = 0.01 s).
"""
import pathlib
import runpy
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sim import autopilot as ap, linearize as lin, mission, trim, visualize  # noqa: E402

D = np.deg2rad
CMDS = [(0, 25.0, 100.0, 0.0), (5, 25.0, 150.0, D(90)), (60, 28.0, 150.0, D(90)),
        (90, 28.0, 120.0, D(45)), (130, 25.0, 120.0, D(45))]


def main():
    for name in ("phase1_plots", "phase2_plots", "phase3_plots", "phase4_plots", "phase5_plots"):
        runpy.run_path(str(ROOT / "scripts" / f"{name}.py"), run_name="__main__")
        # phase4/5 scripts call main() under __main__; phase1 does as well
    t0 = trim.compute_trim(25.0)
    g, _ = ap.design_gains(lin.linearize(t0))
    o = mission.simulate_mission(t0, g, CMDS, 170.0, dt=0.01)
    fig = plt.figure(figsize=(8, 6.5))
    ax = fig.add_subplot(projection="3d")
    visualize.plot_trajectory_3d(ax, o["x"][::10], scale=60.0, n_aircraft=14, z_exaggeration=8.0)
    ax.set_title("Mission trajectory (vertical axis stretched 8x, airframe glyphs enlarged ~20x)")
    fig.tight_layout(); fig.savefig(ROOT / "results" / "trajectory_3d.png", dpi=130); plt.close(fig)
    n = visualize.make_animation(o["x"], o["t"], str(ROOT / "results" / "attitude_animation.gif"), stride=60, fps=20)
    print("animation frames:", n)


if __name__ == "__main__":
    main()
