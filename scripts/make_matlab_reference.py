"""Generate Python reference data for the MATLAB/Simulink validation -> matlab/validation/data/python_reference.mat

Run: python scripts/make_matlab_reference.py     (about 2 minutes)
Contents:
  deriv_*      random states/inputs and the Python 17-state derivative (tests the MATLAB plant equations)
  trim_*, A, B Python trim, linearisation and autopilot gain vector at Va = 25 m/s (tests the MATLAB design chain)
  ol_*         open-loop run (doublets, saturation, steady wind + Dryden gusts), RK4 dt=0.01, NO quaternion renormalisation
  cl_nom_*, cl_wind_*  closed-loop mission runs without / with wind+gusts (no renormalisation), plus *_norm: with renormalisation
"""
import pathlib
import sys

import numpy as np
from scipy.io import savemat

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from params import aerosonde as P  # noqa: E402
from sim import actuators, aircraft, autopilot as ap, linearize as lin, mission, rotations as rot, trim, wind  # noqa: E402

D = np.deg2rad
DT = 0.01
ZONE = {"hold": 0, "climb": 1, "descend": 2}
CMDS = [(0, 25.0, 100.0, 0.0), (5, 25.0, 150.0, D(90)), (60, 28.0, 150.0, D(90)),
        (90, 28.0, 120.0, D(45)), (130, 25.0, 120.0, D(45))]


def gain_vector(g, t):
    """Packing shared with matlab/uav_autopilot.m (1-based indices in MATLAB)."""
    s = g.spec
    return np.array([g.kp_phi, g.kd_phi, g.kp_chi, g.ki_chi, g.kr_yaw, g.kp_theta, g.kd_theta, g.kp_h, g.ki_h,
                     g.kp_V, g.ki_V, g.kp_V2, g.ki_V2, s.phi_max, s.theta_max, s.h_hold, s.tau_wash,
                     t.theta, *t.delta])


def command_array(cmds, n):
    out = np.empty((n + 1, 3))
    i = 0
    for k in range(n + 1):
        while i + 1 < len(cmds) and k * DT >= cmds[i + 1][0]:
            i += 1
        out[k] = cmds[i][1:]
    return out


def derivative_cases(t, n=200, seed=11):
    rng = np.random.default_rng(seed)
    X, U, W, G, XD = [], [], [], [], []
    ac = aircraft.Aircraft(t.x13, t.delta)
    for _ in range(n):
        ang = (rng.uniform(-0.8, 0.8), rng.uniform(-0.6, 0.6), rng.uniform(-3, 3))
        x13 = np.concatenate([rng.normal(size=3) * 100, [rng.uniform(15, 35), rng.normal() * 3, rng.normal() * 4],
                              rot.euler_to_quat(*ang), rng.normal(size=3) * 0.6])
        delta = rng.uniform(actuators.LIMITS_LO, actuators.LIMITS_HI)
        cmd = rng.uniform(1.4 * actuators.LIMITS_LO, 1.4 * actuators.LIMITS_HI)   # beyond limits sometimes
        w, gu = rng.normal(size=3) * 4, rng.normal(size=3) * 1.5
        x17 = np.concatenate([x13, delta])
        X.append(x17); U.append(cmd); W.append(w); G.append(gu)
        XD.append(ac.deriv(x17, cmd, w, gu))
    return dict(deriv_x=np.array(X), deriv_cmd=np.array(U), deriv_wind=np.array(W), deriv_gust=np.array(G), deriv_xdot=np.array(XD))


def open_loop_case(t, T=30.0, seed=3):
    n = int(round(T / DT))
    tt = DT * np.arange(n + 1)
    cmd = np.tile(t.delta, (n + 1, 1))
    cmd[(tt >= 1) & (tt < 2), 0] += 0.05
    cmd[(tt >= 2) & (tt < 3), 0] -= 0.05
    cmd[(tt >= 5) & (tt < 5.5), 1] += 0.04
    cmd[(tt >= 5.5) & (tt < 6), 1] -= 0.04
    cmd[(tt >= 10) & (tt < 11), 2] += 0.05
    cmd[(tt >= 11) & (tt < 12), 2] -= 0.05
    cmd[tt >= 15, 3] += 0.1
    cmd[(tt >= 20) & (tt < 20.4), 1] = 1.2            # beyond the surface limit: exercises saturation
    cmd[(tt >= 22) & (tt < 23), 3] = 1.4
    cmd[(tt >= 24) & (tt < 25), 3] = -0.2
    w = np.array([3.0, -2.0, 0.5])
    g = wind.DrydenGust(DT, 25.0, np.random.default_rng(seed))
    gust = np.array([g.step() for _ in range(n + 1)])
    ac = aircraft.Aircraft(t.x13, t.delta)
    X = np.empty((n + 1, 17))
    X[0] = ac.x
    for k in range(n):
        ac.step(DT, cmd[k], w, gust[k], normalize=False)
        X[k + 1] = ac.x
    return dict(ol_t=tt, ol_cmd=cmd, ol_wind=w, ol_gust=gust, ol_x=X)


def closed_loop_case(prefix, t, g, **kw):
    o = mission.simulate_mission(t, g, CMDS, 170.0, dt=DT, normalize=False, **kw)
    on = mission.simulate_mission(t, g, CMDS, 170.0, dt=DT, normalize=True, **kw)
    n = len(o["t"]) - 1
    return {f"{prefix}_t": o["t"], f"{prefix}_cmd": command_array(CMDS, n), f"{prefix}_gust": o["gust"],
            f"{prefix}_x": o["x"], f"{prefix}_delta_cmd": o["delta_cmd"],
            f"{prefix}_zone": np.array([ZONE[z] for z in o["zone"]]), f"{prefix}_x_norm": on["x"]}


def main():
    out_dir = ROOT / "matlab" / "validation" / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    t = trim.compute_trim(25.0)
    model = lin.linearize(t)
    g, dm = ap.design_gains(model)
    data = dict(trim_x13=t.x13, trim_delta=t.delta, trim_alpha=t.alpha, trim_theta=t.theta, trim_phi=t.phi,
                A=model.A, B=model.B, gains=gain_vector(g, t), dt=DT, Va_trim=25.0,
                design_models=np.array([dm.a_phi1, dm.a_phi2, dm.a_theta1, dm.a_theta2, dm.a_theta3, dm.a_V1, dm.a_V2, dm.a_V3]))
    data.update(derivative_cases(t))
    data.update(open_loop_case(t))
    data.update(closed_loop_case("cl_nom", t, g))
    data.update(closed_loop_case("cl_wind", t, g, wind_ned=(-4.0, 6.0, 0.0), gusts=True, seed=1))
    data["cl_wind_wind"] = np.array([-4.0, 6.0, 0.0])
    savemat(out_dir / "python_reference.mat", data, do_compression=True)
    print("wrote", out_dir / "python_reference.mat")


if __name__ == "__main__":
    main()
