"""Closed-loop mission simulation: full nonlinear plant + actuators + wind/gusts + autopilot."""

import numpy as np

from sim import aircraft, autopilot as ap, wind as windmod


def simulate_mission(trim, gains, commands, T, dt=0.01, wind_ned=(0.0, 0.0, 0.0), gusts=False, seed=0, x0=None):
    """Run the autopilot on the nonlinear plant.

    commands: list of (t_start, Va_c [m/s], h_c [m], chi_c [rad]) steps, sorted by t_start.
    wind_ned: steady wind in NED [m/s]; gusts: add Dryden turbulence (body frame) if True.
    Returns dict of arrays: t, x (N, 17), cmd (N, 3), delta_cmd (N, 4), zone (N,), meas fields Va, h, chi, beta, ...
    """
    x13 = trim.x13.copy() if x0 is None else x0
    plant = aircraft.Aircraft(x13, trim.delta)
    pilot = ap.Autopilot(gains, trim, dt)
    gust = windmod.DrydenGust(dt, trim.Va, np.random.default_rng(seed)) if gusts else None
    n = int(round(T / dt))
    keys = ("Va", "h", "chi", "beta", "phi", "theta", "psi", "p", "q", "r")
    out = {"t": np.arange(n + 1) * dt, "x": np.empty((n + 1, 17)), "cmd": np.empty((n + 1, 3)),
           "delta_cmd": np.empty((n + 1, 4)), "zone": [], **{k: np.empty(n + 1) for k in keys}}
    cmd_i = 0
    w = np.asarray(wind_ned, float)
    g_b = np.zeros(3)
    for k in range(n + 1):
        t = k * dt
        while cmd_i + 1 < len(commands) and t >= commands[cmd_i + 1][0]:
            cmd_i += 1
        c = commands[cmd_i][1:]
        m = ap.measurements(plant.x[:13], w, g_b)
        d = pilot.update(c, m)
        out["x"][k], out["cmd"][k], out["delta_cmd"][k] = plant.x, c, d
        out["zone"].append(pilot.zone)
        for key in keys:
            out[key][k] = getattr(m, key)
        if k < n:
            g_b = gust.step() if gust is not None else g_b
            plant.step(dt, d, w, g_b)
    out["zone"] = np.array(out["zone"])
    return out
