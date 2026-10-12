"""Autopilot by successive loop closure (Beard & McLain Ch. 6 structure), designed from the linear models.

Loops (inner -> outer), all with perfect state feedback (no sensor/estimator models):
    lateral      roll attitude (PD, rate feedback p)  <-  course hold (PI)           delta_a
                 yaw damper (washed-out yaw-rate feedback)                            delta_r
    longitudinal pitch attitude (PD, rate feedback q) <-  altitude hold (PI)         delta_e
                 airspeed hold: throttle (PI) in the hold zone, pitch (PI) in climb/descent zones

Gains come from `design_gains(model)`: second-order design models are read off the numerically
linearised plant (sim.linearize), then gains follow from desired natural frequency / damping, with a
bandwidth separation between successive loops. Uses only numpy -- no Control System Toolbox.
Design finding (see docs/05_autopilot.md): the book's sideslip-hold PI on the rudder *reduces* the Dutch-roll
damping of this airframe in the linear closed-loop analysis, so the rudder channel is a washed-out yaw damper instead.
Sign conventions: all error terms are (command - measurement); controller outputs are deviations added
to the trim feed-forward values.
"""

from dataclasses import dataclass, field

import numpy as np

from params import aerosonde as P
from sim import actuators, linearize as lin, rotations as rot, wind as windmod
from sim.aircraft import forces_moments  # noqa: F401  (re-exported for convenience)


# ------------------------------------------------------------------ design
@dataclass(frozen=True)
class DesignSpec:
    """Design targets. Bandwidth separation W >= 5-10 between a loop and the loop inside it."""
    wn_phi: float = 8.0          # rad/s roll-attitude loop natural frequency (< actuator bandwidth 1/tau = 20)
    zeta_phi: float = 0.9
    W_chi: float = 10.0          # course loop: wn_chi = wn_phi / W_chi
    zeta_chi: float = 1.3
    wn_theta: float = 14.0       # rad/s pitch-attitude loop
    zeta_theta: float = 0.8
    W_h: float = 50.0            # altitude loop: wn_h = wn_theta / W_h
    zeta_h: float = 1.0
    wn_V: float = 0.8            # rad/s airspeed-by-throttle loop
    zeta_V: float = 0.9
    W_V2: float = 15.0           # airspeed-by-pitch loop: wn_V2 = wn_theta / W_V2
    zeta_V2: float = 0.9
    zeta_dr: float = 0.7         # target Dutch-roll damping for the yaw damper
    tau_wash: float = 1.0        # s, washout time constant (passes transients, blocks steady-turn yaw rate)
    phi_max: float = np.deg2rad(35.0)     # bank-angle command limit
    theta_max: float = np.deg2rad(30.0)   # pitch-angle command limit
    h_hold: float = 10.0         # m, half-width of the altitude hold zone


@dataclass(frozen=True)
class Gains:
    kp_phi: float; kd_phi: float
    kp_chi: float; ki_chi: float
    kr_yaw: float                      # yaw damper: delta_r = kr_yaw * r_washout
    kp_theta: float; kd_theta: float; K_theta_DC: float
    kp_h: float; ki_h: float
    kp_V: float; ki_V: float          # airspeed via throttle
    kp_V2: float; ki_V2: float        # airspeed via pitch
    spec: DesignSpec = field(default_factory=DesignSpec)


@dataclass(frozen=True)
class DesignModels:
    """Coefficients of the design transfer functions (all SI)."""
    a_phi1: float; a_phi2: float          # p_dot = -a_phi1 p + a_phi2 delta_a
    a_theta1: float; a_theta2: float; a_theta3: float   # theta_ddot + a1 theta_dot + a2 theta = a3 delta_e
    a_V1: float; a_V2: float; a_V3: float  # Va_dot = -a_V1 (Va-Va*) + a_V2 delta_t - a_V3 theta   (alpha fixed)
    Va: float


def airspeed_coefficients(t):
    """a_V1, a_V2, a_V3 of Va_dot(Va, delta_t, theta) with alpha, beta, rates held at trim, by central differences."""
    alpha = t.alpha

    def va_dot(z):
        Va, th, dt = z
        x = np.array([0, 0, -100.0, Va * np.cos(alpha), 0, Va * np.sin(alpha), 0, th, 0, 0, 0, 0])
        d = t.delta.copy()
        d[3] = dt
        xd = lin.plant_deriv_euler(x, d)
        return np.array([np.cos(alpha) * xd[3] + np.sin(alpha) * xd[5]])

    J = lin.numerical_jacobian(va_dot, [t.Va, t.theta, t.delta[3]], 1e-6)[0]
    return -J[0], J[2], -J[1]


def design_models(model):
    """Read the design-model coefficients off the linearised plant `model` (sim.linearize.LinearModel)."""
    A, B, I = model.A, model.B, lin.IDX
    Va = model.trim.Va
    a_phi1, a_phi2 = -A[I["p"], I["p"]], B[I["p"], lin.IN["da"]]
    a_th1 = -A[I["q"], I["q"]]
    a_th2 = -A[I["q"], I["w"]] * Va            # w = Va * alpha
    a_th3 = B[I["q"], lin.IN["de"]]
    a_V1, a_V2, a_V3 = airspeed_coefficients(model.trim)
    return DesignModels(a_phi1, a_phi2, a_th1, a_th2, a_th3, a_V1, a_V2, a_V3, Va)


def design_yaw_damper(model, zeta_target, k_max=5.0):
    """Smallest gain k with delta_r = k r reaching Dutch-roll damping >= zeta_target on the (v, r) sub-model
    (best available if the target is unreachable). Uses only eigenvalues of a 2x2 matrix."""
    iv, ir = lin.IDX["v"], lin.IDX["r"]
    A2 = model.A[np.ix_([iv, ir], [iv, ir])]
    b2 = model.B[[iv, ir], lin.IN["dr"]]
    best = (0.0, -1.0)
    for k in np.linspace(0.0, k_max, 1001):
        ev = np.linalg.eigvals(A2 + np.outer(b2, [0.0, k]))
        z = -ev[0].real / abs(ev[0]) if abs(ev[0].imag) > 1e-9 else 1.0   # complex pair -> its damping, two real -> overdamped
        if z >= zeta_target:
            return float(k)
        if z > best[1]:
            best = (k, z)
    return float(best[0])


def design_gains(model, spec=DesignSpec()):
    """Gains from the design models and `spec`. Returns (Gains, DesignModels)."""
    m = design_models(model)
    Va = m.Va
    # roll: phi_ddot + (a1 + a2 kd) phi_dot + a2 kp phi = a2 kp phi_c    (kd >= 0: roll damping already large)
    kp_phi = spec.wn_phi**2 / m.a_phi2
    kd_phi = max(0.0, (2 * spec.zeta_phi * spec.wn_phi - m.a_phi1) / m.a_phi2)
    # course: chi_dot = (g/Vg) phi  ->  PI with wn_chi = wn_phi / W
    wn_chi = spec.wn_phi / spec.W_chi
    kp_chi = 2 * spec.zeta_chi * wn_chi * Va / P.gravity
    ki_chi = wn_chi**2 * Va / P.gravity
    # pitch: theta_ddot + (a1 + a3 kd) theta_dot + (a2 + a3 kp) theta = a3 kp theta_c
    kp_th = (spec.wn_theta**2 - m.a_theta2) / m.a_theta3
    kd_th = (2 * spec.zeta_theta * spec.wn_theta - m.a_theta1) / m.a_theta3
    K_dc = m.a_theta3 * kp_th / (m.a_theta2 + m.a_theta3 * kp_th)
    # altitude: h_dot = Va theta  ->  PI around the pitch loop DC gain
    wn_h = spec.wn_theta / spec.W_h
    kp_h = 2 * spec.zeta_h * wn_h / (K_dc * Va)
    ki_h = wn_h**2 / (K_dc * Va)
    # airspeed by throttle: Va_dot = -a1 Va + a2 delta_t
    kp_V = (2 * spec.zeta_V * spec.wn_V - m.a_V1) / m.a_V2
    ki_V = spec.wn_V**2 / m.a_V2
    # airspeed by pitch: Va_dot = -a1 Va - a3 theta, theta = K_dc theta_c  (negative gains: pitch down to speed up)
    wn_V2 = spec.wn_theta / spec.W_V2
    kp_V2 = (m.a_V1 - 2 * spec.zeta_V2 * wn_V2) / (K_dc * m.a_V3)
    ki_V2 = -wn_V2**2 / (K_dc * m.a_V3)
    # yaw damper: scan the gain on the (v, r) sub-model for the target Dutch-roll damping
    kr = design_yaw_damper(model, spec.zeta_dr)
    g = Gains(kp_phi, kd_phi, kp_chi, ki_chi, kr, kp_th, kd_th, K_dc, kp_h, ki_h, kp_V, ki_V, kp_V2, ki_V2, spec)
    return g, m


# ------------------------------------------------------------------ run-time controller
class PI:
    """PI with conditional-integration anti-windup and bumpless re-initialisation."""

    def __init__(self, kp, ki):
        self.kp, self.ki, self.I = kp, ki, 0.0

    def reset(self, output, error, ff=0.0):
        """Choose the integrator so that the output equals `output` for the given error (bumpless transfer)."""
        self.I = output - ff - self.kp * error

    def update(self, error, dt, lo, hi, ff=0.0):
        I_new = self.I + self.ki * error * dt
        u_unsat = ff + self.kp * error + I_new
        if lo <= u_unsat <= hi or (u_unsat > hi and self.ki * error < 0) or (u_unsat < lo and self.ki * error > 0):
            self.I = I_new
        return float(np.clip(ff + self.kp * error + self.I, lo, hi))


def wrap_pi(a):
    """Wrap angle [rad] to (-pi, pi]."""
    return (a + np.pi) % (2 * np.pi) - np.pi


@dataclass
class Measurements:
    """Perfect-state feedback signals (SI, rad). chi = ground course, h = altitude, Va/beta air-relative."""
    phi: float; theta: float; psi: float
    p: float; q: float; r: float
    h: float; Va: float; beta: float; chi: float


def measurements(x13, wind_ned=(0, 0, 0), gust_b=(0, 0, 0)):
    """Build `Measurements` from the 13-state quaternion model state (and the wind, for Va and beta)."""
    R = rot.quat_to_R(x13[6:10])
    phi, theta, psi = rot.quat_to_euler(x13[6:10])
    v_ned = R.T @ x13[3:6]
    Va, _, beta, _ = windmod.airdata(x13[3:6], R, wind_ned, gust_b)
    return Measurements(phi, theta, psi, *x13[10:13], -x13[2], Va, beta, np.arctan2(v_ned[1], v_ned[0]))


class Autopilot:
    """Discrete-time autopilot. update(cmd, meas) -> [da, de, dr, dt] commands (before actuator dynamics).

    cmd = (Va_c [m/s], h_c [m], chi_c [rad]). Zones: |h_c - h| <= h_hold -> hold (altitude by pitch, speed by
    throttle); above/below -> climb (throttle max) / descend (throttle min) with speed by pitch.
    """

    def __init__(self, gains, trim, Ts):
        self.g, self.Ts, self.trim = gains, Ts, trim
        s = gains.spec
        self.chi_pi = PI(gains.kp_chi, gains.ki_chi)
        self.r_wash = 0.0   # washout filter state
        self.h_pi, self.V_thr_pi, self.V_pitch_pi = PI(gains.kp_h, gains.ki_h), PI(gains.kp_V, gains.ki_V), PI(gains.kp_V2, gains.ki_V2)
        self.zone = "hold"
        self.theta_c_prev, self.dt_prev = trim.theta, trim.delta[3]
        self._first = True
        self.limits = actuators.LIMITS_HI
        self.spec = s

    def _zone(self, h_err):
        h = self.spec.h_hold
        if self.zone == "hold":
            if h_err > h:
                return "climb"
            if h_err < -h:
                return "descend"
        elif self.zone == "climb" and h_err < 0.8 * h:      # hysteresis: leave climb when within 80 % of the zone
            return "hold"
        elif self.zone == "descend" and h_err > -0.8 * h:
            return "hold"
        return self.zone

    def update(self, cmd, m):
        Va_c, h_c, chi_c = cmd
        dt, tr, g, s = self.Ts, self.trim, self.g, self.spec
        # ---- lateral
        chi_err = wrap_pi(chi_c - m.chi)
        phi_c = self.chi_pi.update(chi_err, dt, -s.phi_max, s.phi_max)
        d_a = tr.delta[0] + g.kp_phi * (phi_c - m.phi) - g.kd_phi * m.p
        # yaw damper: r_hp = r - x_w, x_w' = (r - x_w)/tau_w (first-order high-pass); exact discrete update
        r_hp = m.r - self.r_wash
        self.r_wash += (m.r - self.r_wash) * (1.0 - np.exp(-dt / s.tau_wash))
        d_r = tr.delta[2] + g.kr_yaw * r_hp
        # ---- longitudinal
        h_err = h_c - m.h
        zone = self._zone(h_err)
        if zone != self.zone and not self._first:
            self._switch_zone(zone, Va_c, h_c, m)
        self.zone = zone
        if zone == "hold":
            theta_c = self.h_pi.update(h_err, dt, -s.theta_max, s.theta_max, ff=tr.theta)
            d_t = self.V_thr_pi.update(Va_c - m.Va, dt, 0.0, 1.0, ff=tr.delta[3])
        else:
            theta_c = self.V_pitch_pi.update(Va_c - m.Va, dt, -s.theta_max, s.theta_max, ff=tr.theta)
            d_t = 1.0 if zone == "climb" else 0.0
        d_e = tr.delta[1] + g.kp_theta * (theta_c - m.theta) - g.kd_theta * m.q
        self.theta_c_prev, self._first = theta_c, False
        return np.clip(np.array([d_a, d_e, d_r, d_t]), actuators.LIMITS_LO, actuators.LIMITS_HI)

    def _switch_zone(self, zone, Va_c, h_c, m):
        """Bumpless transfer: re-initialise the integrators of the newly active loops to the previous pitch/throttle command."""
        tr = self.trim
        if zone == "hold":
            self.h_pi.reset(self.theta_c_prev, h_c - m.h, ff=tr.theta)
            self.V_thr_pi.reset(1.0 if self.zone == "climb" else (0.0 if self.zone == "descend" else tr.delta[3]), Va_c - m.Va, ff=tr.delta[3])
        else:
            self.V_pitch_pi.reset(self.theta_c_prev, Va_c - m.Va, ff=tr.theta)


# ------------------------------------------------------------------ linear closed-loop analysis (hold zone)
_X10 = ["pd", "u", "v", "w", "phi", "theta", "psi", "p", "q", "r"]
Z_NAMES = _X10 + ["d_a", "d_e", "d_r", "d_t", "xi_chi", "x_wash", "xi_h", "xi_V"]


def output_jacobians(model):
    """(4, 12) Jacobians of the outputs [course chi, sideslip beta, altitude h, airspeed Va] about the trim."""
    def outs(x12):
        v_ned = rot.euler_to_R(*x12[6:9]).T @ x12[3:6]
        Va = np.linalg.norm(x12[3:6])
        return np.array([np.arctan2(v_ned[1], v_ned[0]), np.arcsin(x12[4] / Va), -x12[2], Va])
    return lin.numerical_jacobian(outs, model.x_trim, 1e-6)


def closed_loop_hold_zone(model, gains):
    """Linear closed loop in the hold zone (no saturation), including the actuator lags.

    z = [x10 (pd u v w phi theta psi p q r), actuators (da de dr dt), xi_chi, x_wash, xi_h, xi_V]   (see Z_NAMES)
    inputs = [chi_c, h_c, Va_c] (deviations from trim).
    Returns (Acl, Bcl, C); C rows = [chi, h, Va, beta, phi, theta, da, de, dr, dt].
    """
    g, n = gains, 10
    idx = [lin.IDX[s] for s in _X10]
    A, B = model.A[np.ix_(idx, idx)], model.B[idx, :]
    Jchi, Jbeta, Jh, JVa = output_jacobians(model)[:, idx]
    E = np.eye(n)
    phi_r, th_r, p_r, q_r, r_r = (E[_X10.index(k)] for k in ("phi", "theta", "p", "q", "r"))
    N = n + 8
    iD, iX = n, n + 4
    Acl, Bcl = np.zeros((N, N)), np.zeros((N, 3))
    Acl[:n, :n], Acl[:n, iD:iD + 4] = A, B

    Cz, Cr = np.zeros((4, N)), np.zeros((4, 3))          # actuator commands = Cz z + Cr refs
    # course PI -> phi_c -> roll PD
    phi_c_x = -g.kp_chi * Jchi
    Cz[0, :n] = g.kp_phi * (phi_c_x - phi_r) - g.kd_phi * p_r
    Cz[0, iX + 0] = g.kp_phi * g.ki_chi
    Cr[0, 0] = g.kp_phi * g.kp_chi
    # altitude PI -> theta_c -> pitch PD
    th_c_x = -g.kp_h * Jh
    Cz[1, :n] = g.kp_theta * (th_c_x - th_r) - g.kd_theta * q_r
    Cz[1, iX + 2] = g.kp_theta * g.ki_h
    Cr[1, 1] = g.kp_theta * g.kp_h
    # yaw damper: r_hp = r - x_wash
    Cz[2, :n] = g.kr_yaw * r_r
    Cz[2, iX + 1] = -g.kr_yaw
    # throttle PI on airspeed
    Cz[3, :n] = -g.kp_V * JVa
    Cz[3, iX + 3] = g.ki_V
    Cr[3, 2] = g.kp_V
    for k in range(4):                                    # first-order actuators
        Acl[iD + k, :] = Cz[k] / actuators.TAU[k]
        Acl[iD + k, iD + k] -= 1.0 / actuators.TAU[k]
        Bcl[iD + k, :] = Cr[k] / actuators.TAU[k]
    Acl[iX + 0, :n], Bcl[iX + 0, 0] = -Jchi, 1.0          # xi_chi' = chi_c - chi
    Acl[iX + 1, :n], Acl[iX + 1, iX + 1] = r_r / g.spec.tau_wash, -1.0 / g.spec.tau_wash   # x_w' = (r - x_w)/tau_w
    Acl[iX + 2, :n], Bcl[iX + 2, 1] = -Jh, 1.0
    Acl[iX + 3, :n], Bcl[iX + 3, 2] = -JVa, 1.0
    C = np.zeros((10, N))
    C[0, :n], C[1, :n], C[2, :n], C[3, :n], C[4, :n], C[5, :n] = Jchi, Jh, JVa, Jbeta, phi_r, th_r
    for k in range(4):
        C[6 + k, iD + k] = 1.0
    return Acl, Bcl, C
