# Phase 5 — Autopilot (successive loop closure)

Perfect-state feedback (no sensor/estimator models; those are an optional extension). No Control System Toolbox:
design uses numpy eigenvalue/polynomial arithmetic and `scipy.signal` only (closed-form gain formulas), so the
Python reference runs anywhere and the MATLAB port (Phase 7) uses the same formulas in core MATLAB.

## Structure (inner → outer)
```
chi_c ─► course PI ─► phi_c ─► roll PD (rate p) ──────────────► delta_a
r ────► washout ─► yaw damper ────────────────────────────────► delta_r
h_c ──► altitude PI ─► theta_c ─► pitch PD (rate q) ──────────► delta_e      (hold zone)
Va_c ─► airspeed PI ─────────────────────────────────────────► delta_t      (hold zone)
Va_c ─► airspeed PI ─► theta_c (pitch), throttle = 1 or 0  ───► climb / descend zones (|h_c − h| > 10 m)
```
Trim values are feed-forward; controller outputs are deviations. Saturations: surfaces ±45°, throttle 0–1,
bank command ±35°, pitch command ±30°. Integrators use conditional-integration anti-windup and bumpless re-initialisation
on zone changes, and zone switching has hysteresis (leave climb/descend at 80 % of the zone width).
Course error is wrapped to ±π (short way round).

## Gain design (`sim.autopilot.design_gains`)
Design-model coefficients are *read off the numerically linearised plant* (Phase 4), not typed in:
* roll: `ṗ = −a_φ1 p + a_φ2 δa` (a_φ1 = 22.6, a_φ2 = 131) → `kp = ωn²/a_φ2`, `kd = max(0, (2ζωn − a_φ1)/a_φ2)`; roll damping is already large so `kd = 0`.
* course: `χ̇ = (g/V_g) φ` → `kp_χ = 2ζωn_χ V/g`, `ki_χ = ωn_χ² V/g`, `ωn_χ = ωn_φ / W_χ` (W_χ = 10).
* pitch: `θ̈ + a_θ1 θ̇ + a_θ2 θ = a_θ3 δe` (a_θ3 < 0 ⇒ negative gains) → `kp = (ωn² − a_θ2)/a_θ3`, `kd = (2ζωn − a_θ1)/a_θ3`; DC gain `K_θDC = a_θ3 kp/(a_θ2 + a_θ3 kp)`.
* altitude: `ḣ = V θ` → `kp_h = 2ζωn_h/(K_θDC V)`, `ki_h = ωn_h²/(K_θDC V)`, `ωn_h = ωn_θ / W_h`.
* airspeed: `V̇a = −a_V1 Va + a_V2 δt − a_V3 θ` with coefficients by finite differences on the nonlinear plant at fixed α
  (a_V3 = g cosγ exactly) → throttle PI `kp = (2ζωn − a_V1)/a_V2`, `ki = ωn²/a_V2`; pitch PI (negative gains).
* yaw damper: `δr = k_r · r_washout`, k_r = smallest gain giving Dutch-roll ζ ≥ 0.7 on the 2×2 (v, r) sub-model (root scan).

### Design findings (these changed the design; they are covered by tests)
1. **The book's sideslip-hold PI on the rudder was not used.** On this airframe it reduces Dutch-roll damping
   (closed-loop minimum ζ fell with every sideslip gain tried; and the single-DOF design formula with ωn = 1 rad/s
   was unstable through the rudder actuator lag). A washed-out yaw damper (rate feedback) raises the closed-loop minimum
   damping from 0.17 to 0.41 (`test_yaw_damper_improves_dutch_roll_damping`). The rudder sign follows C_nδr < 0 (positive r → positive δr).
2. **Naive bandwidth targets had to be moderated.** With ωn_θ = 12 and W_h = 15–25 the altitude loop overshot by 55 % with 17° pitch for
   a 5 m step because K_θDC is only ≈ 0.5 (short-period stiffness dominates the pitch loop). Final: ωn_θ = 14 rad/s, W_h = 50, ζ_h = 1.0, ζ_χ = 1.3
   (linear step overshoot: course 21 %, altitude 27 %, airspeed 9 %). Targets are verified on the *full 18-state linear
   closed loop including actuator lags*, not on the 2nd-order design models.
3. The open-loop spiral instability (+0.089 s⁻¹) is stabilised by the course/roll loops (all closed-loop poles < −0.15 s⁻¹).

## Verification (tests/test_autopilot.py, 20 tests)
* design coefficients vs analytic aero expressions; gain signs follow plant gain signs; design algebra reproduces ωn, ζ and the loop separations;
* linear closed loop (`closed_loop_hold_zone`: plant + actuators + controller states): stable, min ζ > 0.25, step specs, zero steady-state error;
  gains designed at 25 m/s stay stable on trims from 18 to 30 m/s;
* **the same small steps on the full nonlinear plant match the linear closed-loop prediction within 3 % of the step size** (measured ≈ 1–2 %) — this validates
  the design analysis *and* the implementation (it exposed and fixed a start-up integrator-reset bug);
* PI anti-windup bounded and fast recovery (a naive integrator would take > 40 s), bumpless reset, angle wrapping, course command across ±180°;
* mission (altitude 100→150→120 m, course 0→90→45°, airspeed 25→28→25 m/s): without wind final errors < 0.3 m, 0.3°, 0.2 m/s; with steady wind (−4, 6, 0) m/s and Dryden gusts:
  within 1.5 m, 1.5°, 1 m/s (3 s averages); actuator limits and 40° bank / 32° pitch never exceeded; no zone chattering;
* **course vs heading**: in an 8 m/s crosswind the aircraft holds ground course to 0.2° while its heading crabs by the wind-triangle angle (−18.7°);
* climb at full throttle ≥ 4 m/s (measured ≈ 6 m/s) with airspeed held within 3 m/s using the elevator;
* infeasible altitude command (2 km) saturates without integrator windup.

Figures: `results/phase5_mission.png`, `results/phase5_track.png`.

## Limitations
Perfect state feedback; gains designed at one trim (stable 18–30 m/s but performance not scheduled); 
idle-throttle descents use the extrapolated propeller fit; windy runs start with the trim ground velocity, so a start-up airspeed transient appears in the wind plot.
