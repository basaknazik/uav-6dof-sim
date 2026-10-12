# Phase 3 — Trim

## Definition
A trim is a steady flight condition: body-frame velocity, body rates, Euler angles φ, θ and the controls are
constant, so the aircraft flies a straight line, a constant-rate climb or a helix. Specified by
`(V_a, γ, R, β)` (airspeed, flight-path angle, signed turn radius, sideslip), no wind.

## Unknowns and equations
Steady helical flight ⇒ `φ̇ = θ̇ = 0`, `ψ̇ = V_a cosγ / R`. The Euler-rate relation then fixes the body rates:

    p = −ψ̇ sθ,   q = ψ̇ sφ cθ,   r = ψ̇ cφ cθ

and `(u,v,w) = V_a (cα cβ, sβ, sα cβ)`. Unknowns `z = [α, θ, φ, δa, δe, δr, δt]` (7).
Residuals (7): `u̇ = v̇ = ẇ = 0`, `ṗ = q̇ = ṙ = 0`, `ḣ − V_a sinγ = 0`, evaluated with the *full nonlinear
plant* (`sim.aircraft.plant_deriv`). Solved with `scipy.optimize.least_squares` (trust-region reflective with
actuator bounds), minimising ‖residual‖²; the solver reaches ~1e-15 so it behaves as a root solve.

**Why β is specified, not solved for.** With β free, steady straight flight has a one-parameter family:
wings level, or banked with sideslip (the lateral force from gravity balanced by side force). Fixing β = 0
(coordinated flight) selects the physically sensible member. Without it the 7×8 system is under-determined.

**Prop torque.** The propeller reaction torque is a constant roll moment, so even straight trim has a small aileron/rudder
offset and ≈ −0.03° bank to hold β = 0 (see `test_level_flight_is_nearly_wings_level_and_coordinated`).

## Expected analytical relations (used in tests)
* Level/climb, wings level, β = 0: `θ = α + γ` (from `sinγ = cα cβ sθ − (sβ sφ + sα cβ cφ) cθ`).
* Level coordinated turn: `tanφ ≈ V_a ψ̇ / g` (exact if the side-force terms were zero).
* Slower flight ⇒ larger α, more nose-up elevator; climbing/turning ⇒ more throttle.

## Verification (tests/test_trim.py)
1. Residual < 1e-10 for 8 conditions (level at 18/25/30 m/s, ±5° climb, left/right turn, climbing turn); inputs inside limits.
2. **Independent re-check with the 12-state Euler model** (different attitude parameterisation than the solver):
   accelerations, φ̇, θ̇ = 0; ψ̇ = V_a cosγ/R; ḣ = V_a sinγ; horizontal ground speed = V_a cosγ.
3. **Simulate 10 s from trim** with the full plant (RK4, actuators at trim): velocities and rates stay within 1e-8, altitude rate matches
   V_a sinγ, ground track is a circle of radius |R| (centre located from the initial course).
4. Geometry: γ = θ − α; bank-angle formula; trends vs airspeed/climb/turn; left/right mirror for longitudinal trim.
5. Infeasible requests (V_a = 60 m/s, 10 m radius) raise `RuntimeError`.

## Results
See `results/phase3_trim_map.png`: α, δe, δt versus airspeed and flight-path angle; bank angle and throttle versus turn radius.
Level trim at 25 m/s: α ≈ 2.9°, δe ≈ −0.124 rad, δt ≈ 0.77. Maximum level speed with this model ≈ 33 m/s (throttle saturates).
