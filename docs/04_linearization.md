# Phase 4 — Linearisation & Modal Analysis

## Method
Full nonlinear plant (12-state Euler model, forces/moments from Phase 2) is linearised numerically about a trim
`(x*, u*)` by central differences:

    A_ij = [f_i(x* + h e_j, u*) − f_i(x* − h e_j, u*)] / 2h,    B likewise,   h = 1e-6

giving `δẋ = A δx + B δu` with `δx = x − x*`, `δu = u − u*`, inputs `[δa δe δr δt]` (actual surface positions;
actuator lag is *not* part of the design model, as in the book). Position states enter nothing except the
kinematics; altitude does not affect forces (constant density), so the `pd` column is exactly zero.

Reduced models (book ordering):
* longitudinal `[u w q θ h]`, inputs `[δe δt]` (`h = −pd`);
* lateral `[v p r φ ψ]`, inputs `[δa δr]`.

For modal analysis the pure integrators `h`, `ψ` are dropped (nothing depends on them) → 4-state models.
At level trim the coupling between the two sets is only through the trim bank angle φ0 ≈ −0.03°
(`∂ẇ/∂φ = −g sinφ0`, ...), so the 8-state eigenvalues equal the decoupled 4+4 ones to 1e-4.

## Modes (Va = 25 m/s level trim; full table in `results/phase4_modes.md`)
| Mode | eigenvalue | ω_n [rad/s] | ζ | notes |
|---|---|---|---|---|
| Short period | −4.89 ± 9.87j | 11.0 | 0.44 | pitch-rate/α dominated, speed constant |
| Phugoid | −0.143 ± 0.478j | 0.50 | 0.29 | θ and u dominated, α constant; period 13 s |
| Roll subsidence | −22.4 | — | — | τ = 0.045 s, almost pure roll rate |
| Dutch roll | −1.14 ± 4.66j | 4.79 | 0.24 | β, r, p |
| Spiral | **+0.089** | — | — | **unstable**, doubling time 7.8 s, bank dominated |

The spiral mode is slightly *unstable* in this model (it is the classic, slowly diverging mode of many aircraft;
the autopilot must stabilise it). That is a model result, not a bug: it is consistent with the simulation from trim staying
at trim for >20 s only because the divergence starts from 1e-16 round-off.

Trends with speed (tests): short-period and roll rates ∝ V, phugoid frequency falls with V (opposite to a naive guess),
spiral gets slower and Dutch-roll damping barely changes.

## Verification (tests/test_linearize.py, 25 tests)
1. Jacobian code exact on linear/analytic functions; step-size independence (1e-6 vs 1e-5).
2. Entries against hand derivatives: kinematic (`∂φ̇/∂p = 1`, `∂ψ̇/∂r = cφ/cθ`), gravity (`∂u̇/∂θ = −g cθ`), control
   derivatives from the aerodynamic coefficients (`∂q̇/∂δe = q̄Sc C_mδe/J_y`, `∂ṗ/∂δa = Γ3 q̄Sb C_ℓδa + Γ4 q̄Sb C_nδa`,
   `∂v̇/∂δr = q̄S C_Yδr/m`, `∂u̇/∂δt` from an independent propulsion finite difference).
3. Classic analytical approximations: short period (2-DOF α,q) within 5 % (measured ≈1 % in ω_n, ΔZeta < 0.01); roll
   subsidence `L_p` (single DOF with Γ terms) within 5 % (≈1 %); Dutch roll (β, r) within 12 % (≈5 %); phugoid
   Lanchester `√2 g/V` within 20 % (≈10 %; Lanchester's damping estimate fails here and is not asserted).
4. Eigenvector participation identifies the modes physically (q for short period, θ/u for phugoid, p for roll, φ for spiral).
5. **Linear vs nonlinear**: step inputs on every control, 3 s horizon, error < 6 % of the peak response; and
   *first-order accuracy*: halving the input amplitude shrinks the relative error by > 1.7×. Negative controls: a 5 % corrupted `A` leaves an
   amplitude-independent error floor; a large (0.3 rad) elevator step is badly predicted.

Figures: `phase4_poles.png`, `phase4_step_compare.png`, `phase4_impulse.png`.
