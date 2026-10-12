# Interview prep — questions with model answers

One section per phase, appended as the project grows. Answers are written to be spoken, then
backed by the equation/test that proves it. File references point at this repo.

---

## Phase 1 — Rigid-body kinematics & dynamics

**Q1. Where does the `−ω×V` term in `u̇ = rv − qw + f_x/m` come from? Why does a steady turn need a lateral force `m·r·u` even though the body-frame velocity is constant?**
Newton's law holds in an inertial frame, but we write the velocity in a rotating (body) frame. The
transport theorem gives `dV/dt|_i = dV/dt|_b + ω×V`, so `m(V̇_b + ω×V) = F`. The `ω×V` term is not a
force; it is the rotation of the axes. In a steady turn the body-frame components `u,v,w` are
constant but the *inertial* velocity vector rotates, so a real centripetal force `m·u²/R = m·r·u`
(R = u/r) is needed. Test: `test_steady_turn_is_a_circle_of_radius_u_over_r`; its negative control
removes the force and `v` runs away.

**Q2. Why did the steady-turn test need a pitching moment when Jxz ≠ 0?**
`ω×Jω` for `ω=[0,0,r]` is `[0, −Jxz r², 0]`: with a product of inertia, a pure yaw rate gives angular
momentum with an x-component, and rotating that vector demands a pitch moment. Without it `q̇ = Jxz r²/Jy ≠ 0`
and the aircraft pitches (a real aircraft has the elevator/aero moment to balance this).

**Q3. Why is ‖Jω‖ conserved in torque-free motion but the body-frame vector Jω is not? What is conserved inertially?**
Angular momentum `h` is constant in the *inertial* frame. The body frame rotates, so its components
`Jω` change while the magnitude is fixed. `h_i = R_b→v J ω` must be constant in all three components,
which validates the attitude kinematics too: a wrong rotation (e.g. R instead of Rᵀ) makes it drift
(`test_inertial_angular_momentum_check_detects_wrong_kinematics`). Kinetic energy `½ωᵀJω` is also conserved.

**Q4. Explain gimbal lock physically. Why don't quaternions have it, and why still renormalize?**
At θ = ±90° the roll axis (body x) aligns with the yaw axis (vertical), so φ and ψ rotate about the
same line and only φ±ψ is defined; the rate map has `tan θ` and `1/cos θ` and blows up. It is a flaw of the
3-parameter *representation*, not the aircraft. Quaternions are a 4-parameter, singularity-free cover of
SO(3) (kinematics `ė = ½Ω(ω)e`, bounded). Numerical integration breaks `‖e‖=1` slowly (7e-6 drift in
120 s of fast spin in our test), so we renormalize each step.

**Q5. Why fixed-step RK4 in the sim loop instead of an adaptive solver? What does the convergence plot tell you?**
SIL/HIL runs at a fixed rate with deterministic computation time; adaptive step sizes break real-time
scheduling and make controller sampling inconsistent. RK4 is 4th order: halving dt cuts error ≈16× (we
measured 18.5 and 17.7), so you can choose dt from an error budget. `solve_ivp` (DOP853, 1e-12) is
used only as ground truth in tests.

**Q6. (Bonus) Why do we derive the Γ form instead of just inverting J numerically?**
Closed form with Γ1–Γ8 is cheap, explicit, matches the book and the Simulink port, and exposes the
Jxz coupling between roll and yaw. We verify it against `solve(J, m − ω×Jω)` to machine precision.
