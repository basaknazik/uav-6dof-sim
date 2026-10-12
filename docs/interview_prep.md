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

---

## Phase 2 — Forces, moments, propulsion, actuators, wind

**Q1. Define V_a, α and β precisely. Why must the aerodynamics use air-relative velocity? How does wind enter?**
`V_a^b = v_ground^b − (R_v^b w_steady + w_gust^b)`; `V_a=‖V_a^b‖`, `α=atan2(w_r,u_r)`, `β=asin(v_r/V_a)`.
Lift/drag depend on the flow over the wing, i.e. motion relative to the air mass, not the ground. Wind
only changes `V_a^b`. Verified by the **Galilean-invariance test**: a uniform steady wind must leave
attitude, rates, α, β, V_a unchanged and only shift the ground track by `w·t` (`test_galilean_invariance_steady_wind`),
with a negative control that fails when the wind is ignored.

**Q2. Why blend the lift curve with a sigmoid? What does σ(α) do?**
Linear `C_L0 + C_Lα α` stays valid only until stall; beyond it lift collapses and approaches the flat-plate
value `2 sgn(α) sin²α cos α`. σ is a logistic window equal to ½ at |α|=α0 and ≈0/1 away from it; `C_L=(1−σ)·linear+σ·flat plate`.
It keeps C_L smooth and differentiable (needed for trim/linearisation). M sets the sharpness. Tests: σ(±α0)=½, σ even,
C_L matches each regime, no jumps.

**Q3. Give three sign conventions you verified and why each is physically right.**
`C_mα<0`: positive α produces nose-down moment (static stability). `C_nβ>0`: wind from the right yaws the nose
right, toward the wind (weathercock). `C_mδe<0`: positive elevator → nose-down. Also `C_ℓp, C_mq, C_nr <0` are damping.
Tests apply finite increments and check the sign of the resulting force/moment.

**Q4. What symmetry does the aerodynamic model have and how do you test it?**
Left-right mirror: `β,p,r,δa,δr → −` flips the signs of `f_y,l,n` and leaves `f_x,f_z,m` unchanged. Longitudinal
terms are even, lateral odd in those variables. Test on 100 random states. Break symmetry intentionally with the
prop torque (it only enters `l`), which is why the test targets aero only.

**Q5. Walk through the motor/propeller model. How do you get Ω, thrust and the quadratic?**
`V_in=V_max δt`; motor torque `K_Q((V_in−K_VΩ)/R−i0)`; prop `T=ρn²D⁴C_T(J)`, `Q=ρn²D⁵C_Q(J)`, `J=Va/(nD)`.
Equate torques, substitute the quadratic C_Q(J): `aΩ²+bΩ+c=0` with `a=ρD⁵C_Q0/4π²`, `b=ρD⁴VaC_Q1/2π+K_QK_V/R`,
`c=ρD³Va²C_Q2−K_QV_in/R+K_Qi0`; take the positive root. Checked against an independent root-finder using the raw J definition,
plus efficiency `T·Va/(V·i)<1`. Weakness: C_T(J) is extrapolated at idle (large negative thrust at 25 m/s, δt=0).

**Q6. Why saturate the command before the first-order lag instead of clipping the output?**
Clipping the output after the integrator allows integrator wind-up of the actuator state past its limit (it keeps
"charging"); saturating the command guarantees the state converges to the limit and never exceeds it. Step response
`δ(t)=δ_cmd(1−e^{−t/τ})` checked against the analytic solution.

**Q7. Explain the Dryden gust model and how you verified it.**
Unit white noise through shaping filters: `H_u=σ√(2V/L)/(s+V/L)`, `H_{v,w}=σ√(3V/L)(s+V/(√3L))/(s+V/L)²`, giving the
Dryden spectra of MIL-F-8785C. Verified: output variance (Lyapunov equation) is exactly σ²; autocorrelation equals
the closed forms `σ²e^{−ξ/L}` and `σ²(1−ξ/2L)e^{−ξ/L}`; the discrete implementation (Van Loan exact discretisation) is stationary
for any Ts; long simulation reproduces variance and correlation. L and σ values are flagged unverified.
