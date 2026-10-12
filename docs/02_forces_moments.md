# Phase 2 — Forces & Moments, Propulsion, Actuators, Wind

Frames/units as in Phase 1 (NED inertial, body x-nose/y-right/z-down, SI, rad). All coefficients and
sources: `params/aerosonde.py`. Items **not from the authors' file** are listed in section 9.

## 1. Gravity
Inertial `[0,0,mg]` rotated to body: `f_g = m g R_v^b [0,0,1]ᵀ = m g [−sθ, cθ sφ, cθ cφ]`
(`dynamics.gravity_force_body`).

## 2. Air-relative velocity, angle of attack, sideslip
The aerodynamics depend on motion *relative to the air*, not the ground. With ground velocity `V_g^b = [u v w]`,
steady wind `w_s` (NED) and gust `w_g` (body):

    V_a^b = [u_r v_r w_r]ᵀ = [u v w]ᵀ − ( R_v^b w_s + w_g^b )
    V_a = ‖V_a^b‖,   α = atan2(w_r, u_r),   β = asin(v_r / V_a)

- α: angle between the body x-axis projection onto the xz plane and the relative wind; positive when the
  relative wind comes from *below* the nose (w_r > 0 — the aircraft sinks relative to air).
- β: positive when the relative wind comes from the aircraft's *right* (v_r > 0 — drifting right through the air).
- Wind enters only through `V_a^b`. Steady wind shifts the ground track (crab angle) but leaves α, β, V_a
  unchanged if the aircraft flies with the same air-relative velocity. Gusts are a body-frame disturbance
  (the Dryden model is applied in body axes).
- `V_a → 0` makes α, β undefined; the code floors `V_a` at 1e-3 m/s inside the nondimensional-rate terms only.

## 3. Longitudinal aerodynamics
Dynamic pressure force `q̄S = ½ρV_a²S`. Lift/drag in the stability axes:

    F_lift = q̄S [ C_L(α) + C_Lq (c/2V_a) q + C_Lδe δe ]
    F_drag = q̄S [ C_D(α) + C_Dq (c/2V_a) q + C_Dδe δe ]
    C_D(α) = C_Dp + (C_L0 + C_Lα α)² / (π e AR),   AR = b²/S

**Nonlinear lift with stall blending.** Linear `C_L0 + C_Lα α` is only valid pre-stall; past stall a flat
plate behaves like `2 sgn(α) sin²α cos α`. Blend with a logistic window

    σ(α) = (1 + e^{−M(α−α0)} + e^{M(α+α0)}) / ((1 + e^{−M(α−α0)})(1 + e^{M(α+α0)}))
    C_L(α) = (1 − σ) (C_L0 + C_Lα α) + σ · 2 sgn(α) sin²α cos α

σ ≈ 0 for |α| ≪ α0, σ = ½ exactly at |α| = α0, σ → 1 beyond. (M = 50 sets the transition sharpness.)

Rotation from stability to body axes (rotate by α about y):

    f_x = −cα F_drag + sα F_lift,     f_z = −sα F_drag − cα F_lift

Check: α = 0 gives `f_x = −D, f_z = −L` (lift is up, i.e. −z).

Pitching moment: `m = q̄S c [ C_m0 + C_mα α + C_mq (c/2V_a) q + C_mδe δe ]`.
Signs that matter: `C_mα < 0` (static stability: positive α → nose-down restoring moment),
`C_mq < 0` (pitch damping), `C_mδe < 0` (positive elevator deflection → nose-DOWN moment — and a small
positive `C_Lδe`: lift up).

## 4. Lateral aerodynamics

    f_y = q̄S [ C_Y0 + C_Yβ β + C_Yp (b/2V_a) p + C_Yr (b/2V_a) r + C_Yδa δa + C_Yδr δr ]
    l   = q̄S b [ C_ℓ0 + C_ℓβ β + C_ℓp (b/2V_a) p + C_ℓr (b/2V_a) r + C_ℓδa δa + C_ℓδr δr ]
    n   = q̄S b [ C_n0 + C_nβ β + C_np (b/2V_a) p + C_nr (b/2V_a) r + C_nδa δa + C_nδr δr ]

Signs: `C_Yβ < 0` (side force opposes sideslip); `C_nβ > 0` (weathercock stability: positive β, wind from the
right, yaws the nose right, toward the wind); `C_ℓβ < 0` (dihedral effect: right sideslip rolls the aircraft
left, away); `C_ℓp < 0` (roll damping); `C_nr < 0` (yaw damping); `C_ℓδa > 0` (positive aileron → positive
roll, right wing down).

**Symmetry property** (used as a test): reflecting the aircraft left↔right flips the sign of
`β, p, r, δa, δr` and of `f_y, l, n`, leaving `f_x, f_z, m` unchanged. All lateral terms above are odd in
those variables and the longitudinal ones are even.

## 5. Propulsion (motor + propeller)
The authors' parameter file follows the "addendum by McLain" (electric motor, propeller fit). The file contains
only the coefficients, **not the equations**; the equations below are standard DC-motor + propeller
coefficient definitions, written by me and validated by physical tests (torque balance residual, plausible
thrust). Please treat them as "derived", not "from source".

    V_in = V_max δt
    Motor:  i = (V_in − K_V Ω)/R,   Q_m = K_Q (i − i0)         (Ω = prop speed [rad/s], K_V in V s/rad)
    Prop:   n = Ω/2π,  J = V_a/(nD),   T = ρ n² D⁴ C_T(J),   Q_p = ρ n² D⁵ C_Q(J),   C(J) = C2 J² + C1 J + C0

Steady state `Q_m = Q_p`. Substituting `n, J` gives a polynomial in Ω (no division by Ω, so Ω = 0 is safe):

    a Ω² + b Ω + c = 0
    a = ρD⁵ C_Q0/(4π²),  b = ρD⁴ V_a C_Q1/(2π) + K_Q K_V/R,  c = ρD³ V_a² C_Q2 − K_Q V_in/R + K_Q i0

Positive root `Ω = (−b + √(b²−4ac))/2a`; if `c ≥ 0` the motor cannot overcome no-load current and Ω = 0.
Thrust acts along +x_b; the reaction torque of the prop on the airframe is `−Q_p` about x_b (assumes the prop
spins in the +x_b direction; **sign convention is a modelling assumption**).
Motor dynamics (rotor inertia) are ignored: Ω is the algebraic steady state, the throttle actuator lag below
is the only dynamic.

Limitation: `C_T(J), C_Q(J)` are quadratic fits valid over the J range of the prop data. At low throttle in
fast flight (J ≳ 1.2) the quadratic is extrapolated and gives large negative thrust (−22 N at δt = 0, 25 m/s).
The trim/cruise regime (δt ≈ 0.6) is inside the fit; descents at idle should be interpreted with that in mind.

## 6. Actuators
Per channel: `δ̇ = (sat(δ_cmd) − δ)/τ`. Saturating the command first guarantees the state stays inside limits.
Limits ±45° (surfaces) and 0..1 (throttle); τ = 0.05 s (surfaces), 0.10 s (throttle). Step response:
`δ(t) = δ_cmd (1 − e^{−t/τ})`.

## 7. Wind: steady + Dryden gusts
Steady wind `w_s` in NED. Gust = Dryden turbulence, body frame, from unit white noise through shaping filters
(Va0 = nominal airspeed, L = scale length, σ = intensity):

    H_u(s) = σ_u √(2Va0/L_u) / (s + Va0/L_u)
    H_v(s) = σ_v √(3Va0/L_v) (s + Va0/(√3 L_v)) / (s + Va0/L_v)²        (H_w likewise with σ_w, L_w)

Normalisation check (numerically in tests, via a Lyapunov equation): output variance = σ² for unit-intensity
white noise. The autocorrelation must match Dryden's closed forms (ξ = Va0 τ):
`R_u = σ² e^{−ξ/L}`, `R_v = σ² (1 − ξ/2L) e^{−ξ/L}`.
Discretisation: exact (Van Loan), so statistics don't depend on the sample time.
Values of L, σ: `mavsim_simulink/parameters/wind_parameters.m` of the authors' repo (commit 10363e0, see `params/aerosonde.py`).

## 8. Verification plan (tests/test_forces.py, test_propulsion.py, test_wind.py, test_actuators.py)
| Test | Checks |
|---|---|
| Air data, no wind | `V_a=‖v‖`, `α=atan2(w,u)`, `β=asin(v/V_a)` |
| Wind triangle | tail/head/cross-wind cases against hand-calculated V_a, crab angle; wind-in-NED vs body consistency; gust added in body frame |
| Zero lateral input | β=p=r=δa=δr=0 ⇒ f_y = l = n = 0 (aero) |
| Mirror symmetry | lateral outputs odd, longitudinal even under left↔right reflection (random states) |
| Sign checks | +δe → m<0; +δa → l>0; +δr → n<0 and f_y>0; +α → m<0; +β → n>0, f_y<0, l<0; damping signs |
| Lift model | σ(±α0)=½, σ even; C_L → linear below stall, → flat plate beyond; C_L continuous |
| Rotation stability→body | α=0: f_x=−D, f_z=−L; magnitude `√(fx²+fz²)=√(L²+D²)` for any α |
| Propeller | torque-balance residual ≈ 0 vs. independent `brentq` solve; monotone in throttle; thrust falls with V_a; static thrust positive; cruise thrust ≈ cruise drag at a plausible throttle |
| Actuators | analytic first-order step; saturation; state never exceeds limits; steady state = clipped command |
| Dryden | exact output variance=σ²; simulated variance and autocorrelation vs closed forms; Ts-independence |
| Plant | straight level steady flight needs no NaNs; thrust/drag numbers; Aircraft RK4 vs scipy |

## 9. Items not taken from the authors' parameter file (need your decision/verification)
1. **Dryden L, σ** (200/200/50 m; 1.06/1.06/0.7 m/s) were first typed from memory and then **confirmed against the authors' repo**
   (`mavsim_simulink/parameters/wind_parameters.m`, commit 10363e0, which also gives Va0 = 25 m/s); identical, so no results changed.
   The shaping-filter transfer functions are not in the repo: they follow the standard Dryden form and are verified by tests.
2. **Actuator time constants and limits**: project choices, not Aerosonde data.
3. **Motor/prop equations** (section 5): derived here; coefficient fit validity assumes the standard
   `T = ρn²D⁴C_T(J)` convention.
4. Constant air density (no altitude/temperature model).
