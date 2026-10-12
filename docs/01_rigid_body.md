# Phase 1 — Rigid-Body Kinematics & Dynamics

Status: implemented and tested (Phase 1).
Reference: Beard & McLain (2012), Ch. 3 and Appendix B. Equation numbers are deliberately not
quoted (I cannot check the printed book); verify each result against your copy.

## 0. Assumptions and conventions

| Item | Choice |
|---|---|
| Earth | flat, non-rotating, constant g (no Coriolis/centrifugal Earth terms) |
| Vehicle | rigid, constant mass, xz-plane of symmetry (J_xy = J_yz = 0) |
| Inertial frame `i`/`v` | NED: x North, y East, z Down. Altitude h = −p_d |
| Body frame `b` | x out the nose, y out the right wing, z out the belly, origin at the CG |
| Units | SI; angles in rad |
| Rotation sequence | 3-2-1: yaw ψ about z, pitch θ about y, roll φ about x |
| Notation | `R_a^b` maps components of a vector from frame `a` into frame `b` |

State (12): `p_n, p_e, p_d` (inertial, m), `u, v, w` (body-frame velocity of the CG relative to the
inertial frame, m/s), `φ, θ, ψ` (rad), `p, q, r` (body rates relative to inertial, expressed in `b`, rad/s).
Inputs to the core model: body-frame force `f = [f_x f_y f_z]` (N, **including gravity** — supplied by
the caller) and moment `m = [l m n]` (N·m about the CG).

> Note: `u,v,w` here are the *ground-relative* body velocity. Air-relative velocity (V_a, α, β)
> appears in Phase 2 together with wind. Do not mix them.

## 1. Rotations

Three elementary passive rotations (vehicle → v1 → v2 → body):

    R_v^v1(ψ) = [ cψ  sψ  0 ; −sψ cψ 0 ; 0 0 1 ]
    R_v1^v2(θ) = [ cθ 0 −sθ ; 0 1 0 ; sθ 0 cθ ]
    R_v2^b(φ) = [ 1 0 0 ; 0 cφ sφ ; 0 −sφ cφ ]

    R_v^b = R_v2^b R_v1^v2 R_v^v1 =
      [ cθcψ                 cθsψ                 −sθ   ]
      [ sφsθcψ − cφsψ        sφsθsψ + cφcψ        sφcθ  ]
      [ cφsθcψ + sφsψ        cφsθsψ − sφcψ        cφcθ  ]

`R_b^v = (R_v^b)ᵀ`. Properties (tested): orthogonal, det = +1, so `R⁻¹ = Rᵀ`.
Sanity check: pure yaw ψ=90° maps the inertial North axis `[1 0 0]` to body `[0 −1 0]`
(the aircraft now points East, so North is on its left).

## 2. Translational kinematics

    ṗ_ned = R_b^v · [u v w]ᵀ

Expanded, `ṗ_n = (cθcψ)u + (sφsθcψ − cφsψ)v + (cφsθcψ + sφsψ)w`, etc.
Altitude: `ḣ = −ṗ_d`.

## 3. Translational dynamics — the transport theorem

Newton's law holds in an inertial frame: `m · dV/dt|_i = F`. We want components in the rotating body
frame. For any vector `x`: `dx/dt|_i = dx/dt|_b + ω_{b/i} × x`. Hence

    m ( dV/dt|_b + ω × V ) = F      with V = [u v w]ᵀ, ω = [p q r]ᵀ

    u̇ = r v − q w + f_x/m
    v̇ = p w − r u + f_y/m
    ẇ = q u − p v + f_z/m

The `ω × V` terms are not forces; they appear because the body axes rotate. Physical consequence
(used in a test, §9): flying a steady turn at yaw rate r and speed u requires a lateral body force
`f_y = m r u` just to keep `v̇ = 0` — it is the centripetal force `m u²/R`, R = u/r.

## 4. Rotational kinematics (Euler angles)

Body rates are *not* Euler-angle rates. Resolving `[φ̇ 0 0]`, `[0 θ̇ 0]`, `[0 0 ψ̇]` through the
intermediate frames into the body frame gives
`[p q r]ᵀ = [φ̇ 0 0]ᵀ + R_v2^b(φ)[0 θ̇ 0]ᵀ + R_v2^b(φ)R_v1^v2(θ)[0 0 ψ̇]ᵀ`, i.e.

    p = φ̇ − ψ̇ sθ
    q = θ̇ cφ + ψ̇ sφ cθ
    r = −θ̇ sφ + ψ̇ cφ cθ

Inverting:

    [φ̇]   [ 1   sφ tθ    cφ tθ ] [p]
    [θ̇] = [ 0   cφ      −sφ    ] [q]
    [ψ̇]   [ 0   sφ/cθ   cφ/cθ  ] [r]

### Gimbal lock
At θ = ±90°, `cθ = 0`: `tθ` and `1/cθ` are unbounded. Geometrically the roll axis (body x) lines up
with the yaw axis (vertical), so φ and ψ rotate about the *same* line: only φ ± ψ is defined and the
map (φ,θ,ψ) → R is no longer locally invertible. Numerically: at θ = 89.9°, `1/cθ ≈ 573`, so a
modest 1 rad/s body rate asks for ψ̇ of several hundred rad/s and any integrator step blows up.
This is a property of the *parameterisation*, not of the aircraft. Quaternions avoid it (§6).

## 5. Rotational dynamics

Euler's equation in the body frame (same transport-theorem argument applied to angular momentum
`h = Jω`): `J ω̇ + ω × (Jω) = m`, with the book's inertia matrix (note the minus signs on J_xz):

    J = [ J_x 0 −J_xz ; 0 J_y 0 ; −J_xz 0 J_z ]

    Jω = [ J_x p − J_xz r , J_y q , J_z r − J_xz p ]
    ω × Jω = [ qr(J_z−J_y) − J_xz pq ,
               pr(J_x−J_z) + J_xz(p²−r²) ,
               pq(J_y−J_x) + J_xz qr ]

Invert J (the xz block has determinant Γ = J_x J_z − J_xz²):

    J⁻¹ = [ J_z/Γ 0 J_xz/Γ ; 0 1/J_y 0 ; J_xz/Γ 0 J_x/Γ ]

and `ω̇ = J⁻¹ (m − ω × Jω)` gives the book's Γ form (I re-derived each coefficient; check against the book):

    ṗ = Γ1 pq − Γ2 qr + Γ3 l + Γ4 n
    q̇ = Γ5 pr − Γ6 (p² − r²) + m/J_y
    ṙ = Γ7 pq − Γ1 qr + Γ4 l + Γ8 n

    Γ1 = J_xz (J_x − J_y + J_z)/Γ        Γ5 = (J_z − J_x)/J_y
    Γ2 = (J_z(J_z − J_y) + J_xz²)/Γ      Γ6 = J_xz/J_y
    Γ3 = J_z/Γ                            Γ7 = ((J_x − J_y)J_x + J_xz²)/Γ
    Γ4 = J_xz/Γ                           Γ8 = J_x/Γ

Implementation choice: the code will use the Γ form (matches the book and the later Simulink port) **and**
the test suite will cross-check it against the generic `solve(J, m − ω×Jω)`.

Important: because J_xz ≠ 0, a pure yaw rate `ω=[0,0,r]` is **not** an equilibrium without a moment:
`ω × Jω = [0, −J_xz r², 0]`, so `q̇ = J_xz r²/J_y ≠ 0` unless `m_y = −J_xz r²` is applied. The steady-turn
test accounts for this.

## 6. Quaternions (scalar-first, e = [e₀ e_x e_y e_z]ᵀ, unit norm)

Kinematics (no singularity):

    ė = ½ Ω(ω) e,   Ω = [ 0 −p −q −r ; p 0 r −q ; q −r 0 p ; r q −p 0 ]

Ω is skew-symmetric, so `d(eᵀe)/dt = eᵀΩe = 0`: the continuous-time norm is exactly conserved.
Numerical integration breaks this at truncation-error level, so ‖e‖ drifts.

    R_v^b(e) = [ e₀²+e_x²−e_y²−e_z²   2(e_x e_y + e₀e_z)    2(e_x e_z − e₀e_y) ;
                 2(e_x e_y − e₀e_z)   e₀²−e_x²+e_y²−e_z²   2(e_y e_z + e₀e_x) ;
                 2(e_x e_z + e₀e_y)   2(e_y e_z − e₀e_x)    e₀²−e_x²−e_y²+e_z² ]

    Euler → quaternion (half-angles, c=cos, s=sin):
      e₀ = cψ cθ cφ + sψ sθ sφ
      e_x = cψ cθ sφ − sψ sθ cφ
      e_y = cψ sθ cφ + sψ cθ sφ
      e_z = sψ cθ cφ − cψ sθ sφ
    Quaternion → Euler:
      φ = atan2(2(e₀e_x + e_y e_z), e₀²+e_z²−e_x²−e_y²)
      θ = asin(2(e₀e_y − e_x e_z))
      ψ = atan2(2(e₀e_z + e_x e_y), e₀²+e_x²−e_y²−e_z²)

Drift handling: I will implement *renormalisation after each integration step* (simple, exact for the
norm). An alternative is adding `λ(1−‖e‖²)e` to `ė`. **I am not certain this term appears in the
2012 book — I will not cite it as from the book**; if you want to compare the two, I can add it as an
option and say so.

Sign/handedness risk: quaternion conventions differ between sources (scalar-first vs last, active vs
passive). The test suite checks the quaternion rotation matrix against the Euler one, which pins down
our convention.

The quaternion model has 13 states (position 3, body velocity 3, quaternion 4, rates 3).

## 7. Integration

Fixed-step RK4: `x_{k+1} = x_k + (dt/6)(k₁ + 2k₂ + 2k₃ + k₄)`, local error O(dt⁵), global O(dt⁴).
Fixed step is chosen deliberately: SIL/HIL loops run at a fixed rate with deterministic timing, which
adaptive solvers cannot give. `scipy.solve_ivp` (DOP853, rtol/atol ~1e-12) serves as the reference solution.
The inputs `f, m` are held constant within a step (zero-order hold) — this will matter in Phase 2
when control inputs and wind are updated at a fixed rate.

## 8. Parameters (revised: authors' public repo is the single source)

The printed book is not available to the project owner, so **all parameters come from the authors'
companion repository** (read from the downloaded file, never from memory):

| Item | Value |
|---|---|
| Repo | https://github.com/byu-magicc/mavsim_public |
| Commit | `10363e0eaeb16a61a7463505f8e122a9eb7e2b4e` (merge of PR #18, 2026-06-12) |
| File | `mavsim_python/parameters/aerosonde_parameters.py` |
| mass | 11 kg |
| Jx, Jy, Jz, Jxz | 0.8244, 1.135, 1.759, 0.1204 kg m^2 |
| gravity | 9.81 m/s^2 |

My earlier from-memory guess of m = 13.5 kg was **wrong** (the file says 11 kg); the inertias matched.
Caveat: this is the repo's current set (it accompanies the authors' in-progress second edition) and
could differ from printed 2012 Appendix E. The same file also contains a newer motor model
(KV, KQ, R_motor, i0, V_max, C_Q*, C_T* prop-fit coefficients) - to be reviewed in Phase 2.

## 9. Verification plan

| # | Test | What it proves |
|---|---|---|
| 1 | **Torque-free tumbling**: arbitrary initial ω, f=0, m=0. Check T = ½ωᵀJω conserved; ‖Jω‖ conserved; and **inertial angular momentum `h_i = R_b^v J ω` conserved in all three components** | rotational dynamics + attitude kinematics are mutually consistent (the inertial-vector check also catches kinematics errors the norm cannot) |
| 2 | **Intermediate-axis instability**: principal axes from `eig(J)` (not body axes, since J_xz ≠ 0 makes x,z non-principal); spin about each principal axis with a small perturbation. Growth for the intermediate axis, bounded for min/max | Euler equations correct; tolerances chosen from the eigenvalue spread |
| 3 | **Free fall** (gravity only): `p_d = p_d0 + v_d0 t + ½ g t²`, also with non-zero attitude | gravity-in-body-frame + R + kinematics |
| 4 | **Steady turn / transport term**: `ω=[0,0,r]`, `v=[u,0,0]`, `f=[0, m r u, 0]`, compensating moment `m_y = −J_xz r²` (second variant J_xz=0). Trajectory is a circle of radius u/r: `‖(p_n,p_e) − centre‖ = u/r`, ψ=rt, speed constant, v=w=0 | the `ω×V` term (negative control: set f_y=0 ⇒ v drifts and the path is not that circle) |
| 5 | Quaternion norm ‖e‖≈1 over a long spinning run; report drift with and without renormalisation | integrator/normalisation |
| 6 | Euler vs quaternion trajectories agree away from θ=±90° | two parameterisations equivalent |
| 7 | Rotation-matrix properties and Euler↔quat round trip | convention is pinned down |
| 8 | Gimbal lock: Euler-rate matrix norm → ∞ as θ→90°, quaternion path stays finite | the explanation in §4 |
| 9 | Γ-form `ω̇` equals `solve(J, m − ω×Jω)` for random ω, m | the Γ algebra in §5 |
| 10 | RK4 convergence order ≈ 4 (error ratio ≈16 per halving of dt); RK4 vs DOP853 on a manoeuvring case | integrator correctness |

## 10. Toolbox-independence (for later phases)

Control System Toolbox is not required anywhere in this project. Python uses numpy/scipy; the MATLAB
port will use core `eig`/Simulink `linmod`, with optional Control System Toolbox calls guarded by
`license('test','Control_Toolbox')`. Phase 5 gains come from the book's closed-form formulas, so they
are computed identically either way.

## Next step
After you confirm this note (and the parameter values), I write `params/aerosonde.py`, `sim/rotations.py`,
`sim/dynamics.py`, `sim/integrators.py`, then the tests above.
