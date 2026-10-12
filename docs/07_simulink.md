# Phase 7 — MATLAB/Simulink port and cross-validation

MATLAB R2026a with Simulink (no Aerospace Blockset, no Control System Toolbox, no Optimization Toolbox; see docs/00_setup.md).
Everything below uses only core MATLAB and basic Simulink blocks.

## What was built (`matlab/`)
| File | Role |
|---|---|
| `uav_params.m` | **Generated** from `params/aerosonde.py` by `scripts/export_matlab_params.py` (bit-identical values; checked by `tests/test_matlab_params.py`) |
| `uav_quat_to_R/euler_to_R/euler_to_quat/quat_to_euler.m` | attitude utilities |
| `uav_airdata.m`, `uav_aero.m`, `uav_prop.m`, `uav_gravity.m` | air data, aerodynamics (stall-blended lift), motor/propeller, gravity |
| `uav_rigid_body.m` | 13-state equations of motion (transport theorem, quaternion kinematics, Γ-form rotational dynamics) |
| `uav_forces_moments.m`, `uav_plant_deriv.m`, `uav_plant_deriv_euler.m` | assembled plant (17 states incl. actuators; 12-state Euler version for trim/linearisation) |
| `uav_trim.m`, `uav_linearize.m`, `uav_design_gains.m` | trim (damped Gauss-Newton), central-difference linearisation, autopilot gain design — all core MATLAB |
| `uav_autopilot.m` | the autopilot (persistent PI states, zones, bumpless transfer, yaw damper) |
| `build_uav_models.m` | builds `uav_open_loop.slx` and `uav_closed_loop.slx` **from a script** (no `.slx` is committed: the script is the diffable source) |
| `validation/validate_functions.m`, `validation/run_validation.m` | the cross-validation against Python |

Model structure (pictures: `results/simulink_closed_loop.png`, `results/simulink_plant.png`):
```
Mission commands ─► Autopilot (MATLAB Function, Ts = 0.01 s) ─► UAV Plant ─► x17 ─► ZOH ─┐
                         ▲                                         ▲                       │
                         └─────────────── ZOH(air data) ◄──────────┘  steady wind, gust    │
UAV Plant = Actuators (saturation + 1/τ + integrator) → Air data → Aerodynamics / Propulsion / Gravity
            → Σ force, Σ moment → Rigid-body 6-DOF → Integrator (13 states)
```
Equation code lives in plain `uav_*.m` files that the MATLAB Function blocks call, so the same code is unit-testable in MATLAB, diffable in git and
code-generation compatible (`%#codegen`). Fixed-step `ode4` (classical RK4) at 0.01 s, the same integrator and step as the Python reference.

## Cross-validation (`results/phase7_validation.md`, `.json`, `.png`)
Same inputs, same step, **no quaternion renormalisation on the Python side** (Simulink integrates the quaternion without renormalising; the Python simulator can
switch it off with `normalize=False`).

| check | result |
|---|---|
| 17-state derivative at 200 random states/inputs (incl. over-limit commands) | max abs error 8.5e-14 |
| trim at 25 m/s (core-MATLAB Gauss-Newton vs scipy `least_squares`) | α, θ, controls agree to 1e-16 |
| linearisation A / B (central differences) | 1.8e-9 / 8.1e-11; eigenvalues 6.8e-10 |
| autopilot gains designed in MATLAB vs Python | relative 4.4e-10 |
| Simulink open loop (doublets, saturation, steady wind, Dryden gusts), 30 s | max abs state error **6.8e-13** |
| Simulink closed loop, no wind, 170 s mission | **2.7e-12**; autopilot commands 8e-15; zone sequence identical |
| Simulink closed loop, wind + gusts, 170 s mission | **5.0e-12**; commands 1.2e-14; zone sequence identical |
| versus Python *with* renormalisation | 2.1e-5 (no wind), 6.5e-4 (gusts), the effect of the one thing Simulink does not do |

Errors at the 1e-12 level are floating-point round-off of the same algorithm: the two implementations are numerically the same model. The remaining difference
between Simulink and the production Python simulator (6.5e-4 in the gust case) is entirely the quaternion renormalisation (drift of ‖e‖ is 1e-7 over 170 s).

### Pitfalls found (worth knowing for the interview)
* **Input timing**: RK4 evaluates 4 stages per step; to reproduce a zero-order-hold input, From Workspace blocks must be discrete (sample time = Ts) so the value is constant across the stages.
  Time stamps are shifted by −1e-9 so the row at the exact hit time `k·Ts` is selected despite floating-point round-off of the time vector.
* **The loop is sample-and-hold**: the autopilot reads the state through a ZOH at Ts and its output is held through the RK4 stages: exactly the Python loop `measure → control → step`.
* **Model-workspace constants** (`x13_0`, `delta_0`, `wind_vec`) must be overridden with `setVariable(..., 'Workspace', model)`; the default (base workspace) is silently shadowed
  (first validation run showed a 255 m error before this was fixed).
* **To Workspace** limits logging to 1000 points by default and stores vector signals with time in the *last* dimension.
* Python's mission loop used to draw the gust *after* the measurement, giving the autopilot a one-sample-old gust; it was reordered so controller and plant see the same sample (both implementations now agree).

## How this setup extends to SIL and HIL
* **MIL → SIL**: the controller is already isolated in `uav_autopilot.m` as a pure function of (command, state, air data, gains) with persistent state and `%#codegen`
  compatible code, called at a fixed 100 Hz. With MATLAB Coder / Embedded Coder (not installed here) it generates C that is compiled into a library or an S-function
  and run inside the same Simulink plant in place of the MATLAB Function block. The *same validation harness* then compares
  the SIL output with the MIL output (should match to round-off, apart from single- vs double-precision if the target is `float`).
  Because we already have a bit-level Python reference, regression after every change is automatic (`RUN_MATLAB_TESTS=1 pytest`).
* **SIL → PIL/HIL**: the plant stays in Simulink and runs in real time on a target computer (Simulink Real-Time / Speedgoat or a PC with real-time scheduling).
  Requirements this model already satisfies: fixed-step solver (`ode4`, 10 ms), no variable-step or zero-crossing dependent behaviour, no algebraic loops, discrete controller with
  explicit sample time. Still to do for HIL: replace the perfect-state feedback with sensor models (IMU, GPS, baro, pitot) and the real autopilot's I/O: servo PWM
  to actuators (our actuator model becomes the real servos or a servo-load emulator), serial/CAN messages for the sensor data, and a timing budget check (worst-case
  execution time of plant + I/O inside the 10 ms frame). Actuator lag, saturation and rate limits are already part of the plant, so the controller sees realistic command conditioning.
* **Fault injection / SIL-HIL parity**: because commands, wind and gust are logged time series fed through From Workspace blocks, the identical scenarios can be replayed in
  MIL, SIL and HIL, and results compared with the same metrics used here.

## How to run
```bash
python scripts/export_matlab_params.py          # only after changing params/aerosonde.py
python scripts/make_matlab_reference.py         # Python reference (about 2 min)
matlab -batch "addpath('matlab/validation'); run_validation"     # builds the models, simulates, writes results/phase7_*
RUN_MATLAB_TESTS=1 python -m pytest tests/test_matlab_validation.py   # the same, with assertions
```
