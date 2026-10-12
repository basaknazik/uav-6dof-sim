# 6-DOF Fixed-Wing UAV Simulation & Autopilot

A from-scratch nonlinear 6-DOF simulation of the Aerosonde small UAV, with trim,
linearization/modal analysis, a successive-loop-closure autopilot, and a
MATLAB/Simulink port validated against the Python reference implementation.

Notation, frames and parameters follow Beard & McLain, *Small Unmanned Aircraft:
Theory and Practice* (Princeton, 2012).

> Work in progress — built phase by phase. Each phase has a derivation note in
> [`docs/`](docs/) and physics-based verification tests in [`tests/`](tests/).

## Status

| Phase | Topic | Status |
|---|---|---|
| 0 | Setup, tooling, test harness | done |
| 1 | Rigid-body kinematics & dynamics (Euler + quaternion, RK4) | planned |
| 2 | Forces & moments: aero, propulsion, actuators, wind/Dryden | planned |
| 3 | Trim (level, climb, coordinated turn) | planned |
| 4 | Linearization & modal analysis | planned |
| 5 | Autopilot (successive loop closure) | planned |
| 6 | Visualization & polish | planned |
| 7 | MATLAB/Simulink port + cross-validation | planned |

## Repository layout
```
docs/      derivation note per phase
params/    aerosonde.py — every parameter, with its source
sim/       dynamics, forces, wind, integrator, trim, linearize, autopilot
tests/     physics verification tests (pytest)
matlab/    Simulink model built programmatically from .m scripts (Phase 7)
results/   generated plots and animations
```

## How to run
```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Windows (use .venv/bin on Linux/macOS)
.venv/Scripts/python -m pytest
```
