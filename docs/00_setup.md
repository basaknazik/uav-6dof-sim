# Phase 0 — Environment & Tooling

## Toolchain (checked 2026-10-12, Windows 11)

| Tool | Version | Notes |
|---|---|---|
| Python | 3.14.6 | Reference implementation; venv in `.venv/` |
| numpy / scipy / matplotlib | 2.5.4 / 1.18.1 / 3.11.2 | pinned in `requirements.txt` |
| pytest | 9.1.1 | config in `pyproject.toml` |
| MATLAB | R2026a (26.1.0.3203278) | `matlab -batch "version"` works from the CLI |
| git | 2.54.0 | |

### Installed MATLAB products
MATLAB, Simulink, Simscape, Simscape Multibody, Simscape Driveline,
Simscape Electrical, Simscape Fluids, Simscape Battery.

**Not installed** (matters for Phase 7):
- *Aerospace Blockset / Aerospace Toolbox* — no ready-made 6DOF, atmosphere or
  Dryden blocks. We build the equations of motion from basic Simulink blocks
  (MATLAB Function blocks / integrators). This is arguably better for the
  portfolio: every equation is visible, and it maps 1:1 to the Python code.
- *Control System Toolbox* — no `ss`, `eig`-based `damp`, `linmod` helpers
  beyond core MATLAB. Core `linmod` ships with Simulink; `eig` is core MATLAB,
  so linear analysis in MATLAB is still possible.

## Conventions used throughout
- Frames: inertial NED (flat, non-rotating Earth), body frame x-nose / y-right-wing / z-down.
- Units: SI everywhere (m, s, kg, rad). Degrees only at plotting/IO boundaries.
- Notation: Beard & McLain (2012).

## How to run
```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Windows
.venv/Scripts/python -m pytest
```
