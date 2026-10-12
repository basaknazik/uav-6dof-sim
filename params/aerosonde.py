"""Aerosonde UAV parameter set -- the single source of truth for all numbers.

SOURCE (the only one used in this project)
------------------------------------------
Authors' public companion repository for Beard & McLain, *Small Unmanned
Aircraft: Theory and Practice* (Princeton, 2012):

    repo   : https://github.com/byu-magicc/mavsim_public   (BYU MAGICC lab)
    commit : 10363e0eaeb16a61a7463505f8e122a9eb7e2b4e
             (Merge pull request #18 from byu-magicc/Winter_2026, 2026-06-12)
    file   : mavsim_python/parameters/aerosonde_parameters.py
             (lines under "Physical Parameters": mass, Jx, Jy, Jz, Jxz, gravity)

The same mass/inertia values also appear in
legacy_mavsim_python/parameters/aerosonde_parameters.py at that commit.
Values were read from the downloaded file on 2026-10-12, not typed from memory.

CAVEAT: this is the repo's *current* parameter set, which accompanies the
authors' in-progress second-edition text. It may differ from the printed 2012
Appendix E (the printed book is not available to the project owner, so it
could not be cross-checked). Anything below marked ``NOT FROM SOURCE`` is
derived or a project choice and is explained at its definition.

Rules: SI units; every value carries a citation; parameters are added phase by
phase (aerodynamics, propulsion, actuators come in Phase 2).
"""

# --- Physical parameters (mavsim_python/parameters/aerosonde_parameters.py) ---
mass = 11.0         # kg      "mass = 11. #kg"
Jx = 0.8244         # kg m^2  roll  moment of inertia
Jy = 1.135          # kg m^2  pitch moment of inertia
Jz = 1.759          # kg m^2  yaw   moment of inertia
Jxz = 0.1204        # kg m^2  product of inertia (sign convention: J[0,2] = -Jxz)
gravity = 9.81      # m/s^2   "gravity = 9.81"
