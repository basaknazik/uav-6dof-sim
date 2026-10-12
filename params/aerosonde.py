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

# --- Geometry / atmosphere (same upstream file, "Physical Parameters") ---
S_wing = 0.55       # m^2     wing planform area
b = 2.8956          # m       wingspan
c = 0.18994         # m       mean chord
S_prop = 0.2027     # m^2     propeller disc area (not used by the model below)
rho = 1.2682        # kg/m^3  air density (constant: no altitude dependence)
e_oswald = 0.9      # -       Oswald efficiency ("e = 0.9" upstream; renamed: e is the quaternion elsewhere)
AR = b**2 / S_wing  # -       aspect ratio (derived; upstream computes it the same way)

# --- Longitudinal coefficients (upstream "Longitudinal Coefficients") ---
C_L_0 = 0.23
C_D_0 = 0.0424
C_m_0 = 0.0135
C_L_alpha = 5.61
C_D_alpha = 0.132
C_m_alpha = -2.74
C_L_q = 7.95
C_D_q = 0.0
C_m_q = -38.21
C_L_delta_e = 0.13
C_D_delta_e = 0.0135
C_m_delta_e = -0.99
M_blend = 50.0      # stall-blending sharpness ("M = 50.0" upstream; renamed, M is a moment elsewhere)
alpha0 = 0.47       # rad     stall-blending angle of attack
epsilon = 0.16      # -       upstream value; unused by the model in this repo (listed for completeness)
C_D_p = 0.043       # -       parasitic drag coefficient

# --- Lateral coefficients (upstream "Lateral Coefficients") ---
C_Y_0 = 0.0
C_ell_0 = 0.0
C_n_0 = 0.0
C_Y_beta = -0.98
C_ell_beta = -0.13
C_n_beta = 0.073
C_Y_p = 0.0
C_ell_p = -0.51
C_n_p = 0.069
C_Y_r = 0.0
C_ell_r = 0.25
C_n_r = -0.095
C_Y_delta_a = 0.075
C_ell_delta_a = 0.17
C_n_delta_a = -0.011
C_Y_delta_r = 0.19
C_ell_delta_r = 0.0024
C_n_delta_r = -0.069

# --- Propeller / motor (upstream "Propeller thrust / torque parameters (see addendum by McLain)") ---
D_prop = 20 * 0.0254                 # m      prop diameter ("20*(0.0254)")
KV_rpm_per_volt = 145.0              # RPM/V  motor speed constant
KV = (1.0 / KV_rpm_per_volt) * 60.0 / (2.0 * 3.141592653589793)  # V s/rad, back-emf constant (upstream formula)
KQ = KV                              # N m/A  torque constant (upstream sets KQ = KV)
R_motor = 0.042                      # ohm
i0 = 1.5                             # A      no-load current
ncells = 12.0
V_max = 3.7 * ncells                 # V      battery voltage at full throttle
C_Q2, C_Q1, C_Q0 = -0.01664, 0.004970, 0.005230   # prop torque-coefficient fit C_Q(J)
C_T2, C_T1, C_T0 = -0.1079, -0.06044, 0.09357      # prop thrust-coefficient fit C_T(J)

# ---------------------------------------------------------------------------
# NOT FROM SOURCE (project choices) -- the upstream repo does not contain these.
# ---------------------------------------------------------------------------
# Actuators (Phase 2 requirement: saturation + first-order dynamics).
# Time constants are generic small-UAV servo/ESC values chosen by the project, NOT Aerosonde data.
tau_servo = 0.05                     # s   aileron/elevator/rudder first-order lag
tau_throttle = 0.10                  # s   throttle (ESC/motor command) first-order lag
delta_surface_max = 0.7853981633974483   # rad  +-45 deg surface deflection limit (project choice)
delta_t_min, delta_t_max = 0.0, 1.0      # -    throttle limits

# Dryden gust model parameters: SOURCE = authors' repo, same commit as above,
#   mavsim_simulink/parameters/wind_parameters.m  (WIND.L_u/L_v/L_w = 200/200/50, WIND.sigma_u/v/w = 1.06/1.06/0.7,
#   WIND.Va0 = 25). (The Python/MATLAB wind_simulation files in that repo are student skeletons without values.)
# The shaping-filter transfer functions themselves are not in the repo; they are implemented from the standard Dryden
# form and verified against its closed-form variance and autocorrelation in tests/test_wind_dryden.py.
dryden_Lu, dryden_Lv, dryden_Lw = 200.0, 200.0, 50.0       # m   turbulence scale lengths
dryden_sigma_u, dryden_sigma_v, dryden_sigma_w = 1.06, 1.06, 0.7   # m/s  turbulence intensities
