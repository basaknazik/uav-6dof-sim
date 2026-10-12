# Phase 7 validation: MATLAB/Simulink vs Python

MATLAB 26.1.0.3203278 (R2026a). Python reference: `scripts/make_matlab_reference.py` (RK4, dt = 0.01 s, no quaternion renormalisation).

## Plain-MATLAB functions vs Python (identical inputs)
| check | max error |
|---|---|
| 17-state derivative, 200 random states/inputs (abs) | 8.53e-14 |
| trim at 25 m/s: alpha / theta / controls (abs) | 2.1e-17 / 1.4e-17 / 1.1e-16 |
| linearisation A / B (abs) | 1.8e-09 / 8.1e-11 |
| eigenvalues (abs) | 6.8e-10 |
| autopilot gains designed in MATLAB vs Python (rel) | 4.4e-10 |

## Simulink models vs Python time histories (ode4 = RK4, Ts = 0.01 s)
| run | steps | max abs state error | position | velocity | quaternion | rates | actuators | quat norm drift |
|---|---|---|---|---|---|---|---|---|
| open loop (doublets + saturation + wind + gusts) | 3000 | 6.82e-13 | 6.8e-13 m | 1.4e-14 m/s | 1.3e-15 | 1.3e-15 rad/s | 5.6e-17 | 1.5e-08 |
| closed loop, no wind | 17000 | 2.73e-12 | 2.7e-12 m | 2.2e-14 m/s | 3.9e-15 | 3.3e-14 rad/s | 5.9e-15 | 8.8e-09 |
| closed loop, wind + Dryden gusts | 17000 | 5.00e-12 | 5.0e-12 m | 8.4e-14 m/s | 6.3e-15 | 4.3e-14 rad/s | 7.4e-15 | 1.3e-07 |

closed_loop_nominal: autopilot command error 8.0e-15, zone mismatches 0, difference to Python WITH renormalisation 2.1e-05 (the quaternion renormalisation Simulink does not do).

closed_loop_wind: autopilot command error 1.2e-14, zone mismatches 0, difference to Python WITH renormalisation 6.5e-04 (the quaternion renormalisation Simulink does not do).
