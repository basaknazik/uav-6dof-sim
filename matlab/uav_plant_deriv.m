function xdot = uav_plant_deriv(x17, delta_cmd, wind_ned, gust_b)
%UAV_PLANT_DERIV  Derivative of the 17-state plant [x13; delta_act] (rigid body + forces + first-order actuators).
%   Command saturation, then d(delta)/dt = (sat(cmd) - delta)./tau  (same as sim/actuators.py).
%#codegen
P = uav_params();
lo = [-P.delta_surface_max; -P.delta_surface_max; -P.delta_surface_max; P.delta_t_min];
hi = [ P.delta_surface_max;  P.delta_surface_max;  P.delta_surface_max; P.delta_t_max];
tau = [P.tau_servo; P.tau_servo; P.tau_servo; P.tau_throttle];
x13 = x17(1:13); delta = x17(14:17);
[f_b, m_b] = uav_forces_moments(x13, delta, wind_ned, gust_b);
xdot = [uav_rigid_body(x13, f_b, m_b);
        (min(max(delta_cmd(:), lo), hi) - delta)./tau];
end
