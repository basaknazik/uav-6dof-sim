function [f_b, m_b, ad] = uav_forces_moments(x13, delta, wind_ned, gust_b)
%UAV_FORCES_MOMENTS  Total body force/moment (gravity + aero + propulsion) and air data ad = [Va; alpha; beta].
%   delta = ACTUAL [delta_a delta_e delta_r delta_t].
%#codegen
[Va, alpha, beta] = uav_airdata(x13, wind_ned, gust_b);
[fa, ma] = uav_aero(Va, alpha, beta, x13(11:13), delta(1:3));
[T, Q] = uav_prop(Va, delta(4));
f_b = uav_gravity(x13) + fa + [T; 0; 0];
m_b = ma + [-Q; 0; 0];
ad = [Va; alpha; beta];
end
