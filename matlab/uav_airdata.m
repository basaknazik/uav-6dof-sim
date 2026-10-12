function [Va, alpha, beta] = uav_airdata(x13, wind_ned, gust_b)
%UAV_AIRDATA  Airspeed [m/s], angle of attack and sideslip [rad].
%   x13 = [pn pe pd u v w e0 ex ey ez p q r]; wind_ned [m/s] steady wind in NED; gust_b [m/s] gust in body axes.
%   v_air_b = v_b - (R_v_to_b*wind + gust);  Va = |v_air_b|, alpha = atan2(w_r,u_r), beta = asin(v_r/Va).
%#codegen
R = uav_quat_to_R(x13(7:10));
v_air = x13(4:6) - (R*wind_ned(:) + gust_b(:));
Va = norm(v_air);
alpha = atan2(v_air(3), v_air(1));
beta = asin(max(-1, min(1, v_air(2)/max(Va, 1e-3))));
end
