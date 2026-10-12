function f = uav_gravity(x13)
%UAV_GRAVITY  Gravity force [N] in the body frame: m*g*R_v_to_b*[0;0;1].
%#codegen
P = uav_params();
R = uav_quat_to_R(x13(7:10));
f = P.mass*P.gravity*R(:, 3);
end
