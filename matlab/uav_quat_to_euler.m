function [phi, theta, psi] = uav_quat_to_euler(e)
%UAV_QUAT_TO_EULER  Euler angles [rad] from quaternion e = [e0 ex ey ez].
%#codegen
e0 = e(1); ex = e(2); ey = e(3); ez = e(4);
phi   = atan2(2*(e0*ex + ey*ez), e0^2 + ez^2 - ex^2 - ey^2);
theta = asin(max(-1, min(1, 2*(e0*ey - ex*ez))));
psi   = atan2(2*(e0*ez + ex*ey), e0^2 + ex^2 - ey^2 - ez^2);
end
