function R = uav_quat_to_R(e)
%UAV_QUAT_TO_R  Rotation matrix R_v_to_b (NED -> body) from quaternion e = [e0 ex ey ez] (scalar first).
%#codegen
e0 = e(1); ex = e(2); ey = e(3); ez = e(4);
R = [e0^2+ex^2-ey^2-ez^2, 2*(ex*ey+e0*ez),       2*(ex*ez-e0*ey);
     2*(ex*ey-e0*ez),     e0^2-ex^2+ey^2-ez^2,   2*(ey*ez+e0*ex);
     2*(ex*ez+e0*ey),     2*(ey*ez-e0*ex),       e0^2-ex^2-ey^2+ez^2];
end
