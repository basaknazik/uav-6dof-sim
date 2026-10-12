function R = uav_euler_to_R(phi, theta, psi)
%UAV_EULER_TO_R  Rotation matrix R_v_to_b (NED -> body), 3-2-1 Euler angles [rad].
%#codegen
cphi = cos(phi); sphi = sin(phi); cth = cos(theta); sth = sin(theta); cpsi = cos(psi); spsi = sin(psi);
R = [cth*cpsi,                    cth*spsi,                    -sth;
     sphi*sth*cpsi-cphi*spsi,     sphi*sth*spsi+cphi*cpsi,     sphi*cth;
     cphi*sth*cpsi+sphi*spsi,     cphi*sth*spsi-sphi*cpsi,     cphi*cth];
end
