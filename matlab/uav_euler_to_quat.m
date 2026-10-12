function e = uav_euler_to_quat(phi, theta, psi)
%UAV_EULER_TO_QUAT  Unit quaternion [e0 ex ey ez] (scalar first) from 3-2-1 Euler angles [rad].
%#codegen
cphi = cos(phi/2); sphi = sin(phi/2); cth = cos(theta/2); sth = sin(theta/2); cpsi = cos(psi/2); spsi = sin(psi/2);
e = [cpsi*cth*cphi + spsi*sth*sphi;
     cpsi*cth*sphi - spsi*sth*cphi;
     cpsi*sth*cphi + spsi*cth*sphi;
     spsi*cth*cphi - cpsi*sth*sphi];
end
