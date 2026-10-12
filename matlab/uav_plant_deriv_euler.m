function xdot = uav_plant_deriv_euler(x12, delta)
%UAV_PLANT_DERIV_EULER  12-state Euler-angle plant derivative (no wind) used for trim and linearisation.
%   x12 = [pn pe pd u v w phi theta psi p q r]; delta = actual [da de dr dt].
P = uav_params();
x13 = [x12(1:6); uav_euler_to_quat(x12(7), x12(8), x12(9)); x12(10:12)];
[f_b, m_b] = uav_forces_moments(x13, delta, [0; 0; 0], [0; 0; 0]);
phi = x12(7); th = x12(8);
v = x12(4:6); pqr = x12(10:12); p = pqr(1); q = pqr(2); r = pqr(3);
R = uav_euler_to_R(x12(7), x12(8), x12(9));
S = [1, sin(phi)*tan(th), cos(phi)*tan(th);
     0, cos(phi),        -sin(phi);
     0, sin(phi)/cos(th), cos(phi)/cos(th)];
G  = P.Jx*P.Jz - P.Jxz^2;
g1 = P.Jxz*(P.Jx - P.Jy + P.Jz)/G;  g2 = (P.Jz*(P.Jz - P.Jy) + P.Jxz^2)/G;
g3 = P.Jz/G;  g4 = P.Jxz/G;  g5 = (P.Jz - P.Jx)/P.Jy;  g6 = P.Jxz/P.Jy;
g7 = ((P.Jx - P.Jy)*P.Jx + P.Jxz^2)/G;  g8 = P.Jx/G;
l = m_b(1); m = m_b(2); n = m_b(3);
xdot = [R.'*v;
        -cross(pqr, v) + f_b/P.mass;
        S*pqr;
        g1*p*q - g2*q*r + g3*l + g4*n;
        g5*p*r - g6*(p^2 - r^2) + m/P.Jy;
        g7*p*q - g1*q*r + g4*l + g8*n];
end
