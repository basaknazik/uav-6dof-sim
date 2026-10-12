function xdot = uav_rigid_body(x13, f_b, m_b)
%UAV_RIGID_BODY  13-state rigid-body equations (flat Earth, constant mass).
%   x13 = [pn pe pd u v w e0 ex ey ez p q r] (NED position, body velocity, quaternion, body rates).
%   f_b [N] total body-frame force (incl. gravity), m_b [N m] total moment about the CG.
%   p_dot = R^T v;  v_dot = -w x v + f/m;  e_dot = 0.5*Omega(w) e;  w_dot via Gamma_1..Gamma_8.
%#codegen
P = uav_params();
v = x13(4:6); e = x13(7:10); pqr = x13(11:13);
p = pqr(1); q = pqr(2); r = pqr(3);
R = uav_quat_to_R(e);

G  = P.Jx*P.Jz - P.Jxz^2;
g1 = P.Jxz*(P.Jx - P.Jy + P.Jz)/G;
g2 = (P.Jz*(P.Jz - P.Jy) + P.Jxz^2)/G;
g3 = P.Jz/G;  g4 = P.Jxz/G;
g5 = (P.Jz - P.Jx)/P.Jy;  g6 = P.Jxz/P.Jy;
g7 = ((P.Jx - P.Jy)*P.Jx + P.Jxz^2)/G;
g8 = P.Jx/G;

l = m_b(1); m = m_b(2); n = m_b(3);
pdot = g1*p*q - g2*q*r + g3*l + g4*n;
qdot = g5*p*r - g6*(p^2 - r^2) + m/P.Jy;
rdot = g7*p*q - g1*q*r + g4*l + g8*n;

Om = [0 -p -q -r; p 0 r -q; q -r 0 p; r q -p 0];
xdot = [R.'*v;
        -cross(pqr, v) + f_b(:)/P.mass;
        0.5*Om*e;
        pdot; qdot; rdot];
end
