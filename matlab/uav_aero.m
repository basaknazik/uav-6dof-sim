function [f, m] = uav_aero(Va, alpha, beta, pqr, delta_aer)
%UAV_AERO  Aerodynamic force [N] and moment [N m] in the body frame (no gravity, no propulsion).
%   Va [m/s], alpha/beta [rad], pqr [rad/s], delta_aer = [delta_a delta_e delta_r] [rad].
%   Stall-blended lift, parasitic + induced drag, linear lateral model (Beard & McLain Ch. 4).
%#codegen
P = uav_params();
p = pqr(1); q = pqr(2); r = pqr(3);
da = delta_aer(1); de = delta_aer(2); dr = delta_aer(3);
Va_s = max(Va, 1e-3);
qS = 0.5*P.rho*Va^2*P.S_wing;
cq = P.c/(2*Va_s);  bq = P.b/(2*Va_s);

e1 = exp(-P.M_blend*(alpha - P.alpha0));
e2 = exp( P.M_blend*(alpha + P.alpha0));
sigma = (1 + e1 + e2)/((1 + e1)*(1 + e2));
CL = (1 - sigma)*(P.C_L_0 + P.C_L_alpha*alpha) + sigma*2*sign(alpha)*sin(alpha)^2*cos(alpha);
CD = P.C_D_p + (P.C_L_0 + P.C_L_alpha*alpha)^2/(pi*P.e_oswald*P.AR);
F_lift = qS*(CL + P.C_L_q*cq*q + P.C_L_delta_e*de);
F_drag = qS*(CD + P.C_D_q*cq*q + P.C_D_delta_e*de);

fx = -cos(alpha)*F_drag + sin(alpha)*F_lift;
fz = -sin(alpha)*F_drag - cos(alpha)*F_lift;
fy = qS*(P.C_Y_0 + P.C_Y_beta*beta + P.C_Y_p*bq*p + P.C_Y_r*bq*r + P.C_Y_delta_a*da + P.C_Y_delta_r*dr);
l  = qS*P.b*(P.C_ell_0 + P.C_ell_beta*beta + P.C_ell_p*bq*p + P.C_ell_r*bq*r + P.C_ell_delta_a*da + P.C_ell_delta_r*dr);
m_ = qS*P.c*(P.C_m_0 + P.C_m_alpha*alpha + P.C_m_q*cq*q + P.C_m_delta_e*de);
n  = qS*P.b*(P.C_n_0 + P.C_n_beta*beta + P.C_n_p*bq*p + P.C_n_r*bq*r + P.C_n_delta_a*da + P.C_n_delta_r*dr);
f = [fx; fy; fz];
m = [l; m_; n];
end
