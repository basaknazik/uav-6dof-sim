function [G, dm] = uav_design_gains(T, A, B)
%UAV_DESIGN_GAINS  Autopilot gains from the linearised model, core MATLAB only (port of sim/autopilot.py design_gains).
%   G = gain vector, packed as
%     [kp_phi kd_phi kp_chi ki_chi kr_yaw kp_theta kd_theta kp_h ki_h kp_V ki_V kp_V2 ki_V2 ...
%      phi_max theta_max h_hold tau_wash trim_theta trim_delta(1:4)]       (22 entries)
%   dm = [a_phi1 a_phi2 a_theta1 a_theta2 a_theta3 a_V1 a_V2 a_V3] design-model coefficients.
P = uav_params();
Va = T.Va;
% design spec (same defaults as DesignSpec in sim/autopilot.py)
wn_phi = 8.0;  zeta_phi = 0.9;  W_chi = 10.0;  zeta_chi = 1.3;
wn_theta = 14.0;  zeta_theta = 0.8;  W_h = 50.0;  zeta_h = 1.0;
wn_V = 0.8;  zeta_V = 0.9;  W_V2 = 15.0;  zeta_V2 = 0.9;  zeta_dr = 0.7;  tau_wash = 1.0;
phi_max = deg2rad(35.0);  theta_max = deg2rad(30.0);  h_hold = 10.0;

p = 10; q = 11; v = 5; w = 6; r = 12;                       % state indices (pn pe pd u v w phi theta psi p q r)
a_phi1 = -A(p, p);  a_phi2 = B(p, 1);
a_th1 = -A(q, q);   a_th2 = -A(q, w)*Va;  a_th3 = B(q, 2);
[a_V1, a_V2, a_V3] = airspeed_coefficients(T);
dm = [a_phi1 a_phi2 a_th1 a_th2 a_th3 a_V1 a_V2 a_V3];

kp_phi = wn_phi^2/a_phi2;
kd_phi = max(0, (2*zeta_phi*wn_phi - a_phi1)/a_phi2);
wn_chi = wn_phi/W_chi;
kp_chi = 2*zeta_chi*wn_chi*Va/P.gravity;
ki_chi = wn_chi^2*Va/P.gravity;
kp_th = (wn_theta^2 - a_th2)/a_th3;
kd_th = (2*zeta_theta*wn_theta - a_th1)/a_th3;
K_dc = a_th3*kp_th/(a_th2 + a_th3*kp_th);
wn_h = wn_theta/W_h;
kp_h = 2*zeta_h*wn_h/(K_dc*Va);
ki_h = wn_h^2/(K_dc*Va);
kp_V = (2*zeta_V*wn_V - a_V1)/a_V2;
ki_V = wn_V^2/a_V2;
wn_V2 = wn_theta/W_V2;
kp_V2 = (a_V1 - 2*zeta_V2*wn_V2)/(K_dc*a_V3);
ki_V2 = -wn_V2^2/(K_dc*a_V3);
kr = yaw_damper_gain(A, B, zeta_dr);

G = [kp_phi kd_phi kp_chi ki_chi kr kp_th kd_th kp_h ki_h kp_V ki_V kp_V2 ki_V2 ...
     phi_max theta_max h_hold tau_wash T.theta T.delta(:).'].';
end

function [a_V1, a_V2, a_V3] = airspeed_coefficients(T)
% Va_dot(Va, theta, delta_t) with alpha, beta, rates fixed at trim; central differences on the nonlinear plant.
f = @(z) va_dot(z, T);
z0 = [T.Va; T.theta; T.delta(4)];
J = zeros(1, 3);
for j = 1:3
    d = zeros(3, 1); d(j) = 1e-6;
    J(j) = (f(z0 + d) - f(z0 - d))/2e-6;
end
a_V1 = -J(1); a_V2 = J(3); a_V3 = -J(2);
end

function y = va_dot(z, T)
Va = z(1); th = z(2); dt = z(3); alpha = T.alpha;
x12 = [0; 0; -100; Va*cos(alpha); 0; Va*sin(alpha); 0; th; 0; 0; 0; 0];
d = T.delta; d(4) = dt;
xd = uav_plant_deriv_euler(x12, d);
y = cos(alpha)*xd(4) + sin(alpha)*xd(6);
end

function k = yaw_damper_gain(A, B, zeta_target)
% smallest k with delta_r = k*r reaching Dutch-roll damping >= zeta_target on the (v, r) sub-model
A2 = A([5 12], [5 12]);  b2 = B([5 12], 3);
best_k = 0; best_z = -1;
for kk = linspace(0, 5, 1001)
    ev = eig(A2 + b2*[0 kk]);
    [~, imax] = max(imag(ev));                 % member with positive imaginary part (or first of a real pair)
    if abs(imag(ev(imax))) > 1e-9
        z = -real(ev(imax))/abs(ev(imax));
    else
        z = 1.0;
    end
    if z >= zeta_target
        k = kk; return;
    end
    if z > best_z, best_z = z; best_k = kk; end
end
k = best_k;
end
