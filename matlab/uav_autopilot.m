function [delta_cmd, zone] = uav_autopilot(cmd, x17, ad, G)
%UAV_AUTOPILOT  Successive-loop-closure autopilot, one call per controller sample (port of sim/autopilot.py).
%   cmd = [Va_c; h_c; chi_c]  (m/s, m, rad)     x17 = [x13; actuators]      ad = [Va; alpha; beta]
%   G   = 23-vector: the 22-entry gain vector of uav_design_gains followed by the sample time Ts [s].
%   Returns delta_cmd = [da de dr dt] (before actuator dynamics) and zone (0 hold, 1 climb, 2 descend).
%   Persistent state: PI integrators, yaw-damper washout, zone, previous pitch command. Use `clear uav_autopilot`
%   to reset between runs in plain MATLAB (Simulink resets it at the start of every simulation).
%#codegen
persistent I_chi I_h I_Vt I_Vp r_wash zone_s theta_c_prev first
if isempty(first)
    I_chi = 0; I_h = 0; I_Vt = 0; I_Vp = 0; r_wash = 0; zone_s = 0; theta_c_prev = G(18); first = true;
end
P = uav_params();
kp_phi = G(1); kd_phi = G(2); kp_chi = G(3); ki_chi = G(4); kr = G(5); kp_th = G(6); kd_th = G(7);
kp_h = G(8); ki_h = G(9); kp_Vt = G(10); ki_Vt = G(11); kp_Vp = G(12); ki_Vp = G(13);
phi_max = G(14); theta_max = G(15); h_hold = G(16); tau_w = G(17);
tr_theta = G(18); tr_da = G(19); tr_de = G(20); tr_dr = G(21); tr_dt = G(22); dt = G(23);

Va_c = cmd(1); h_c = cmd(2); chi_c = cmd(3);

% ---- measurements (perfect state feedback)
x13 = x17(1:13);
[phi, theta, ~] = uav_quat_to_euler(x13(7:10));
R = uav_quat_to_R(x13(7:10));
v_ned = R.'*x13(4:6);
chi = atan2(v_ned(2), v_ned(1));
pb = x13(11); qb = x13(12); rb = x13(13);
h = -x13(3);
Va = ad(1);

% ---- lateral: course PI -> roll command -> roll PD
chi_err = mod(chi_c - chi + pi, 2*pi) - pi;
[phi_c, I_chi] = pi_update(chi_err, dt, -phi_max, phi_max, 0, kp_chi, ki_chi, I_chi);
d_a = tr_da + kp_phi*(phi_c - phi) - kd_phi*pb;
% yaw damper with washout
r_hp = rb - r_wash;
r_wash = r_wash + (rb - r_wash)*(1 - exp(-dt/tau_w));
d_r = tr_dr + kr*r_hp;

% ---- longitudinal
h_err = h_c - h;
new_zone = zone_s;
if zone_s == 0
    if h_err > h_hold
        new_zone = 1;
    elseif h_err < -h_hold
        new_zone = 2;
    end
elseif zone_s == 1 && h_err < 0.8*h_hold
    new_zone = 0;
elseif zone_s == 2 && h_err > -0.8*h_hold
    new_zone = 0;
end
if new_zone ~= zone_s && ~first                          % bumpless transfer
    if new_zone == 0
        I_h = theta_c_prev - tr_theta - kp_h*h_err;
        if zone_s == 1
            out_t = 1.0;
        elseif zone_s == 2
            out_t = 0.0;
        else
            out_t = tr_dt;
        end
        I_Vt = out_t - tr_dt - kp_Vt*(Va_c - Va);
    else
        I_Vp = theta_c_prev - tr_theta - kp_Vp*(Va_c - Va);
    end
end
zone_s = new_zone;
if zone_s == 0
    [theta_c, I_h] = pi_update(h_err, dt, -theta_max, theta_max, tr_theta, kp_h, ki_h, I_h);
    [d_t, I_Vt] = pi_update(Va_c - Va, dt, 0, 1, tr_dt, kp_Vt, ki_Vt, I_Vt);
else
    [theta_c, I_Vp] = pi_update(Va_c - Va, dt, -theta_max, theta_max, tr_theta, kp_Vp, ki_Vp, I_Vp);
    if zone_s == 1
        d_t = 1.0;
    else
        d_t = 0.0;
    end
end
d_e = tr_de + kp_th*(theta_c - theta) - kd_th*qb;
theta_c_prev = theta_c;
first = false;

lim = P.delta_surface_max;
delta_cmd = [min(max(d_a, -lim), lim); min(max(d_e, -lim), lim); min(max(d_r, -lim), lim); min(max(d_t, P.delta_t_min), P.delta_t_max)];
zone = zone_s;
end

function [u, I] = pi_update(e, dt, lo, hi, ff, kp, ki, I)
% PI with conditional-integration anti-windup (same logic as sim.autopilot.PI.update)
I_new = I + ki*e*dt;
u_unsat = ff + kp*e + I_new;
if (u_unsat >= lo && u_unsat <= hi) || (u_unsat > hi && ki*e < 0) || (u_unsat < lo && ki*e > 0)
    I = I_new;
end
u = min(max(ff + kp*e + I, lo), hi);
end
