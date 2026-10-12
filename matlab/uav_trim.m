function T = uav_trim(Va, gamma, R_turn, beta, h)
%UAV_TRIM  Steady-flight trim (core MATLAB only: damped Gauss-Newton, no Optimization Toolbox).
%   T = uav_trim(Va, gamma, R_turn, beta, h): airspeed [m/s], flight-path angle [rad], signed turn radius [m]
%   (Inf = straight), sideslip [rad], altitude [m]. Defaults: gamma = 0, R_turn = Inf, beta = 0, h = 100.
%   Unknowns z = [alpha theta phi da de dr dt]; residuals u_dot,v_dot,w_dot,p_dot,q_dot,r_dot = 0 and
%   h_dot - Va*sin(gamma) = 0 (same formulation as sim/trim.py).
if nargin < 2, gamma = 0; end
if nargin < 3, R_turn = Inf; end
if nargin < 4, beta = 0; end
if nargin < 5, h = 100; end
P = uav_params();
psi_dot = 0; if isfinite(R_turn), psi_dot = Va*cos(gamma)/R_turn; end
z = [0.05; 0.05 + gamma; atan(Va*psi_dot/P.gravity); 0; 0; 0; 0.5];
res = @(z) trim_residual(z, Va, gamma, psi_dot, beta);
for it = 1:100
    r = res(z);
    if norm(r, inf) < 1e-13, break; end
    J = zeros(7, 7);
    for j = 1:7
        dz = zeros(7, 1); dz(j) = 1e-7;
        J(:, j) = (res(z + dz) - res(z - dz))/2e-7;
    end
    step = -J\r;
    t = 1;
    while norm(res(z + t*step)) > norm(r) && t > 1e-6
        t = t/2;
    end
    z = z + t*step;
end
if norm(res(z), inf) > 1e-9
    error('uav_trim:noConvergence', 'trim did not converge (Va=%g, gamma=%g, R=%g)', Va, gamma, R_turn);
end
x13 = trim_state(Va, z(1), beta, z(2), z(3), psi_dot, h);
T = struct('Va', Va, 'gamma', gamma, 'R', R_turn, 'beta', beta, 'alpha', z(1), 'theta', z(2), 'phi', z(3), ...
           'delta', z(4:7), 'x13', x13, 'psi_dot', psi_dot, 'residual', norm(res(z), inf));
end

function x13 = trim_state(Va, alpha, beta, theta, phi, psi_dot, h)
u = Va*cos(alpha)*cos(beta); v = Va*sin(beta); w = Va*sin(alpha)*cos(beta);
p = -psi_dot*sin(theta); q = psi_dot*sin(phi)*cos(theta); r = psi_dot*cos(phi)*cos(theta);
x13 = [0; 0; -h; u; v; w; uav_euler_to_quat(phi, theta, 0); p; q; r];
end

function r = trim_residual(z, Va, gamma, psi_dot, beta)
x13 = trim_state(Va, z(1), beta, z(2), z(3), psi_dot, 100);
xd = uav_plant_deriv([x13; z(4:7)], z(4:7), [0; 0; 0], [0; 0; 0]);   % actuators at their commanded value
r = [xd(4:6); xd(11:13); -xd(3) - Va*sin(gamma)];
end
