function [A, B, x12, delta] = uav_linearize(T)
%UAV_LINEARIZE  Central-difference linearisation of the 12-state Euler plant about the trim T (from uav_trim).
%   d(dx)/dt = A dx + B du,  x = [pn pe pd u v w phi theta psi p q r],  u = [da de dr dt]  (step 1e-6).
x13 = T.x13;
[phi, th, psi] = uav_quat_to_euler(x13(7:10));
x12 = [x13(1:6); phi; th; psi; x13(11:13)];
delta = T.delta;
A = zeros(12, 12); B = zeros(12, 4);
h = 1e-6;
for j = 1:12
    d = zeros(12, 1); d(j) = h;
    A(:, j) = (uav_plant_deriv_euler(x12 + d, delta) - uav_plant_deriv_euler(x12 - d, delta))/(2*h);
end
for j = 1:4
    d = zeros(4, 1); d(j) = h;
    B(:, j) = (uav_plant_deriv_euler(x12, delta + d) - uav_plant_deriv_euler(x12, delta - d))/(2*h);
end
end
