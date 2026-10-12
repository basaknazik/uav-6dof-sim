function M = validate_functions(ref)
%VALIDATE_FUNCTIONS  Compare the plain-MATLAB model/trim/linearisation/gain-design against the Python reference.
%   ref = struct loaded from matlab/validation/data/python_reference.mat. Returns a struct of error metrics.
M = struct();

% ---- A: 17-state derivative at 200 random states/inputs (same inputs, same equations)
n = size(ref.deriv_x, 1);
err = zeros(n, 17); scale = zeros(n, 17);
for k = 1:n
    xd = uav_plant_deriv(ref.deriv_x(k, :).', ref.deriv_cmd(k, :).', ref.deriv_wind(k, :).', ref.deriv_gust(k, :).');
    err(k, :) = abs(xd.' - ref.deriv_xdot(k, :));
    scale(k, :) = abs(ref.deriv_xdot(k, :));
end
M.deriv_max_abs_err = max(err(:));
M.deriv_max_rel_err = max(err(:)./max(scale(:), 1));    % relative to max(|value|, 1)
M.deriv_cases = n;

% ---- B: trim, linearisation, gains at Va = 25 m/s (core MATLAB solver vs scipy least_squares)
T = uav_trim(ref.Va_trim);
M.trim_alpha_err = abs(T.alpha - ref.trim_alpha);
M.trim_theta_err = abs(T.theta - ref.trim_theta);
M.trim_delta_max_err = max(abs(T.delta(:) - ref.trim_delta(:)));
M.trim_residual = T.residual;
[A, B] = uav_linearize(T);
M.A_max_abs_err = max(abs(A - ref.A), [], 'all');
M.B_max_abs_err = max(abs(B - ref.B), [], 'all');
M.A_max_rel_err = max(abs(A - ref.A)./max(abs(ref.A), 1), [], 'all');
[G, dm] = uav_design_gains(T, A, B);
g_ref = ref.gains(:);
M.gains_max_rel_err = max(abs(G - g_ref)./max(abs(g_ref), 1e-12));
M.design_models_max_rel_err = max(abs(dm(:) - ref.design_models(:))./abs(ref.design_models(:)));
M.gains = G(:).';
M.gains_python = g_ref(:).';

% ---- C: eigenvalues of the linearised model (modal analysis in core MATLAB: eig)
idx_lon = [4 6 8 11]; idx_lat = [5 7 10 12];     % u w theta q ; v phi p r
ev_lon = sort(eig(A(idx_lon, idx_lon))); ev_lat = sort(eig(A(idx_lat, idx_lat)));
ev_lon_py = sort(eig(ref.A(idx_lon, idx_lon))); ev_lat_py = sort(eig(ref.A(idx_lat, idx_lat)));
M.eig_max_abs_err = max([abs(ev_lon - ev_lon_py); abs(ev_lat - ev_lat_py)]);
M.eig_lon_re = real(ev_lon(:)).'; M.eig_lon_im = imag(ev_lon(:)).';
M.eig_lat_re = real(ev_lat(:)).'; M.eig_lat_im = imag(ev_lat(:)).';
end
