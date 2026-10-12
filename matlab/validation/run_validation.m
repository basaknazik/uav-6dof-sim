function M = run_validation(skipSimulink)
%RUN_VALIDATION  Cross-validate the MATLAB/Simulink port against the Python reference.
%   Run from anywhere:  matlab -batch "run('matlab/validation/run_validation.m')"   (or call run_validation).
%   Needs matlab/validation/data/python_reference.mat  (python scripts/make_matlab_reference.py).
%   Writes results/phase7_validation.{json,md,png}. Errors out if a tolerance is exceeded.
if nargin < 1, skipSimulink = false; end
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fileparts(here)); addpath(here);
resDir = fullfile(root, 'results'); if ~isfolder(resDir), mkdir(resDir); end
ref = load(fullfile(here, 'data', 'python_reference.mat'));

M = struct();
M.functions = validate_functions(ref);                       % plain-MATLAB equations, trim, linearisation, gains
if ~skipSimulink
    build_uav_models(fileparts(here));                       % writes uav_*.slx next to the build script
    export_diagrams(resDir);                                  % block-diagram pictures for the README
    M.open_loop = run_open_loop(ref);
    M.closed_loop_nominal = run_closed_loop(ref, 'cl_nom');
    M.closed_loop_wind = run_closed_loop(ref, 'cl_wind');
end
M.matlab_version = version;
write_reports(M, resDir);
check_tolerances(M);
end

% ---------------------------------------------------------------------------------------------
function export_diagrams(resDir)
load_system('uav_closed_loop');
print('-suav_closed_loop', '-dpng', '-r110', fullfile(resDir, 'simulink_closed_loop.png'));
print('-suav_closed_loop/UAV Plant', '-dpng', '-r110', fullfile(resDir, 'simulink_plant.png'));
close_system('uav_closed_loop', 0);
end

function out = ts_matrix(t, data)
% time-stamped data for From Workspace blocks: stamps slightly BEFORE the sample instants so that, at the exact
% discrete hit time k*Ts, the block returns row k (immune to floating-point round-off in the time vector)
stamps = t(:) - 1e-9; stamps(1) = 0;
out = [stamps, data];
end

function r = run_open_loop(ref)
mdl = 'uav_open_loop';
n = numel(ref.ol_t);
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StopTime', num2str(ref.ol_t(end)));
in = in.setVariable('x13_0', ref.trim_x13(:), 'Workspace', mdl);
in = in.setVariable('delta_0', ref.trim_delta(:), 'Workspace', mdl);
in = in.setVariable('wind_vec', ref.ol_wind(:), 'Workspace', mdl);
in = in.setVariable('ol_cmd_ts', ts_matrix(ref.ol_t, ref.ol_cmd));
in = in.setVariable('gust_ts', ts_matrix(ref.ol_t, ref.ol_gust));
tic; out = sim(in); r.sim_seconds = toc;
X = to_rows(out.x17_log, n);
r.steps = n - 1;
r = compare_states(r, X, ref.ol_x, ref.ol_t, 'ol');
end

function r = run_closed_loop(ref, prefix)
mdl = 'uav_closed_loop';
t = ref.([prefix '_t']); n = numel(t);
dt = ref.dt;
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StopTime', num2str(t(end)));
in = in.setVariable('x13_0', ref.trim_x13(:), 'Workspace', mdl);
in = in.setVariable('delta_0', ref.trim_delta(:), 'Workspace', mdl);
if isfield(ref, [prefix '_wind']), w = ref.([prefix '_wind']); else, w = [0; 0; 0]; end
in = in.setVariable('wind_vec', w(:), 'Workspace', mdl);
in = in.setVariable('cl_cmd_ts', ts_matrix(t, ref.([prefix '_cmd'])));
in = in.setVariable('gust_ts', ts_matrix(t, ref.([prefix '_gust'])));
in = in.setVariable('ap_G', [ref.gains(:); dt]);
tic; out = sim(in); r.sim_seconds = toc;
X = to_rows(out.x17_log, n);
r.steps = n - 1;
r = compare_states(r, X, ref.([prefix '_x']), t, prefix);
Xn = ref.([prefix '_x_norm']);
r.max_abs_state_err_vs_python_renormalised = max(abs(X - Xn), [], 'all');
dcmd = to_rows(out.dcmd_log, n); zone = to_rows(out.zone_log, n);
r.delta_cmd_max_abs_err = max(abs(dcmd - ref.([prefix '_delta_cmd'])), [], 'all');
r.zone_mismatches = sum(zone ~= ref.([prefix '_zone'])(:));
r.final_h_err_m = abs(-X(end, 3) - ref.([prefix '_cmd'])(end, 2));
r.quat_norm_drift = max(abs(vecnorm(X(:, 7:10), 2, 2) - 1));
end

function A = to_rows(a, nT)
% To Workspace (Array format) stores vector signals with time in the LAST dimension; return an nT-by-width matrix
if ndims(a) == 3
    A = reshape(permute(a, [3 1 2]), size(a, 3), []);
elseif size(a, 1) == nT
    A = a;
else
    A = a.';
end
assert(size(A, 1) == nT, 'logged %d samples, expected %d', size(A, 1), nT);
end

function r = compare_states(r, X, Xref, t, tag)
E = abs(X - Xref);
groups = struct('position', 1:3, 'body_velocity', 4:6, 'quaternion', 7:10, 'body_rates', 11:13, 'actuators', 14:17);
r.max_abs_state_err = max(E, [], 'all');
fn = fieldnames(groups);
for i = 1:numel(fn)
    r.(['err_' fn{i}]) = max(E(:, groups.(fn{i})), [], 'all');
end
r.rms_state_err = sqrt(mean(E.^2, 'all'));
r.quat_norm_drift = max(abs(vecnorm(X(:, 7:10), 2, 2) - 1));
assignin('base', ['val_' tag '_X'], X); assignin('base', ['val_' tag '_Xref'], Xref); assignin('base', ['val_' tag '_t'], t);
end

function write_reports(M, resDir)
fid = fopen(fullfile(resDir, 'phase7_validation.json'), 'w'); fprintf(fid, '%s', jsonencode(M, 'PrettyPrint', true)); fclose(fid);
f = M.functions;
L = {};
L{end+1} = '# Phase 7 validation: MATLAB/Simulink vs Python';
L{end+1} = '';
L{end+1} = sprintf('MATLAB %s. Python reference: `scripts/make_matlab_reference.py` (RK4, dt = 0.01 s, no quaternion renormalisation).', M.matlab_version);
L{end+1} = '';
L{end+1} = '## Plain-MATLAB functions vs Python (identical inputs)';
L{end+1} = '| check | max error |'; L{end+1} = '|---|---|';
L{end+1} = sprintf('| 17-state derivative, %d random states/inputs (abs) | %.2e |', f.deriv_cases, f.deriv_max_abs_err);
L{end+1} = sprintf('| trim at 25 m/s: alpha / theta / controls (abs) | %.1e / %.1e / %.1e |', f.trim_alpha_err, f.trim_theta_err, f.trim_delta_max_err);
L{end+1} = sprintf('| linearisation A / B (abs) | %.1e / %.1e |', f.A_max_abs_err, f.B_max_abs_err);
L{end+1} = sprintf('| eigenvalues (abs) | %.1e |', f.eig_max_abs_err);
L{end+1} = sprintf('| autopilot gains designed in MATLAB vs Python (rel) | %.1e |', f.gains_max_rel_err);
if isfield(M, 'open_loop')
    L{end+1} = ''; L{end+1} = '## Simulink models vs Python time histories (ode4 = RK4, Ts = 0.01 s)';
    L{end+1} = '| run | steps | max abs state error | position | velocity | quaternion | rates | actuators | quat norm drift |';
    L{end+1} = '|---|---|---|---|---|---|---|---|---|';
    runs = {'open_loop', 'open loop (doublets + saturation + wind + gusts)'; ...
            'closed_loop_nominal', 'closed loop, no wind'; 'closed_loop_wind', 'closed loop, wind + Dryden gusts'};
    for i = 1:size(runs, 1)
        r = M.(runs{i, 1});
        L{end+1} = sprintf('| %s | %d | %.2e | %.1e m | %.1e m/s | %.1e | %.1e rad/s | %.1e | %.1e |', runs{i, 2}, r.steps, ...
            r.max_abs_state_err, r.err_position, r.err_body_velocity, r.err_quaternion, r.err_body_rates, r.err_actuators, r.quat_norm_drift);
    end
    for nm = {'closed_loop_nominal', 'closed_loop_wind'}
        r = M.(nm{1});
        L{end+1} = '';
        L{end+1} = sprintf('%s: autopilot command error %.1e, zone mismatches %d, difference to Python WITH renormalisation %.1e (the quaternion renormalisation Simulink does not do).', ...
            nm{1}, r.delta_cmd_max_abs_err, r.zone_mismatches, r.max_abs_state_err_vs_python_renormalised);
    end
end
fid = fopen(fullfile(resDir, 'phase7_validation.md'), 'w'); fprintf(fid, '%s\n', L{:}); fclose(fid);

if isfield(M, 'open_loop')
    fig = figure('Visible', 'off', 'Position', [100 100 1100 700]);
    tags = {'ol', 'cl_nom', 'cl_wind'}; titles = {'open loop', 'closed loop, no wind', 'closed loop, wind + gusts'};
    for i = 1:3
        X = evalin('base', ['val_' tags{i} '_X']); Xr = evalin('base', ['val_' tags{i} '_Xref']); t = evalin('base', ['val_' tags{i} '_t']);
        subplot(2, 3, i); plot(t, -X(:, 3), 'b', t, -Xr(:, 3), 'r--'); grid on; title(titles{i}); ylabel('altitude [m]'); legend('Simulink', 'Python', 'Location', 'best');
        subplot(2, 3, 3 + i); semilogy(t, max(abs(X - Xr), 1e-17) , 'LineWidth', 0.5); grid on; xlabel('t [s]'); ylabel('|error| (all 17 states)');
    end
    saveas(fig, fullfile(resDir, 'phase7_validation.png')); close(fig);
end
end

function check_tolerances(M)
f = M.functions;
assert(f.deriv_max_abs_err < 1e-9, 'derivative mismatch %.2e', f.deriv_max_abs_err);
assert(f.trim_delta_max_err < 1e-8 && f.A_max_abs_err < 1e-6 && f.gains_max_rel_err < 1e-6, 'trim/linearisation/gain mismatch');
if isfield(M, 'open_loop')
    assert(M.open_loop.max_abs_state_err < 1e-6, 'open-loop mismatch %.2e', M.open_loop.max_abs_state_err);
    assert(M.closed_loop_nominal.max_abs_state_err < 1e-5, 'closed-loop (nominal) mismatch %.2e', M.closed_loop_nominal.max_abs_state_err);
    assert(M.closed_loop_wind.max_abs_state_err < 1e-5, 'closed-loop (wind) mismatch %.2e', M.closed_loop_wind.max_abs_state_err);
    assert(M.closed_loop_nominal.zone_mismatches == 0 && M.closed_loop_wind.zone_mismatches == 0, 'zone logic mismatch');
end
end
