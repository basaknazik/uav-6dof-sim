function build_uav_models(outDir)
%BUILD_UAV_MODELS  Build the Simulink models programmatically (reproducible, diffable; no .slx is committed).
%   build_uav_models(outDir) writes
%     uav_open_loop.slx    plant only, driven by logged command / wind / gust time series
%     uav_closed_loop.slx  plant + autopilot (MATLAB Function block, 100 Hz) driven by [Va_c h_c chi_c] and gusts
%   Structure of the plant subsystem ("UAV Plant"):
%     Actuators (saturation + first-order lags) -> Air data -> Aerodynamics / Propulsion / Gravity -> summed force and
%     moment -> Rigid-body 6-DOF equations (13 states) -> Integrator.   Fixed-step ode4 (RK4), Ts = 0.01 s, like the Python RK4.
%   Equation code lives in plain .m files (uav_*.m) called from MATLAB Function blocks, so it is unit-testable and diffable.
if nargin < 1, outDir = fileparts(mfilename('fullpath')); end
P = uav_params();
build_model('uav_open_loop', false, outDir, P);
build_model('uav_closed_loop', true, outDir, P);
end

function build_model(name, closed, outDir, P)
if bdIsLoaded(name), close_system(name, 0); end
if isfile(fullfile(outDir, [name '.slx'])), delete(fullfile(outDir, [name '.slx'])); end   % avoid name shadowing
new_system(name);
set_param(name, 'SolverType', 'Fixed-step', 'Solver', 'ode4', 'FixedStep', '0.01', 'StopTime', '30', ...
    'SaveFormat', 'Dataset', 'SignalLogging', 'off', 'ReturnWorkspaceOutputs', 'on');
% constants live in the model workspace (overridden per run through Simulink.SimulationInput.setVariable)
mw = get_param(name, 'ModelWorkspace');
assignin(mw, 'actuator_lo', [-P.delta_surface_max; -P.delta_surface_max; -P.delta_surface_max; P.delta_t_min]);
assignin(mw, 'actuator_hi', [ P.delta_surface_max;  P.delta_surface_max;  P.delta_surface_max; P.delta_t_max]);
assignin(mw, 'actuator_tau', [P.tau_servo; P.tau_servo; P.tau_servo; P.tau_throttle]);
assignin(mw, 'x13_0', [0; 0; -100; 25; 0; 0; 1; 0; 0; 0; 0; 0; 0]);
assignin(mw, 'delta_0', [0; 0; 0; 0.5]);
assignin(mw, 'wind_vec', [0; 0; 0]);

plant = [name '/UAV Plant'];
add_block('built-in/Subsystem', plant);
build_plant(plant);

% ---- gust input (always from workspace, sampled at the controller rate)
add_block('built-in/FromWorkspace', [name '/Gust (body, from workspace)'], 'VariableName', 'gust_ts', ...
    'SampleTime', '0.01', 'Interpolate', 'off', 'OutputAfterFinalValue', 'Holding final value');
add_block('built-in/Constant', [name '/Steady wind NED'], 'Value', 'wind_vec');
add_line(name, 'Steady wind NED/1', 'UAV Plant/2');
add_line(name, 'Gust (body, from workspace)/1', 'UAV Plant/3');

add_block('built-in/ToWorkspace', [name '/Log x17'], 'VariableName', 'x17_log', 'SaveFormat', 'Array', 'SampleTime', '-1', 'MaxDataPoints', 'inf');
add_block('built-in/ToWorkspace', [name '/Log air data'], 'VariableName', 'ad_log', 'SaveFormat', 'Array', 'SampleTime', '-1', 'MaxDataPoints', 'inf');
add_line(name, 'UAV Plant/1', 'Log x17/1');
add_line(name, 'UAV Plant/2', 'Log air data/1');

if ~closed
    add_block('built-in/FromWorkspace', [name '/Surface and throttle commands'], 'VariableName', 'ol_cmd_ts', ...
        'SampleTime', '0.01', 'Interpolate', 'off', 'OutputAfterFinalValue', 'Holding final value');
    add_line(name, 'Surface and throttle commands/1', 'UAV Plant/1');
else
    add_block('built-in/FromWorkspace', [name '/Mission commands [Va_c h_c chi_c]'], 'VariableName', 'cl_cmd_ts', ...
        'SampleTime', '0.01', 'Interpolate', 'off', 'OutputAfterFinalValue', 'Holding final value');
    add_block('built-in/ZeroOrderHold', [name '/ZOH state'], 'SampleTime', '0.01');
    add_block('built-in/ZeroOrderHold', [name '/ZOH air data'], 'SampleTime', '0.01');
    add_block('built-in/Constant', [name '/Autopilot gains + Ts'], 'Value', 'ap_G');
    ap = [name '/Autopilot'];
    add_mfb(ap, sprintf(['function [delta_cmd, zone] = fcn(cmd, x17, ad, G)\n' ...
        '%% successive-loop-closure autopilot, executed at every 0.01 s sample\n' ...
        '[delta_cmd, zone] = uav_autopilot(cmd, x17, ad, G);\n']));
    add_block('built-in/Terminator', [name '/unused']);
    add_block('built-in/ToWorkspace', [name '/Log delta cmd'], 'VariableName', 'dcmd_log', 'SaveFormat', 'Array', 'SampleTime', '-1', 'MaxDataPoints', 'inf');
    add_block('built-in/ToWorkspace', [name '/Log zone'], 'VariableName', 'zone_log', 'SaveFormat', 'Array', 'SampleTime', '-1', 'MaxDataPoints', 'inf');
    add_assign = @(a, b) add_line(name, a, b);
    add_assign('Mission commands [Va_c h_c chi_c]/1', 'Autopilot/1');
    add_assign('UAV Plant/1', 'ZOH state/1');
    add_assign('ZOH state/1', 'Autopilot/2');
    add_assign('UAV Plant/2', 'ZOH air data/1');
    add_assign('ZOH air data/1', 'Autopilot/3');
    add_assign('Autopilot gains + Ts/1', 'Autopilot/4');
    add_assign('Autopilot/1', 'UAV Plant/1');
    add_assign('Autopilot/1', 'Log delta cmd/1');
    add_assign('Autopilot/2', 'Log zone/1');
    delete_block([name '/unused']);
end
Simulink.BlockDiagram.arrangeSystem(name);
Simulink.BlockDiagram.arrangeSystem(plant);
save_system(name, fullfile(outDir, [name '.slx']));
close_system(name, 0);
end

function build_plant(sys)
% ---- ports
add_block('built-in/Inport', [sys '/delta_cmd'], 'Port', '1');
add_block('built-in/Inport', [sys '/wind_ned'], 'Port', '2');
add_block('built-in/Inport', [sys '/gust_body'], 'Port', '3');
add_block('built-in/Outport', [sys '/x17 (state + actuators)'], 'Port', '1');
add_block('built-in/Outport', [sys '/air data [Va alpha beta]'], 'Port', '2');

% ---- actuators: saturate command, first-order lag  d(delta)/dt = (sat(cmd) - delta)/tau
add_block('simulink/Discontinuities/Saturation', [sys '/Command saturation'], 'UpperLimit', 'actuator_hi', 'LowerLimit', 'actuator_lo');
add_block('built-in/Sum', [sys '/Servo error'], 'Inputs', '+-');
add_block('built-in/Gain', [sys '/1 over tau'], 'Gain', '1./actuator_tau', 'Multiplication', 'Element-wise(K.*u)');
add_block('built-in/Integrator', [sys '/Actuator states'], 'InitialCondition', 'delta_0');
add_line(sys, 'delta_cmd/1', 'Command saturation/1');
add_line(sys, 'Command saturation/1', 'Servo error/1');
add_line(sys, 'Actuator states/1', 'Servo error/2');
add_line(sys, 'Servo error/1', '1 over tau/1');
add_line(sys, '1 over tau/1', 'Actuator states/1');

% ---- model blocks (code in matlab/uav_*.m)
add_mfb([sys '/Air data'], sprintf(['function ad = fcn(x13, wind, gust)\n[Va, alpha, beta] = uav_airdata(x13, wind, gust);\nad = [Va; alpha; beta];\n']));
add_mfb([sys '/Aerodynamics'], sprintf(['function [f_aero, m_aero] = fcn(x13, ad, delta)\n[f_aero, m_aero] = uav_aero(ad(1), ad(2), ad(3), x13(11:13), delta(1:3));\n']));
add_mfb([sys '/Propulsion'], sprintf(['function [f_prop, m_prop] = fcn(ad, delta)\n[T, Q] = uav_prop(ad(1), delta(4));\nf_prop = [T; 0; 0];\nm_prop = [-Q; 0; 0];\n']));
add_mfb([sys '/Gravity'], sprintf(['function f_g = fcn(x13)\nf_g = uav_gravity(x13);\n']));
add_mfb([sys '/Rigid body 6-DOF'], sprintf(['function xdot = fcn(x13, f_b, m_b)\nxdot = uav_rigid_body(x13, f_b, m_b);\n']));
add_block('built-in/Sum', [sys '/Total force'], 'Inputs', '+++');
add_block('built-in/Sum', [sys '/Total moment'], 'Inputs', '++');
add_block('built-in/Integrator', [sys '/State integrator (13)'], 'InitialCondition', 'x13_0');
add_block('built-in/Mux', [sys '/Mux x17'], 'Inputs', '2');

% ---- wiring
add_line(sys, 'State integrator (13)/1', 'Air data/1');
add_line(sys, 'wind_ned/1', 'Air data/2');
add_line(sys, 'gust_body/1', 'Air data/3');
add_line(sys, 'State integrator (13)/1', 'Aerodynamics/1');
add_line(sys, 'Air data/1', 'Aerodynamics/2');
add_line(sys, 'Actuator states/1', 'Aerodynamics/3');
add_line(sys, 'Air data/1', 'Propulsion/1');
add_line(sys, 'Actuator states/1', 'Propulsion/2');
add_line(sys, 'State integrator (13)/1', 'Gravity/1');
add_line(sys, 'Aerodynamics/1', 'Total force/1');
add_line(sys, 'Propulsion/1', 'Total force/2');
add_line(sys, 'Gravity/1', 'Total force/3');
add_line(sys, 'Aerodynamics/2', 'Total moment/1');
add_line(sys, 'Propulsion/2', 'Total moment/2');
add_line(sys, 'State integrator (13)/1', 'Rigid body 6-DOF/1');
add_line(sys, 'Total force/1', 'Rigid body 6-DOF/2');
add_line(sys, 'Total moment/1', 'Rigid body 6-DOF/3');
add_line(sys, 'Rigid body 6-DOF/1', 'State integrator (13)/1');
add_line(sys, 'State integrator (13)/1', 'Mux x17/1');
add_line(sys, 'Actuator states/1', 'Mux x17/2');
add_line(sys, 'Mux x17/1', 'x17 (state + actuators)/1');
add_line(sys, 'Air data/1', 'air data [Va alpha beta]/1');
end

function add_mfb(path, code)
% add a MATLAB Function block and set its script
add_block('simulink/User-Defined Functions/MATLAB Function', path);
rt = sfroot;
chart = rt.find('-isa', 'Stateflow.EMChart', 'Path', path);
chart.Script = code;
end
