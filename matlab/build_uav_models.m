function build_uav_models(outDir)
%BUILD_UAV_MODELS  Build the Simulink models programmatically (reproducible, diffable; no .slx is committed).
%   build_uav_models(outDir) writes
%     uav_open_loop.slx    plant only, driven by logged command / wind / gust time series
%     uav_closed_loop.slx  plant + autopilot (MATLAB Function block, 100 Hz) driven by [Va_c h_c chi_c] and gusts
%   Structure of the plant subsystem ("UAV Plant"):
%     Actuators (saturation + first-order lags) -> Air data -> Aerodynamics / Propulsion / Gravity -> summed force and
%     moment -> Rigid-body 6-DOF equations (13 states) -> Integrator.   Fixed-step ode4 (RK4), Ts = 0.01 s, like the Python RK4.
%   Equation code lives in plain .m files (uav_*.m) called from MATLAB Function blocks, so it is unit-testable and diffable.
%   Layout: signal flow left to right, feedback below the forward path, Goto/From tags for the long state feedback.
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
set_param(plant, 'Position', [560 100 700 360]);

fw = {'SampleTime', '0.01', 'Interpolate', 'off', 'OutputAfterFinalValue', 'Holding final value'};
add_blk(name, 'built-in/FromWorkspace', 'Gust (body)', [40 330 200 360], 'VariableName', 'gust_ts', fw{:});
add_blk(name, 'built-in/Constant', 'Steady wind NED', [40 250 200 280], 'Value', 'wind_vec');
add_blk(name, 'built-in/ToWorkspace', 'Log x17', [820 110 920 140], 'VariableName', 'x17_log', 'SaveFormat', 'Array', 'SampleTime', '-1', 'MaxDataPoints', 'inf');
add_blk(name, 'built-in/ToWorkspace', 'Log air data', [820 190 920 220], 'VariableName', 'ad_log', 'SaveFormat', 'Array', 'SampleTime', '-1', 'MaxDataPoints', 'inf');
add_line(name, 'Steady wind NED/1', 'UAV Plant/2', 'autorouting', 'smart');
add_line(name, 'Gust (body)/1', 'UAV Plant/3', 'autorouting', 'smart');
add_line(name, 'UAV Plant/1', 'Log x17/1', 'autorouting', 'smart');
add_line(name, 'UAV Plant/2', 'Log air data/1', 'autorouting', 'smart');

if ~closed
    add_blk(name, 'built-in/FromWorkspace', 'Surface and throttle commands', [40 120 240 150], 'VariableName', 'ol_cmd_ts', fw{:});
    add_line(name, 'Surface and throttle commands/1', 'UAV Plant/1', 'autorouting', 'smart');
else
    add_blk(name, 'built-in/FromWorkspace', 'Mission commands [Va_c h_c chi_c]', [40 40 260 70], 'VariableName', 'cl_cmd_ts', fw{:});
    add_blk(name, 'built-in/Constant', 'Autopilot gains + Ts', [40 170 260 200], 'Value', 'ap_G');
    ap = [name '/Autopilot'];
    add_mfb(ap, sprintf(['function [delta_cmd, zone] = fcn(cmd, x17, ad, G)\n' ...
        '%% successive-loop-closure autopilot, executed at every 0.01 s sample\n' ...
        '[delta_cmd, zone] = uav_autopilot(cmd, x17, ad, G);\n']));
    set_param(ap, 'Position', [340 60 480 200]);
    % feedback path (below the forward path, flowing right to left)
    add_blk(name, 'built-in/ZeroOrderHold', 'ZOH state', [420 420 480 450], 'SampleTime', '0.01', 'Orientation', 'left');
    add_blk(name, 'built-in/ZeroOrderHold', 'ZOH air data', [420 490 480 520], 'SampleTime', '0.01', 'Orientation', 'left');
    add_blk(name, 'built-in/ToWorkspace', 'Log delta cmd', [820 40 920 70], 'VariableName', 'dcmd_log', 'SaveFormat', 'Array', 'SampleTime', '-1', 'MaxDataPoints', 'inf');
    add_blk(name, 'built-in/ToWorkspace', 'Log zone', [820 270 920 300], 'VariableName', 'zone_log', 'SaveFormat', 'Array', 'SampleTime', '-1', 'MaxDataPoints', 'inf');
    add_line(name, 'Mission commands [Va_c h_c chi_c]/1', 'Autopilot/1', 'autorouting', 'smart');
    add_line(name, 'Autopilot gains + Ts/1', 'Autopilot/4', 'autorouting', 'smart');
    add_line(name, 'Autopilot/1', 'UAV Plant/1', 'autorouting', 'smart');
    add_line(name, 'Autopilot/1', 'Log delta cmd/1', 'autorouting', 'smart');
    add_line(name, 'Autopilot/2', 'Log zone/1', 'autorouting', 'smart');
    add_line(name, 'UAV Plant/1', 'ZOH state/1', 'autorouting', 'smart');
    add_line(name, 'ZOH state/1', 'Autopilot/2', 'autorouting', 'smart');
    add_line(name, 'UAV Plant/2', 'ZOH air data/1', 'autorouting', 'smart');
    add_line(name, 'ZOH air data/1', 'Autopilot/3', 'autorouting', 'smart');
end
save_system(name, fullfile(outDir, [name '.slx']));
close_system(name, 0);
end

function build_plant(sys)
% Left-to-right flow:  inputs | actuators | air data | forces and moments | sum | rigid body | integrator | outputs.
% Long feedback signals use Goto/From tags (x13, delta_act, ad) instead of wires that cross the diagram.
add_blk(sys, 'built-in/Inport', 'delta_cmd', [30 63 60 77], 'Port', '1');
add_blk(sys, 'built-in/Inport', 'wind_ned', [30 253 60 267], 'Port', '2');
add_blk(sys, 'built-in/Inport', 'gust_body', [30 303 60 317], 'Port', '3');
add_blk(sys, 'built-in/Outport', 'x17 (state + actuators)', [1240 303 1270 317], 'Port', '1');
add_blk(sys, 'built-in/Outport', 'air data [Va alpha beta]', [1240 383 1270 397], 'Port', '2');

% ---- actuators: saturate command, first-order lag  d(delta)/dt = (sat(cmd) - delta)/tau
add_blk(sys, 'simulink/Discontinuities/Saturation', 'Command saturation', [100 55 150 85], 'UpperLimit', 'actuator_hi', 'LowerLimit', 'actuator_lo');
add_blk(sys, 'built-in/Sum', 'Servo error', [190 55 210 85], 'Inputs', '+-', 'ShowName', 'off');
add_blk(sys, 'built-in/Gain', '1 over tau', [250 55 330 85], 'Gain', '1./actuator_tau', 'Multiplication', 'Element-wise(K.*u)');
add_blk(sys, 'built-in/Integrator', 'Actuator states', [370 55 420 85], 'InitialCondition', 'delta_0');
add_blk(sys, 'built-in/Goto', 'tag delta_act', [460 62 520 78], 'GotoTag', 'delta_act', 'TagVisibility', 'local');
add_blk(sys, 'built-in/From', 'delta_act (feedback)', [150 120 210 136], 'GotoTag', 'delta_act');
add_line(sys, 'delta_cmd/1', 'Command saturation/1', 'autorouting', 'smart');
add_line(sys, 'Command saturation/1', 'Servo error/1', 'autorouting', 'smart');
add_line(sys, 'Servo error/1', '1 over tau/1', 'autorouting', 'smart');
add_line(sys, '1 over tau/1', 'Actuator states/1', 'autorouting', 'smart');
add_line(sys, 'Actuator states/1', 'tag delta_act/1', 'autorouting', 'smart');
add_line(sys, 'delta_act (feedback)/1', 'Servo error/2', 'autorouting', 'smart');

% ---- air data
add_blk(sys, 'built-in/From', 'x13 (to air data)', [130 213 190 227], 'GotoTag', 'x13');
add_mfb([sys '/Air data'], sprintf(['function ad = fcn(x13, wind, gust)\n[Va, alpha, beta] = uav_airdata(x13, wind, gust);\nad = [Va; alpha; beta];\n']));
set_param([sys '/Air data'], 'Position', [260 210 360 330]);
add_blk(sys, 'built-in/Goto', 'tag ad', [400 262 450 278], 'GotoTag', 'ad', 'TagVisibility', 'local');
add_line(sys, 'x13 (to air data)/1', 'Air data/1', 'autorouting', 'smart');
add_line(sys, 'wind_ned/1', 'Air data/2', 'autorouting', 'smart');
add_line(sys, 'gust_body/1', 'Air data/3', 'autorouting', 'smart');
add_line(sys, 'Air data/1', 'tag ad/1', 'autorouting', 'smart');
add_blk(sys, 'built-in/From', 'ad (output)', [1100 383 1160 397], 'GotoTag', 'ad');
add_line(sys, 'ad (output)/1', 'air data [Va alpha beta]/1', 'autorouting', 'smart');

% ---- forces and moments
add_blk(sys, 'built-in/From', 'x13 (aero)', [490 403 550 417], 'GotoTag', 'x13');
add_blk(sys, 'built-in/From', 'ad (aero)', [490 433 550 447], 'GotoTag', 'ad');
add_blk(sys, 'built-in/From', 'delta_act (aero)', [490 463 550 477], 'GotoTag', 'delta_act');
add_mfb([sys '/Aerodynamics'], sprintf(['function [f_aero, m_aero] = fcn(x13, ad, delta)\n[f_aero, m_aero] = uav_aero(ad(1), ad(2), ad(3), x13(11:13), delta(1:3));\n']));
set_param([sys '/Aerodynamics'], 'Position', [620 390 740 490]);
add_blk(sys, 'built-in/From', 'ad (prop)', [490 253 550 267], 'GotoTag', 'ad');
add_blk(sys, 'built-in/From', 'delta_act (prop)', [490 293 550 307], 'GotoTag', 'delta_act');
add_mfb([sys '/Propulsion'], sprintf(['function [f_prop, m_prop] = fcn(ad, delta)\n[T, Q] = uav_prop(ad(1), delta(4));\nf_prop = [T; 0; 0];\nm_prop = [-Q; 0; 0];\n']));
set_param([sys '/Propulsion'], 'Position', [620 240 740 320]);
add_blk(sys, 'built-in/From', 'x13 (gravity)', [490 153 550 167], 'GotoTag', 'x13');
add_mfb([sys '/Gravity'], sprintf(['function f_g = fcn(x13)\nf_g = uav_gravity(x13);\n']));
set_param([sys '/Gravity'], 'Position', [620 140 740 180]);
add_line(sys, 'x13 (gravity)/1', 'Gravity/1', 'autorouting', 'smart');
add_line(sys, 'x13 (aero)/1', 'Aerodynamics/1', 'autorouting', 'smart');
add_line(sys, 'ad (aero)/1', 'Aerodynamics/2', 'autorouting', 'smart');
add_line(sys, 'delta_act (aero)/1', 'Aerodynamics/3', 'autorouting', 'smart');
add_line(sys, 'ad (prop)/1', 'Propulsion/1', 'autorouting', 'smart');
add_line(sys, 'delta_act (prop)/1', 'Propulsion/2', 'autorouting', 'smart');

add_blk(sys, 'built-in/Sum', 'Total force', [800 140 830 190], 'Inputs', '+++', 'IconShape', 'round', 'ShowName', 'off');
add_blk(sys, 'built-in/Sum', 'Total moment', [800 330 830 380], 'Inputs', '++', 'IconShape', 'round', 'ShowName', 'off');
add_text(sys, 'sum F [N]', [838 128 890 144]);
add_text(sys, 'sum M [N m]', [838 318 900 334]);
add_line(sys, 'Gravity/1', 'Total force/1', 'autorouting', 'smart');
add_line(sys, 'Aerodynamics/1', 'Total force/2', 'autorouting', 'smart');
add_line(sys, 'Propulsion/1', 'Total force/3', 'autorouting', 'smart');
add_line(sys, 'Aerodynamics/2', 'Total moment/1', 'autorouting', 'smart');
add_line(sys, 'Propulsion/2', 'Total moment/2', 'autorouting', 'smart');

% ---- rigid body and state integrator
add_blk(sys, 'built-in/From', 'x13 (rigid body)', [820 83 880 97], 'GotoTag', 'x13');
add_mfb([sys '/Rigid body 6-DOF'], sprintf(['function xdot = fcn(x13, f_b, m_b)\nxdot = uav_rigid_body(x13, f_b, m_b);\n']));
set_param([sys '/Rigid body 6-DOF'], 'Position', [900 120 1020 220]);
add_blk(sys, 'built-in/Integrator', 'State integrator (13)', [1060 155 1110 185], 'InitialCondition', 'x13_0');
add_blk(sys, 'built-in/Goto', 'tag x13', [1060 205 1120 221], 'GotoTag', 'x13', 'TagVisibility', 'local');
add_blk(sys, 'built-in/Mux', 'Mux x17', [1060 285 1065 335], 'Inputs', '2');
add_blk(sys, 'built-in/From', 'x13 (output)', [960 288 1020 302], 'GotoTag', 'x13');
add_blk(sys, 'built-in/From', 'delta_act (output)', [960 318 1020 332], 'GotoTag', 'delta_act');
add_line(sys, 'x13 (rigid body)/1', 'Rigid body 6-DOF/1', 'autorouting', 'smart');
add_line(sys, 'Total force/1', 'Rigid body 6-DOF/2', 'autorouting', 'smart');
add_line(sys, 'Total moment/1', 'Rigid body 6-DOF/3', 'autorouting', 'smart');
add_line(sys, 'Rigid body 6-DOF/1', 'State integrator (13)/1', 'autorouting', 'smart');
add_line(sys, 'State integrator (13)/1', 'tag x13/1', 'autorouting', 'smart');
add_line(sys, 'x13 (output)/1', 'Mux x17/1', 'autorouting', 'smart');
add_line(sys, 'delta_act (output)/1', 'Mux x17/2', 'autorouting', 'smart');
add_line(sys, 'Mux x17/1', 'x17 (state + actuators)/1', 'autorouting', 'smart');
end

function add_text(sys, str, pos)
% free text annotation (used instead of block names that would collide with wires)
a = Simulink.Annotation([sys '/' str]);
a.Position = pos;
end

function add_blk(sys, lib, name, pos, varargin)
% add a block with a given name and position [left top right bottom], then optional parameter/value pairs
path = [sys '/' name];
add_block(lib, path);
set_param(path, 'Position', pos);
for k = 1:2:numel(varargin)
    set_param(path, varargin{k}, varargin{k+1});
end
end

function add_mfb(path, code)
% add a MATLAB Function block and set its script
add_block('simulink/User-Defined Functions/MATLAB Function', path);
rt = sfroot;
chart = rt.find('-isa', 'Stateflow.EMChart', 'Path', path);
chart.Script = code;
end
