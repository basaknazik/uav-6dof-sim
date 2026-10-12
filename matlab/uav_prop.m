function [T, Q] = uav_prop(Va, delta_t)
%UAV_PROP  Propeller thrust T [N] (along +x body) and aerodynamic torque Q [N m] from the motor/prop torque balance.
%   V_in = V_max*delta_t;  motor torque K_Q((V_in - K_V w)/R - i0) = prop torque  ->  a w^2 + b w + c = 0.
%#codegen
P = uav_params();
dt = min(max(delta_t, P.delta_t_min), P.delta_t_max);
V_in = P.V_max*dt;
D = P.D_prop;
a = P.rho*D^5*P.C_Q0/(4*pi^2);
b = P.rho*D^4*Va*P.C_Q1/(2*pi) + P.KQ*P.KV/P.R_motor;
c = P.rho*D^3*Va^2*P.C_Q2 - P.KQ*V_in/P.R_motor + P.KQ*P.i0;
if c >= 0
    w = 0;
else
    w = (-b + sqrt(b*b - 4*a*c))/(2*a);
end
T = P.rho*D^2*Va^2*P.C_T2 + P.rho*D^3*Va*P.C_T1*w/(2*pi) + P.rho*D^4*P.C_T0*w^2/(4*pi^2);
Q = P.rho*D^3*Va^2*P.C_Q2 + P.rho*D^4*Va*P.C_Q1*w/(2*pi) + P.rho*D^5*P.C_Q0*w^2/(4*pi^2);
end
