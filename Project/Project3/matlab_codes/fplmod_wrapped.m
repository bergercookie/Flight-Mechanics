function [xdot_ode] = fplmod_wrapped(tval, xstate)
% The purpose of FPLMOD_WRAPPED is to pass valid input to fplmod.m and
% filter its output so that *only the derivatives of the ode states are
% returned*.

% INPUT:
%     tval: time instance
%     xstate: u, w, q, theta, x, h, m, de, dp
% OUTPUT:
%     xdot_ode: udot, wdot, qdot, thetadot, xdot, hdot, mdot

% The function can also add to global lists the output variables gathered
% from fplmod, so that they may later be processed by the caller.

%% Variable initialization
global de dp
global alfa_list V_list nz_list % ouput variables needed
global de_list dp_list
global output_i

%% Wrapping fplmod

% xstate should be compatible with the ode standard - same "type" as xdot_ode
% add the control variables 
x = vertcat(xstate, de(tval), dp(tval));
[xdot]=fplmod(tval,x);


%% Aditional Processing
state_range = 1:7;
output_range = length(state_range)+1:length(xdot);

% extract the output vars
output_vars = xdot(output_range);
V = output_vars(1);
alfa = output_vars(2);
% vcal = output_vars(3);
% rmach = output_vars(4);
% CL = output_vars(5);
nz = output_vars(6);

% pass them into the global lists
V_list(output_i) = V;
alfa_list(output_i) = alfa;
nz_list(output_i) = nz;

de_list(output_i) = de(tval);
dp_list(output_i) = dp(tval);

% update the arrays index
output_i = output_i + 1;
%% Return the derivatives - ode45 compatible - must be a column vector
xdot_ode = xdot(state_range)';
end