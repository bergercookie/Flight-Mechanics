function [out] = find_a(T, m, gamma)
% FIND_A computes a solving the nonlinear algebraic equation
% Tsin(a+e_t) + L(a) - mg*cos(gamma) = 0
% a is returned in radians [!]

global tepsr g sref
cla = evalin('caller', 'cla');
ca0 = evalin('caller', 'ca0');
qdyn = evalin('caller', 'qdyn');

init_estimation = to_rad(3); %3 degrees initial aoa estimation
L = @(a) qdyn*sref*cla*(to_degrees(a)-ca0);
eqn_a = @(a) T*sin(a + tepsr) + L(a) - m*g*cos(gamma);
a_found = fzero(eqn_a, init_estimation);

out = a_found;
end