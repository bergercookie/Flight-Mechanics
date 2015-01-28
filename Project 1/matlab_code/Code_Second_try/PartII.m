%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% nkoukis, January 2014              %
% Flight Mechanics                   %
% Project Work I                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PartII

clearvars

%% Variable Definitions
global sref tepsr xtpref xtpfpl rmfpl rmfull
global rmfuel iebk
global g

% parameters_partI
fuel_percent = 1;
gamma = 0;
gamma_rad = to_rad(gamma);

% Initialize data structures

initfpl35
initrm6

g = 9.81;

iebk=1;            % Engine at full afterburner
rmfuel=fuel_percent*rmfull;
mass=rmfpl+rmfuel;  % Total mass of aircraft [kg]

%% End-of-Flight Requirements
fuel_req = 0.3 * rmfuel;
alt_req = 11; %[km]
mach_req = 1.5;

%% Flight Simulation

t_s = 0;
t_f = 10; %how much will the duration last
t_span = [t_s, t_f];

v0 = 100; %initial velocity
h0 = 100; %initial height

y0 = [v0; h0; mass]; % initial estimation of v, h, m 
[t_sim, y_out] = ode45(@rhs_eqns, t_span, y0);

v = y_out(:, 1); h = y_out(:, 2) / 1000; rem_fuel = y_out(:, 3) - rmfpl; %h [km]


%% horizontal position of airplane
%todo

%% Plotting the results
figure(3); clf()
nplots = 3;

% Velocity plot
subplot(311)
hold on
grid on
title('Airplane velocity [m/s]');
% velocity limit plot
plot(t_sim, v, 'b')


% Altitude plot
subplot(312)
hold on
grid on
title('Altitude [km]');
plot(t_span, [alt_req, alt_req], 'g', 'LineWidth', 2)
plot(t_sim, h, 'm')

% Fuel remaining plot
subplot(313)
hold on
grid on
title('Fuel Remaining[kg]');
plot(t_span, [fuel_req, fuel_req], 'r', 'LineWidth', 2)
plot(t_sim, rem_fuel, 'g')

%% Printing Results

%final mach number
h_final = h(length(h));
[rho,aspeed,temp,press] = stdatm(h_final); %given altitude compute these properties
v_final = v(length(v));
alt_final = h(length(h));
mach_final = v(length(v)) / aspeed;
fuel_final = rem_fuel(length(rem_fuel));

disp('***At the end of the simulation***');
disp(sprintf('\t v = %.2f m/s,\n\t mach = %.2f\n\t Alt = %.2f km\n\t Fuel Rem = %.2f Required\n', v_final, mach_final, alt_final, fuel_final/fuel_req));

