%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% nkoukis, March 2015                %
% Flight Mechanics                   %
% Project - Part III                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project Work III - Part III Manuevers
%
% CONTENTS
% todo

% Remember that this model is for low speed but high alfa

clc;
clearvars

%% variable definitions - initializations
% globals
global xcg; xcg = 10.2;
global de dp
global V_list alfa_list nz_list output_i de_list dp_list
global path_to_inputs

% paths
path_to_plots = '../../Report/Drawings/MatlabFigures/';
path_to_logs = '../logfiles/';
path_to_inputs = '../inputfiles/';




% locals

newton_iter = 10;
converg_lim = 10^-4; % if under this bound break newton, already converged
l_step = 10^(-5); % linearization step
float_tol = 10^(-15); % tolerance for defining differnet floats
max_size = 10000; % max number of elements in array


% plotting
fontsize = 9;
linewidth = 1.2;

color_list;



%% Initialize global data structures

initfpl35
initrm6

%% Intial state
% h_init = 11000; % in meters for fplmod and for time integration
% h_init = 1000;
h_init = 10000;
% V_kmh = 500; 
V_kmh = 700;
V =  V_kmh/3.6; % in m/s -> this one is used.
theta_init = 0.0;
fuel_rem = 500;

%% Read list of values for de, dp
de_total = csvread(strcat(path_to_inputs, 'de_inputs_loop.csv'));
% de_total = csvread(strcat(path_to_inputs, 'de_inputs_cobra.csv'));
timeforde = de_total(:, 1);
de_values = de_total(:, 2).*pi/180;

dp_total = csvread(strcat(path_to_inputs, 'dp_inputs_loop.csv'));
% dp_total = csvread(strcat(path_to_inputs, 'dp_inputs_cobra.csv'));
timefordp = dp_total(:, 1);
dp_values = dp_total(:, 2);

% de_init = -0.0544410758220454;
de_init = 0;
dp_init = 0.288038080583802;

%% Trim the state
h_init_km = h_init/1000;
[rho,aspeed,temp,press] = stdatm(h_init_km);
machset = V/aspeed;


% gamma = 0, psi = 0, so..
% beta = 0 that's why I can write these formulas
u_init = V*cos(theta_init);
w_init = V*sin(theta_init);


x=[
    u_init              % u (m/s)
    w_init              % w (m/s)
    0.0                 % q (rad/s)
    theta_init          % theta (rad)
    0.0                 % Distance (m)
    h_init                % alt (m)
    fuel_rem            % fuel (kg)
    de_init             % de (rad)
    dp_init             % dp
    ];

[xdot]=fplmod(0,x);

% Example uses illustrated below

%% Linearization
J = linearize_mdl(x,l_step);

%% Trim loop

% Only u, v, theta, de, dp change
ivar=[1 2 4 8 9];   % Selects variables to change during trim iteration
% udot, wdot, qdot, altdot, rmach
ifun=[1 2 3 6 11];  % Selects functions that should be zero in trim
xtrim=[u_init, w_init, theta_init, de_init, dp_init]'; % Initial guess

%% Iterate until convergence
% Newton-Raphson method used

for iter=1:newton_iter
    
    x(ivar)=xtrim;
    [xdot]=fplmod(0,x);
    
    J = linearize_mdl(x,l_step);
    
    ftrim=xdot(ifun)';
    ftrim(5)=ftrim(5)-machset;
    Jtrim=J(ifun,ivar);
    xtrim=xtrim-Jtrim\ftrim; % actual change of xtrim
    fprintf('|f| %e\n',norm(ftrim));
    
    % if I am already below the convergence limit set,
    % exit the loop
    if norm(ftrim) <= converg_lim
        break;
    end
end % End of Newton iteration
disp('***Done***');



%% Define the inputs
% initial input
de_init = x(end-1);
dp_init = x(end);

% replace the zeros with the trim values
de_values(de_values == 0) = de_init;
dp_values(dp_values == 0) = dp_init;

% define the de, dp as interpolations of the given data
de = @(t) interp1(timeforde, de_values, t);
dp = @(t) interp1(timefordp, dp_values, t);

%% Time integration scheme  - setup
% x is the initial state
x_init = x(1:end-2); % exclude de, dp from the state

% time interval for the integration
t_start = timeforde(1); % [s]
t_end = timeforde(end); % [s]
t_interval = [t_start, t_end];

V_list = zeros(1, max_size);
alfa_list = zeros(1, max_size);
nz_list = zeros(1, max_size);
de_list = zeros(1, max_size);
dp_list = zeros(1, max_size);

output_i = 1;

%% Time integration scheme  - execution

% call ode45 with the fplmod_wrapped function for evaluating the rhs of the
% dot equations
[time,y_ode] = ode45(@fplmod_wrapped, t_interval, x_init);
% ode45(@fplmod_wrapped, t_interval, x_init)
%% Postprocessing

% state variables - extraction
u_list = y_ode(:, 1);
w_list = y_ode(:, 2);
q_list = y_ode(:, 3);
theta_list = y_ode(:, 4);
x_list = y_ode(:, 5);
h_list = y_ode(:, 6);
m_list = y_ode(:, 7);

% output variables - filter out the zeros
V_list = V_list(1:output_i - 1);
alfa_list = alfa_list(1:output_i - 1);
nz_list = nz_list(1:output_i - 1);
de_list = de_list(1: output_i - 1);
dp_list = dp_list(1: output_i - 1);

time_output = linspace(t_start, t_end, length(alfa_list));
%% Plotting the results
figure(1); clf()

% Distance
subplot(6,2,1)
hold on; grid on;
plot(time, x_list, 'Color', colors_array(1, :),...
    'LineWidth', linewidth);
title('Distance Covered [m]', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);

% Altitude
subplot(6,2,2)
hold on; grid on;
plot(time, h_list, 'Color', colors_array(3, :),...
    'LineWidth', linewidth);
title('Altitude [m]', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);

% Theta angle
subplot(6,2,3)
hold on; grid on;
plot(time, theta_list.*180/pi, 'Color', colors_array(5, :),...
    'LineWidth', linewidth);
title('Attitude angle [deg]', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);

% Alfa angle
subplot(6,2,4)
hold on; grid on;
plot(time_output, alfa_list.*180/pi, 'Color', colors_array(7, :),...
    'LineWidth', linewidth);
title('Angle of Attack [deg]', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);

% Velocity
subplot(6,2,5)
hold on; grid on;
plot(time_output, V_list, 'Color', colors_array(9, :),...
    'LineWidth', linewidth);
title('Velocity [m/s]', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);

% Mass
subplot(6,2,6)
hold on; grid on;
plot(time, m_list, 'Color', colors_array(11, :),...
    'LineWidth', linewidth);
title('Mass [kg]', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);

% load factor
subplot(6,2,7)
hold on; grid on;
plot(time_output, nz_list, 'Color', colors_array(12, :),...
    'LineWidth', linewidth);
title('Load Factor', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);

% Elevator setting
subplot(6,2,9:10)
hold on; grid on;
plot(time_output, de_list*180/pi, 'Color', colors_array(13, :),...
    'LineWidth', linewidth);
ylabel('Elev. Setting [deg]', 'FontSize', fontsize);
title('Inputs', 'FontSize', fontsize+2, 'Color', 'red');


% Thrust level
subplot(6,2,11:12)
hold on; grid on;
plot(time_output, dp_list, 'Color', colors_array(15, :),...
    'LineWidth', linewidth);
ylabel('Thrust level', 'FontSize', fontsize);
xlabel('Time [s]', 'FontSize', fontsize);


% title for states_outptuts
p=mtit('States & Outputs',...
    'fontsize',fontsize+2,'color',[1 0 0]);

% save the figure
print(gcf, '-depsc2',[path_to_plots,'cobra_status']);
saveas(gcf, [path_to_plots, 'cobra_status'], 'fig');

figure(2); clf()
hold on; grid on;
title('Aircraft trajectory', 'FontSize', fontsize+2, 'Color', 'red');
xlabel('Distance [m]', 'FontSize', fontsize);
ylabel('Altitude [m]', 'FontSize', fontsize);
plot(x_list, h_list, 'Color', colors_array(15, :),...
    'LineWidth', linewidth);

% save the figure
print(gcf, '-depsc2',[path_to_plots,'cobra1']);
saveas(gcf, [path_to_plots, 'cobra1'], 'fig');