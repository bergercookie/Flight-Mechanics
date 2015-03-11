%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% nkoukis, January 2014              %
% Flight Mechanics                   %
% Project Work I                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PartII
% standard version of the optimization part
% **FIRST** run part I to have the optical results of the SEP Graph
%
% An interactive version of the algorithm also exists in PartII_2.m

clearvars -except qdyn_lim2 alfa_lim2 alt_range mach_range energy_2dmat powers_2dmat

%% Variable Definitions
% globals
global sref tepsr xtpref xtpfpl rmfpl rmfull
global rmfuel iebk
global g
global gammavec

% various parameters
path_to_plots = '../../../Report/Drawings/MatlabFigures/';
linewidth = 1.5;
fontsize = 9;

%% General Configuration
plot_energy_lines = 1;

% parameters_partI
fuel_percent = 1;
mach15alt11 = 0; %go for the mach = 1.5, alt = 11km question

% Initialize data structures

initfpl35
initrm6

g = 9.81;

iebk=1;            % Engine at full afterburner
rmfuel=fuel_percent*rmfull;
mass=rmfpl+rmfuel;  % Total mass of aircraft [kg]


%% gamma func times & values
% more handy to be here than in the function
t1 = 33; val1 = 0;
t2 = 45; val2 = 0.56;
t3 = 130; val3 = 0.19;
t4 = 300; val4 = 0;
t5 = 350; val5 = -0.2;
t6 = 370; val6 = 0.;
t7 = 450; val7 = 0.02;
t8 = 580; val8 = 0.01;
t9 = 650; val9 = -0.2;
t10 = 670; val10 = 0;

gammavec = [0, 0;
    t1, val1;
    t2, val2;
    t3, val3;
    t4, val4;
    t5, val5;
    t6, val6;
    t7, val7;
    t8, val8;
    t9, val9;
    t10, val10;
    10000, 0];

%% Flight Simulation

t_s = 0;
t_f = 600; %how much will the duration last
t_span = [t_s, t_f];

v0 = 100; %initial velocity
h0 = 100; %initial height

y0 = [v0; h0; mass]; % initial estimation of v, h, m 
[t_sim, y_out] = ode45(@rhs_eqns, t_span, y0);

v = y_out(:, 1); h = y_out(:, 2) / 1000; rem_fuel = y_out(:, 3) - rmfpl; %h [km]

%% PostProcessing
rho = zeros(size(h)); aspeed = zeros(size(h));
for i = 1:length(h)
    [rho(i), aspeed(i)] = stdatm(h(i));
end
mach = v ./ aspeed;

%final mach number
h_final = h(length(h));
[rho_f,aspeed_f,temp_f,press_f] = stdatm(h_final); %given altitude compute these properties
v_final = v(length(v));
alt_final = h(length(h));
mach_final = v(length(v)) / aspeed_f;
fuel_final = rem_fuel(length(rem_fuel));

%% End-of-Flight Requirements
fuel_req = 0.3 * rmfuel;
alt_req_1 = 11; %[km]
mach_req_1 = 1.5;
vspeed_req_1 = mach_req_1 * aspeed;

%% horizontal position of airplane
t_vals = linspace(t_s, t_f, length(v));
gamma_vals = arrayfun(@(t) cos(gamma_fun(t)), t_vals');
xe_f = trapz(t_vals,v.*gamma_vals);

%% Plotting the results
figure(3); clf()

% Velocity plot
subplot(3,4,1:3)
hold on
grid on
title('Airplane velocity [m/s]', 'FontSize', fontsize);
plot(t_sim, vspeed_req_1, 'g', 'LineWidth', linewidth);
plot(t_sim, v, 'b')
h_legend = legend('Mach = 1.5', 'Velocity');
set(h_legend, 'FontSize', fontsize);


% Altitude plot
subplot(3,4,5:7)
hold on
grid on
title('Altitude [km]', 'FontSize', fontsize);
plot(t_span, [alt_req_1, alt_req_1], 'g', 'LineWidth', linewidth)
plot(t_sim, h, 'm')

% Fuel remaining plot
subplot(3,4,9:11)
hold on
grid on
title('Fuel Remaining[kg]', 'FontSize', fontsize);
plot(t_span, [fuel_req, fuel_req], 'r', 'LineWidth', linewidth)
plot(t_sim, rem_fuel, 'g')

% Control Variable Data Plotting
subplot(3,4,[4 8 12])
hold on
grid on
title('\gamma (t)', 'FontSize', fontsize);
t_range = linspace(t_span(1), t_span(2), 1000);
plot(t_range, arrayfun(@gamma_fun, t_range));

figure(2);
if mach15alt11
    plot(mach_req_1, alt_req_1, 'x', 'LineWidth', linewidth);
end

if plot_energy_lines && exist('mach_range', 'var')
    contour(mach_range, alt_range, energy_2dmat, 10);
end

curve_in_sep = plot(mach, h, 'm', 'LineWidth', linewidth);
if exist('qdyn_lim2', 'var') && exist('alfa_lim2', 'var')
    h_legend = legend([alfa_lim2, qdyn_lim2, curve_in_sep], {'Angle of attack < 15 degrees', ...
        'Max Dynamic Pressure', ...
        'Airplane Trajectory'});
    set(h_legend, 'FontSize', fontsize);
    
    text1 = sprintf('\t v = %.2f m/s, mach = %.2f\n\t Alt = %.2f km\n\t Fuel Rem = %.2f m_{req}\n\t Xe = %.2f km\n Sim.Time = %.2f s', ...
    v_final, mach_final, alt_final, fuel_final/fuel_req, xe_f*10^(-3), t_f);

    annotation('textbox',...
    [0.13 0.80 0.13 0.12],...
    'String', text1,...
    'FontSize',fontsize,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[1 1 0],...
    'LineWidth',linewidth,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
end


