%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% nkoukis, January 2014              %
% Flight Mechanics                   %
% Project Work I                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PartII
% alternative implementation of the optimization procedure
% run PartI first, then run Part II and follow the instructions 

clearvars -except qdyn_lim2 alfa_lim2 alt_range mach_range energy_2dmat powers_2dmat

%% Variable Definitions
global sref tepsr xtpref xtpfpl rmfpl rmfull
global rmfuel iebk
global g
global gamma_f gamma_s  t_s t_f 

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

% varius configurations
gamma_f = 0;
t_f = 0;
v0 = 100; %initial velocity
h0 = 100; %initial height
mass0 = mass;
y0 = [v0; h0; mass0]; % initial estimation of v, h, m 

t_init = 0;
%% End-of-Flight Requirements
[rho_init,aspeed_init,temp_init,press_init] = stdatm(h0/1000);
fuel_req = 0.3 * rmfuel;
alt_req_1 = 11; %[km]
mach_req_1 = 1.5;
vspeed_req_1 = mach_req_1 * aspeed_init;

%% Initial Plotting 
figure(3); clf()
subplot(3,4,1:3)
hold on
grid on
title('Airplane velocity [m/s]');

subplot(3,4,5:7)
hold on
grid on
title('Altitude [km]');
subplot(3,4,9:11)
hold on
grid on
title('Fuel Remaining[kg]');


subplot(3,4,[4 8 12])
hold on
grid on
title('\gamma (t)');

figure(2);
hold all;
xlim([0 2]);
ylim([0 16]);

if mach15alt11
    plot(mach_req_1, alt_req_1, 'x', 'LineWidth', 3);
end

if plot_energy_lines && exist('mach_range', 'var')
    contour(mach_range, alt_range, energy_2dmat, 10);
end

annotation('textbox',...
    [0.13 0.74 0.19 0.19],...
    'String', 'Simulation has not began yet. Put figure in full screen',...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[1 1 0],...
    'LineWidth',2,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
    
%% Flight Simulation - Interactive way

% Flag initialization
stop = 0;

% Total matrices configuration
max_els = 10000;
total_cnt = 1;
v_total = zeros(max_els, 1);
h_total = zeros(max_els, 1);
rem_fuel_total = zeros(max_els, 1);

% write to file descriptor to keep the results
fid = fopen('output.txt', 'w');



%% loop
while stop == 0
    % times
    t_diff = input(sprintf('How long should I simulate for:\t'));
    while t_diff == 0
        t_diff = input('How long should I simulate for:\t');
    end
    
    t_s = t_f;
    t_f = t_s + t_diff;
    t_span = [t_s, t_f];
    gamma_s = gamma_f;
    
    % input gamma
    gamma_slope = input('Specify the slope of gamma function for the simulated time');
    gamma_f = gamma_slope * (t_f - t_s) + gamma_s;
    
    % od45e running
    [t_sim, y_out] = ode45(@rhs_eqns, t_span, y0);
    
    v = y_out(:, 1); h = y_out(:, 2) / 1000; rem_fuel = y_out(:, 3) - rmfpl; %h [km]
    
    % 'Initial conditions update
    v0 = v(length(v));
    h0 = h(length(h))*1000; %want it in meters
    mass0_all = y_out(:,3);
    mass0 = mass0_all(length(mass0_all));
    y0 = [v0, h0, mass0];
    
    %update mach - v_req
    [rho_init,aspeed_init,temp_init,press_init] = stdatm(h0/1000);
    vspeed_req_1 = mach_req_1 * aspeed_init;
    aspeed = zeros(size(h));
    for i = 1:length(h)
        [rho_dump, aspeed(i)] = stdatm(h(i));
    end
    mach = v ./ aspeed;
    
    % write time data to file descriptor
    fprintf(fid, '%f\t%f\n', t_f, gamma_fun_2(t_f));
    
    %% Plotting the results
    figure(3);
    
    % Velocity plot
    subplot(3,4,1:3)
    hold on
    h1 = plot(t_sim, vspeed_req_1, 'g', 'LineWidth', 2);
    h2 = plot(t_sim, v, 'b');
    legend(h1, 'Mach = 1.5');
    legend(h2, 'Velocity');
    
    % Altitude plot
    subplot(3,4,5:7)
    hold on
    plot(t_span, [alt_req_1, alt_req_1], 'g', 'LineWidth', 2)
    plot(t_sim, h, 'm')
    
    % Fuel remaining plot
    subplot(3,4,9:11)
    hold on
    plot(t_span, [fuel_req, fuel_req], 'r', 'LineWidth', 2)
    plot(t_sim, rem_fuel, 'g')
    
    % Control Variable Data Plotting
    subplot(3,4,[4 8 12])
    hold on
    t_range = linspace(t_span(1), t_span(2), 1000);
    plot(t_range, arrayfun(@gamma_fun_2, t_range));
    
    figure(2);
    curve_in_sep = plot(mach, h, 'm', 'LineWidth', 1.5);
    text1 = sprintf('\t v = %.2f m/s, mach = %.2f\n\t Alt = %.2f km\n\t Fuel Rem = %.2f m_{req}\n\t\n Sim.Time = %.2f s', ...
        v(length(v)), mach(length(mach)), h(length(h)), rem_fuel(length(rem_fuel))/fuel_req, t_init + t_f);
    
    annotation('textbox',...
        [0.13 0.74 0.19 0.19],...
        'String', text1,...
        'FontSize',14,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 0],...
        'LineWidth',2,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0.84 0.16 0]);
    
    %% update total matrices
    v_total(total_cnt:total_cnt + length(v) - 1) = v;
    h_total(total_cnt:total_cnt + length(v) - 1) = h;
    rem_fuel_total(total_cnt:total_cnt + length(v) - 1) = rem_fuel;
    
    % update total counter
    total_cnt = total_cnt + length(v);
    
    
    % Another iteration
    stop = input('Stop the simulation (1 / 0)');

end


% close file descriptor
fclose(fid);
%% Post Plotting
figure(2)
if exist('qdyn_lim2', 'var') && exist('alfa_lim2', 'var')
    legend([alfa_lim2, qdyn_lim2, curve_in_sep], {'Angle of attack < 15 degrees', ...
        'Max Dynamic Pressure', ...
        'Airplane Trajectory'});
end

%% PostProcessing
aspeed = zeros(size(h_total));
for i = 1:length(h_total)
    [rho_dump, aspeed(i)] = stdatm(h_total(i));
end
mach = v_total ./ aspeed;

%final mach number
h_final = h_total(length(h_total));
[rho_f,aspeed_f,temp_f,press_f] = stdatm(h_final); %given altitude compute these properties
v_final = v_total(length(v_total));
alt_final = h_total(length(h_total));
mach_final = v_total(length(v_total)) / aspeed_f;
fuel_final = rem_fuel_total(length(rem_fuel_total));


%% horizontal position of airplane
t_vals = linspace(t_init, t_f, length(v_total));
gamma_vals = arrayfun(@(t) cos(gamma_fun(t)), t_vals');
xe_f = trapz(t_vals,v.*gamma_vals);


%% Textbox reloaded
if exist('qdyn_lim2', 'var') && exist('alfa_lim2', 'var')
    legend([alfa_lim2, qdyn_lim2, curve_in_sep], {'Angle of attack < 15 degrees', ...
        'Max Dynamic Pressure', ...
        'Airplane Trajectory'});
    
    text1 = sprintf('Simulation is over.\n\t v = %.2f m/s, mach = %.2f\n\t Alt = %.2f km\n\t Fuel Rem = %.2f m_{req}\n\t Xe = %.2f km\n Sim.Time = %.2f s', ...
        v_final, mach_final, alt_final, fuel_final/fuel_req, xe_f*10^(-3), t_f);
    

end