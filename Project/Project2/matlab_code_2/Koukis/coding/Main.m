%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% Flight Mechanics                   %
% Group D, February 2015             %
% Project Work II                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flight Mechanics, Project work 2. Main file
%
% Main file for computing the Clp Clb coefficients for the aircraft model.
% Reads the data out of the experimental_data folder for each one of the
% teams. (Do not change the location of that folder).


%% Initializing stuff
clearvars;
clc;

%% Supplementary functions

pattern1 = 'a0'; pattern2 = 'a00';
a0exists = @(str) strfind2(str, pattern1, pattern2);

%% Variable definitions
global path_to_experimental suffix path_to_plot graph_type
global c_mech sref span Ix

% Estimate torsional stiffness of the flexible sting
%% General aircraft data
L=0.81;         % Length in meters
d=0.01;         % Diameter in meters
E=206*1000^3;   % Modulus of elasticity for this steel type
nu=0.3;         % Poisson ratio
G=E/(2*(1+nu)); % Shear modulus
K=pi/2*(d/2)^4; % 
span_real= 9.4;
k=(G*K/L);      % Estimated stiffness
kexp=98.2951;   % From an experiment
sref_real=50.0;            % Reference wing area in m*m
% Ix=13050;   % kgm2 % with pilot & fuel

% Scaling factors
span_scale = 1/14.7;
sref_scale = 1/14.7^2;
% model data
span = span_real*span_scale; 
sref = sref_real*sref_scale;

%% Varius variable configuration
path_to_experimental = '../experimental_data/';
path_to_plot = '../../../../Report/Drawings/MatlabFigures/'; 
graph_type = 'fig';
suffix = '.log';
max_len = 1000; % intial size of array with unknown number of elements at start
per_plots = 10; % How many logfiles per plot.
% how many points does it take to consider this particular angle
lower_calc_limit_clb = 15; 

linewidth = 1.5;
fontsize = 9;

%% Ground Vibration Test - Calculate the C_mech coefficient
% load the data - get meaningfull range
[time, vibration] = load_data('D', 'a00_u0');
% turning volts into p_values
vibration = v_to_p(vibration);
% calculate the coefficients of the approximation curve
[omega, damping, phase, amplitude] = approx_data(time, vibration);

% construct the general equation out of the coefficients
y_fun = @(t) amplitude.* exp(damping*t).* sin(omega*t + phase);
y_fit = amplitude.* exp(damping*time).* sin(omega*time + phase);

% see model for ODE and thus for the damping coefficient
c_mech = -damping;

% X - Moment of Inertia
Ix = kexp/omega^2;

%% Initial plotting of the curve
figure(1); clf()
set(gcf,'paperunits','centimeters') 
plot(time,vibration,'b',time,y_fit,'r');
h_legend = legend('exp data','fit curve');
set(h_legend,'FontSize',fontsize);
xlabel('time (s)', 'FontSize', fontsize);
ylabel('Angular Velocity - P [deg/sec]')

%% Calculation of Clp

folders = dir(path_to_experimental);
folder_offset = 3; % start from 4th element of folds

% arrays for storing the v-n values
v_vals_clp = zeros(1, max_len);
n_vals = zeros(1, max_len);
a0_cnt = 1;

% arrays for storing v-wmega values
v_vals_clb = zeros(1, max_len);
wm_vals = zeros(1, max_len);
a_vals = zeros(1, max_len);
anon0_cnt = 1;

plot_cnt = 1;

% for every group
for i = folder_offset+1:numel(folders)
    grp_name = folders(i).name
    letter = grp_name(end); % get the group letter
    the_path = [path_to_experimental, grp_name, '/'];
    logfiles = dir(the_path);
    log_offset = 3;
    
    % for every a
    logs_range = log_offset+1:numel(logfiles);
    for i_log = logs_range
        name = logfiles(i_log).name(1:end-4); % without the .log suffix
        
        a = parse_angle(name);
        disp(sprintf('Loading file: %s', [the_path, name, suffix]));
        
        % What's the velocity
        vel = parse_vel(name);
        
        [time, vibration] = load_data(letter, name);
        % turning volts into p_values
        vibration = v_to_p(vibration);
        % calculate the coefficients of the approximation curve
        [omega, damping, phase, amplitude] = approx_data(time, vibration);
        
        % construct the general equation out of the coefficients
        y_fit = amplitude.* exp(damping*time).* sin(omega*time + phase);
        
        % plot data to verify the curve follows the data
%         if ~rem(i_log, per_plots)
%             figure(1); clf()
%             plot(time,vibration,'b',time,y_fit,'r')
%             set(gcf,'paperunits','centimeters')
%             h_legend = legend('exp data','fit curve');
%             set(h_legend,'FontSize',fontsize+2);
%             xlabel('time (s)', 'FontSize', fontsize+2);
%             ylabel('Angular Velocity - P [deg/sec]','FontSize', fontsize+2);
%             title(sprintf('Approximation Results - aoa = %.2f [deg], airspeed = %.2f [m/s]',...
%                 a, vel), 'FontSize', fontsize+4);
%             plot_name = sprintf('datafit_clp%d', plot_cnt);
%             saveas(gcf, strcat(path_to_plot, plot_name), graph_type);
%             print(gcf, '-depsc2',[path_to_plot,plot_name]);
%             plot_cnt = plot_cnt + 1;
%         end
        % gather values for Clp
        if a == 0
            % store velocity, damping
            v_vals_clp(a0_cnt) = vel;
            n_vals(a0_cnt) = damping;
            
            % update the a0 counter, a0/a00 logfile is found
            a0_cnt = a0_cnt + 1;
        
        % gather values for Clb
        else
            v_vals_clb(anon0_cnt) = vel;
            wm_vals(anon0_cnt) = omega;
            a_vals(anon0_cnt) = a;
            anon0_cnt = anon0_cnt + 1;
        end
        
    end
end

%% Processing of Clp
v_vals_clp = v_vals_clp(1:find(n_vals,1,'last'));
n_vals = n_vals(1:find(n_vals,1,'last'));

% sort the arrays with regards to v_vals_clp
[v_vals_clp, sort_index] = sort(v_vals_clp);
n_vals = n_vals(sort_index);


% initial estimation of coefficients
coeff0= [0, 0.05];

%% Implementation of the least squares method
coeff = lsqcurvefit(@linear_model,coeff0, v_vals_clp, n_vals);
% extract the coefficients
y0 = coeff(1); slope = coeff(2);
% I got the magnitude of the slope, the sign is negative
slope = -slope;

damp_fit = @(v) y0 + slope *v;
v_vals_tmp = linspace(v_vals_clp(1), max(v_vals_clp), 1000);
aprox_vals = arrayfun(damp_fit, v_vals_tmp);

%% Get the Clp coefficient
clp = extract_clp(slope);

str = sprintf('Clp as calculated from the a = 0 experimental data:\n\tClp = %.3f', clp);
print_result(str)

%% Plotting the damping fit curve
figure(2); clf()
hold on;
plot(v_vals_clp, n_vals, 'o');
plot(v_vals_tmp, aprox_vals, 'r-');
set(gcf,'paperunits','centimeters') 
xlabel('Velocity [m/s]', 'FontSize', fontsize);
ylabel('Damping Coefficient - n', 'FontSize', fontsize);
h_legend = legend('Experimental Data', 'Approximation Curve');
set(h_legend,'FontSize',fontsize);
title('Velocity - Damping Curve', 'FontSize', fontsize+2);
% saveas(gcf, strcat(path_to_plot, 'v_n_graph'), graph_type);
% print(gcf, '-depsc2',[path_to_plot, 'v_n_graph']);


%% Processing for Clb
v_vals_clb = v_vals_clb(1:find(wm_vals,1,'last'));
a_vals = a_vals(1:find(wm_vals,1,'last'));
wm_vals = wm_vals(1:find(wm_vals, 1, 'last'));

% take the curve with regards to wmega^2
wm_vals_2 = wm_vals.^2;
% initial estimation of coefficients
%see quadratic_model, x(1)*t^2+x(2)*t+x(3)
coeff0= [-1, -1, 0.05]; 


[a_vals_unique, a_vals_frequency] = count_unique(a_vals);

% consider only the alphas that appear more than lower_calc_limit_clb times
alphas2consider = a_vals_unique(a_vals_frequency > lower_calc_limit_clb);
clb_list = zeros(1, length(alphas2consider));

figure(3); clf()
eachline = 3; 
if length(alphas2consider) <= eachline
    nplots_y = length(alphas2consider);
else
    nplots_y = eachline;
end
nplots_x = ceil(length(alphas2consider)/eachline);

% assemble the results string 
results_strs = cell(length(alphas2consider), 1);
results_strs{1} = sprintf('**Calculated Clb values for different values of aoa**:\n');

for a_i = 1:length(alphas2consider)
    a = alphas2consider(a_i);
    % how many of these values are there?
    a__equals_a_vals = (a == a_vals);
    
    v_vals_specific = v_vals_clb(a__equals_a_vals);
    wm_vals_2_specific = wm_vals_2(a__equals_a_vals);
    
    % filter out 'wrong' points using mean value + standard deviation
    valid_range = [mean(wm_vals_2_specific) - 3*std(wm_vals_2_specific), ...
        mean(wm_vals_2_specific) + 3*std(wm_vals_2_specific)];
    
    wm_vals_2_specific = wm_vals_2_specific(...
        wm_vals_2_specific > valid_range(1) &...
        wm_vals_2_specific < valid_range(2));
    v_vals_specific = v_vals_specific(...
        wm_vals_2_specific > valid_range(1) &...
        wm_vals_2_specific < valid_range(2));
    
    %% Implementation of the least squares method
    coeff = lsqcurvefit(@quadratic_model,coeff0, v_vals_specific, wm_vals_2_specific);
    % extract the coefficients
    y2 = coeff(1); y1 = coeff(2); y0 = coeff(3);
    
    damp_fit = @(v) y0 + y1 *v + y2*v^2;
    v_vals_tmp = linspace(v_vals_specific(1), max(v_vals_specific), 1000);
    aprox_vals = arrayfun(damp_fit, v_vals_tmp);
    
    
    %% Plotting the fit curves altogether
    subplot(sprintf('%d%d%d', nplots_x, nplots_y, a_i));
    hold on; grid on;
    h1 = plot(v_vals_specific, wm_vals_2_specific, 'o');
    h2 = plot(v_vals_tmp, aprox_vals, 'r-');
    xlabel('Velocity [m/s]', 'FontSize', fontsize);
    ylabel('Wn^2', 'FontSize', fontsize);
    title(sprintf('aoa = %.2f [deg]', a), 'FontSize', fontsize);


    
    %% Get the Clp coefficient
    clb = extract_clb(y2, clp, a);
    
    % store current clb
    clb_list(a_i) = clb;
    
    % update the results string
    update_str = sprintf('aoa = %2.2f deg:\tClb = %.4f\n', a, clb);
    results_strs{a_i+1} = update_str;
end
% saveas(gcf, strcat(path_to_plot, 'v_wmega2_graph'), graph_type);
% print(gcf, '-depsc2',[path_to_plot, 'v_wmega2_graph']);


% print the results for Clb
for one_str_i = 1:length(results_strs)
    disp(results_strs{one_str_i})
end