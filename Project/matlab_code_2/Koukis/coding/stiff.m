%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% Flight Mechanics                   %
% Group D, February 2015             %
% Project Work II                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing stuff
clearvars;
clc;

%% Supplementary functions

pattern1 = 'a0'; pattern2 = 'a00';
a0exists = @(str) strfind2(str, pattern1, pattern2);

%% Variable definitions
global path_to_experimental suffix
global c_mech sref span Ix

% Estimate torsional stiffness of the flexible sting
%% General aircraft data
L=0.81;         % Length in meters
d=0.01;         % Diameter in meters
E=206*1000^3;   % Modulus of elasticity for this steel type
nu=0.3;         % Poisson ratio
G=E/(2*(1+nu)); % Shear modulus
K=pi/2*(d/2)^4; % 
span= 9.4;
k=(G*K/L);      % Estimated stiffness
kexp=98.2951;   % From an experiment
sref=50.0;            % Reference wing area in m*m
Ix=13050;   % kgm2 % with pilot & fuel

%% Varius variable configuration
path_to_experimental = '../experimental_data/';
suffix = '.log';
max_len = 1000; % intial size of array with unknown number of elements at start
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

%% Initial plotting of the curve
figure(1); clf()
plot(time,vibration,'b',time,y_fit,'r');
h_legend = legend('exp data','fit curve');
set(h_legend,'FontSize',12);
xlabel('time (s)', 'FontSize', 12);
ylabel('Angular Velocity - P [deg/sec]')

% What's this?? *todo*
% % Mass moment of inertia
% Ix=k/(omega^2+damping^2);

%% Calculation of Clp

folds = dir(path_to_experimental);
fold_offset = 3; % start from 4th element of folds

% arrays for storing the v-n values
v_vals = zeros(1, max_len);
n_vals = zeros(1, max_len);
a0_cnt = 1;

% for every group
for i = fold_offset+1:numel(folds)
    g_name = folds(i).name
    letter = g_name(end); % get the group letter
    the_path = [path_to_experimental, g_name, '/'];
    logfiles = dir(the_path);
    log_offset = 2;
    
    % for every a
    for i_log = log_offset+1:numel(logfiles)
        name = logfiles(i_log).name(1:end-4); % without the .log suffix
        
        a = parse_angle(name);
        % for aoa = 0
        if a == 0
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
            figure(1); clf()
            plot(time,vibration,'b',time,y_fit,'r')
            h_legend = legend('exp data','fit curve');
            set(h_legend,'FontSize',12);
            xlabel('time (s)', 'FontSize', 12);
            ylabel('Angular Velocity - P [deg/sec]')
            title(sprintf('Approximation Results - aoa = %.2f, airspeed = %.2f',...
                a, vel), 'FontSize', 12);
            
            % store velocity, damping
            v_vals(a0_cnt) = vel;
            n_vals(a0_cnt) = damping;
            
            % update the a0 counter, a0/a00 logfile is found
            a0_cnt = a0_cnt + 1;
        end
    end
end
v_vals = v_vals(1:find(n_vals,1,'last'));
n_vals = n_vals(1:find(n_vals,1,'last'));
    
% sort the arrays with regards to v_vals
[v_vals, sort_index] = sort(v_vals);
n_vals = n_vals(sort_index);

% initial estimation of coefficients
coeff0= [0, 0.05];

%% Implementation of the least squares method
coeff = lsqcurvefit(@linear_model,coeff0, v_vals, n_vals);
% extract the coefficients
y0 = coeff(1); slope = coeff(2);
% I got the magnitude of the slope, the sign is negative
slope = -slope;

damp_fit = @(v) y0 + slope *v;
v_vals_tmp = linspace(v_vals(1), max(v_vals), 1000);
aprox_vals = arrayfun(damp_fit, v_vals_tmp);

%% Plotting the damping fit curve
figure(2); clf()
hold on;
plot(v_vals, n_vals, 'o');
plot(v_vals_tmp, aprox_vals, 'r-');
xlabel('Velocity', 'FontSize', 12);
ylabel('Damping Coefficient', 'FontSize', 12);
h_legend = legend('Experimental Data', 'Approximation Curve');
set(h_legend,'FontSize',12);
title('Velocity - Damping Curve', 'FontSize', 12);

%% Get the Clp coefficient
clp = extract_clp(slope);

str = sprintf('Clp as calculated from the a = 0 experimental data:\n\tClp = %.3f', clp);
print_result(str)
