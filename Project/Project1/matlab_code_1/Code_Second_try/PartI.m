 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% nkoukis, January 2014              %
% Flight Mechanics                   %
% Project Work I                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This is the implementation of part I of the 1st project work in Flight
% Mechanics. The purpose of the module is to use TexSep.m to calculate the
% excess thrust and power of the Draken model and make visual graphs out of
% the results. 

clearvars

path_to_plots = '../../../Report/Drawings/MatlabFigures/';
linewidth = 1.5;
fontsize = 9;

%% Global variables definitions
global sref tepsr xtpref xtpfpl rmfpl rmfull qmax
global rmfuel iebk
global g
global qdyn_max qdyn_index qdyn_mach qdyn_alt qdyn_tex

%% parameters_partI
fuel_percent = 0.3;
gamma = 0;
gamma_rad = to_rad(gamma);

%% Initialize data structures

initfpl35
initrm6

g=9.81;

iebk=1;            % Engine at full afterburner
rmfuel=fuel_percent*rmfull; % Fuel tank is 60% full
mass=rmfpl+rmfuel;  % Total mass of aircraft [kg]

%% Basic Structures Initialization
% Store everything in an array of structs
% Each struct shall contain the current altitude and 3 arrays
% corresponding to the excess thrusts [Tex], the specific excess powers [Sep] 
% and the aoa angles calculated for the given mach nums array
alt_range = 0:1:17;
mach_range = linspace(0, 2, 30); % same for all curves

% limits configuration 
alfa_index = 1; qdyn_index = 1;
max_alfa =15;
v_max = 1350 / 3.6; %velocity limit
rho_0 = stdatm(0);
qdyn_max = 0.5 * rho_0 * v_max^2;


% alatitude, thrusts, powers
alt_curve = struct('alt', 0,...
    'thrusts', zeros(length(mach_range), 1), ...
    'powers', zeros(length(mach_range), 1), ...
    'alfas', zeros(length(mach_range), 1), ...
    'eq_energy', zeros(length(mach_range), 1));

% replicate style into an array
curves = repmat(alt_curve, length(alt_range), 1);

% usefull arrays in printing - initialization
colors_array = {'g', 'r', 'b', 'c', 'k', 'y'};
styles_array = {'-', '--',  '-.', ':', 'o-'};
plot_config = cell(length(alt_range), 1);

% equal energy lines

%% Call TexSep - Calculate Tex, Sep, alfa, ..
for i = 1:length(alt_range)
    temp_array_tex = zeros(1, length(mach_range)); %temporary array - storing Tex for cur altitude
    temp_array_sep = zeros(1, length(mach_range)); %temporary array - storing Sep for cur altitude
    temp_array_aoa = zeros(1, length(mach_range)); %temporary array - storing aoa for cur altitude
    temp_array_energy = zeros(1, length(mach_range)); %temporary array - storing eq_energy for cur altitude

    curves(i).alt = alt_range(i);
    for j = 1:length(mach_range)
        %[Tex, Sep] = TexSep(alt_range(i), mach_range(j), mass, gamma_rad);
        [Tex, Sep, alfa_val, energy] = TexSep(alt_range(i), mach_range(j), mass, gamma_rad);
        temp_array_tex(j) = Tex; temp_array_sep(j) = Sep; temp_array_aoa(j) = to_degrees(alfa_val);
        temp_array_energy(j) = energy;
        % store the limits - alfa
        if abs(alfa_val < to_rad(max_alfa)) < 0.1
            alfa_alt(alfa_index) = alt_range(i);
            alfa_mach(alfa_index) = mach_range(j);
            alfa_tex(alfa_index) = Tex;
            alfa_index = alfa_index + 1;
        end
    end
    curves(i).thrusts = temp_array_tex;
    curves(i).powers = temp_array_sep;
    curves(i).alfas = temp_array_aoa;
    curves(i).eq_energy = temp_array_energy;

    % plotting configuration
    color_pos = rem(i-1, length(colors_array)-1)+1; the_color = colors_array{color_pos};
    style_pos = max(1, ceil(i/length(colors_array))); the_style = styles_array{style_pos};
    plot_config{i} = strcat(the_color, the_style);
end


%% Plotting all the curves
figure(1); clf();
hold on
h1 = plot(mach_range, curves(1).thrusts, ...
    plot_config{1}, 'LineWidth', linewidth); % first and last bolted
for i = 2:length(curves)-1
    plot(mach_range, curves(i).thrusts, plot_config{i});
end
h2 = plot(mach_range, curves(length(curves)).thrusts, plot_config{length(curves)}, 'LineWidth', linewidth);

grid on
xlabel('Mach', 'FontSize', fontsize);
ylabel('Excess Thrust [N]', 'FontSize', fontsize);
set(gca,'Ylim',[-2000, 0.6*10^5]);

powers_2dmat = zeros(length(curves), length(mach_range));
energy_2dmat = zeros(length(curves), length(mach_range));
for i = 1:length(curves)
    powers_2dmat(i, :) = curves(i).powers;
    energy_2dmat(i, :) = curves(i).eq_energy;
end
figure(2); clf();
% You have to use the f***** contour!
hold on
% plotting the contour + label for each curve
[C,h] = contour(mach_range, alt_range, powers_2dmat, 0:5:200);
clabel(C, h);
ylabel('Altitude [km]', 'FontSize', fontsize);
xlabel('Mach', 'FontSize', fontsize);
grid on

%% Envelope limits calculations

%need to filter out the data to take the curve
alfa_mach2 = zeros(size(alfa_mach)); alfa_alt2 = zeros(size(alfa_alt));
qdyn_mach2 = zeros(size(qdyn_mach)); qdyn_alt2 = zeros(size(qdyn_alt));

new_i = 1;
for i = 2:length(alfa_alt)
    if alfa_alt(i) ~= alfa_alt(i-1)
        alfa_alt2(new_i) = alfa_alt(i-1); %take the previous
        alfa_mach2(new_i) = alfa_mach(i-1);
        new_i = new_i + 1;
    end
end

new_i = 1;
for i = 2:length(qdyn_alt)
    if qdyn_alt(i) ~= qdyn_alt(i-1)
        qdyn_mach2(new_i) = qdyn_mach(i); %take that one! I want the curve to the left
        qdyn_alt2(new_i) = qdyn_alt(i);
        new_i = new_i + 1;
    end
end


%% curve fitting through limits
% filter out the last zeros - alfa
 %the limit should be over 0.2 for mach - robustness reason
alfa_mach2_filt = alfa_mach2(alfa_mach2 ~= 0 & alfa_mach2 > 0.2);
alfa_alt2_filt = alfa_alt2(alfa_mach2 ~= 0 & alfa_mach2 > 0.2);
% polynomial fitting through data
poly_order = 2; x_range_alfa = linspace(0.1, 0.65, 100);
p_coefs = polyfit(alfa_mach2_filt, alfa_alt2_filt, poly_order);
f2_alfa = polyval(p_coefs, x_range_alfa);

% filter out the last zeros - qdyn
qdyn_mach2_filt = qdyn_mach2(qdyn_mach2 ~= 0);
qdyn_alt2_filt = qdyn_alt2(qdyn_mach2 ~= 0);
% polynomial fitting through data
poly_order = 2; x_range_qdyn = linspace(1.1, 2, 100);
p_coefs = polyfit(qdyn_mach2_filt, qdyn_alt2_filt, poly_order);
f2_qdyn = polyval(p_coefs, x_range_qdyn);

%% alfa - qdyn limits plotting - Saving the Plots
figure(1);% Tex
hold on
title('Excess Thrust Graph', 'FontSize', fontsize+2)
alfa_lim1 = plot(alfa_mach, alfa_tex, 'rx', 'LineWidth', linewidth);
qdyn_lim1 = plot(qdyn_mach, qdyn_tex, 'bx', 'LineWidth', linewidth);
h1_text = strcat('altitude = ', num2str(curves(1).alt), 'km');
h2_text = strcat('altitude = ', num2str(curves(length(curves)).alt), 'km');
h3_text = 'Angle of attack < 15 degrees';
h4_text = 'Max Dynamic Pressure';
h_legend = legend([h1, h2, alfa_lim1, qdyn_lim1], {h1_text, h2_text, h3_text, h4_text});
set(h_legend,'FontSize',fontsize);
% % save the plot
% print(gcf, '-depsc2',[path_to_plots,'Tex_Graph']);
% saveas(gcf, [path_to_plots, 'Tex_Graph'], 'fig');

figure(2); % SEP
hold on
title(sprintf('Specific Excess Power Graph'), 'FontSize', fontsize+2);
alfa_lim2 = plot(x_range_alfa, f2_alfa, 'r', 'LineWidth', linewidth);
% qdyn_lim2 = plot(qdyn_mach2_filt, qdyn_alt2_filt, 'bx', 'LineWidth', linewidth);
qdyn_lim2 = plot(x_range_qdyn, f2_qdyn, 'b', 'LineWidth', linewidth);
h_legend = legend([alfa_lim2, qdyn_lim2], {'Angle of attack < 15 degrees', ...
    'Max Dynamic Pressure'});
set(h_legend,'FontSize',fontsize);
limsy=get(gca,'YLim');
set(gca,'Ylim',[0, limsy(2)])
% % save the plot
% print(gcf, '-depsc2',[path_to_plots,'SEP_Graph']);
% saveas(gcf, [path_to_plots, 'SEP_Graph'], 'fig');
