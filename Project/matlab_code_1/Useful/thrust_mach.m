%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% nkoukis, January 2014              %
% Flight Mechanics                   %
% Project Work I                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global sref tepsr xtpref xtpfpl rmfpl rmfull qmax
global e_t_rad rmfuel alfa alfadeg iebk
global g

% parameters_partI
fuel_percent = 0.3;
gamma = 0;
gamma_rad = gamma * pi / 180;

% Initialize data structures

initfpl35
initrm6

g=9.81;

iebk=1;            % Engine at full afterburner
rmfuel=fuel_percent*rmfull; % Fuel tank is 60% full
mass=rmfpl+rmfuel;  % Total mass of aircraft [kg]

% Use fixed angle of attack or insert code to find equilibrium value

alfa=0.0873;            % Radians
alfadeg=alfa*180/pi;    % Degrees

e_t = 5; % Thrust angle for the Draken
e_t_rad = 5* pi / 180;

% Store everything in an array of structs
% Each struct shall contain the current altitude and one array
% corresponding to the excess thrusts calculated for the given mach nums
% array
alt_range = 0:1:16;
mach_range = linspace(0, 2, 30);

alt_curve = struct('alt', 0, 'thrusts', zeros(length(mach_range), 1));
curves = repmat(alt_curve, length(alt_range), 1);

% usefull arrays in printing - initialization
colors_array = {'r', 'g', 'b', 'c', 'k', 'y'};
styles_array = {'-', '--',  '-.', ':', 'o-'};
plot_config = cell(length(alt_range), 1);

for i = 1:length(alt_range)
    temp_array = zeros(1, length(mach_range));
    curves(i).alt = alt_range(i);
    for j = 1:length(mach_range)
        excess_thrust = texcess35(alt_range(i), mach_range(j), mass, gamma_rad);
        temp_array(j) = excess_thrust;
    end
    curves(i).thrusts = temp_array;
    
    % plotting configuration
    color_pos = rem(i-1, length(colors_array)-1)+1; the_color = colors_array{color_pos};
    style_pos = max(1, ceil(i/length(colors_array))); the_style = styles_array{style_pos};
    plot_config{i} = strcat(the_color, the_style);
end


% Plotting all the curves
clf()
hold on
h1 = plot(mach_range, curves(1).thrusts, plot_config{1}, 'LineWidth', 2); % first and last bolted
for i = 2:length(curves)-1
    plot(mach_range, curves(i).thrusts, plot_config{i});
end
h2 = plot(mach_range, curves(length(curves)).thrusts, plot_config{length(curves)}, 'LineWidth', 2);

h1_text = strcat('altitude = ', num2str(curves(1).alt), 'km');
h2_text = strcat('altitude = ', num2str(curves(length(curves)).alt), 'km');
legend([h1, h2], {h1_text, h2_text});
grid on
ylabel('Excess Thrust [N]');
xlabel('Mach');
limsy=get(gca,'YLim');
set(gca,'Ylim',[0, limsy(2)]);

% What about the envelope limits?