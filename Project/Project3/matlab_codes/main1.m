%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology %
% School of Engineering Sciences     %
% nkoukis, March 2015                %
% Flight Mechanics                   %
% Project - Part III                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project Work III - Part I Equilibrium
% 
% CONTENTS
%     1a. Level Flight Trim Condition
%         Altitudes (0, 5000, 10000)m
%         Mach (0.1 - 0.7)
%     1b. Center of Gravity Influence
%         Altitude (5000m)
%     1c. Elevator per g
%         Altitudes (0, 5000, 10000)m
%         Mach (0.5)
%
% Remember that this model is for low speed but high alfa

clc;
clearvars -except states_array states_str

%% variable definitions - initializations

% Which part of the project to run
part1a = 0; part1b = 1; part1c = 0;

% paths
path_to_plots = '../../Report/Drawings/MatlabFigures/';
path_to_logs = '../logfiles/';

% globals 
global xcg;

% locals

mach_list = 0.1:0.05:0.7;
altm_list = [0, 5000, 10000]; % in meters
newton_iter = 10;
converg_lim = 10^-4; % if under this bound break newton, already converged
l_step = 10^(-5); % linearization step

fuel_rem = 500;
de_init = -0.0544410758220454;
dp_init = 0.288038080583802;

float_tol = 10^(-15); % tolerance for defining differnet floats

% plotting
fontsize = 12;
linewidth = 1.2;




%% Initialize global data structures

initfpl35
initrm6

% run it if explicitly told, or if part1b or c has to run and the states array
% isn't set
if part1a || ((part1b || part1c) && ~exist('states_array', 'var') && ~exist('states_str', 'var'))
    
    disp('******* PART Ia *********');
    % array for storing the equilibrium points
    states_str = struct('u', 0, 'v', 0, 'theta', 0, ...
        'de', 0, 'dp', 0,...
        'norm', 0, 'aspeed', 0, ...
        'nz', 0);
    % replicate style into an array
    states_array = repmat(states_str, length(altm_list), length(mach_list));
    
    %% Call model for trimmed point
    xcg=10.0;
    
    
    for altm_i = 1:length(altm_list)
        altm = altm_list(altm_i);
        % get the properties at the specific altitude
        altkm = altm/1000;
        [rho,aspeed,temp,press] = stdatm(altkm);
        for mach_i = 1:length(mach_list)
            machset = mach_list(mach_i);
            V = machset*aspeed; % get the magnitude of velocity
            
            % initial point
            theta_init = 0.0;
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
                altm                % alt (m)
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
            xtrim=[u_init, w_init, theta_init, 0, 0]'; % Initial guess
            
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
            
            % store the current state in the struct
            states_array(altm_i, mach_i).u = xtrim(1);
            states_array(altm_i, mach_i).w = xtrim(2);
            states_array(altm_i, mach_i).theta = xtrim(3);
            states_array(altm_i, mach_i).de = xtrim(4);
            states_array(altm_i, mach_i).dp = xtrim(5);
            states_array(altm_i, mach_i).aspeed = sqrt(xtrim(1)^2+xtrim(2)^2);
            states_array(altm_i, mach_i).norm = norm(ftrim);

        end
    end
    
    %% Plotting the figures
    
    % three curves per plot, one for each altitude
    figure(1); clf()
    % initialization of subplots, titles
    subplot(3,1,1);
    hold on; grid on;
    ylabel('Angle of Attack [deg]', 'FontSize', fontsize);
    xlabel('Airspeed [m/s]', 'FontSize', fontsize);
    subplot(3,1,2);
    hold on; grid on;
    ylabel('Elevator Setting [deg]', 'FontSize', fontsize);
    xlabel('Airspeed [m/s]', 'FontSize', fontsize);
    subplot(3,1,3)
    hold on; grid on;
    ylabel('Thrust Level', 'FontSize', fontsize);
    xlabel('Airspeed [m/s]', 'FontSize', fontsize);
    
    
    legend_strings = {length(altm_list), 1};
    legend_plots = [];
    colors = {'r', 'b', 'm'};
    for altm_i = 1:length(altm_list)
        
        altm = altm_list(altm_i);
        
        % extract data - convert to degrees for the graphs
        aspeed_list = extractfield(states_array(altm_i, :), 'aspeed');
        theta_list = extractfield(states_array(altm_i, :), 'theta')*180/pi;
        de_list = extractfield(states_array(altm_i, :), 'de')*180/pi;
        dp_list = extractfield(states_array(altm_i, :), 'dp');
        norm_list = extractfield(states_array(altm_i, :), 'norm');
        
        % filter out the ones for which newton didn't converge
        converged = norm_list < converg_lim;
        aspeed_list = aspeed_list(converged);
        theta_list = theta_list(converged);
        de_list = de_list(converged);
        dp_list = dp_list(converged);
        
        subplot(3,1,1);
        h_plot = plot(aspeed_list, theta_list, colors{altm_i}, 'LineWidth', linewidth);
        legend_plots =  [legend_plots, h_plot];
        legend_strings{altm_i} = sprintf('alt = %d m', altm);
        subplot(3,1,2);
        plot(aspeed_list, de_list, colors{altm_i}, 'LineWidth', linewidth);
        subplot(3,1,3);
        plot(aspeed_list, dp_list, colors{altm_i}, 'LineWidth', linewidth);
    end
    
    % axis limits
    subplot(3,1,1)
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[limsy(1), 40]);
    subplot(3,1,2)
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[-10, limsy(2)]);
    subplot(3,1,3)

    
    p=mtit('Level Flight Trim Conditions',...
        'fontsize',fontsize+4,'color',[1 0 0]);
    % refine title using its handle <p.th>
    set(p.th,'edgecolor',.5*[1 1 1]);
    h_legend = legend(legend_plots, legend_strings);
    set(h_legend, 'FontSize', fontsize);
    
    % save the plots
    print(gcf, '-depsc2',[path_to_plots,'equilibrium_conditions']);
    saveas(gcf, [path_to_plots, 'equilibrium_conditions'], 'fig');
end

if part1b
    
    disp('******* PART Ib *********');
    
    %% Computations
    xcg_list = [10.0, 10.30];
    % xcg = 10.0 - has already been caclculated
    % pick an altitude
    altm_i = 2; % 5000m
    altm = altm_list(altm_i);
    % extract data - convert to degrees for the graphs
    aspeed_list = extractfield(states_array(altm_i, :), 'aspeed');
    theta_list = extractfield(states_array(altm_i, :), 'theta')*180/pi;
    de_list = extractfield(states_array(altm_i, :), 'de')*180/pi;
    dp_list = extractfield(states_array(altm_i, :), 'dp');
    norm_list = extractfield(states_array(altm_i, :), 'norm');
    % filter out the ones for which newton didn't converge
    converged = norm_list < converg_lim;
    aspeed_list = aspeed_list(converged);
    theta_list = theta_list(converged);
    de_list = de_list(converged);
    dp_list = dp_list(converged);
    de_101 = de_list;
    aspeed_101 = aspeed_list;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % xcg = 10.3                                                     % 
    % **WARNING**                                                    %
    % Change every modification you make to the converging procedure %
    % above.                                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run the procedure again, find a de_list of values corresponding to
    % the mach_numbers given
    %find list of airspeeds
    xcg = xcg_list(2);
    
    % array for storing the equilibrium points
    % replicate style into an array
    states_xcg103 = repmat(states_str, length(xcg_list) - 1, length(mach_list));

    altkm = altm/1000;
    [rho,aspeed,temp,press] = stdatm(altkm);
    for mach_i = 1:length(mach_list)
        machset = mach_list(mach_i);
        V = machset*aspeed; % get the magnitude of velocity
        
        % initial point
        theta_init = 0.0;
        % gamma = 0, psi = 0, so..
        % beta = 0 that's why I can write these formulas
        u_init = V*cos(theta_init);
        w_init = V*sin(theta_init);
        
        
        x=[
            u_init    % u (m/s)
            w_init    % w (m/s)
            0.0                 % q (rad/s)
            theta_init          % theta (rad)
            0.0                 % Distance (m)
            altm                % alt (m)
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
        % udot, vdot, qdot, altdot, rmach
        ifun=[1 2 3 6 11];  % Selects functions that should be zero in trim
        xtrim=[u_init, w_init, theta_init, 0, 0]'; % Initial guess
        
        %% Iterate until convergence
        % Newton-Raphson method used
        
        for iter=1:newton_iter
            
            x(ivar)=xtrim;
            [xdot]=fplmod(0,x);
            
            J = linearize_mdl(x,l_step);
            
            ftrim=xdot(ifun)';
            ftrim(5)=ftrim(5)-machset;
            Jtrim=J(ifun,ivar);
            xtrim=xtrim-Jtrim\ftrim;
            fprintf('|f| %e\n',norm(ftrim));
            
            
            % if I am already below the convergence limit set,
            % exit the loop
            if norm(ftrim) <= converg_lim
                break;
            end
        end % End of Newton iteration
        disp('***Done***');
        
        % store the current state in the struct
        states_xcg103(mach_i).u = xtrim(1);
        states_xcg103(mach_i).w = xtrim(2);
        states_xcg103(mach_i).theta = xtrim(3);
        states_xcg103(mach_i).de = xtrim(4);
        states_xcg103(mach_i).dp = xtrim(5);
        states_xcg103(mach_i).aspeed = sqrt(xtrim(1)^2+xtrim(2)^2);
        states_xcg103(mach_i).norm = norm(ftrim);
    end
    % extract data - convert to degrees for the graphs
    aspeed_list = extractfield(states_xcg103, 'aspeed');
    theta_list = extractfield(states_xcg103, 'theta')*180/pi;
    de_list = extractfield(states_xcg103, 'de')*180/pi;
    dp_list = extractfield(states_xcg103, 'dp');
    norm_list = extractfield(states_xcg103, 'norm');
    
    % filter out the ones for which newton didn't converge
    converged = norm_list < converg_lim;
    aspeed_list = aspeed_list(converged);
    theta_list = theta_list(converged);
    de_list = de_list(converged);
    dp_list = dp_list(converged);
    
    de_103 = de_list;
    aspeed_103 = aspeed_list;
    
    %% Plotting
    figure(2); 
    % title, labels ..
    hold on; grid on;
    ylabel('Elevator Setting [deg]', 'FontSize', fontsize);
    xlabel('Airspeed [m/s]', 'FontSize', fontsize);
    title(sprintf('Influence of xcg - Altitude = %dm', altm), 'FontSize', fontsize+4);
    
    % plot for xcg = 10.0
    h_plot10 = plot(aspeed_101, de_101, 'r', 'LineWidth', linewidth);
    h_plot103 = plot(aspeed_103, de_103, 'b', 'LineWidth', linewidth);
    
    h_legend = legend([h_plot10, h_plot103], ...
        {sprintf('xcg = %.2fm', xcg_list(1)), sprintf('xcg = %.2fm', xcg_list(2))});
    set(h_legend, 'Fontsize', fontsize);
    
    % save the plots
    print(gcf, '-depsc2',[path_to_plots,'xcg_investigation']);
    saveas(gcf, [path_to_plots, 'xcg_investigation'], 'fig');
end

if part1c
    fid = fopen(strcat(path_to_logs, 'elevator_per_g.log'), 'w');
    fprintf(fid, '******* PART Ic *********\n\n');
    fprintf(fid, '***Elevator per g***\n');
    
    % values for dp
    de_list = -[2, 4, 6, 8, 10]*pi/180;
    
    % arrays for storing Dde and nz values
    de_per_alt = repmat(de_list, length(altm_list), 1);
    

    % find the elevator per g for mach = 0.5 and altitudes = 0, 5000, 10000
    
    n_vals = zeros(1, length(de_list));
    norm_vals = zeros(1, length(de_list));
    elevator_per_g = cell(length(altm_list), 1);

    for altm_i = 1:length(altm_list)
        machset = 0.5; % mach is constant at part 1c
        
        altm = altm_list(altm_i);
        % get the properties at the specific altitude
        altkm = altm/1000;
        [rho,aspeed,temp,press] = stdatm(altkm);
        V = machset*aspeed; % get the magnitude of velocity
        
        % Do not start the elevator settign from a predefined value but
        % rather vary it so that later you can take the mean of the values
        % gathered.
        for de_i = 1:length(de_list)
            de = de_list(de_i);
            % initial point
            theta_init = 0.0;
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
                altm                % alt (m)
                fuel_rem            % fuel (kg)
                de                  % de (rad)
                dp_init             % dp
                ];
            
            [xdot]=fplmod(0,x);
            
            % Example uses illustrated below
            
            %% Linearization
            J = linearize_mdl(x,l_step);
            
            %% Trim loop
            
            % Only u, v, theta, de, dp change
            ivar=[1 2 4 8 9];   % Selects variables to change during trim iteration
            % udot, vdot, qdot, altdot, rmach
            ifun=[1 2 3 6 11];  % Selects functions that should be zero in trim
            xtrim=[u_init, w_init, theta_init, de, 0]'; % Initial guess
            
            %% Iterate until convergence
            % Newton-Raphson method used
            
            for iter=1:newton_iter
                
                x(ivar)=xtrim;
                [xdot]=fplmod(0,x);
                
                J = linearize_mdl(x,l_step);
                
                ftrim=xdot(ifun)';
                ftrim(5)=ftrim(5)-machset;
                Jtrim=J(ifun,ivar);
                xtrim=xtrim-Jtrim\ftrim;
                fprintf('|f| %e\n',norm(ftrim));
                
                
                % if I am already below the convergence limit set,
                % exit the loop
                if norm(ftrim) <= converg_lim
                    break;
                end
            end % End of Newton iteration
            disp('***Done***');
            
            % store the current state in the arrays
            n_vals(de_i) = xdot(13);
            norm_vals(de_i) = norm(ftrim);
            

        end
    
        % filter out the ones for which newton didn't converge
        converged = norm_vals < converg_lim;
        norm_vals = norm_vals(converged);
        n_vals = n_vals(converged);
        
        de_per_alt_conv = de_per_alt(altm_i, :);
        
        Delta_de = diff(de_per_alt_conv);
        elevator_per_g{altm_i} = Delta_de./(n_vals(2:end)-1);

    end
    %% Printing the results
    for altm_i = 1:length(altm_list)
        altm = altm_list(altm_i);
        fprintf(fid, '*Altitude = %.1f\n', altm);
        fprintf(fid, 'Elevator per g = %.3f +- %.3f\n', mean(elevator_per_g{altm_i}),...
            std(elevator_per_g{altm_i}));
    end
    fclose(fid);

end