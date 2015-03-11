% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KTH, Royal Institute of Technology   %
% School of Engineering Sciences       %
% nkoukis, March 2015                  %
% Flight Mechanics                     %
% Project - Part III                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project Work III - Part II Stability

clc;

%% variable definitions - initializations

% globals 
global xcg; xcg = 10.1;
% global path_to_plots path_to_logs

% paths
path_to_plots = '../../Report/Drawings/MatlabFigures/';
path_to_logs = '../logfiles/';


% locals

mach_list = 0.1:0.07:0.7;
altm_list = [0, 5000, 10000]; % in meters
newton_iter = 10;
converg_lim = 10^-4; % if under this bound break newton, already converged
l_step = 10^(-5); % linearization step

fuel_rem = 500;
de_init = -0.0544410758220454;
dp_init = 0.288038080583802;

float_tol = 10^(-15); % tolerance for defining differnet floats

% plotting configuration

% inserts a color list in RGB form to work with
color_list;

fontsize = 9;
linewidth = 1.5;

% the longitudinal stability is determined by the following equations (see
% flpmod.m' return arguements)
% u, w, q, theta
stab_indices = [1, 2, 3, 4];


%% Initialize global data structures

initfpl35
initrm6

%% Iterations scheme
figure(1); clf()
% file discriptor for the Eigenvalue analysis output
fid = fopen(strcat(path_to_logs, 'eigenvalue_analysis.log'), 'w');

for altm_i = 1:length(altm_list)
    % set the current subplot
    subplot(length(altm_list), 1, altm_i);
    hold on;
    
    altm = altm_list(altm_i);
    fprintf('**Computing for altitude = %d m**\n', altm);
    % get the properties at the specific altitude
    altkm = altm/1000;
    [rho,aspeed,temp,press] = stdatm(altkm);
    
    % legend information storing
    h_list = zeros(1, length(mach_list));
    h_strs = cell(1, length(mach_list));
    
    % index for converged data
    converged_i = 1;
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
        
        % Only u, v, theta, de, dp change
        ivar=[1 2 4 8 9];   % Selects variables to change during trim iteration
        % udot, wdot, qdot, altdot, rmach
        ifun=[1 2 3 6 11];  % Selects functions that should be zero in trim
        xtrim=[u_init, w_init, theta_init, 0, 0]'; % Initial guess
        
        %% Iterate until convergence
        % Newton-Raphson method used
        has_converged = 0; % set the convergence flag
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
                has_converged = 1;
                break;
            end
        end % End of Newton iteration
        disp('***Done***');
        % what about the non converged
        
        % Plot only if Newton method has converged
        if has_converged == 1;
            % extract the needed part
            J_longstab = J(stab_indices, stab_indices);
            [eigvectors, eigvalues] = eig(J_longstab);
            eigvalues = diag(eigvalues);
            %         eigvs = eig(Jtrim);
            
            % plot the results
            hold on; % just to make sure %todo - does it have any effect?
            h = plot(eigvalues, 'x', 'LineWidth', linewidth, ...
                'Color', colors_array(2*converged_i, :));
            h_string = sprintf('V = %.1f', V);
            
            h_list(converged_i) = h;
            h_strs{converged_i} = h_string;

            
            % print the results only for the highest and the lowest mach
            if converged_i == 1 || mach_i == length(mach_list)
                print_eigs_data(altm, machset, eigvectors, eigvalues, ...
                    fid);
            end
            
            converged_i = converged_i + 1;
        end
    end
    
    %% Finishing the figures off;
    
    % plot the legend - only for the ones that have converged
    h_strs_conv = h_strs(1:converged_i - 1);
    h_list_conv = h_list(1:converged_i - 1);
    h_legend = legend(h_list_conv, h_strs_conv);
    set(h_legend, 'FontSize', fontsize-2, 'Location', 'West');
    % axis, titles..
    hold on;
    xlabel('Real Axis', 'FontSize', fontsize);
    ylabel('Imaginary Axis', 'FontSize', fontsize);
    title(sprintf('Altitude = %d m', altm_list(altm_i)), 'FontSize', fontsize+2);
    grid on;
    
    % plo the vertical 0 line - investigate stability
    limsy=get(gca,'YLim');
    plot([0, 0], [limsy(1), limsy(2)], 'k', 'LineWidth', linewidth)
    % set the x-limits so that we can see the stability region clearlly
    limsx = get(gca, 'Xlim');
    set(gca, 'Xlim', [limsx(1)-0.7, limsx(2)]); % legend position problems
    if limsx(2) <= 0
        set(gca,'Xlim',[limsx(1), 0.1]);
    end
end
fclose(fid);

% save the figure
print(gcf, '-depsc2',[path_to_plots,'linear_stability']);
saveas(gcf, [path_to_plots, 'linear_stability'], 'fig');