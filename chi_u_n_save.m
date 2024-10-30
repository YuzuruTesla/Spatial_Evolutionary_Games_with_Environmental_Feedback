% Script generates Figure 6(c), which illustrates the corresponding decrease in the environmental metric n 
% as chi_u increases from 0.05 to 0.1. 
% 
% To ensure the reliability of our results, we loop over 50 distinct random seeds and calculate the average outcome. 
% The environmental metrics corresponding to different values of chi_u are stored in a file named average_n.csv.
% For the average environmental metric calculation, both the right and left Riemann sums are computed and averaged for each seed.
% 
% Additionally, a horizontal line at n = 0.5 is included to represent uniform n, and a vertical line at chi_u = 0.0646 is added 
% to indicate the critical value, chi_u^*.


clc; clear; close all;

L = 10;
maxt = 2000000;

x = linspace(0, L, 101); 
t = linspace(0, maxt, 2000000); 

% Environmental motion sensitivity (chi_u) values
chi_u_values = 0.0500:0.001:0.1;

% Initialize the final average n values over seeds
average_n_final_over_seeds = zeros(size(chi_u_values));

% Number of seeds
num_seeds = 50;

% To store all n_final values for each seed and each chi_u
n_final_all_seeds = zeros(num_seeds, length(chi_u_values), length(x));

% Loop over seeds
for seed = 1:num_seeds
    
    rng(seed, "twister");
    
    % Loop over chi_u values
    for i = 1:length(chi_u_values)
        chi_u = chi_u_values(i);
        
        m = 0; % Symmetry parameter

        sol = pdepe(m, @(x,t,u,DuDx)pdeFunction(x, t, u, DuDx, chi_u), @initialConditions, @boundaryConditions, x, t);

        % Extract solutions
        u_1 = sol(:,:,1);
        v_1 = sol(:,:,2);
        n_1 = sol(:,:,3);

        % Get the solution at t = maxt
        n_final = n_1(end, :);

        % Store n_final for this seed and chi_u
        n_final_all_seeds(seed, i, :) = n_final;
        
        dx = x(2) - x(1);

        % Compute Left Riemann sum (excluding the last element)
        total_n_final_left = sum(n_final(1:end-1)) * dx;

        % Compute Right Riemann sum (excluding the first element)
        total_n_final_right = sum(n_final(2:end)) * dx;

        % Compute Average Riemann sum
        total_n_final_avg = (total_n_final_left + total_n_final_right) / 2;

        % Accumulate average n for this chi_u
        average_n_final_over_seeds(i) = average_n_final_over_seeds(i) + total_n_final_avg / L;
    end
end

% Compute the average n values across all seeds
average_n_final_over_seeds = average_n_final_over_seeds / num_seeds;

% The environmental metrics corresponding to different values of chi_u are stored in a file named average_n.csv.
csvwrite('average_n_im.csv', [chi_u_values' average_n_final_over_seeds']);

% Plot for n (Figure 6(c))
figure;
hold on;
dotSize = 80;
h_scatter_n = scatter(chi_u_values, average_n_final_over_seeds, 'filled', 'DisplayName', 'Environmental metric $n$');
h_line_factor = yline(0.5, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $n$');
h_line_critical = xline(0.0646, '--b','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$');
ylim([0.15, 0.65]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
ylabel('Environmental Metric $n$', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
legend([h_scatter_n, h_line_factor, h_line_critical], ...
    {'Environmental metric $n$', 'Uniform $n$', 'Critical $\chi_u^*$'}, ...
    'Location', 'best', 'fontsize', 17, 'fontname', 'arial','Interpreter', 'latex');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;
hold off;


% PDE function definition
function [c, f, s] = pdeFunction(x, t, u, DuDx, chi_u)
    % Parameters
    D_u = 0.01;
    D_n = 0.01;
    D_v = 0.01;
    chi_v = 0.02;
    r = 1;
    a = 0.5;
    e_L = 0.2;
    e_H = 0.5;
    epsilon = 10;
    kappa = 0.0001;
    
    % Game parameters for pi_L and pi_H
    R_0 = 0.35;
    S_0 = 0.3;
    T_0 = 0.5;
    P_0 = 0.2;
    R_1 = 0.6;
    S_1 = 0.35;
    T_1 = 0.6;
    P_1 = 0.3;
    
    % Solution variables
    u1 = u(1);
    u2 = u(2);
    u3 = u(3);
    
    % Spatial derivatives
    Du1Dx = DuDx(1);
    Du2Dx = DuDx(2);
    Du3Dx = DuDx(3);
    
    % Compute f_L and f_H
    if u1 + u2 == 0
        frac = 0;
    else
        frac = u1 / (u1 + u2);
    end
    
    f_L = (1 - u3) * (R_0 * frac + S_0 * (1 - frac)) + u3 * (R_1 * frac + S_1 * (1 - frac));
    f_H = (1 - u3) * (T_0 * frac + P_0 * (1 - frac)) + u3 * (T_1 * frac + P_1 * (1 - frac));
    
    % Define the PDEs
     c = [1; 1; 1];
    
    f = [D_u * Du1Dx - chi_u * u1 * Du3Dx;
         D_v * Du2Dx - chi_v * u2 * Du3Dx;
         D_n * Du3Dx];
    
    s = [u1 * (f_L - kappa * (u1 + u2)) / epsilon;
         u2 * (f_H - kappa * (u1 + u2)) / epsilon;
         (r - a * (e_L * u3 + e_H * (1 - u3))) * (u1 / (u1 + u2) - u3)];
end

% Boundary conditions function definition
function [pl, ql, pr, qr] = boundaryConditions(xl, ul, xr, ur, t)
    % Left boundary conditions (zero flux)
    pl = [0; 0; 0];
    ql = [1; 1; 1];
    
    % Right boundary conditions (zero flux)
    pr = [0; 0; 0];
    qr = [1; 1; 1];
end

% Initial conditions function definition
function u0 = initialConditions(x)
    u0 = [2000 + 0.2 * randn - 0.1; 2000 + 0.2 * randn - 0.1; 0.5 + 0.002 * randn - 0.001];
end