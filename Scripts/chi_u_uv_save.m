% Script generates Figure 6(a), which illustrates the corresponding decrease in density of L-type u and increase 
% in density of H-type v as chi_u increases from 0.05 to 0.1. 
% 
% To ensure the reliability of our results, we loop over 50 distinct random seeds and calculate the average outcome. 
% The density of individuals corresponding to different values of chi_u are stored in a file named average_uv.csv.
% For the average density of individuals calculation, both the right and left Riemann sums are computed and averaged for each seed.
% 
% Additionally, a horizontal line at u = v = 2000 is included to represent uniform n, and a vertical line at chi_u = 0.0646 is added 
% to indicate the critical value, chi_u^*.


clc; clear; close all;

L = 10;
maxt = 2000000;

x = linspace(0, L, 101);
t = linspace(0, maxt, 2000000);

% Parameters
R_0 = 0.35;
S_0 = 0.3;
T_0 = 0.5;
P_0 = 0.2;
R_1 = 0.6;
S_1 = 0.35;
T_1 = 0.6;
P_1 = 0.3;

chi_u_values = 0.0500:0.001:0.1;

% Number of seeds
num_seeds = 50;

% Initialize arrays to store final u and v values for each seed and chi_u
average_u_final_over_seeds = zeros(size(chi_u_values));
average_v_final_over_seeds = zeros(size(chi_u_values));

% To store all u_final and v_final values for each seed and chi_u
u_final_all_seeds = zeros(num_seeds, length(chi_u_values), length(x));
v_final_all_seeds = zeros(num_seeds, length(chi_u_values), length(x));

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

        % Get the solution at t = maxt
        u_final = u_1(end, :);
        v_final = v_1(end, :);
        
        dx = x(2) - x(1);

        % Compute Left Riemann sums (excluding the last element)
        total_u_final_left = sum(u_final(1:end-1)) * dx;
        total_v_final_left = sum(v_final(1:end-1)) * dx;

        % Compute Right Riemann sums (excluding the first element)
        total_u_final_right = sum(u_final(2:end)) * dx;
        total_v_final_right = sum(v_final(2:end)) * dx;

        % Compute Average Riemann sums
        total_u_final_avg = (total_u_final_left + total_u_final_right) / 2;
        total_v_final_avg = (total_v_final_left + total_v_final_right) / 2;

        % Accumulate average values over seeds
        average_u_final_over_seeds(i) = average_u_final_over_seeds(i) + total_u_final_avg / L;
        average_v_final_over_seeds(i) = average_v_final_over_seeds(i) + total_v_final_avg / L;
    end
end

% Compute the average u and v values across all seeds
average_u_final_over_seeds = average_u_final_over_seeds / num_seeds;
average_v_final_over_seeds = average_v_final_over_seeds / num_seeds;


csvwrite('average_uv.csv', [chi_u_values' average_u_final_over_seeds' average_v_final_over_seeds']);

% Plot for u and v (Figure 6(a))
figure;
hold on;
dotSize = 80;
h_scatter_u = scatter(chi_u_values, average_u_final_over_seeds, dotSize, "red", 'filled', 'DisplayName', 'L-types $u$');
dotSize = 50;
h_scatter_v = scatter(chi_u_values, average_v_final_over_seeds, dotSize, "blue",'filled', 'DisplayName', 'H-types $v$');
h_line_populations = yline(2000, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $u$ \& $v$');
h_line_critical = xline(0.0646, '--b','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$');
ylim([500, 3000]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
ylabel('Density of Individuals', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
legend([h_scatter_u, h_scatter_v, h_line_populations, h_line_critical], ...
    {'L-types $u$', 'H-types $v$', 'Uniform $u$ \& $v$', 'Critical $\chi_u^*$'}, ...
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
    
    % Compute pi_L and pi_H
    if u1 + u2 == 0
        frac = 0;
    else
        frac = u1 / (u1 + u2);
    end
    
    pi_L = (1 - u3) * (R_0 * frac + S_0 * (1 - frac)) + u3 * (R_1 * frac + S_1 * (1 - frac));
    pi_H = (1 - u3) * (T_0 * frac + P_0 * (1 - frac)) + u3 * (T_1 * frac + P_1 * (1 - frac));
    
    % Define the PDEs
     c = [1; 1; 1];
    
    f = [D_u * Du1Dx - chi_u * u1 * Du3Dx;
         D_v * Du2Dx - chi_v * u2 * Du3Dx;
         D_n * Du3Dx];
    
    s = [u1 * (pi_L - kappa * (u1 + u2)) / epsilon;
         u2 * (pi_H - kappa * (u1 + u2)) / epsilon;
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
