% Script generates Figure 6(a), which illustrates the corresponding decrease in the payoff f_L and f_H 
% as chi_u increases from 0.05 to 0.1. 
% 
% To ensure the reliability of our results, we loop over 50 distinct random seeds and calculate the average outcome. 
% The average payoff corresponding to different values of chi_u are stored in a file named average_f_LH.csv.
% For the average payoffs calculation, both the right and left Riemann sums are computed and averaged for each seed.
% 
% Additionally, a horizontal line at f_L = f_H = 0.5 is included to represent uniform payoff, and a vertical line at chi_u = 0.0646 is added 
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

% Initialize arrays to store final f_L and _H values for each chi_u
average_f_L_over_seeds = zeros(size(chi_u_values));
average_f_H_over_seeds = zeros(size(chi_u_values));

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
        u_final = u_1(end, :);
        v_final = v_1(end, :);
        n_final = n_1(end, :);

        % Compute payoffs f_L and f_H at t = maxt
        frac_final = u_final ./ (u_final + v_final);
        frac_final(isnan(frac_final)) = 0; % Handle division by zero

        f_L_final = (1 - n_final) .* (R_0 .* frac_final + S_0 .* (1 - frac_final)) + n_final .* (R_1 .* frac_final + S_1 .* (1 - frac_final));
        f_H_final = (1 - n_final) .* (T_0 .* frac_final + P_0 .* (1 - frac_final)) + n_final .* (T_1 .* frac_final + P_1 .* (1 - frac_final));

        dx = x(2) - x(1);

        % Compute Left Riemann sum (excluding the last element)
        total_f_L_left = sum(f_L_final(1:end-1) .* u_final(1:end-1)) * dx;
        total_f_H_left = sum(f_H_final(1:end-1) .* v_final(1:end-1)) * dx;
        total_u_left = sum(u_final(1:end-1)) * dx;
        total_v_left = sum(v_final(1:end-1)) * dx;

        % Compute Right Riemann sum (excluding the first element)
        total_f_L_right = sum(f_L_final(2:end) .* u_final(2:end)) * dx;
        total_f_H_right = sum(f_H_final(2:end) .* v_final(2:end)) * dx;
        total_u_right = sum(u_final(2:end)) * dx;
        total_v_right = sum(v_final(2:end)) * dx;

        % Compute Average Riemann sums
        total_f_L_avg = (total_f_L_left + total_f_L_right) / 2;
        total_f_H_avg = (total_f_H_left + total_f_H_right) / 2;
        total_u_avg = (total_u_left + total_u_right) / 2;
        total_v_avg = (total_v_left + total_v_right) / 2;

        % Compute the average payoffs
        average_f_L_final = total_f_L_avg / total_u_avg;
        average_f_H_final = total_f_H_avg / total_v_avg;

        % Accumulate average values over seeds
        average_f_L_over_seeds(i) = average_f_L_over_seeds(i) + average_f_L_final;
        average_f_H_over_seeds(i) = average_f_H_over_seeds(i) + average_f_H_final;
    end
end

% Compute the average payoffs across all seeds
average_f_L_over_seeds = average_f_L_over_seeds / num_seeds;
average_f_H_over_seeds = average_f_H_over_seeds / num_seeds;

% The average payoffs corresponding to different values of chi_u are stored in a file named average_f_LH.csv.
csvwrite('average_f_LH.csv', [chi_u_values' average_f_L_over_seeds' average_f_H_over_seeds']);

% Plot for n (Figure 6(a))
figure;
hold on;
dotSize = 80;
h_scatter_L = scatter(chi_u_values, average_f_L_over_seeds, dotSize, "red", 'filled', 'DisplayName', 'Average payoff $f_L$');
dotSize = 50;
h_scatter_H = scatter(chi_u_values, average_f_H_over_seeds, dotSize,"blue", 'filled', 'DisplayName', 'Average payoff $f_H$');
h_line_payoffs = yline(0.4, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $f_L$ \& $f_H$', 'Interpreter', 'latex');
h_line_critical = xline(0.0646, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$', 'Interpreter', 'latex');
ylim([0.25, 0.5]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Average Payoff', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
legend([h_scatter_L, h_scatter_H, h_line_payoffs, h_line_critical], ...
    {'Average payoff $f_L$', 'Average payoff $f_H$', 'Uniform $f_L$ \& $f_H$', 'Critical $\chi_u^*$'}, ...
    'Location', 'best', 'fontsize', 17, 'fontname', 'arial', 'Interpreter', 'latex');
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
