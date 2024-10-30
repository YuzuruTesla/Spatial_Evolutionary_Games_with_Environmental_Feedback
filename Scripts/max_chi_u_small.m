% This script generates Figure 7(d), which depict the spatial maxima of payoffs in the small range, 
% and the environmental metric as the sensitivity of environmental motion, chi_u, increases from 0.064 to 0.067.
%
% To ensure the reliability of the results, we perform simulations using 50 distinct random seeds and calculate the average outcomes. 
% 
% The spatial maxima of the environmental metric for different values of chi_u are saved in the file average_max_f_small.csv.


clc; clear; close all;

L = 10;
maxt = 2000000;
x = linspace(0, L, 101); 
t = linspace(0, maxt, 2000000); 

% Parameters for payoffs
R_0 = 0.35;
S_0 = 0.3;
T_0 = 0.5;
P_0 = 0.2;
R_1 = 0.6;
S_1 = 0.35;
T_1 = 0.6;
P_1 = 0.3;

chi_u_values = 0.064:0.0001:0.067;

% Number of seeds
num_seeds = 50;

% Initialize arrays to store maximum values for different chi_u and seeds
max_f_L = zeros(num_seeds, length(chi_u_values));
max_f_H = zeros(num_seeds, length(chi_u_values));

% Loop over seeds
for seed = 1:num_seeds
    rng(seed, "twister"); 
    
    for i = 1:length(chi_u_values)
        chi_u = chi_u_values(i);

        m = 0;
        sol = pdepe(m, @(x, t, u, DuDx) pdeFunction(x, t, u, DuDx, chi_u), ...
                    @initialConditions, @boundaryConditions, x, t);
        
        u = sol(:,:,1);
        v = sol(:,:,2);
        n = sol(:,:,3);
        
        % Get the solution at t = maxt
        u_final = u(end, :);
        v_final = v(end, :);
        n_final = n(end, :);
        
        % Compute final payoffs f_L and f_H
        frac_final = u_final ./ (u_final + v_final);
        frac_final(isnan(frac_final)) = 0;
        
        f_L_final = (1 - n_final) .* (R_0 .* frac_final + S_0 .* (1 - frac_final)) + ...
                      n_final .* (R_1 .* frac_final + S_1 .* (1 - frac_final));
        f_H_final = (1 - n_final) .* (T_0 .* frac_final + P_0 .* (1 - frac_final)) + ...
                      n_final .* (T_1 .* frac_final + P_1 .* (1 - frac_final));

        % Store maximum values for this seed and chi_u
        max_f_L(seed, i) = max(f_L_final);
        max_f_H(seed, i) = max(f_H_final);
    end
end

% Compute the average across all seeds for each chi_u value
avg_f_L = mean(max_f_L, 1);
avg_f_H = mean(max_f_H, 1);

% Store them as files
csvwrite('average_max_f_small.csv', [chi_u_values' avg_f_L' avg_f_H']);

% Plot for max f_L and f_H (Figure 7(a))
figure;
hold on;
dotSize = 80;
scatter(chi_u_values, avg_f_L, dotSize, 'r', 'filled', 'DisplayName', 'Maximal payoff $f_L$');
dotSize = 50;
scatter(chi_u_values, avg_f_H, dotSize, 'b', 'filled', 'DisplayName', 'Maximal payoff $f_H$');
h_line_payoffs = yline(0.4, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $f_L$ \& $f_H$', 'Interpreter', 'latex');
h_line_critical = xline(0.0646, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$', 'Interpreter', 'latex');
% ylim([0.35, 0.44]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'Interpreter', 'latex', 'FontSize', 23, 'fontname', 'arial');
ylabel('Spatial Maximum of Payoff', 'Interpreter', 'latex', 'FontSize', 23, 'fontname', 'arial');
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 17, 'fontname', 'arial');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23; 
ax.YLabel.FontSize = 23;
grid on;



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

    % Payoff parameters
    R_0 = 0.35;
    S_0 = 0.3;
    T_0 = 0.5;
    P_0 = 0.2;
    R_1 = 0.6;
    S_1 = 0.35;
    T_1 = 0.6;
    P_1 = 0.3;

    % Variables
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

% Boundary conditions function
function [pl, ql, pr, qr] = boundaryConditions(xl, ul, xr, ur, t)
    pl = [0; 0; 0];
    ql = [1; 1; 1];
    pr = [0; 0; 0];
    qr = [1; 1; 1];
end

% Initial conditions function
function u0 = initialConditions(x)
    u0 = [2000 + 0.2 * randn - 0.1;
          2000 + 0.2 * randn - 0.1;
          0.5 + 0.002 * randn - 0.001];
end
