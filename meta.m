% This script generates Figure 8(a), which depicts the oscillations in average payoffs over time, from t = 0 to t = 1200. 
% To ensure accuracy in the calculations, both left and right Riemann sums are computed, and their average is taken.
% 
% For the choice of environmental sensitivity, chi_u, we set chi_u = 0.1, which is significantly above the critical threshold.

clc; clear; close all;

seed = 7;
rng(seed, "twister");

L = 10;
maxt = 1.2e3;

x = linspace(0, L, 100); 
t = linspace(0, maxt, 1200);

m = 0; % Symmetry parameter

sol = pdepe(m, @pdeFunction, @initialConditions, @boundaryConditions, x, t);

% Extract solutions
u = sol(:,:,1);
v = sol(:,:,2);
n = sol(:,:,3);

% Parameters for payoff calculation
R_0 = 0.35; S_0 = 0.3; T_0 = 0.5; P_0 = 0.2;
R_1 = 0.6; S_1 = 0.35; T_1 = 0.6; P_1 = 0.3;

% Initialize arrays to store average payoffs over time
average_payoff_f_L = zeros(length(t), 1);
average_payoff_f_H = zeros(length(t), 1);

dx = x(2) - x(1);

% Helper function to compute Riemann sums
function [total_f_L, total_f_H, total_u, total_v] = computeRiemannSums(f_L, f_H, u, v, dx)
    total_f_L = sum(f_L .* u) * dx;
    total_f_H = sum(f_H .* v) * dx;
    total_u = sum(u) * dx;
    total_v = sum(v) * dx;
end

% Compute payoffs \pi_L and \pi_H at each time step
for i = 1:length(t)
    frac = u(i,:) ./ (u(i,:) + v(i,:));
    frac(isnan(frac)) = 0; % Handle division by zero
    
    f_L = (1 - n(i,:)) .* (R_0 .* frac + S_0 .* (1 - frac)) + n(i,:) .* (R_1 .* frac + S_1 .* (1 - frac));
    f_H = (1 - n(i,:)) .* (T_0 .* frac + P_0 .* (1 - frac)) + n(i,:) .* (T_1 .* frac + P_1 .* (1 - frac));
    
    % Calculate left and right Riemann sums
    [total_f_L_left, total_f_H_left, total_u_left, total_v_left] = computeRiemannSums(f_L(1:end-1), f_H(1:end-1), u(i,1:end-1), v(i,1:end-1), dx);
    [total_f_L_right, total_f_H_right, total_u_right, total_v_right] = computeRiemannSums(f_L(2:end), f_H(2:end), u(i,2:end), v(i,2:end), dx);

    % Compute average Riemann sums
    total_f_L_avg = (total_f_L_left + total_f_L_right) / 2;
    total_f_H_avg = (total_f_H_left + total_f_H_right) / 2;
    total_u_avg = (total_u_left + total_u_right) / 2;
    total_v_avg = (total_v_left + total_v_right) / 2;

    % Compute the average payoffs
    average_payoff_f_L(i) = total_f_L_avg / total_u_avg;
    average_payoff_f_H(i) = total_f_H_avg / total_v_avg;
end

figure;
plot(t, average_payoff_f_L, 'LineWidth', 4, 'DisplayName', 'Average payoff $f_L$');
hold on;
plot(t, average_payoff_f_H, 'LineWidth', 4, 'DisplayName', 'Average payoff $f_H$');
xlabel('Time t','fontsize',20,'fontname','arial', 'Interpreter', 'latex');
ylabel('Average Payoff','fontsize',20,'fontname','arial', 'Interpreter', 'latex');
legend('show', 'Location', 'best', 'fontsize', 17, 'fontname','arial', 'Interpreter', 'latex');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;

% PDE function definition
function [c, f, s] = pdeFunction(x, t, u, DuDx)
    % Parameters
    D_u = 0.01;
    D_n = 0.01;
    D_v = 0.01;
    chi_v = 0.02;
    chi_u = 0.1;
    r = 1;
    a = 0.5;
    e_L = 0.2;
    e_H = 0.5;
    epsilon = 10;
    kappa = 0.0001;
    
    % Game parameters for f_L and f_H
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
