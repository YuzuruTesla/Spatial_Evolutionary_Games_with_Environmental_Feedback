% This script generates Figure 5 (a), (b), and (c), illustrating the spatial patterns of populations u & v, 
% payoffs f_L & f_H, and the environmental metric n. Additionally, it computes the average Riemann sum at t = 10^6. 
% To improve accuracy, both left and right Riemann sums are calculated, and their average is taken.
%
% The sensitivity of environmental motion for L-type individuals, chi_u, is set significantly above the critical value 
% chi_u = 0.1. At this level of chi_u, multiple wavenumbers are present.


clc; clear; close all;

seed = 4;

rng(seed,"twister");

L = 10;
maxt = 1000000;

x = linspace(0, L, 101); 
t = linspace(0, maxt, 1000000); 

m = 0; % Symmetry parameter

sol = pdepe(m, @pdeFunction, @initialConditions, @boundaryConditions, x, t);

% Extract solutions
u = sol(:,:,1);
v = sol(:,:,2);
n = sol(:,:,3);

% Get the solution at t = maxt
u_final = u(end, :);
v_final = v(end, :);
n_final = n(end, :);

% Get the solution at t = 0
u_initial = u(1, :);
v_initial = v(1, :);
n_initial = n(1, :);

%Parameters
R_0 = 0.35;
S_0 = 0.3;
T_0 = 0.5;
P_0 = 0.2;
R_1 = 0.6;
S_1 = 0.35;
T_1 = 0.6;
P_1 = 0.3;


% Compute payoffs f_L and f_H at t = maxt
frac_final = u_final ./ (u_final + v_final);
frac_final(isnan(frac_final)) = 0; % Handle division by zero

f_L_final = (1 - n_final) .* (R_0 .* frac_final + S_0 .* (1 - frac_final)) + n_final .* (R_1 .* frac_final + S_1 .* (1 - frac_final));
f_H_final = (1 - n_final) .* (T_0 .* frac_final + P_0 .* (1 - frac_final)) + n_final .* (T_1 .* frac_final + P_1 .* (1 - frac_final));

% Compute the left and right Riemann sums of payoffs at t=maxt
dx = x(2) - x(1);

% Left Riemann sums (excluding the last element)
total_payoff_f_L_final_left = sum(f_L_final(1:end-1) .* u_final(1:end-1)) * dx;
total_payoff_f_H_final_left = sum(f_H_final(1:end-1) .* v_final(1:end-1)) * dx;
total_u_final_left = sum(u_final(1:end-1)) * dx;
total_v_final_left = sum(v_final(1:end-1)) * dx;

% Right Riemann sums (excluding the first element)
total_payoff_f_L_final_right = sum(f_L_final(2:end) .* u_final(2:end)) * dx;
total_payoff_f_H_final_right = sum(f_H_final(2:end) .* v_final(2:end)) * dx;
total_u_final_right = sum(u_final(2:end)) * dx;
total_v_final_right = sum(v_final(2:end)) * dx;

% Average Riemann sum
total_payoff_f_L_final_avg = (total_payoff_f_L_final_left + total_payoff_f_L_final_right) / 2;
total_payoff_f_H_final_avg = (total_payoff_f_H_final_left + total_payoff_f_H_final_right) / 2;
total_u_final_avg = (total_u_final_left + total_u_final_right) / 2;
total_v_final_avg = (total_v_final_left + total_v_final_right) / 2;

average_payoff_f_L_final = total_payoff_f_L_final_avg / total_u_final_avg;
average_payoff_f_H_final = total_payoff_f_H_final_avg / total_v_final_avg;

fprintf('Average Riemann sum of f_L at t=maxt: %f\n', average_payoff_f_L_final);
fprintf('Average Riemann sum of f_H at t=maxt: %f\n', average_payoff_f_H_final);


% Compute the left and right Riemann sums of u, v, and n at t=maxt
% Left Riemann sums (excluding the last element)
total_u_final_left = sum(u_final(1:end-1)) * dx;
total_v_final_left = sum(v_final(1:end-1)) * dx;
total_n_final_left = sum(n_final(1:end-1)) * dx;

% Right Riemann sums (excluding the first element)
total_u_final_right = sum(u_final(2:end)) * dx;
total_v_final_right = sum(v_final(2:end)) * dx;
total_n_final_right = sum(n_final(2:end)) * dx;

% Average Riemann sums
total_u_final_avg = (total_u_final_left + total_u_final_right) / 2;
total_v_final_avg = (total_v_final_left + total_v_final_right) / 2;
total_n_final_avg = (total_n_final_left + total_n_final_right) / 2;

average_u_final = total_u_final_avg / L;
average_v_final = total_v_final_avg / L;
average_n_final = total_n_final_avg / L;

fprintf('Average Riemann sum of u at t=maxt: %f\n', average_u_final);
fprintf('Average Riemann sum of v at t=maxt: %f\n', average_v_final);
fprintf('Average Riemann sum of n at t=maxt: %f\n', average_n_final);


% Plot f_L and f_H 
figure;
plot(x, f_L_final, 'r', 'LineWidth', 4, 'DisplayName', 'Patterned $f_L$');
hold on;
plot(x, f_H_final, 'b', 'LineWidth', 4, 'DisplayName', 'Patterned $f_H$');
yline(0.4, 'k--', 'LineWidth', 4, 'DisplayName', 'Uniform $f_L$ \& $f_H$');
xlabel('Position $x$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Payoff', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex'); 
ylim([0.22, 0.5]);
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;


% Plot n
figure;
plot(x, n_final, 'LineWidth', 4, 'DisplayName', 'Patterned $n$');
hold on;
yline(0.5, '--k', 'LineWidth', 4, 'DisplayName', 'Uniform $n$');
xlabel('Position $x$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Environmental Metric $n$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylim([0.05, 0.65]); 
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;


% Plot u and v
figure;
plot(x, u_final, 'r-', 'LineWidth', 4, 'DisplayName', 'Patterned $u$');
hold on;
plot(x, v_final, 'b-', 'LineWidth', 4, 'DisplayName', 'Patterned $v$');
yline(2000, 'k--', 'LineWidth', 4, 'DisplayName', 'Uniform $u$ \& $v$');
xlabel('Position $x$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Density of Individuals', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylim([0, 4000]); 
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
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
    u0 = [2000 + 0.2 * randn - 0.1; 2000 + 0.2 * randn - 0.1; 0.5 + 0.002 * randn - 0.001]; % Initial conditions for u, v, and n
end
