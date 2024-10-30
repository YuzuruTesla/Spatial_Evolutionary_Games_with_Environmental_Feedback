% This script generates Figure 4 (a), (b), (c), and (d), which show heatmaps of u, v, f_L, and f_H. 
% Since the spatial patterns are derived from the q, p, n system, the heatmaps serve not only to verify 
% that the u, v, n system aligns with the q, p, n system under simulation, but also to depict the time of spatial pattern formation.
% 
% Key parameter: The sensitivity of environmental motion for L-type individuals, chi_u, is set slightly above the critical value, chi_u = 0.0647. 
% At this chi_u, only one wavenumber, m = 7, is present.

clc; clear; close all;

seed = 4;

rng(seed,"twister");

L = 10;
maxt = 100000;

% Parameters for pi_L and pi_H
R_0 = 0.35;
S_0 = 0.3;
T_0 = 0.5;
P_0 = 0.2;
R_1 = 0.6;
S_1 = 0.35;
T_1 = 0.6;
P_1 = 0.3;

x = linspace(0, L, 101); 
t = linspace(0, maxt, 100000);

m = 0; % Symmetry parameter

sol = pdepe(m, @pdeFunction, @initialConditions, @boundaryConditions, x, t);

% Extract solutions
u = sol(:,:,1);
v = sol(:,:,2);
n = sol(:,:,3);

% Preallocate arrays for pi_L and pi_H
pi_L = zeros(length(t), length(x));
pi_H = zeros(length(t), length(x));

% Calculate pi_L and pi_H at each time step
for i = 1:length(t)
    frac = u(i, :) ./ (u(i, :) + v(i, :));
    frac(isnan(frac)) = 0; % Handle division by zero

    pi_L(i, :) = (1 - n(i, :)) .* (R_0 .* frac + S_0 .* (1 - frac)) + n(i, :) .* (R_1 .* frac + S_1 .* (1 - frac));
    pi_H(i, :) = (1 - n(i, :)) .* (T_0 .* frac + P_0 .* (1 - frac)) + n(i, :) .* (T_1 .* frac + P_1 .* (1 - frac));
end

% Plot heat map for u (Figure 3(a))
figure;
imagesc(x, t, u);
set(gca, 'YDir', 'normal');
xlabel('Position $x$','fontsize',23,'fontname','arial', 'Interpreter', 'latex')
ylabel('Time $t$','fontsize',23,'fontname','arial','Interpreter', 'latex')
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
colorbar;

% Plot heat map for v (Figure 3(b))
figure;
imagesc(x, t, v);
set(gca, 'YDir', 'normal');
xlabel('Position $x$','fontsize',23,'fontname','arial', 'Interpreter', 'latex')
ylabel('Time $t$','fontsize',23,'fontname','arial','Interpreter', 'latex')
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
colorbar;

% Plot heat map for pi_L (Figure 3(c))
figure;
imagesc(x, t, pi_L);
set(gca, 'YDir', 'normal');
xlabel('Position $x$','fontsize',23,'fontname','arial', 'Interpreter', 'latex')
ylabel('Time $t$','fontsize',23,'fontname','arial','Interpreter', 'latex')
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
colorbar;

% Plot heat map for pi_H (Figure 3(d))
figure;
imagesc(x, t, pi_H);
set(gca, 'YDir', 'normal');
xlabel('Position $x$','fontsize',23,'fontname','arial', 'Interpreter', 'latex')
ylabel('Time $t$','fontsize',23,'fontname','arial','Interpreter', 'latex')
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
colorbar;


% PDE function definition
function [c, f, s] = pdeFunction(x, t, u, DuDx)
    % Parameters
    D_u = 0.01;
    D_n = 0.01;
    D_v = 0.01;
    chi_v = 0.02;
    chi_u = 0.0647;
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
