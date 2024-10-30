% Script used to generate the Figure 1 showing the critical value of chi_u^*(m) 
% as a function of the perturbation mode m. The plot highlights the minimum 
% value of chi_u^*(m) occurring at m^* = 6.78 with \chi_u^* = 0.0646. 


% Define initial conditions and parameters for the simulation
p0 = 0.5;               % Initial condition of p
n0 = 0.5;               % Initial condition of n
L = 10;                 % Length of the spatial domain
Du = 0.01;              % Diffusion coefficient of u
Dn = 0.01;              % Diffusion coefficient of n
chi_v = 0.02;           % Sensitivity of environment-driven motion for H-type individuals
r = 1;
a = 0.5; 
epsilon = 10;
eL = 0.2;
eH = 0.5;

%Game parameter settings
R0 = 0.35;
S0 = 0.3; 
T0 = 0.5;
P0 = 0.2;
R1 = 0.6;
S1 = 0.35;
T1 = 0.6;
P1 = 0.3;


c32 = r - a * (eL * n0 + eH * (1 - n0));
kp = @(p, n) (1 - n) * (R0 - S0 - T0 + P0) + n * (R1 - S1 - T1 + P1); 
kn = @(p, n) p * (R1 - R0 - T1 + T0) + (1 - p) * (S1 - S0 - P1 + P0); 
c22 = p0 * (1 - p0) * kp(p0, n0) / epsilon; 
c23 = p0 * (1 - p0) * kn(p0, n0) / epsilon; 
c33 = -c32;

% Calculate parameters b, c, and d
b = Du * Dn / (c32 * p0 * (1 - p0));
c = (c22 * c33 - c32 * c23) / (c32 * p0 * (1 - p0));
d = chi_v - (c22 * Dn + Du * c33) / (c32 * p0 * (1 - p0));


% Define the function for critical sensitivity $\chi_u^*(m)$ as a function of wavenumber m
chi_u_star = @(m) b * (m * pi / L)^2 + c / (m * pi / L)^2 + d; 

% Define the range of wavenumber m
m_values = linspace(1, 60, 10000); % Generate 10,000 points between m = 1 and m = 60
chi_u_values = arrayfun(chi_u_star, m_values); % Evaluate the critical sensitivity for each m value

figure;
plot(m_values, chi_u_values, 'LineWidth', 5); 
xlabel('Wavenumber $m$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex'); 
ylabel('Critical Sensitivity $\chi_u^*(m)$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylim([0, 0.25]);
set(gca, 'FontSize', 16); 
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;


