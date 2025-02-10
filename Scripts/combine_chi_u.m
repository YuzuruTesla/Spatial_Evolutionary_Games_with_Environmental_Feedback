% Script used to generate the Figure 9 which plots two thresholds strengths of 
% environmental-driven motion of as a function of the varying diffusion coefficient $D_v$.


% Define initial conditions and parameters for the simulation
p_0 = 0.5;          % Initial condition of p
n_0 = 0.5;          % Initial condition of n
q_0 = 4000;         % Initial condition of q
L = 10;             % Length of the spatial domain
D_u = 0.01;         % Diffusion coefficient of u
D_n = 0.01;         % Diffusion coefficient of n
chi_v = 0.02;       % Sensitivity of environment-driven motion
r = 1;
a = 0.5;
e_L = 0.2;
e_H = 0.5;
R_0 = 0.35;
S_0 = 0.3;
T_0 = 0.5;
P_0 = 0.2;
R_1 = 0.6;
S_1 = 0.35;
T_1 = 0.6;
P_1 = 0.3;
epsilon = 10;
kappa = 0.0001;

% Calculate payoffs
f_L = (1 - n_0) * (R_0 * p_0 + S_0 * (1 - p_0)) + n_0 * (R_1 * p_0 + S_1 * (1 - p_0)); 
f_H = (1 - n_0) * (T_0 * p_0 + P_0 * (1 - p_0)) + n_0 * (T_1 * p_0 + P_1 * (1 - p_0));

f_L_p = (1 - n_0) * (R_0 - S_0) + n_0 * (R_1 - S_1);
f_H_p = (1 - n_0) * (T_0 - P_0) + n_0 * (T_1 - P_1);
f_L_n = p_0 * (R_1 - R_0) + (1 - p_0) * (S_1 - S_0);
f_H_n = p_0 * (T_1 - T_0) + (1 - p_0) * (P_1 - P_0);
k_p = f_L_p - f_H_p; 
k_n = f_L_n - f_H_n;

% Calculate coefficients
c11 = ((1 - p_0) * f_H + p_0 * f_L - 2 * kappa * q_0) / epsilon;  
c12 = (q_0 * (f_H_p + p_0 * k_p)) / epsilon;
c13 = (q_0 * (f_H_n + p_0 * k_n)) / epsilon;
c22 = (p_0 * (1 - p_0) * k_p) / epsilon;
c23 = (p_0 * (1 - p_0) * k_n) / epsilon;
c32 = r - a * (e_L * n_0 + e_H * (1 - n_0));
c33 = -(r - a * (e_L * n_0 + e_H * (1 - n_0)));


ChiU1 = @(m, D_v) ...
    (((c32 * chi_v * (D_u - D_v) * (m * pi / L)^4 * (1 - p_0)^2 * p_0) + ...
     c23 * c32 * (c11 - (m * pi / L)^2 * (D_v * (1 - p_0) + D_u * p_0))) + ...
    (-c32 * chi_v * (m * pi / L)^2 * (1 - p_0) * p_0 * (c11 - (m * pi / L)^2 * (D_v * (1 - p_0) + D_u * p_0))) + ...
    (-c13 * c32 * (D_v - D_u) * (m * pi / L)^2 * (1 - p_0) * p_0 / q_0) + ...
    (-(c33 - (m * pi / L)^2 * D_n) * (c11 - (m * pi / L)^2 * (D_v * (1 - p_0) + D_u * p_0)) * ...
     (c22 - (m * pi / L)^2 * (D_u * (1 - p_0) + D_v * p_0))) + ...
    ((c33 - (m * pi / L)^2 * D_n) * ((D_v - D_u) * (m * pi / L)^2 * (1 - p_0) * p_0 * ...
     (c12 - (m * pi / L)^2 * (D_u - D_v) * q_0)) / q_0)) / ...
    (-c32 * (m * pi / L)^4 * (D_u - D_v) * (1 - p_0) * p_0^2 - ...
     c32 * (m * pi / L)^2 * (1 - p_0) * p_0 * (c11 - (m * pi / L)^2 * (D_v * (1 - p_0) + D_u * p_0)));


ChiU2 = @(m, D_v) ...
    (-2 * c13 * c32 * (D_u - D_v) * L^4 * m^2 * (-1 + p_0) * p_0 * pi^2 + ...
     2 * c32 * chi_v * (D_u - D_v) * L^2 * m^4 * (-1 + p_0)^2 * p_0 * pi^4 * q_0 + ...
     2 * c23 * c32 * L^4 * (c11 * L^2 + c22 * L^2 + c33 * L^2 - D_n * m^2 * pi^2 - ...
        D_u * m^2 * pi^2 - D_v * m^2 * pi^2) * q_0 + ...
     2 * c32 * chi_v * L^2 * m^2 * (-1 + p_0) * p_0 * pi^2 * (c11 * L^2 + c22 * L^2 + ...
        c33 * L^2 - D_n * m^2 * pi^2 - D_u * m^2 * pi^2 - D_v * m^2 * pi^2) * q_0 + ...
     (c33 * L^2 - D_n * m^2 * pi^2)^2 * (c11 * L^2 + c22 * L^2 + c33 * L^2 - ...
        D_n * m^2 * pi^2 - D_u * m^2 * pi^2 - D_v * m^2 * pi^2) * q_0 - ...
     (c11 * L^2 + c22 * L^2 + c33 * L^2 - D_n * m^2 * pi^2 - D_u * m^2 * pi^2 - ...
        D_v * m^2 * pi^2)^3 * q_0 + ...
     2 * c23 * c32 * L^4 * (c11 * L^2 + m^2 * (D_v * (-1 + p_0) - D_u * p_0) * pi^2) * q_0 - ...
     2 * c32 * chi_v * L^2 * m^2 * (1 - p_0) * p_0 * pi^2 * (c11 * L^2 + ...
        m^2 * (D_v * (-1 + p_0) - D_u * p_0) * pi^2) * q_0 + ...
     (c11 * L^2 + c22 * L^2 + c33 * L^2 - D_n * m^2 * pi^2 - D_u * m^2 * pi^2 - ...
        D_v * m^2 * pi^2) * (c11 * L^2 + m^2 * (D_v * (-1 + p_0) - D_u * p_0) * pi^2)^2 * q_0 + ...
     (c11 * L^2 + c22 * L^2 + c33 * L^2 - D_n * m^2 * pi^2 - D_u * m^2 * pi^2 - ...
        D_v * m^2 * pi^2) * (c22 * L^2 + m^2 * (D_u * (-1 + p_0) - D_v * p_0) * pi^2)^2 * q_0 - ...
     2 * (D_u - D_v) * m^2 * (1 - p_0) * p_0 * pi^2 * (c11 * L^2 + c22 * L^2 + ...
        c33 * L^2 - D_n * m^2 * pi^2 - D_u * m^2 * pi^2 - D_v * m^2 * pi^2) * ...
        (c12 * L^2 + (-D_u + D_v) * m^2 * pi^2 * q_0) + ...
     2 * (c33 * L^2 - D_n * m^2 * pi^2) * (c12 * (D_u - D_v) * L^2 * m^2 * (-1 + ...
        p_0) * p_0 * pi^2 + (-c11 * c22 * L^4 + c22 * L^2 * m^2 * (D_v + D_u * p_0 - D_v * p_0) * pi^2 + ...
        c11 * L^2 * m^2 * (D_u - D_u * p_0 + D_v * p_0) * pi^2 - D_u * D_v * m^4 * pi^4) * q_0)) / ...
     (2 * c32 * L^2 * m^2 * (-1 + p_0) * p_0 * pi^2 * (2 * c11 * L^2 + c22 * L^2 + ...
        c33 * L^2 - D_n * m^2 * pi^2 - D_u * m^2 * pi^2 - 2 * D_v * m^2 * pi^2) * q_0);

% Define range for Dv and m values
DvValues = 0:0.00001:0.04;
mValues = 0:1:100;

minChiUStar2 = zeros(size(DvValues));
for i = 1:length(DvValues)
    D_v = DvValues(i);
    chiUValues2 = arrayfun(@(m) ChiU2(m, D_v), mValues);
    minChiUStar2(i) = min(real(chiUValues2));
end

minChiUStar1 = zeros(size(DvValues)); 
for i = 1:length(DvValues)
    D_v = DvValues(i);
    chiUValues1 = arrayfun(@(m) ChiU1(m, D_v), mValues);
    minChiUStar1(i) = min(real(chiUValues1));
end 


figure;
plot(DvValues, minChiUStar1, 'b', 'LineWidth', 5, 'DisplayName', '$\chi^{\textnormal{II}}_{u}\left(D_v \right)$'); 
hold on;
plot(DvValues, minChiUStar2, 'r', 'LineWidth', 5, 'DisplayName', '$\chi^{\textnormal{III}}_{u}\left(D_v \right)$');
xlabel('Diffusivity of H-type Individuals $D_v$', 'fontsize', 20, 'Interpreter', 'latex');
ylabel('Minimum Sensitivity $\chi_u^*(D_v)$', 'fontsize', 20, 'Interpreter', 'latex');
ylim([0, 0.25]);
legend('show', 'Location', 'best', 'fontsize', 17, 'Interpreter', 'latex', 'fontname', 'arial');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;
