% Script used to generate Figure 2 which plots the critical value of $\chi_u^*(m)$ as a function of the varying 
% diffusion coefficient $D_v$.

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

% Define range of D_v values (varying diffusion coefficients)
D_v_values = 0:0.00001:0.04;
m_values = 0:100;

% Initialize matrix to store minimum chi_u* for each D_v
min_chi_u_star = zeros(length(D_v_values), 1);

% Loop over each D_v value
for i = 1:length(D_v_values)
    D_v = D_v_values(i);
    chi_u_star_values = zeros(length(m_values), 1);  % Store chi_u* for each m
    
    % Loop over each m value
    for j = 1:length(m_values)
        m = m_values(j);
        
        term1 = (c32 * chi_v * (D_u - D_v) * (m * pi / L)^4 * (1 - p_0)^2 * p_0) + ...
                c23 * c32 * (c11 - (m * pi / L)^2 * (D_v * (1 - p_0) + D_u * p_0));
        term2 = - c32 * chi_v * (m * pi / L)^2 * (1 - p_0) * p_0 * (c11 - (m * pi / L)^2 * (D_v * (1 - p_0) + D_u * p_0));
        term3 = - (c13 * c32 * (D_v - D_u) * (m * pi / L)^2 * (1 - p_0) * p_0) / q_0;
        term4 = - (c33 - (m * pi / L)^2 * D_n) * (c11 - (m * pi / L)^2 * (D_v * (1 - p_0) + D_u * p_0)) * ...
                 (c22 - (m * pi / L)^2 * (D_u * (1 - p_0) + D_v * p_0));
        term5 = (c33 - (m * pi / L)^2 * D_n) * ((D_v - D_u) * (m * pi / L)^2 * (1 - p_0) * p_0 * (c12 - (m * pi / L)^2 * (D_u - D_v) * q_0)) / q_0;
        
        % Numerator of chi_u_star
        chi_u_num = term1 + term2 + term3 + term4 +term5;
        
        % Denominator of chi_u_star
        chi_u_den = -c32 * (m * pi / L)^4 * (D_u - D_v) * (1 - p_0) * p_0^2 - ...
                    c32 * (m * pi / L)^2 * (1 - p_0) * p_0 * (c11 - (m * pi / L)^2 * (D_v * (1 - p_0) + D_u * p_0));
        
        chi_u_star = chi_u_num / chi_u_den;
        
        chi_u_star_values(j) = chi_u_star;  % Store the calculated value
    end
    
    % Find the minimum chi_u_star for the current D_v
    min_chi_u_star(i) = min(chi_u_star_values);
   
end

% Figure 2(a)
figure;
plot(D_v_values, min_chi_u_star, 'LineWidth', 5); 
xlabel('Diffusivity of H-type Individuals $D_v$', 'fontsize', 20, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Minimum Sensitivity $\chi_u^*(D_v)$', 'fontsize', 20, 'fontname', 'arial', 'Interpreter', 'latex');
ylim([0, 0.25]);
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;


% Figure 2(b)
figure;
plot(D_v_values, min_chi_u_star, 'LineWidth', 5); 
xlabel('Diffusivity of H-type Individuals $D_v$', 'fontsize', 20, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Minimum Sensitivity $\chi_u^*(D_v)$', 'fontsize', 20, 'fontname', 'arial', 'Interpreter', 'latex');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;
