%Script used to generate Figure 7(b) based on file average_max_uv.csv

% Load the data from average_max_uv.csv
data = load('average_max_uv.csv'); 

chi_u_values = data(:, 1);  % First column: Sensitivity of Environmental Motion
avg_u = data(:, 2);         % Second column: Spatial Max Density of L-type Individuals
avg_v = data(:, 3);         % Third column: Spatial Max Density of H-type Individuals

% Plot for max u and v (Figure 7(b))
figure;
hold on;
dotSize = 80;
scatter(chi_u_values, avg_u, dotSize, 'r', 'filled', 'DisplayName', 'Maximal L-type $u$');
dotSize = 50;
scatter(chi_u_values, avg_v, dotSize, 'b', 'filled', 'DisplayName', 'Maximal H-type $v$');
h_line_populations = yline(2000, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $u$ \& $v$');
h_line_critical = xline(0.0646, '--b','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$');
ylim([1900, 3400]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'Interpreter', 'latex', 'FontSize', 23, 'fontname', 'arial');
ylabel('Spatial Max Density of Individuals', 'Interpreter', 'latex', 'FontSize', 23, 'fontname', 'arial');
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 17, 'fontname', 'arial');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;