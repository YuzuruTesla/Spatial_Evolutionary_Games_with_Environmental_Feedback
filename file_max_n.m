%Script used to generate Figure 7(c) based on file average_max_n.csv

% Load the data from average_max_n.csv
data = load('average_max_n.csv'); 

chi_u_values = data(:, 1);  % First column: Sensitivity of Environmental Motion
avg_n = data(:, 2);         % Second column: Spatial Max of Environmental Metric


% Plot for max n (Figure 7(c))
figure;
hold on;
dotSize = 80;
scatter(chi_u_values, avg_n, dotSize, 'filled', 'DisplayName', 'Maximal $n$');
h_line_factor = yline(0.5, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $n$');
h_line_critical = xline(0.0646, '--b','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$');
ylim([0.34, 0.6]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Spatial Max of Environmental Metric', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
legend('show', 'Location', 'best', 'fontsize', 17, 'fontname', 'arial', 'Interpreter', 'latex');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;  
ax.YLabel.FontSize = 23;  
grid on;