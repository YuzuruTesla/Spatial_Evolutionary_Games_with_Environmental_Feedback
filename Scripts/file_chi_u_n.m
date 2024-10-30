%Script used to generate Figure 6(c) based on file average_n.csv

% Load the data from average_n.csv
data = load('average_n.csv'); 

chi_u = data(:, 1);  % First column: Sensitivity of Environmental Motion
n = data(:, 2);      % Second column: Environmental Metric

% Plot for n (Figure 6(c))
figure;
hold on;
dotSize = 80;
h_scatter_n = scatter(chi_u, n, dotSize, 'filled', 'DisplayName', 'Environmental metric $n$');
h_line_factor = yline(0.5, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $n$');
h_line_critical = xline(0.0646, '--b','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$');
ylim([0.15, 0.65]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
ylabel('Environmental Metric $n$', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
legend([h_scatter_n, h_line_factor, h_line_critical], ...
    {'Environmental metric $n$', 'Uniform $n$', 'Critical $\chi_u^*$'}, ...
    'Location', 'best', 'fontsize', 17, 'fontname', 'arial','Interpreter', 'latex');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;
hold off;
