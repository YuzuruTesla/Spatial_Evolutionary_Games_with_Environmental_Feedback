%Script used to generate Figure 6(a) based on file average_n.csv

% Load the data from your file
data = load('average_uv.csv'); 

chi_u = data(:, 1);  % First column: Sensitivity of Environmental Motion
u = data(:, 2);      % Second column: Average Density of u
v = data(:, 3);      % Third column: Average Density of v

% Plot for u and v (Figure 6(b))
figure;
hold on;
dotSize = 80;
h_scatter_u = scatter(chi_u, u, dotSize, "red", 'filled', 'DisplayName', 'L-types $u$');
dotSize = 50;
h_scatter_v = scatter(chi_u, v, dotSize, "blue",'filled', 'DisplayName', 'H-types $v$');
h_line_populations = yline(2000, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $u$ \& $v$');
h_line_critical = xline(0.0646, '--b','Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$');
ylim([500, 3000]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
ylabel('Density of Individuals', 'fontsize', 23, 'fontname', 'arial','Interpreter', 'latex');
legend([h_scatter_u, h_scatter_v, h_line_populations, h_line_critical], ...
    {'L-types $u$', 'H-types $v$', 'Uniform $u$ \& $v$', 'Critical $\chi_u^*$'}, ...
    'Location', 'best', 'fontsize', 17, 'fontname', 'arial','Interpreter', 'latex');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;
hold off;
