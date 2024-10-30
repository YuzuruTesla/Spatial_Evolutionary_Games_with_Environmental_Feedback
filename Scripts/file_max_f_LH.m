%Script used to generate Figure 7(a) based on file average_max_f_LH.csv

% Load the data from average_max_f_LH.csv
data = load('average_max_f_LH.csv'); 

chi_u_values = data(:, 1);  % First column: Sensitivity of Environmental Motion
avg_f_L = data(:, 2);       % Second column: Spatial Maximum of Payoff for L-type
avg_f_H = data(:, 3);       % Third column: Spatial Maximum of Payoff for H-type

% Plot for max f_L and f_H (Figure 7(a))
figure;
hold on;
dotSize = 80;
scatter(chi_u_values, avg_f_L, dotSize, 'r', 'filled', 'DisplayName', 'Maximal payoff $f_L$');
dotSize = 50;
scatter(chi_u_values, avg_f_H, dotSize, 'b', 'filled', 'DisplayName', 'Maximal payoff $f_H$');
h_line_payoffs = yline(0.4, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $f_L$ \& $f_H$', 'Interpreter', 'latex');
h_line_critical = xline(0.0646, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$', 'Interpreter', 'latex');
ylim([0.35, 0.44]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'Interpreter', 'latex', 'FontSize', 23, 'fontname', 'arial');
ylabel('Spatial Maximum of Payoff', 'Interpreter', 'latex', 'FontSize', 23, 'fontname', 'arial');
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 17, 'fontname', 'arial');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23; 
ax.YLabel.FontSize = 23;
grid on;
