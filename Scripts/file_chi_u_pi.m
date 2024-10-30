%Script used to generate Figure 6(a) based on file average_f_LH.csv
data = load('average_f_LH.csv'); 

% Assuming the data is stored as a matrix with three columns
chi_u = data(:, 1); % First column: Sensitivity of Environmental Motion
f_L = data(:, 2);   % Second column: Average payoff f_L
f_H = data(:, 3);   % Third column: Average payoff f_H

% Plot for n (Figure 6(a))
figure;
hold on;
dotSize = 80;
h_scatter_L = scatter(chi_u, f_L, dotSize, "red", 'filled', 'DisplayName', 'Average payoff $f_L$');
dotSize = 50;
h_scatter_H = scatter(chi_u, f_H, dotSize,"blue", 'filled', 'DisplayName', 'Average payoff $f_H$');
h_line_payoffs = yline(0.4, '--k', 'LineWidth', 3, 'DisplayName', 'Uniform $f_L$ \& $f_H$', 'Interpreter', 'latex');
h_line_critical = xline(0.0646, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'DisplayName', 'Critical $\chi_u^*$', 'Interpreter', 'latex');
ylim([0.25, 0.5]);
xlabel('Sensitivity of Environmental Motion $\chi_u$', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
ylabel('Average Payoff', 'fontsize', 23, 'fontname', 'arial', 'Interpreter', 'latex');
legend([h_scatter_L, h_scatter_H, h_line_payoffs, h_line_critical], ...
    {'Average payoff $f_L$', 'Average payoff $f_H$', 'Uniform $f_L$ \& $f_H$', 'Critical $\chi_u^*$'}, ...
    'Location', 'best', 'fontsize', 17, 'fontname', 'arial', 'Interpreter', 'latex');
set(gca, 'FontSize', 16);
ax = gca;
ax.XLabel.FontSize = 23;
ax.YLabel.FontSize = 23;
grid on;
hold off;
