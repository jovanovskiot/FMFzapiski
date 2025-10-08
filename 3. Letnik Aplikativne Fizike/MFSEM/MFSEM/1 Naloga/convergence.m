clc; clear; close all;

n_terms = 15;  % Number of terms in the asymptotic series
x_values = [8];  % Different x-values to analyze

colors = ['r', 'g', 'b', 'm'];  % Different colors for each x
hold on;

for i = 1:length(x_values)
    x = x_values(i);
    terms = zeros(n_terms, 1);  % Store individual terms
    sum_series = zeros(n_terms, 1);  % Store partial sums

    % Compute the asymptotic series for erf(x)
    for n = 1:n_terms
        term = ((-1)^n) * (factorial(2*n) / (2^n * factorial(n) * (2*n + 1))) * x^(-2*n-1);
        terms(n) = term;
        sum_series(n) = sum(terms(1:n));  % Compute the partial sum
    end

    % Plot the convergence for this x
    plot(1:n_terms, sum_series, '-o', 'Color', colors(i), 'LineWidth', 1.5, 'DisplayName', ['x = ', num2str(x)]);
end

% Formatting the plot
xlabel('Število členov (n)');
ylabel('Približek erf(z)');
title('Konvergenca asimptotske vrste za erf(x)');
legend('show');
grid on;
hold off;

