% Parameters
intervals = {[-1, 1], [-2, 2], [-3, 3]};
mu_values = [0, 0.5, 1, 1.5, 2]; % Importance sampling means
N_values = round(logspace(4, 6, 5)); % MC sample sizes: 1e4 to 1e6
n_trap_values = round(logspace(2, 4, 5)); % Trapezoidal points: 100 to 10000

% Define manual PDF and CDF (to avoid needing the Statistics Toolbox)
normpdf_manual = @(x) (1/sqrt(2*pi)) * exp(-x.^2/2); % PDF for N(0,1)
normcdf_manual = @(x) 0.5 * (1 + erf(x / sqrt(2))); % CDF for N(0,1)

% Preallocate results
exact_vals = zeros(size(intervals));
errors_trap = cell(size(intervals));
times_trap = cell(size(intervals));
errors_mc = cell(size(intervals));
times_mc = cell(size(intervals));
errors_is = cell(size(intervals));
times_is = cell(size(intervals));

for k = 1:length(intervals)
    a = intervals{k}(1);
    b = intervals{k}(2);
    exact_vals(k) = normcdf_manual(b) - normcdf_manual(a); % Exact integral
    
    % Trapezoidal method
    errors_trap{k} = zeros(size(n_trap_values));
    times_trap{k} = zeros(size(n_trap_values));
    f = @(x) normpdf_manual(x); % Use manual PDF
    for i = 1:length(n_trap_values)
        n = n_trap_values(i);
        x = linspace(a, b, n);
        y = f(x);
        h = (b - a)/(n - 1);
        tic;
        I = h/2 * (y(1) + 2*sum(y(2:end-1)) + y(end));
        times_trap{k}(i) = toc;
        errors_trap{k}(i) = abs(I - exact_vals(k));
    end
    
    % Standard Monte Carlo
    errors_mc{k} = zeros(size(N_values));
    times_mc{k} = zeros(size(N_values));
    for i = 1:length(N_values)
        N = N_values(i);
        tic;
        x_uniform = a + (b - a)*rand(N, 1);
        I = (b - a) * mean(normpdf_manual(x_uniform)); % Manual PDF
        times_mc{k}(i) = toc;
        errors_mc{k}(i) = abs(I - exact_vals(k));
    end
    
    % Importance Sampling
    errors_is{k} = zeros(length(mu_values), length(N_values));
    times_is{k} = zeros(length(mu_values), length(N_values));
    for m = 1:length(mu_values)
        mu = mu_values(m);
        for i = 1:length(N_values)
            N = N_values(i);
            tic;
            x_samples = mu + randn(N, 1); % Sample from N(mu, 1)
            in_interval = (x_samples >= a) & (x_samples <= b);
            % Manual PDF calculations:
            f = normpdf_manual(x_samples);
            g = (1/sqrt(2*pi)) * exp(-(x_samples - mu).^2 / 2); % PDF N(mu,1)
            f_over_g = f ./ g;
            I = mean(in_interval .* f_over_g);
            times_is{k}(m, i) = toc;
            errors_is{k}(m, i) = abs(I - exact_vals(k));
        end
    end
end


% Plotting results for one interval (e.g., [-1, 1])
k = 1; % Change index for other intervals
a = intervals{k}(1);
b = intervals{k}(2);

% Error vs N for Trapezoidal and MC
figure;
loglog(n_trap_values, errors_trap{k}, 'b-o', 'DisplayName', 'Trapezoidal');
hold on;
loglog(N_values, errors_mc{k}, 'r-s', 'DisplayName', 'Monte Carlo');
xlabel('Number of points/samples');
ylabel('Absolute error');
legend;
title(['Error Comparison: [', num2str(a), ', ', num2str(b), ']']);
grid on;

% Time vs N for Trapezoidal and MC
figure;
loglog(n_trap_values, times_trap{k}, 'b-o', 'DisplayName', 'Trapezoidal');
hold on;
loglog(N_values, times_mc{k}, 'r-s', 'DisplayName', 'Monte Carlo');
xlabel('Number of points/samples');
ylabel('Time (s)');
legend;
title(['Time Comparison: [', num2str(a), ', ', num2str(b), ']']);
grid on;

% Error vs mu for Importance Sampling (fixed N)
figure;
N_idx = 3; % Example: N = 1e5
semilogy(mu_values, errors_is{k}(:, N_idx), 'b-o');
xlabel('\mu for Importance Sampling');
ylabel('Absolute error');
title(['Error vs \mu (N=1e5): [', num2str(a), ', ', num2str(b), ']']);
grid on;