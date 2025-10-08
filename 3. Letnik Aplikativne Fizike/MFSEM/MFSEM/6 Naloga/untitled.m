% =========================================================================
% MATEMATIČNO-FIZIKALNI SEMINAR 2024/25
% 6. naloga: Enačbe hoda
% MATLAB Code Implementation
% =========================================================================
clear; close all; clc;

fprintf('MATEMATIČNO-FIZIKALNI SEMINAR 2024/25 - Naloga 6\n');
fprintf('==================================================\n\n');

% =========================================================================
% Task 1: du/dx = a*exp(-2x) - u^4, u(0) = 0
% =========================================================================
fprintf('Task 1: Solving du/dx = a*exp(-2x) - u^4\n');
fprintf('----------------------------------------\n');

% Define the ODE function for Task 1
odefun1 = @(x, u, a) a * exp(-2*x) - u.^4;

% Parameters for Task 1
a_val = 10;
u0 = 0;
x_span = [0, 5]; % Integration interval [x_start, x_end]

% --- Part 1a: Method Comparison for a = 10 ---
fprintf('\nPart 1a: Comparing methods for a = %d\n', a_val);

step_sizes = [0.5, 0.1, 0.05, 0.01]; % Different step sizes to test
method_names = {'Euler', 'Midpoint', 'RK4', 'ode45 (RKF)', 'ode113 (Adams)'};
results = struct(); % Store results for comparison

% --- Reference Solution (using ode45 with very tight tolerance) ---
options_ref = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[x_ref, u_ref] = ode45(@(x, u) odefun1(x, u, a_val), x_span, u0, options_ref);
fprintf('Generated reference solution using ode45 with high accuracy.\n');
%%
figure('Name', 'Task 1: Method Comparison (a=10)');
plot(x_ref, u_ref, 'k-', 'LineWidth', 2, 'DisplayName', 'Reference (ode45)');
hold on;
grid on;
xlabel('x');
ylabel('u(x)');
ylim([0 6])
xlim([0 0.5])
title(sprintf('Comparison of Methods for du/dx = %d*exp(-2x) - u^4, u(0)=0', a_val));
colors = lines(length(step_sizes) * 3 + 2); % Colors for plots
line_styles = {'-', '--', ':', '-.'};
marker_styles = {'o', 's', 'd', '^'};
plot_idx = 1;
%%
final_errors = zeros(length(method_names) - 2, length(step_sizes)); % Store final point error

for i = 1:length(step_sizes)
    h = step_sizes(i);
    fprintf('Processing step size h = %.3f...\n', h);
    x_vec = x_span(1):h:x_span(2);
    n_steps = length(x_vec);

    % --- Euler Method ---
    u_euler = zeros(1, n_steps);
    u_euler(1) = u0;
    for k = 1:(n_steps - 1)
        u_euler(k+1) = u_euler(k) + h * odefun1(x_vec(k), u_euler(k), a_val);
    end
    results(i).h = h;
    results(i).x_euler = x_vec;
    results(i).u_euler = u_euler;
    plot(x_vec, u_euler, 'LineStyle', line_styles{1}, 'Color', colors(plot_idx,:), ...
         'Marker', 'none', 'DisplayName', sprintf('Euler (h=%.2f)', h));
    plot_idx = plot_idx + 1;
    u_ref_interp = interp1(x_ref, u_ref, x_vec(end));
    final_errors(1, i) = abs(u_euler(end) - u_ref_interp);


    % --- Midpoint Method ---
    u_midpoint = zeros(1, n_steps);
    u_midpoint(1) = u0;
    for k = 1:(n_steps - 1)
        k1 = h * odefun1(x_vec(k), u_midpoint(k), a_val);
        k2 = h * odefun1(x_vec(k) + h/2, u_midpoint(k) + k1/2, a_val);
        u_midpoint(k+1) = u_midpoint(k) + k2;
    end
    results(i).x_midpoint = x_vec;
    results(i).u_midpoint = u_midpoint;
    plot(x_vec, u_midpoint, 'LineStyle', line_styles{2}, 'Color', colors(plot_idx,:), ...
         'Marker', 'none', 'DisplayName', sprintf('Midpoint (h=%.2f)', h));
     plot_idx = plot_idx + 1;
    final_errors(2, i) = abs(u_midpoint(end) - u_ref_interp);

    % --- Runge-Kutta 4 (RK4) Method ---
    u_rk4 = zeros(1, n_steps);
    u_rk4(1) = u0;
    for k = 1:(n_steps - 1)
        yk = u_rk4(k);
        xk = x_vec(k);
        k1 = odefun1(xk,       yk,         a_val);
        k2 = odefun1(xk + h/2, yk + h*k1/2, a_val);
        k3 = odefun1(xk + h/2, yk + h*k2/2, a_val);
        k4 = odefun1(xk + h,   yk + h*k3,   a_val);
        u_rk4(k+1) = yk + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    results(i).x_rk4 = x_vec;
    results(i).u_rk4 = u_rk4;
    plot(x_vec, u_rk4, 'LineStyle', line_styles{3}, 'Color', colors(plot_idx,:), ...
         'Marker', 'none', 'DisplayName', sprintf('RK4 (h=%.2f)', h));
    plot_idx = plot_idx + 1;
    final_errors(3, i) = abs(u_rk4(end) - u_ref_interp);

end

% --- MATLAB built-in solvers (for comparison) ---
options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); % Standard tolerances

% % ode45 (Runge-Kutta-Fehlberg)
% [x_ode45, u_ode45] = ode45(@(x, u) odefun1(x, u, a_val), x_span, u0, options_ode);
% plot(x_ode45, u_ode45, '--', 'Color', colors(plot_idx,:), 'LineWidth', 1.5, 'DisplayName', 'ode45 (Built-in)');
% plot_idx = plot_idx + 1;
% 
% % ode113 (Adams-Bashforth-Moulton PECE)
% [x_ode113, u_ode113] = ode113(@(x, u) odefun1(x, u, a_val), x_span, u0, options_ode);
% plot(x_ode113, u_ode113, ':', 'Color', colors(plot_idx,:), 'LineWidth', 1.5, 'DisplayName', 'ode113 (Built-in)');

legend('show', 'Location', 'best');
hold off;

% Display final point errors in a table format
fprintf('\nError at x = %.1f (compared to high-accuracy ode45):\n', x_span(2));
fprintf('%-10s |', 'h');
fprintf(' %-12s |', method_names{1:3});
fprintf('\n%s\n', repmat('-', 1, 10 + 15*3));
for i = 1:length(step_sizes)
    fprintf('%-10.3f |', step_sizes(i));
    fprintf(' %-12.4e |', final_errors(:, i));
    fprintf('\n');
end
fprintf('Note: Lower error indicates higher accuracy.\n');

% --- Accuracy vs Step Size Plot ---
figure('Name', 'Task 1: Accuracy vs Step Size (a=10)');
loglog(step_sizes, final_errors(1,:), 'o-', 'DisplayName', 'Euler Error');
hold on;
loglog(step_sizes, final_errors(2,:), 's-', 'DisplayName', 'Midpoint Error');
loglog(step_sizes, final_errors(3,:), 'd-', 'DisplayName', 'RK4 Error');
% Add lines representing theoretical order O(h), O(h^2), O(h^4)
% loglog(step_sizes, (step_sizes/step_sizes(1)).^1 * final_errors(1,1), 'k--', 'DisplayName', 'O(h)');
% loglog(step_sizes, (step_sizes/step_sizes(1)).^2 * final_errors(2,1), 'k:', 'DisplayName', 'O(h^2)');
% loglog(step_sizes, (step_sizes/step_sizes(1)).^4 * final_errors(3,1), 'k-.', 'DisplayName', 'O(h^4)');
grid on;
xlabel('Step size (h)');
ylabel('Absolute Error at x_{end}');
title('Error Convergence for Different Methods (a=10)');
legend('show', 'Location', 'best');
hold off;

fprintf('\nObservations (a=10):\n');
fprintf('- Euler method is the least accurate (converges as O(h)).\n');
fprintf('- Midpoint method is better (converges as O(h^2) locally, error O(h^2) globally for fixed interval).\n');
fprintf('- RK4 method is significantly more accurate (converges as O(h^4) locally, error O(h^4) globally).\n');
fprintf('- MATLAB''s ode45 and ode113 provide accurate solutions efficiently using adaptive steps.\n');
fprintf('- Stability: For this problem, all methods seem stable for the chosen step sizes.\n');

% --- Part 1b: Family of Solutions for a = 2 to 18 ---
fprintf('\nPart 1b: Calculating family of solutions for a = 2:4:18\n');

a_values = 2:4:18;
chosen_method = 'ode45'; % Choose a reliable method (e.g., ode45 or RK4 with small h)
fprintf('Using %s for calculating the family of solutions.\n', chosen_method);

figure('Name', 'Task 1: Family of Solutions');
hold on;
grid on;
max_temps = zeros(length(a_values), 1);
x_at_max_temps = zeros(length(a_values), 1);
colors_a = jet(length(a_values)); % Different colors for different 'a'

for i = 1:length(a_values)
    a = a_values(i);
    odefun_a = @(x, u) odefun1(x, u, a);

    switch chosen_method
        case 'ode45'
            options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
            [x_sol, u_sol] = ode45(odefun_a, x_span, u0, options);
        case 'RK4'
             h_chosen = 0.01; % Use a small step size for RK4
             x_sol = x_span(1):h_chosen:x_span(2);
             n_steps = length(x_sol);
             u_sol = zeros(1, n_steps);
             u_sol(1) = u0;
             for k = 1:(n_steps - 1)
                yk = u_sol(k);
                xk = x_sol(k);
                k1 = odefun_a(xk,       yk);
                k2 = odefun_a(xk + h_chosen/2, yk + h_chosen*k1/2);
                k3 = odefun_a(xk + h_chosen/2, yk + h_chosen*k2/2);
                k4 = odefun_a(xk + h_chosen,   yk + h_chosen*k3);
                u_sol(k+1) = yk + (h_chosen/6) * (k1 + 2*k2 + 2*k3 + k4);
             end
             u_sol = u_sol'; % Transpose for consistency
             x_sol = x_sol';
        otherwise
            error('Chosen method not implemented for family calculation.');
    end

    plot(x_sol, u_sol, 'Color', colors_a(i,:), 'DisplayName', sprintf('a = %d', a));

    % Find maximum temperature and its location
    [max_u, idx_max] = max(u_sol);
    max_temps(i) = max_u;
    x_at_max_temps(i) = x_sol(idx_max);
end

xlabel('x');
ylabel('u(x)');
title('Family of Solutions for du/dx = a*exp(-2x) - u^4');
legend('show', 'Location', 'best');
hold off;

% Display max temperatures
fprintf('\nMaximum Temperature (u_max) and Location (x_at_max) for different ''a'':\n');
fprintf('  a   |  u_max   | x_at_max \n');
fprintf('------|----------|----------\n');
for i = 1:length(a_values)
    fprintf(' %4d | %8.4f | %8.4f \n', a_values(i), max_temps(i), x_at_max_temps(i));
end

fprintf('\nPremislek (Thought Experiment):\n');
fprintf('To accurately find the maximum temperature and its time:\n');
fprintf('1. Use a high-order adaptive step-size method (like ode45 or ode113) with tight tolerances.\n');
fprintf('2. Locate the maximum in the resulting dense output or use event detection.\n');
fprintf('3. Alternatively, find the root of the derivative du/dx = a*exp(-2x) - u^4 = 0.\n');
fprintf('   This means finding x where a*exp(-2x) = u(x)^4. This requires solving the ODE simultaneously\n');
fprintf('   or interpolating the solution u(x) and then finding the root.\n');
fprintf('   Using event detection capabilities of solvers like ode45 is often the most robust way.\n');


% =========================================================================
% Task 2: dy/dt = y^2 + 2t^2, y(0) = 1
% =========================================================================
fprintf('\n\nTask 2: Solving dy/dt = y^2 + 2*t^2\n');
fprintf('------------------------------------\n');

% Define the ODE function for Task 2
odefun2 = @(t, y) y.^2 + 2*t.^2;

% Parameters for Task 2
y0_task2 = 1;
t_span_task2 = [0, 0.95]; % Integrate close to the point of divergence (~0.9)
h_values_task2 = [0.1, 0.01, 0.001];

fprintf('Comparing methods for dy/dt = y^2 + 2*t^2, y(0)=1\n');
fprintf('Solution diverges near t = 0.9.\n');

figure('Name', 'Task 2: Method Comparison');
hold on;
grid on;
xlabel('t');
ylabel('y(t)');
title('Comparison of Methods for dy/dt = y^2 + 2t^2, y(0)=1');
colors_t2 = lines(length(h_values_task2) * 3 + 2);
plot_idx_t2 = 1;

% --- Reference Solution (using ode45 with tight tolerance) ---
options_ref_t2 = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
% Need to handle potential warnings/errors near singularity
try
    [t_ref_t2, y_ref_t2] = ode45(odefun2, t_span_task2, y0_task2, options_ref_t2);
    plot(t_ref_t2, y_ref_t2, 'k-', 'LineWidth', 2, 'DisplayName', 'Reference (ode45)');
catch ME
    warning('ode45 failed near singularity, plotting partial result if available.');
    if exist('t_ref_t2', 'var') && ~isempty(t_ref_t2)
         plot(t_ref_t2, y_ref_t2, 'k-', 'LineWidth', 2, 'DisplayName', 'Reference (ode45)');
    end
    fprintf('Solver Warning/Error: %s\n', ME.message);
end


for i = 1:length(h_values_task2)
    h = h_values_task2(i);
    t_vec = t_span_task2(1):h:t_span_task2(2);
    n_steps = length(t_vec);
     if t_vec(end) < t_span_task2(2) % Ensure the span is covered if h doesn't divide evenly
         t_vec = [t_vec, t_span_task2(2)];
         n_steps = n_steps + 1;
     end

    fprintf('Processing step size h = %.3f...\n', h);

    % --- Euler Method ---
    y_euler = zeros(1, n_steps);
    y_euler(1) = y0_task2;
    for k = 1:(n_steps - 1)
        % Check for potential overflow before calculation
        f_eval = odefun2(t_vec(k), y_euler(k));
        if isinf(f_eval) || isnan(f_eval) || abs(y_euler(k)) > 1e10 % Stop if diverging uncontrollably
             y_euler(k+1:end) = NaN; % Mark rest as NaN
             fprintf('Euler diverged at t=%.4f for h=%.3f\n', t_vec(k), h);
             break;
        end
        y_euler(k+1) = y_euler(k) + h * f_eval;
    end
    plot(t_vec, y_euler, '--', 'Color', colors_t2(plot_idx_t2,:), 'DisplayName', sprintf('Euler (h=%.3f)', h));
    plot_idx_t2 = plot_idx_t2 + 1;


    % --- Midpoint Method ---
    y_midpoint = zeros(1, n_steps);
    y_midpoint(1) = y0_task2;
     for k = 1:(n_steps - 1)
        yk = y_midpoint(k);
        tk = t_vec(k);
        if abs(yk) > 1e10 % Stop if diverging
             y_midpoint(k+1:end) = NaN;
             fprintf('Midpoint diverged at t=%.4f for h=%.3f\n', tk, h);
             break;
        end
        k1_m = h * odefun2(tk, yk);
        k2_m = h * odefun2(tk + h/2, yk + k1_m/2);
         if isinf(k1_m) || isnan(k1_m) || isinf(k2_m) || isnan(k2_m)
             y_midpoint(k+1:end) = NaN;
             fprintf('Midpoint diverged (k val) at t=%.4f for h=%.3f\n', tk, h);
             break;
         end
        y_midpoint(k+1) = yk + k2_m;
    end
    plot(t_vec, y_midpoint, '-.', 'Color', colors_t2(plot_idx_t2,:), 'DisplayName', sprintf('Midpoint (h=%.3f)', h));
    plot_idx_t2 = plot_idx_t2 + 1;

    % --- Runge-Kutta 4 (RK4) Method ---
    y_rk4 = zeros(1, n_steps);
    y_rk4(1) = y0_task2;
    for k = 1:(n_steps - 1)
        yk = y_rk4(k);
        tk = t_vec(k);
         if abs(yk) > 1e10 % Stop if diverging
             y_rk4(k+1:end) = NaN;
             fprintf('RK4 diverged at t=%.4f for h=%.3f\n', tk, h);
             break;
         end
        k1 = odefun2(tk,       yk);
        k2 = odefun2(tk + h/2, yk + h*k1/2);
        k3 = odefun2(tk + h/2, yk + h*k2/2);
        k4 = odefun2(tk + h,   yk + h*k3);
         if isinf(k1) || isnan(k1) || isinf(k2) || isnan(k2) || isinf(k3) || isnan(k3) || isinf(k4) || isnan(k4)
            y_rk4(k+1:end) = NaN;
            fprintf('RK4 diverged (k val) at t=%.4f for h=%.3f\n', tk, h);
            break;
         end
        y_rk4(k+1) = yk + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    plot(t_vec, y_rk4, ':', 'Color', colors_t2(plot_idx_t2,:), 'DisplayName', sprintf('RK4 (h=%.3f)', h));
    plot_idx_t2 = plot_idx_t2 + 1;

end

% --- MATLAB built-in solvers (example) ---
% Using ode45 and ode113 for comparison with default settings
% try
%     [t_ode45_t2, y_ode45_t2] = ode45(odefun2, t_span_task2, y0_task2);
%     plot(t_ode45_t2, y_ode45_t2, '-', 'Color', colors_t2(plot_idx_t2,:), 'LineWidth', 1.5, 'DisplayName', 'ode45 (Built-in)');
% catch ME
%     warning('ode45 failed near singularity in comparison run.');
%      if exist('t_ode45_t2', 'var') && ~isempty(t_ode45_t2)
%          plot(t_ode45_t2, y_ode45_t2, '-', 'Color', colors_t2(plot_idx_t2,:), 'LineWidth', 1.5, 'DisplayName', 'ode45 (Built-in)');
%     end
% end
% plot_idx_t2 = plot_idx_t2 + 1;
% 
% try
%     [t_ode113_t2, y_ode113_t2] = ode113(odefun2, t_span_task2, y0_task2);
%     plot(t_ode113_t2, y_ode113_t2, '--', 'Color', colors_t2(plot_idx_t2,:), 'LineWidth', 1.5, 'DisplayName', 'ode113 (Built-in)');
% catch ME
%      warning('ode113 failed near singularity in comparison run.');
%      if exist('t_ode113_t2', 'var') && ~isempty(t_ode113_t2)
%          plot(t_ode113_t2, y_ode113_t2, '--', 'Color', colors_t2(plot_idx_t2,:), 'LineWidth', 1.5, 'DisplayName', 'ode113 (Built-in)');
%      end
% end


% Adjust plot limits if necessary due to divergence
max_y_plot = 50; % Set a reasonable upper limit for y-axis
current_ylim = ylim;
if current_ylim(2) > max_y_plot
    ylim([current_ylim(1), max_y_plot]);
    fprintf('Adjusted y-axis limit for better visualization due to divergence.\n')
end


legend('show', 'Location', 'NorthWest');
hold off;

fprintf('\nObservations (Task 2):\n');
fprintf('- The solution clearly diverges rapidly as t approaches ~0.9.\n');
fprintf('- With larger step sizes (e.g., h=0.1), simpler methods like Euler significantly deviate earlier and might overshoot the singularity.\n');
fprintf('- Smaller step sizes (h=0.01, h=0.001) allow fixed-step methods to follow the true diverging solution more closely before failing.\n');
fprintf('- RK4 generally performs better than Euler and Midpoint for the same step size.\n');
fprintf('- Adaptive step-size methods (ode45, ode113) automatically reduce step size near the singularity, providing a more accurate path until they might eventually fail or stop if the singularity is within the interval.\n');
fprintf('- Be cautious when interpreting results near a singularity, especially with fixed-step methods.\n');

fprintf('\nEnd of simulation.\n');

% =========================================================================
% Task 2: Error Divergence Plots
% =========================================================================
fprintf('\n\nTask 2: Plotting Error Divergence vs Time\n');
fprintf('-------------------------------------------\n');

t_span_error_plot = [0, 0.9]; % Focus on the range leading to divergence

% --- Recalculate High-Accuracy Reference Solution for Error Plot ---
% (Ensure it covers the specific range [0, 0.9] densely)
options_ref_t2_err = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
t_ref_err = []; y_ref_err = []; % Initialize
fprintf('Calculating high-accuracy reference solution for error plots over [%.1f, %.1f]...\n', t_span_error_plot);
try
    % Use ode45 to get dense output within the interval
    sol_ref = ode45(odefun2, t_span_error_plot, y0_task2, options_ref_t2_err);
    % Evaluate the solution at many points for smooth interpolation later
    t_ref_err = linspace(t_span_error_plot(1), t_span_error_plot(2), 1000)'; % Dense time points
    y_ref_err = deval(sol_ref, t_ref_err)';
    fprintf('Reference solution calculated successfully.\n');
catch ME
    warning('Reference ode45 failed near singularity even for error plot. Error: %s', ME.message);
     % If it failed, we might still have partial data from sol_ref if it exists
     if exist('sol_ref', 'var') && isstruct(sol_ref) && ~isempty(sol_ref.x)
        t_ref_err = sol_ref.x';
        y_ref_err = sol_ref.y';
        fprintf('Using partial reference solution up to t=%.4f.\n', t_ref_err(end));
     else
        fprintf('Could not obtain reference solution for error plot.\n');
        t_ref_err = []; % Ensure it's empty if failed completely
     end
end

% --- Loop through step sizes again for error calculation ---
if ~isempty(t_ref_err) % Only proceed if we have a reference solution
    for i = 1:length(h_values_task2)
        h = h_values_task2(i);
        t_vec = t_span_error_plot(1):h:t_span_error_plot(2);
        % Adjust t_vec if h doesn't divide evenly, only up to the end of span
         if isempty(t_vec) % Handle case where h is larger than span
             continue;
         end
         if t_vec(end) < t_span_error_plot(2)
             t_vec = [t_vec, t_span_error_plot(2)];
         end
         % Make sure t_vec does not exceed the calculated reference range end
         t_vec = t_vec(t_vec <= t_ref_err(end));
         if length(t_vec) < 2 % Need at least two points
             fprintf('Skipping error plot for h=%.3f, too few points in range.\n',h);
             continue;
         end
        n_steps = length(t_vec);

        fprintf('Calculating errors for h = %.3f...\n', h);

        % --- Recalculate Numerical Solutions (Euler, Midpoint, RK4) ---
        % (Copy the loops from the previous section, ensuring they use t_vec)

        % Euler
        y_euler_err = zeros(1, n_steps); y_euler_err(1) = y0_task2;
        for k = 1:(n_steps - 1)
            current_h = t_vec(k+1) - t_vec(k); f_eval = odefun2(t_vec(k), y_euler_err(k));
            if abs(y_euler_err(k)) > 1e10 || isnan(f_eval) || isinf(f_eval)
                 y_euler_err(k+1:end) = NaN; break;
            end; y_euler_err(k+1) = y_euler_err(k) + current_h * f_eval;
        end

        % Midpoint
        y_midpoint_err = zeros(1, n_steps); y_midpoint_err(1) = y0_task2;
         for k = 1:(n_steps - 1)
            current_h = t_vec(k+1) - t_vec(k); yk = y_midpoint_err(k); tk = t_vec(k);
            if abs(yk) > 1e10 || isnan(yk) || isinf(yk); y_midpoint_err(k+1:end) = NaN; break; end
            k1_m = current_h * odefun2(tk, yk); y_mid = yk + k1_m/2;
            if abs(y_mid) > 1e10 || isnan(k1_m) || isinf(k1_m) || isnan(y_mid) || isinf(y_mid); y_midpoint_err(k+1:end) = NaN; break; end
            k2_m = current_h * odefun2(tk + current_h/2, y_mid);
            if isnan(k2_m) || isinf(k2_m); y_midpoint_err(k+1:end) = NaN; break; end
            y_midpoint_err(k+1) = yk + k2_m;
        end

        % RK4
        y_rk4_err = zeros(1, n_steps); y_rk4_err(1) = y0_task2;
        for k = 1:(n_steps - 1)
            current_h = t_vec(k+1) - t_vec(k); yk = y_rk4_err(k); tk = t_vec(k);
            if abs(yk) > 1e10 || isnan(yk) || isinf(yk); y_rk4_err(k+1:end) = NaN; break; end
            k1 = odefun2(tk, yk); k2 = odefun2(tk + current_h/2, yk + current_h*k1/2);
            k3 = odefun2(tk + current_h/2, yk + current_h*k2/2); k4 = odefun2(tk + current_h, yk + current_h*k3);
            if isnan(k1)||isinf(k1)||isnan(k2)||isinf(k2)||isnan(k3)||isinf(k3)||isnan(k4)||isinf(k4)||abs(k1)>1e10||abs(k2)>1e10||abs(k3)>1e10||abs(k4)>1e10
                y_rk4_err(k+1:end) = NaN; break;
            end; y_rk4_err(k+1) = yk + (current_h/6) * (k1 + 2*k2 + 2*k3 + k4);
        end

        % --- Interpolate Reference Solution onto t_vec ---
        % Use pchip for shape-preserving interpolation, good for steep gradients
        y_ref_interp = interp1(t_ref_err, y_ref_err, t_vec, 'pchip');

        % --- Calculate Absolute Errors ---
        error_euler = abs(y_euler_err - y_ref_interp);
        error_midpoint = abs(y_midpoint_err - y_ref_interp);
        error_rk4 = abs(y_rk4_err - y_ref_interp);

        % --- Create Error Plot for this h ---
        figure('Name', sprintf('Task 2: Error Divergence (h=%.3f)', h));
        semilogy(t_vec, error_euler, 'ro-', 'MarkerSize', 4, 'DisplayName', 'Euler Error');
        hold on;
        semilogy(t_vec, error_midpoint, 'gs--', 'MarkerSize', 4, 'DisplayName', 'Midpoint Error');
        semilogy(t_vec, error_rk4, 'bd:', 'MarkerSize', 4, 'DisplayName', 'RK4 Error');
        grid on;
        xlabel('Time (t)');
        ylabel('Absolute Error |y_{num} - y_{ref}|');
        title(sprintf('Task 2: Growth of Absolute Error (h=%.3f)', h));
        legend('show', 'Location', 'NorthWest'); % Usually errors start small
        hold off;
         % Optional: Adjust y-limits if needed, e.g., ylim([1e-10 1e2])
    end
else
    fprintf('Skipping error divergence plots because reference solution calculation failed.\n');
end