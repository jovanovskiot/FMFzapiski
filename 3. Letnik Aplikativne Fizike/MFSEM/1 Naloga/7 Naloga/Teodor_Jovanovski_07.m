% MATEMATIČNO-FIZIKALNI SEMINAR 2024/25
% 7. naloga: Newtonov zakon - Simulation of a Simple Pendulum

clear; clc; close all;

% --- Parameters ---
x0 = 1.0;       % Initial angle (rad)
v0 = 0.0;       % Initial angular velocity (rad/s)
Tmax = 60;      % Total simulation time (s) -> ~10 periods
% Try different step sizes: h = 0.1; h = 0.01; h = 0.001;
h = 0.01;       % Time step (s)

% --- Setup ---
N = round(Tmax / h); % Number of steps
t = (0:N)*h;         % Time vector

% Define the acceleration function (f(x) in the problem text, here a(x))
acceleration = @(x) -sin(x);

% Define the energy function (Hamiltonian)
energy = @(x, v) 1 - cos(x) + 0.5 * v.^2;

% --- Initialization ---
% Arrays to store results [time, angle, velocity, energy]
results_euler = zeros(N+1, 4);
results_verlet = zeros(N+1, 4);
results_rk4 = zeros(N+1, 4);

% Set initial conditions
results_euler(1,:) = [t(1), x0, v0, energy(x0, v0)];
results_verlet(1,:) = [t(1), x0, v0, energy(x0, v0)];
results_rk4(1,:) = [t(1), x0, v0, energy(x0, v0)];

% --- Simulation Loop ---

% Temporary variables for RK4
x_rk4 = x0;
v_rk4 = v0;

for n = 1:N
    % Current state for each method
    t_n = results_euler(n, 1);
    x_euler_n = results_euler(n, 2);
    v_euler_n = results_euler(n, 3);

    x_verlet_n = results_verlet(n, 2);
    v_verlet_n = results_verlet(n, 3);

    % --- 1. Forward Euler Method ---
    a_euler_n = acceleration(x_euler_n);
    x_euler_next = x_euler_n + h * v_euler_n;
    v_euler_next = v_euler_n + h * a_euler_n;
    results_euler(n+1, :) = [t(n+1), x_euler_next, v_euler_next, energy(x_euler_next, v_euler_next)];

    % --- 2. Velocity Verlet Method ---
    % (Corresponds to the y_n+1, v_n+1 formulas in the text)
    a_verlet_n = acceleration(x_verlet_n);
    % Calculate next position
    x_verlet_next = x_verlet_n + h * v_verlet_n + 0.5 * h^2 * a_verlet_n;
    % Calculate acceleration at the *next* position
    a_verlet_next = acceleration(x_verlet_next);
    % Calculate next velocity
    v_verlet_next = v_verlet_n + 0.5 * h * (a_verlet_n + a_verlet_next);
    results_verlet(n+1, :) = [t(n+1), x_verlet_next, v_verlet_next, energy(x_verlet_next, v_verlet_next)];

    % --- 3. Runge-Kutta 4 (RK4) Method ---
    % Define the system derivative function for RK4: dy/dt = f(t, y)
    % where y = [x; v], so dy/dt = [v; acceleration(x)]
    f_rk4 = @(t_rk, y_rk) [y_rk(2); acceleration(y_rk(1))];

    y_rk4_n = [x_rk4; v_rk4]; % State vector at step n

    k1 = h * f_rk4(t_n, y_rk4_n);
    k2 = h * f_rk4(t_n + h/2, y_rk4_n + k1/2);
    k3 = h * f_rk4(t_n + h/2, y_rk4_n + k2/2);
    k4 = h * f_rk4(t_n + h, y_rk4_n + k3);

    y_rk4_next = y_rk4_n + (k1 + 2*k2 + 2*k3 + k4) / 6;

    x_rk4 = y_rk4_next(1); % Update x for the next iteration's start
    v_rk4 = y_rk4_next(2); % Update v for the next iteration's start
    results_rk4(n+1, :) = [t(n+1), x_rk4, v_rk4, energy(x_rk4, v_rk4)];

end

% --- Analysis and Plotting ---

figure('Position', [100, 100, 1000, 800]); % Create a larger figure window

% 1. Plot Angle (Position) vs Time
subplot(2, 1, 1); % Create subplot (2 rows, 1 column, plot 1)
plot(results_euler(:,1), results_euler(:,2), 'r:', 'LineWidth', 1.5); hold on;
plot(results_verlet(:,1), results_verlet(:,2), 'b--', 'LineWidth', 1.5);
plot(results_rk4(:,1), results_rk4(:,2), 'k-', 'LineWidth', 1);
hold off;
title(['Pendulum Angle vs Time (h = ', num2str(h), ')']);
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Euler', 'Velocity Verlet', 'RK4', 'Location', 'northeast');
grid on;
ylim([-1.1*x0, 1.1*x0]); % Adjust y-axis limits slightly beyond initial amplitude

% 2. Plot Energy vs Time
subplot(2, 1, 2); % Create subplot (2 rows, 1 column, plot 2)
plot(results_euler(:,1), results_euler(:,4), 'r:', 'LineWidth', 1.5); hold on;
plot(results_verlet(:,1), results_verlet(:,4), 'b--', 'LineWidth', 1.5);
plot(results_rk4(:,1), results_rk4(:,4), 'k-', 'LineWidth', 1);
hold off;
title(['Energy vs Time (h = ', num2str(h), ')']);
xlabel('Time (s)');
ylabel('Energy (Arbitrary Units)');
legend('Euler', 'Velocity Verlet', 'RK4', 'Location', 'southeast');
grid on;
% Zoom in on energy variation if needed
initial_energy = results_rk4(1,4);
ylim([initial_energy - 0.1*abs(initial_energy), initial_energy + 0.1*abs(initial_energy)]); % Adjust y-axis limit for energy plot if needed
% ylim([min(results_verlet(:,4))*0.99, max(results_verlet(:,4))*1.01]); % Alternative zoom focused on Verlet


% --- Discussion based on plots (add this as comments or separate text) ---

fprintf('\n--- Analysis for h = %f ---\n', h);

% Accuracy:
% - Compare the curves in the Angle vs Time plot. For small h (e.g., 0.001),
%   Verlet and RK4 should be very close. Euler will likely diverge significantly
%   even for h=0.01 over long times.
% - To check for 3 decimal places: Run with h, then h/2, then h/4. Compare
%   the final angle x(Tmax) or the angle at a specific time (e.g., t=10s).
%   When the values agree to 3-4 decimal places, the smaller h is sufficient.
%   RK4 will likely achieve this accuracy with a larger h than Verlet or Euler.
%   Example check:
final_x_euler = results_euler(end, 2);
final_x_verlet = results_verlet(end, 2);
final_x_rk4 = results_rk4(end, 2);
fprintf('Final angle (x) at T=%.1f s:\n', Tmax);
fprintf('  Euler:  %.6f\n', final_x_euler);
fprintf('  Verlet: %.6f\n', final_x_verlet);
fprintf('  RK4:    %.6f\n', final_x_rk4);


% Stability (Amplitude):
% - Observe the Angle vs Time plot over the full Tmax.
% - Euler: Amplitude often increases unrealistically (unstable).
% - Verlet: Amplitude should remain very constant (stable).
% - RK4: Amplitude should also be very stable, potentially even better than
%   Verlet for a given *large* h, but less guaranteed structurally.
fprintf('\nStability (Amplitude):\n');
fprintf('  Euler: Check plot for amplitude increase.\n');
fprintf('  Verlet: Check plot for constant amplitude.\n');
fprintf('  RK4: Check plot for constant amplitude.\n');


% Energy Conservation:
% - Observe the Energy vs Time plot.
initial_energy = results_rk4(1,4);
final_energy_euler = results_euler(end, 4);
final_energy_verlet = results_verlet(end, 4);
final_energy_rk4 = results_rk4(end, 4);
fprintf('\nEnergy Conservation (Initial E = %.6f):\n', initial_energy);
fprintf('  Euler Final E:  %.6f (Change: %.2e)\n', final_energy_euler, final_energy_euler - initial_energy);
fprintf('  Verlet Final E: %.6f (Change: %.2e)\n', final_energy_verlet, final_energy_verlet - initial_energy);
fprintf('  RK4 Final E:    %.6f (Change: %.2e)\n', final_energy_rk4, final_energy_rk4 - initial_energy);
% - Euler: Energy likely drifts significantly (usually increases).
% - Verlet: Energy should oscillate very slightly around the initial value
%   but show no long-term drift (excellent conservation).
% - RK4: Energy should be well-conserved for small h, possibly with a tiny
%   drift over very long times, but generally much better than Euler.

% --- How to Use ---
% 1. Run the script with h = 0.1. Observe the plots and printed analysis.
% 2. Change h to 0.01 and run again. Note the improvements, especially for Euler and Verlet.
% 3. Change h to 0.001 and run again. Verlet and RK4 solutions should now be visually identical and highly accurate. Euler might still show issues over 60s.
% 4. To find the step size for 3 decimal places: Compare the final angle (or angle at a specific time) for h=0.01 and h=0.001 using RK4 as a reference. If they match to 3-4 decimals, h=0.01 might be sufficient for RK4. Repeat comparison for Verlet. Euler will likely require a much smaller h.