% MATEMATIČNO-FIZIKALNI SEMINAR 2024/25
% 7. naloga: Newtonov zakon - Simulation of a Simple Pendulum
% Includes comparison with ode45 reference solution

clear; clc; close all;

% --- Parameters ---
x0 = 1.0;       % Initial angle (rad)
v0 = 0.0;       % Initial angular velocity (rad/s)
Tmax = 60;      % Total simulation time (s) -> ~10 periods
% Try different step sizes: h = 0.1; h = 0.01; h = 0.001;
h = 0.001;       % Time step for fixed-step methods (s)

% --- Setup ---
N = round(Tmax / h); % Number of steps for fixed-step methods
t = linspace(0, Tmax, N+1); % Time vector for fixed-step methods

% Define the acceleration function
acceleration = @(x) -sin(x);

% Define the energy function (Hamiltonian)
energy = @(x, v) 1 - cos(x) + 0.5 * v.^2;

% --- ODE Function for ode45 ---
% Input: t (time), y (state vector [x; v])
% Output: dy/dt (derivative vector [v; acceleration(x)])
pendulum_ode = @(t, y) [y(2); acceleration(y(1))];

% --- Generate Reference Solution using ode45 ---
fprintf('Calculating reference solution using ode45...\n');
y0_ode = [x0; v0]; % Initial condition vector for ode45
options_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12); % Set high accuracy
[t_ref, y_ref] = ode45(pendulum_ode, t, y0_ode, options_ode); % Solve ODE
% y_ref(:,1) contains angle (x), y_ref(:,2) contains velocity (v)
% t_ref will be the same as t because we specified output times
x_ref = y_ref(:,1);
v_ref = y_ref(:,2);
energy_ref = energy(x_ref, v_ref);
fprintf('Reference solution calculated.\n\n');

% --- Initialization for Fixed-Step Methods ---
results_euler = zeros(N+1, 4);
results_verlet = zeros(N+1, 4);
results_rk4 = zeros(N+1, 4);

% Set initial conditions
results_euler(1,:) = [t(1), x0, v0, energy(x0, v0)];
results_verlet(1,:) = [t(1), x0, v0, energy(x0, v0)];
results_rk4(1,:) = [t(1), x0, v0, energy(x0, v0)];

% --- Simulation Loop (Fixed-Step Methods) ---
fprintf('Running fixed-step simulations...\n');

% Temporary variables for RK4
x_rk4 = x0;
v_rk4 = v0;

for n = 1:N
    % Current state for each method
    t_n = results_euler(n, 1); % same as t(n)
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
    a_verlet_n = acceleration(x_verlet_n);
    x_verlet_next = x_verlet_n + h * v_verlet_n + 0.5 * h^2 * a_verlet_n;
    a_verlet_next = acceleration(x_verlet_next);
    v_verlet_next = v_verlet_n + 0.5 * h * (a_verlet_n + a_verlet_next);
    results_verlet(n+1, :) = [t(n+1), x_verlet_next, v_verlet_next, energy(x_verlet_next, v_verlet_next)];

    % --- 3. Runge-Kutta 4 (RK4) Method ---
    f_rk4 = @(t_rk, y_rk) [y_rk(2); acceleration(y_rk(1))]; % Use the same form as ode45
    y_rk4_n = [x_rk4; v_rk4];
    k1 = h * f_rk4(t_n, y_rk4_n);
    k2 = h * f_rk4(t_n + h/2, y_rk4_n + k1/2);
    k3 = h * f_rk4(t_n + h/2, y_rk4_n + k2/2);
    k4 = h * f_rk4(t_n + h, y_rk4_n + k3);
    y_rk4_next = y_rk4_n + (k1 + 2*k2 + 2*k3 + k4) / 6;
    x_rk4 = y_rk4_next(1);
    v_rk4 = y_rk4_next(2);
    results_rk4(n+1, :) = [t(n+1), x_rk4, v_rk4, energy(x_rk4, v_rk4)];

end
fprintf('Fixed-step simulations finished.\n');

% --- Analysis and Plotting ---

figure('Position', [100, 100, 1200, 900]); % Slightly wider figure

% 1. Plot Angle (Position) vs Time
subplot(2, 2, 1); % Top-left
plot(results_euler(:,1), results_euler(:,2), 'r:', 'LineWidth', 1.5); hold on;
plot(results_verlet(:,1), results_verlet(:,2), 'b--', 'LineWidth', 1.5);
plot(results_rk4(:,1), results_rk4(:,2), 'k-', 'LineWidth', 1);
plot(t_ref, x_ref, 'g-', 'LineWidth', 2); % Add reference solution
hold off;
title(['Pendulum Angle vs Time (h = ', num2str(h), ')']);
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Euler', 'Velocity Verlet', 'RK4', 'Reference (ode45)', 'Location', 'northeast');
grid on;
ylim([-1.1*x0, 1.1*x0]);

% 2. Plot Energy vs Time
subplot(2, 2, 2); % Top-right
plot(results_euler(:,1), results_euler(:,4), 'r:', 'LineWidth', 1.5); hold on;
plot(results_verlet(:,1), results_verlet(:,4), 'b--', 'LineWidth', 1.5);
plot(results_rk4(:,1), results_rk4(:,4), 'k-', 'LineWidth', 1);
plot(t_ref, energy_ref, 'g-', 'LineWidth', 2); % Add reference energy
hold off;
title(['Energy vs Time (h = ', num2str(h), ')']);
xlabel('Time (s)');
ylabel('Energy (Arbitrary Units)');
legend('Euler', 'Velocity Verlet', 'RK4', 'Reference (ode45)', 'Location', 'southeast');
grid on;
initial_energy = energy_ref(1); % Use reference initial energy
max_energy_dev = max(abs(energy_ref - initial_energy));
ylim_margin = max(0.05 * abs(initial_energy), 5 * max_energy_dev); % Sensible margin
ylim([initial_energy - ylim_margin, initial_energy + ylim_margin]);

% 3. Plot Angle Error vs Time
subplot(2, 2, 3); % Bottom-left
error_euler = results_euler(:,2) - x_ref;
error_verlet = results_verlet(:,2) - x_ref;
error_rk4 = results_rk4(:,2) - x_ref;
plot(t, error_euler, 'r:', 'LineWidth', 1); hold on;
plot(t, error_verlet, 'b--', 'LineWidth', 1);
plot(t, error_rk4, 'k-', 'LineWidth', 1);
hold off;
title(['Angle Error vs Time (Difference from ode45) (h = ', num2str(h), ')']);
xlabel('Time (s)');
ylabel('Angle Error (rad)');
legend('Euler Error', 'Verlet Error', 'RK4 Error', 'Location', 'northeast');
grid on;

% 4. Placeholder for potential future plot or keep empty
subplot(2, 2, 4);
axis off; % Turn off axes for the empty subplot

% --- Error Calculation and Display ---
fprintf('\n--- Analysis for h = %f ---\n', h);

% Final Angle Error
final_x_euler = results_euler(end, 2);
final_x_verlet = results_verlet(end, 2);
final_x_rk4 = results_rk4(end, 2);
final_x_ref = x_ref(end);

fprintf('Final angle (x) at T=%.1f s:\n', Tmax);
fprintf('  Reference (ode45): %.8f\n', final_x_ref);
fprintf('  Euler:  %.8f (Error: %.2e)\n', final_x_euler, final_x_euler - final_x_ref);
fprintf('  Verlet: %.8f (Error: %.2e)\n', final_x_verlet, final_x_verlet - final_x_ref);
fprintf('  RK4:    %.8f (Error: %.2e)\n', final_x_rk4, final_x_rk4 - final_x_ref);

% Root Mean Square Error (RMSE) for Angle
rmse_euler = sqrt(mean(error_euler.^2));
rmse_verlet = sqrt(mean(error_verlet.^2));
rmse_rk4 = sqrt(mean(error_rk4.^2));

fprintf('\nRoot Mean Square Error (RMSE) in Angle over %.1f s:\n', Tmax);
fprintf('  Euler:  %.3e\n', rmse_euler);
fprintf('  Verlet: %.3e\n', rmse_verlet);
fprintf('  RK4:    %.3e\n', rmse_rk4);

% Energy Conservation Summary
final_energy_euler = results_euler(end, 4);
final_energy_verlet = results_verlet(end, 4);
final_energy_rk4 = results_rk4(end, 4);
final_energy_ref = energy_ref(end);

fprintf('\nEnergy Conservation (Initial E = %.8f):\n', initial_energy);
fprintf('  Ref (ode45) Final E: %.8f (Change: %.2e)\n', final_energy_ref, final_energy_ref - initial_energy);
fprintf('  Euler Final E:       %.8f (Change: %.2e)\n', final_energy_euler, final_energy_euler - initial_energy);
fprintf('  Verlet Final E:      %.8f (Change: %.2e)\n', final_energy_verlet, final_energy_verlet - initial_energy);
fprintf('  RK4 Final E:         %.8f (Change: %.2e)\n', final_energy_rk4, final_energy_rk4 - initial_energy);

% --- How to Use ---
% 1. Run with different h values (0.1, 0.01, 0.001).
% 2. Observe plots:
%    - Angle Plot: How close are the methods to the green reference line?
%    - Energy Plot: How flat is the energy curve? Compare to the very flat green reference.
%    - Error Plot: Shows how the error accumulates over time for each method.
% 3. Check Command Window output:
%    - Final Angle Error: Quantifies accuracy at the end point.
%    - RMSE: Gives an average measure of accuracy over the whole simulation.
%    - Energy Change: Quantifies energy conservation.
% 4. Accuracy Check (3 decimal places): Look at the RMSE or Final Error. For RK4 with h=0.01, the error might be around 1e-6 or better. For Verlet, it might be 1e-4 to 1e-5. You need the error to be significantly smaller than 0.001 (e.g., < 1e-4) to trust 3 decimal places. Find the `h` for Verlet and RK4 that achieve this. Euler will likely need a much smaller `h`.