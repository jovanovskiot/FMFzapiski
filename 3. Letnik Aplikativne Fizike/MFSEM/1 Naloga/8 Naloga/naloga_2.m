% MATEMATIČNO-FIZIKALNI SEMINAR 2024/25
% 8. naloga: Robni problem lastnih vrednosti - Glavna Naloga

clear;
clc;
close all;

% --- Parameters ---
hbar = 1; % Reduced Planck's constant (can be set to 1 for simplicity)
m = 0.5;  % Mass of the particle (setting 2m=1 simplifies Schrödinger eq. to -d^2psi/dx^2 = E psi)
a = 1;    % Width of the infinite potential well (from x=0 to x=a)

N = 100;  % Number of discretization intervals for FDM (N+1 points total)
num_states_to_find = 4; % Number of lowest energy states to find and compare

% --- Analytical Solution ---
E_analytical = @(n) (n * pi * hbar / a)^2 / (2 * m);
psi_analytical = @(x, n) sqrt(2/a) * sin(n * pi * x / a);

fprintf('--- Analytical Solutions (for n=1 to %d) ---\n', num_states_to_find);
E_an_vals = zeros(num_states_to_find, 1);
for n = 1:num_states_to_find
    E_an_vals(n) = E_analytical(n);
    fprintf('E_%d = %.6f\n', n, E_an_vals(n));
end
fprintf('\n');

% --- Finite Difference Method (FDM) ---
fprintf('--- Finite Difference Method (N=%d intervals) ---\n', N);

h = a / N; % Step size
x_internal = linspace(h, a-h, N-1)';
main_diag = 2 * ones(N-1, 1);
off_diag  = -1 * ones(N-2, 1);
A_fdm_scaled = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
A_fdm = (1/h^2) * A_fdm_scaled;

[eigenvectors_fdm_internal, eigenvalues_fdm_raw] = eig(A_fdm);
E_fdm_vals_num = diag(eigenvalues_fdm_raw);
[E_fdm_vals_num_sorted, sort_idx_fdm] = sort(E_fdm_vals_num);
eigenvectors_fdm_internal_sorted = eigenvectors_fdm_internal(:, sort_idx_fdm);

E_fdm_vals_phys = E_fdm_vals_num_sorted * (hbar^2 / (2*m));

psi_fdm_all = cell(num_states_to_find, 1);
x_fdm_full = linspace(0, a, N+1)';

for i = 1:num_states_to_find
    psi_temp_internal = eigenvectors_fdm_internal_sorted(:, i);
    psi_full = [0; psi_temp_internal; 0];
    norm_const = sqrt(trapz(x_fdm_full, psi_full.^2));
    psi_normalized = psi_full / norm_const;
    
    idx_significant = find(abs(psi_normalized) > 1e-2*max(abs(psi_normalized)), 1, 'first');
    if ~isempty(idx_significant) && psi_normalized(idx_significant) < 0
        psi_normalized = -psi_normalized;
    end
    
    psi_fdm_all{i} = psi_normalized;
    fprintf('FDM E_%d = %.6f (Error: %.2e)\n', i, E_fdm_vals_phys(i), E_fdm_vals_phys(i) - E_an_vals(i));
end
fprintf('\n');


% --- Shooting Method ---
fprintf('--- Shooting Method ---\n');
ode_system = @(x, y, E_num_guess) [y(2); -E_num_guess * y(1)];
y0 = [0; 1];
x_span = [0, a];

% Define default ODE and fzero options for the main loop
options_ode_default = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
options_fzero_default = optimset('TolX', 1e-9);

% Modify objective_func to accept odeset options
objective_func_shooting = @(E_num_guess, current_ode_options) ...
    shoot_and_get_boundary_val(ode_system, x_span, y0, E_num_guess, current_ode_options);

E_shoot_vals_phys = zeros(num_states_to_find, 1);
psi_shoot_all = cell(num_states_to_find, 1);
x_shoot_all = cell(num_states_to_find, 1);

for n = 1:num_states_to_find
    E_num_an_n = (n*pi/a)^2;
    E_low_guess = ((n-0.5)*pi/a)^2;
    E_high_guess = ((n+0.5)*pi/a)^2;
    if n == 1, E_low_guess = max(0.1, ((n-0.5)*pi/a)^2); end

    % Pass the default ODE options to the objective function for fzero
    [E_num_found, fval, exitflag] = fzero(@(E_g) objective_func_shooting(E_g, options_ode_default), ...
                                          [E_low_guess, E_high_guess], options_fzero_default);

    if exitflag > 0
        E_shoot_vals_phys(n) = E_num_found * (hbar^2 / (2*m));
        [x_sol, y_sol] = ode45(@(x,y) ode_system(x,y,E_num_found), x_span, y0, options_ode_default);
        psi_temp = y_sol(:,1);
        norm_const_shoot = sqrt(trapz(x_sol, psi_temp.^2));
        psi_normalized_shoot = psi_temp / norm_const_shoot;
        
        idx_significant_shoot = find(abs(psi_normalized_shoot) > 1e-2*max(abs(psi_normalized_shoot)), 1, 'first');
        if ~isempty(idx_significant_shoot) && psi_normalized_shoot(idx_significant_shoot) < 0
            psi_normalized_shoot = -psi_normalized_shoot;
        end

        psi_shoot_all{n} = psi_normalized_shoot;
        x_shoot_all{n} = x_sol;
        fprintf('Shooting E_%d = %.6f (Error: %.2e), fzero fval: %.2e\n', n, E_shoot_vals_phys(n), E_shoot_vals_phys(n) - E_an_vals(n), fval);
    else
        fprintf('Shooting method failed to find E_%d (fzero exitflag: %d)\n', n, exitflag);
        E_shoot_vals_phys(n) = NaN;
    end
end
fprintf('\n');

% --- Comparison Table ---
fprintf('--- Summary of Energies (E_phys) ---\n');
fprintf(' n | E_Analytical | E_FDM (N=%d) | Err_FDM   | E_Shooting | Err_Shoot |\n', N);
fprintf('---|--------------|--------------|-----------|------------|-----------|\n');
for n = 1:num_states_to_find
    err_fdm = E_fdm_vals_phys(n) - E_an_vals(n);
    err_shoot = NaN; % Default if shooting failed for this state
    if n <= length(E_shoot_vals_phys) && ~isnan(E_shoot_vals_phys(n))
         err_shoot = E_shoot_vals_phys(n) - E_an_vals(n);
    end
    fprintf('%2d | %12.6f | %12.6f | %9.2e | %10.6f | %9.2e |\n', ...
            n, E_an_vals(n), E_fdm_vals_phys(n), err_fdm, E_shoot_vals_phys(n), err_shoot);
end
fprintf('\n');

% --- Plotting Wavefunctions ---
figure('Name', 'Wavefunctions for Infinite Potential Well', 'Position', [100, 100, 1000, 700]);
x_plot_fine = linspace(0, a, 500)';
for n = 1:min(num_states_to_find, 4)
    subplot(2, 2, n);
    hold on;
    plot(x_plot_fine, psi_analytical(x_plot_fine, n), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
    if n <= length(psi_fdm_all) && ~isempty(psi_fdm_all{n})
        plot(x_fdm_full, psi_fdm_all{n}, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('FDM (N=%d)', N));
    end
    if n <= length(psi_shoot_all) && ~isempty(x_shoot_all{n}) && ~any(isnan(psi_shoot_all{n}))
        plot(x_shoot_all{n}, psi_shoot_all{n}, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Shooting');
    end
    hold off;
    title(sprintf('\\psi_%d(x), E_{an} = %.3f', n, E_an_vals(n)));
    xlabel('x (position)'); ylabel('\psi(x)'); xlim([0, a]); grid on;
    legend('show', 'Location', 'best');
end
sgtitle(sprintf('Wavefunctions and Energies (hbar=%.1f, m=%.1f, a=%.1f)', hbar, m, a));

% --- Effect of N on FDM Accuracy (Example for ground state E_1) ---
fprintf('\n--- Effect of N on FDM Accuracy for E_1 ---\n');
N_values = [20, 50, 100, 200, 500];
E1_an = E_analytical(1);
fprintf('  N  | E1_FDM     | Error_E1_FDM | Rel. Error |\n');
fprintf('-----|------------|--------------|------------|\n');
for N_test = N_values
    h_test = a / N_test;
    main_diag_test = 2 * ones(N_test-1, 1); off_diag_test  = -1 * ones(N_test-2, 1);
    A_fdm_scaled_test = diag(main_diag_test) + diag(off_diag_test, 1) + diag(off_diag_test, -1);
    A_fdm_test = (1/h_test^2) * A_fdm_scaled_test;
    eigenvalues_fdm_raw_test = eig(A_fdm_test);
    E_fdm_vals_num_sorted_test = sort(eigenvalues_fdm_raw_test);
    E1_fdm_phys_test = E_fdm_vals_num_sorted_test(1) * (hbar^2 / (2*m));
    error_E1 = E1_fdm_phys_test - E1_an;
    rel_error_E1 = error_E1 / E1_an;
    fprintf('%4d | %10.6f | %12.4e | %10.2e |\n', N_test, E1_fdm_phys_test, error_E1, rel_error_E1);
end

% --- NEW: Effect of ODE Tolerance on Shooting Method Accuracy for E_1 ---
fprintf('\n--- Effect of ODE Tolerance on Shooting Method Accuracy for E_1 ---\n');
RelTol_values = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12];
E1_an_shoot_test = E_analytical(1); % Ground state analytical energy
fprintf(' RelTol  | E1_Shoot   | Error_E1_Shoot| Rel. Error | ODE Steps|\n');
fprintf('---------|------------|---------------|------------|----------|\n');

options_fzero_tight = optimset('TolX', 1e-12); % Keep fzero tolerance tight

for RelTol_test = RelTol_values
    current_options_ode_test = odeset('RelTol', RelTol_test, 'AbsTol', RelTol_test*1e-2); % Link AbsTol to RelTol

    % Define objective function for fzero for this specific test
    obj_func_shoot_test = @(E_num_g) shoot_and_get_boundary_val(ode_system, x_span, y0, E_num_g, current_options_ode_test);

    % Find E_num for ground state (n=1)
    E_num_an_1 = (1*pi/a)^2;
    E_low_guess_1 = max(0.1, ((1-0.5)*pi/a)^2);
    E_high_guess_1 = ((1+0.5)*pi/a)^2;

    [E1_num_found_test, ~, exitflag_test] = fzero(obj_func_shoot_test, ...
        [E_low_guess_1, E_high_guess_1], options_fzero_tight);

    if exitflag_test > 0
        E1_shoot_phys_test = E1_num_found_test * (hbar^2 / (2*m));
        error_E1_shoot = E1_shoot_phys_test - E1_an_shoot_test;
        rel_error_E1_shoot = error_E1_shoot / E1_an_shoot_test;

        % Get number of steps for final solution
        [x_sol_test, ~] = ode45(@(x,y) ode_system(x,y,E1_num_found_test), x_span, y0, current_options_ode_test);
        num_ode_steps = length(x_sol_test);

        fprintf('%8.1e | %10.6f | %13.4e | %10.2e | %8d |\n', ...
            RelTol_test, E1_shoot_phys_test, error_E1_shoot, rel_error_E1_shoot, num_ode_steps);
    else
        fprintf('%8.1e | Failed to converge                                  |\n', RelTol_test);
    end
end

% --- NEW: Visual Comparison of Relative Errors ---
figure('Name', 'Relative Errors of Numerical Methods', 'Position', [200, 200, 800, 600]);
rel_err_fdm = zeros(num_states_to_find, 1);
rel_err_shoot = zeros(num_states_to_find, 1);

for n = 1:num_states_to_find
    rel_err_fdm(n) = (E_fdm_vals_phys(n) - E_an_vals(n)) / E_an_vals(n);
    if n <= length(E_shoot_vals_phys) && ~isnan(E_shoot_vals_phys(n))
        rel_err_shoot(n) = (E_shoot_vals_phys(n) - E_an_vals(n)) / E_an_vals(n);
    else
        rel_err_shoot(n) = NaN; % Mark as NaN if shooting failed
    end
end

bar_data = [abs(rel_err_fdm), abs(rel_err_shoot)]; % Plot absolute values of relative errors for log scale
bar_labels = cell(num_states_to_find, 1);
for i=1:num_states_to_find; bar_labels{i} = sprintf('E_%d',i); end

b = bar(bar_data);
set(gca, 'XTickLabel', bar_labels);
set(gca, 'YScale', 'log'); % Use log scale for y-axis
ylabel('Absolute Relative Error |(E_{num} - E_{an}) / E_{an}|');
xlabel('Energy State');
title(sprintf('Comparison of Relative Errors (FDM N=%d, Shooting RelTol=%.0e)', N, options_ode_default.RelTol));
legend('FDM', 'Shooting', 'Location', 'best');
grid on;
ylim([1e-10, 1]); % Adjust y-axis limits if necessary

% Make sure all figures are drawn
drawnow;


% --- Helper function for Shooting Method's fzero ---
% Modified to accept odeset options
function boundary_val = shoot_and_get_boundary_val(ode_system_handle, x_span_f, y0_f, E_num_guess_f, current_options_ode)
    % Use the passed odeset options
    [~, y_sol_f] = ode45(@(x,y) ode_system_handle(x,y,E_num_guess_f), x_span_f, y0_f, current_options_ode);
    boundary_val = y_sol_f(end, 1); % Value of psi at x=a
end