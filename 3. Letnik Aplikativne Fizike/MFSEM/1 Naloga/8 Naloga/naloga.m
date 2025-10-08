% MATEMATIČNO-FIZIKALNI SEMINAR 2024/25
% 8. naloga: Robni problem lastnih vrednosti - Glavna Naloga

clear;
clc;
close all;

% --- Parameters ---
hbar = 1; % Reduced Planck's constant (can be set to 1 for simplicity)
m = 0.5;  % Mass of the particle (setting 2m=1 simplifies Schrödinger eq. to -d^2psi/dx^2 = E psi)
% Note: If hbar=1 and 2m=1 (i.e. m=0.5), then hbar^2/(2*m) = 1.
% The Schrödinger equation inside the well (V=0) becomes:
% - (hbar^2 / (2*m)) * d^2psi/dx^2 = E_phys * psi
% If hbar^2/(2*m) = 1, then -d^2psi/dx^2 = E_phys * psi.
% Let E_num be the eigenvalue from -d^2psi/dx^2 = E_num * psi. Then E_phys = E_num.

a = 1;    % Width of the infinite potential well (from x=0 to x=a)

N = 100;  % Number of discretiza_valtion intervals for FDM (N+1 points total)
          % The problem statement mentions "N točk (xi = 0 + a/N · i)"
          % If i goes from 0 to N, this is N+1 points.
          % If we solve for N-1 internal points, then N is the number of intervals.
          % Let's use N as the number of intervals, so N-1 internal points.

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
x_internal = linspace(h, a-h, N-1)'; % Internal grid points (N-1 of them)
% For N-1 internal points, the matrix A will be (N-1)x(N-1)

% Construct the matrix A for -d^2psi/dx^2 = E_num * psi
% The discretized form is (-psi_{i-1} + 2*psi_i - psi_{i+1}) / h^2 = E_num * psi_i
% So, A_ij * psi_j = E_num * psi_i, where A is (1/h^2) * tridiagonal(-1, 2, -1)
main_diag = 2 * ones(N-1, 1);
off_diag  = -1 * ones(N-2, 1);
A_fdm_scaled = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
A_fdm = (1/h^2) * A_fdm_scaled;

% Solve the eigenvalue problem A * psi = E_num * psi
[eigenvectors_fdm_internal, eigenvalues_fdm_raw] = eig(A_fdm);
E_fdm_vals_num = diag(eigenvalues_fdm_raw); % These are E_num

% Sort eigenvalues and corresponding eigenvectors
[E_fdm_vals_num_sorted, sort_idx_fdm] = sort(E_fdm_vals_num);
eigenvectors_fdm_internal_sorted = eigenvectors_fdm_internal(:, sort_idx_fdm);

E_fdm_vals_phys = E_fdm_vals_num_sorted * (hbar^2 / (2*m)); % Convert E_num to E_phys
                                                         % If hbar^2/(2m)=1, then E_phys = E_num

psi_fdm_all = cell(num_states_to_find, 1);
x_fdm_full = linspace(0, a, N+1)'; % Full grid including boundaries

for i = 1:num_states_to_find
    % Add boundary conditions (psi=0 at x=0 and x=a)
    psi_temp_internal = eigenvectors_fdm_internal_sorted(:, i);
    psi_full = [0; psi_temp_internal; 0];

    % Normalize: integral(abs(psi)^2 dx) = 1
    % Using trapezoidal rule for integration
    norm_const = sqrt(trapz(x_fdm_full, psi_full.^2));
    psi_normalized = psi_full / norm_const;

    % Ensure the first lobe is positive for consistent comparison with sin(n*pi*x/a)
    % Find the first non-zero segment to check sign
    first_lobe_idx = find(abs(psi_normalized) > 1e-6, 1, 'first');
    if ~isempty(first_lobe_idx) && psi_normalized(first_lobe_idx+1) < 0 % Check point after zero
         % Or check midpoint of first expected lobe
        check_idx = round((N+1)/(2*i)); % Approximate index for first lobe peak
        if check_idx > 1 && check_idx < N+1 && psi_normalized(check_idx) < 0
             psi_normalized = -psi_normalized;
        elseif psi_normalized(2) < 0 && i==1 % Simpler check for ground state
             psi_normalized = -psi_normalized;
        end
    end
    % A more robust check: if the integral of the first half is negative
    mid_point_index = floor(length(psi_normalized)/2);
    if trapz(x_fdm_full(1:mid_point_index), psi_normalized(1:mid_point_index)) < 0 && i==1 % More robust for n=1
         psi_normalized = -psi_normalized;
    end
    % Check if the first significant value is negative
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

% Define the system of ODEs for d^2psi/dx^2 = -E_num * psi
% y(1) = psi, y(2) = dpsi/dx
% dy(1)/dx = y(2)
% dy(2)/dx = -E_num * y(1)
% E_num corresponds to (2*m*E_phys)/hbar^2
ode_system = @(x, y, E_num_guess) [y(2); -E_num_guess * y(1)];

% Initial conditions: psi(0) = 0, dpsi/dx(0) = 1 (arbitrary non-zero slope)
y0 = [0; 1];
x_span = [0, a]; % Integration interval

% Objective function: we want psi(a) = 0
% This function will take E_num_guess and return psi(a)
objective_func = @(E_num_guess) ...
    shoot_and_get_boundary_val(ode_system, x_span, y0, E_num_guess);

E_shoot_vals_phys = zeros(num_states_to_find, 1);
psi_shoot_all = cell(num_states_to_find, 1);
x_shoot_all = cell(num_states_to_find, 1);
options_ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10); % ODE solver options
options_fzero = optimset('TolX', 1e-9); % fzero options

for n = 1:num_states_to_find
    % Estimate initial search interval for E_num using analytical E_phys
    % E_num_analytical = (n*pi/a)^2
    E_num_an_n = (n*pi/a)^2;

    % Bracket the root for fzero
    % For n=1, E_num should be around (pi/a)^2.
    % For n=2, E_num should be around (2*pi/a)^2.
    % Search range: [(n-0.5)^2*(pi/a)^2, (n+0.5)^2*(pi/a)^2]
    E_low_guess = ((n-0.5)*pi/a)^2;
    E_high_guess = ((n+0.5)*pi/a)^2;
    if n == 1 % Ensure lower bound is > 0
        E_low_guess = max(0.1, ((n-0.5)*pi/a)^2);
    end

    % Find the E_num that makes psi(a) = 0
    [E_num_found, fval, exitflag] = fzero(objective_func, [E_low_guess, E_high_guess], options_fzero);

    if exitflag > 0
        E_shoot_vals_phys(n) = E_num_found * (hbar^2 / (2*m)); % Convert E_num to E_phys

        % Solve ODE one more time with the found E_num to get the wavefunction
        [x_sol, y_sol] = ode45(@(x,y) ode_system(x,y,E_num_found), x_span, y0, options_ode);
        psi_temp = y_sol(:,1);

        % Normalize
        norm_const_shoot = sqrt(trapz(x_sol, psi_temp.^2));
        psi_normalized_shoot = psi_temp / norm_const_shoot;

        % Ensure positive first lobe
        if psi_normalized_shoot(find(abs(psi_normalized_shoot) > 1e-6, 1, 'first')+1) < 0
             psi_normalized_shoot = -psi_normalized_shoot;
        end
        % Check if the first significant value is negative
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
        psi_shoot_all{n} = NaN;
        x_shoot_all{n} = NaN;
    end
end
fprintf('\n');

% --- Comparison Table ---
fprintf('--- Summary of Energies (E_phys) ---\n');
fprintf(' n | E_Analytical | E_FDM (N=%d) | Err_FDM   | E_Shooting | Err_Shoot |\n', N);
fprintf('---|--------------|--------------|-----------|------------|-----------|\n');
for n = 1:num_states_to_find
    err_fdm = E_fdm_vals_phys(n) - E_an_vals(n);
    err_shoot = E_shoot_vals_phys(n) - E_an_vals(n);
    fprintf('%2d | %12.6f | %12.6f | %9.2e | %10.6f | %9.2e |\n', ...
            n, E_an_vals(n), E_fdm_vals_phys(n), err_fdm, E_shoot_vals_phys(n), err_shoot);
end
fprintf('\n');

%% --- Plotting Wavefunctions ---
figure('Name', 'Wavefunctions for Infinite Potential Well', 'Position', [100, 100, 1000, 700]);
x_plot_fine = linspace(0, a, 500)'; % Fine grid for analytical solution

for n = 1:min(num_states_to_find, 4) % Plot up to 4 states
    subplot(2, 2, n);
    hold on;

    % Analytical
    plot(x_plot_fine, psi_analytical(x_plot_fine, n), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Analytical');

    % FDM
    if n <= length(psi_fdm_all) && ~isempty(psi_fdm_all{n})
        plot(x_fdm_full, psi_fdm_all{n}, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('FDM (N=%d)', N));
    end

    % Shooting
    if n <= length(psi_shoot_all) && ~isempty(x_shoot_all{n}) && ~any(isnan(psi_shoot_all{n}))
        plot(x_shoot_all{n}, psi_shoot_all{n}, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Shooting');
    end

    hold off;
    title(sprintf('\\psi_%d(x), E_{an} = %.3f', n, E_an_vals(n)));
    xlabel('x (position)');
    ylabel('\psi(x)');
    xlim([0, a]);
    grid on;
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
    main_diag_test = 2 * ones(N_test-1, 1);
    off_diag_test  = -1 * ones(N_test-2, 1);
    A_fdm_scaled_test = diag(main_diag_test) + diag(off_diag_test, 1) + diag(off_diag_test, -1);
    A_fdm_test = (1/h_test^2) * A_fdm_scaled_test;
    eigenvalues_fdm_raw_test = eig(A_fdm_test);
    E_fdm_vals_num_sorted_test = sort(eigenvalues_fdm_raw_test);
    E1_fdm_phys_test = E_fdm_vals_num_sorted_test(1) * (hbar^2 / (2*m));
    error_E1 = E1_fdm_phys_test - E1_an;
    rel_error_E1 = error_E1 / E1_an;
    fprintf('%4d | %10.6f | %12.4e | %10.2e |\n', N_test, E1_fdm_phys_test, error_E1, rel_error_E1);
end
% Note: The error in FDM eigenvalues typically scales as O(h^2) for this scheme.


% --- Helper function for Shooting Method's fzero ---
function boundary_val = shoot_and_get_boundary_val(ode_system_handle, x_span_f, y0_f, E_num_guess_f)
    options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9); % Slightly less stringent for fzero calls
    [~, y_sol_f] = ode45(@(x,y) ode_system_handle(x,y,E_num_guess_f), x_span_f, y0_f, options);
    boundary_val = y_sol_f(end, 1); % Value of psi at x=a
end