% Clear workspace and command window
clear; clc; close all;

% Define Matrix A1
A1 = [
    1000.     0.1      0.01     0.001    0.0001   0.00001  0.000001;
    0.1      1000.     0.1      0.01     0.001    0.0001   0.00001;
    0.01     0.1      1000.     0.1      0.01     0.001    0.0001;
    0.001    0.01     0.1      1000.     0.1      0.01     0.001;
    0.0001   0.001    0.01     0.1      1000.     0.1      0.01;
    0.00001  0.0001   0.001    0.01     0.1      1000.     0.1;
    0.000001 0.00001  0.0001   0.001    0.01     0.1      1000.
];

% Define Matrix A2
A2 = [
    1e6      0.1      0.01     0.001    0.0001   0.00001  0.000001;
    0.1      1e6      0.1      0.01     0.001    0.0001   0.00001;
    0.01     0.1      1e6      0.1      0.01     0.001    0.0001;
    0.001    0.01     0.1      1e6      0.1      0.01     0.001;
    0.0001   0.001    0.01     0.1      1e6      0.1      0.01;
    0.00001  0.0001   0.001    0.01     0.1      1e6      0.1;
    0.000001 0.00001  0.0001   0.001    0.01     0.1      1e6
];

% Set tolerance and max iterations for iterative methods
tol = 1e-9; % Increased precision slightly
max_iter = 1000000; % Reduced max_iter - 1e7 is too extreme for plotting/waiting
max_plot_iter = 50000; % Limit iterations shown in plots for clarity

fprintf('MATEMATIČNO-FIZIKALNI SEMINAR 2024/25\n');
fprintf('5. naloga: Lastne vrednosti simetričnega tenzorja\n\n');

%% Reference Solution using built-in eig

fprintf('--- Reference Solution using built-in eig ---\n');
n1 = size(A1, 1);
n2 = size(A2, 1);

[V_ref1, D_ref1] = eig(A1, 'vector'); % Use 'vector' option for direct eigenvalues
[lambda_ref1_sorted, sort_idx1] = sort(D_ref1, 'descend');
[V_ref1_full, ~] = eig(A1); % Get eigenvectors separately if needed later
V_ref1_sorted = V_ref1_full(:, sort_idx1);

fprintf('Matrix A1:\n');
fprintf('Eigenvalues (sorted):\n');
fprintf('%15.8f\n', lambda_ref1_sorted);

[V_ref2, D_ref2] = eig(A2, 'vector');
[lambda_ref2_sorted, sort_idx2] = sort(D_ref2, 'descend');
[V_ref2_full, ~] = eig(A2);
V_ref2_sorted = V_ref2_full(:, sort_idx2);

fprintf('\nMatrix A2:\n');
fprintf('Eigenvalues (sorted):\n');
fprintf('%15.8f\n', lambda_ref2_sorted);
fprintf('----------------------------------------------\n\n');

%% MODIFIED Iterative QR Function with History Tracking

function [lambda, V, iter, diag_history, off_diag_history] = iterative_qr_tracked(A, max_iter, tol)
    [m, n] = size(A);
    if m ~= n
        error('Matrix must be square.');
    end
    V = eye(n); % Accumulate eigenvectors
    Ak = A;

    % History tracking (pre-allocate a reasonable amount initially)
    prealloc_size = min(max_iter, 5000); % Adjust preallocation size
    diag_history = zeros(n, prealloc_size);
    off_diag_history = zeros(1, prealloc_size);
    history_count = 0;

    converged = false;
    for iter = 1:max_iter
        [Q, R] = qr(Ak);
        Ak = R * Q;
        V = V * Q; % Accumulate orthogonal transformations

        % --- Store History (dynamically grow if needed) ---
        history_count = history_count + 1;
        if history_count > size(diag_history, 2)
             % Double the size if preallocation runs out
             diag_history = [diag_history, zeros(n, size(diag_history, 2))];
             off_diag_history = [off_diag_history, zeros(1, size(off_diag_history, 2))];
        end
        diag_history(:, history_count) = diag(Ak);

        % Calculate max absolute off-diagonal element for convergence check
        if n > 1
             off_diag_indices = ~eye(n); % Logical matrix, true for off-diagonals
             max_off_diag = max(abs(Ak(off_diag_indices)));
             off_diag_history(history_count) = max_off_diag; % Store the max value
        else
            max_off_diag = 0; % Scalar matrix is already diagonal
            off_diag_history(history_count) = 0;
        end

        % --- CONVERGENCE CHECK ---
        if max_off_diag < tol
            converged = true;
            break; % Exit loop if converged
        end
    end

    iter = history_count; % Actual number of iterations performed

    if ~converged
        warning('Iterative QR did not converge within %d iterations based on max off-diagonal tolerance (%.2e). Max off-diag reached: %.2e', max_iter, tol, max_off_diag);
    end

    % Trim history to actual size
    diag_history = diag_history(:, 1:iter);
    off_diag_history = off_diag_history(1:iter);

    lambda = diag(Ak);
end

%% Run QR Iterations and Collect Data

results = struct(); % Structure to store results

fprintf('--- Iterative QR Method Calculations ---\n');

% --- Matrix A1 ---
fprintf('Running Matrix A1 (Full Matrix)...\n');
tic;
[lambda_qr1, V_qr1, iter_qr1, diag_hist_qr1, off_diag_hist_qr1] = iterative_qr_tracked(A1, max_iter, tol);
time_qr1 = toc;
[lambda_qr1_sorted, sort_idx_qr1] = sort(lambda_qr1, 'descend');
V_qr1_sorted = V_qr1(:, sort_idx_qr1);
diff_eigvals1 = norm(lambda_ref1_sorted - lambda_qr1_sorted); % Norm of difference vector
results.A1_Full = struct('lambda', lambda_qr1_sorted, 'V', V_qr1_sorted, 'iter', iter_qr1, 'time', time_qr1, 'error_norm', diff_eigvals1, 'diag_hist', diag_hist_qr1, 'off_diag_hist', off_diag_hist_qr1);
fprintf('Done. Iterations: %d, Time: %.4f s, Error Norm: %.2e\n', iter_qr1, time_qr1, diff_eigvals1);

% --- Matrix A1 - Tridiagonal ---
fprintf('Running Matrix A1 (Tridiagonalized First)...\n');
tic;
[Q_hess1, T1] = hess(A1); % For symmetric matrices, hess produces tridiagonal Q'*A*Q = T
[lambda_qr_t1, V_qr_t1_local, iter_qr_t1, diag_hist_t1, off_diag_hist_t1] = iterative_qr_tracked(T1, max_iter, tol);
V_qr_t1 = Q_hess1 * V_qr_t1_local; % Transform eigenvectors back: V_A = Q_hess * V_T
time_qr_t1 = toc;
[lambda_qr_t1_sorted, sort_idx_qr_t1] = sort(lambda_qr_t1, 'descend');
V_qr_t1_sorted = V_qr_t1(:, sort_idx_qr_t1);
diff_eigvals_t1 = norm(lambda_ref1_sorted - lambda_qr_t1_sorted);
results.A1_Tri = struct('lambda', lambda_qr_t1_sorted, 'V', V_qr_t1_sorted, 'iter', iter_qr_t1, 'time', time_qr_t1, 'error_norm', diff_eigvals_t1, 'diag_hist', diag_hist_t1, 'off_diag_hist', off_diag_hist_t1);
fprintf('Done. Iterations: %d, Time: %.4f s, Error Norm: %.2e\n', iter_qr_t1, time_qr_t1, diff_eigvals_t1);


% --- Matrix A2 ---
fprintf('Running Matrix A2 (Full Matrix)...\n');
tic;
[lambda_qr2, V_qr2, iter_qr2, diag_hist_qr2, off_diag_hist_qr2] = iterative_qr_tracked(A2, max_iter, tol);
time_qr2 = toc;
[lambda_qr2_sorted, sort_idx_qr2] = sort(lambda_qr2, 'descend');
V_qr2_sorted = V_qr2(:, sort_idx_qr2);
diff_eigvals2 = norm(lambda_ref2_sorted - lambda_qr2_sorted);
results.A2_Full = struct('lambda', lambda_qr2_sorted, 'V', V_qr2_sorted, 'iter', iter_qr2, 'time', time_qr2, 'error_norm', diff_eigvals2, 'diag_hist', diag_hist_qr2, 'off_diag_hist', off_diag_hist_qr2);
fprintf('Done. Iterations: %d, Time: %.4f s, Error Norm: %.2e\n', iter_qr2, time_qr2, diff_eigvals2);

% --- Matrix A2 - Tridiagonal ---
fprintf('Running Matrix A2 (Tridiagonalized First)...\n');
tic;
[Q_hess2, T2] = hess(A2); % Tridiagonal form
[lambda_qr_t2, V_qr_t2_local, iter_qr_t2, diag_hist_t2, off_diag_hist_t2] = iterative_qr_tracked(T2, max_iter, tol);
V_qr_t2 = Q_hess2 * V_qr_t2_local; % Transform eigenvectors back
time_qr_t2 = toc;
[lambda_qr_t2_sorted, sort_idx_qr_t2] = sort(lambda_qr_t2, 'descend');
V_qr_t2_sorted = V_qr_t2(:, sort_idx_qr_t2);
diff_eigvals_t2 = norm(lambda_ref2_sorted - lambda_qr_t2_sorted);
results.A2_Tri = struct('lambda', lambda_qr_t2_sorted, 'V', V_qr_t2_sorted, 'iter', iter_qr_t2, 'time', time_qr_t2, 'error_norm', diff_eigvals_t2, 'diag_hist', diag_hist_t2, 'off_diag_hist', off_diag_hist_t2);
fprintf('Done. Iterations: %d, Time: %.4f s, Error Norm: %.2e\n', iter_qr_t2, time_qr_t2, diff_eigvals_t2);

fprintf('----------------------------------------------\n\n');

%% Visualizations

disp('Generating Visualizations...');

% ----- Visualization 1: Eigenvalue Convergence over Iterations -----
figure('Name', 'Eigenvalue Convergence during QR Iteration');
sgtitle('Eigenvalue Convergence during QR Iteration');

% A1 Full
subplot(2, 2, 1);
iters_to_plot = min(results.A1_Full.iter, max_plot_iter);
plot(1:iters_to_plot, results.A1_Full.diag_hist(:, 1:iters_to_plot)');
title(sprintf('A1 Full (First %d of %d Iter)', iters_to_plot, results.A1_Full.iter));
xlabel('Iteration'); ylabel('Eigenvalue Estimate'); grid on;
legend(arrayfun(@(x) sprintf('\\lambda_%d',x), 1:n1, 'UniformOutput', false), 'Location', 'best');

% A1 Tridiagonal
subplot(2, 2, 2);
iters_to_plot = min(results.A1_Tri.iter, max_plot_iter);
plot(1:iters_to_plot, results.A1_Tri.diag_hist(:, 1:iters_to_plot)');
title(sprintf('A1 Tridiagonal (First %d of %d Iter)', iters_to_plot, results.A1_Tri.iter));
xlabel('Iteration'); ylabel('Eigenvalue Estimate'); grid on;
legend(arrayfun(@(x) sprintf('\\lambda_%d',x), 1:n1, 'UniformOutput', false), 'Location', 'best');

% A2 Full
subplot(2, 2, 3);
iters_to_plot = min(results.A2_Full.iter, max_plot_iter);
plot(1:iters_to_plot, results.A2_Full.diag_hist(:, 1:iters_to_plot)');
title(sprintf('A2 Full (First %d of %d Iter)', iters_to_plot, results.A2_Full.iter));
xlabel('Iteration'); ylabel('Eigenvalue Estimate'); grid on;
legend(arrayfun(@(x) sprintf('\\lambda_%d',x), 1:n2, 'UniformOutput', false), 'Location', 'best');

% A2 Tridiagonal
subplot(2, 2, 4);
iters_to_plot = min(results.A2_Tri.iter, max_plot_iter);
plot(1:iters_to_plot, results.A2_Tri.diag_hist(:, 1:iters_to_plot)');
title(sprintf('A2 Tridiagonal (First %d of %d Iter)', iters_to_plot, results.A2_Tri.iter));
xlabel('Iteration'); ylabel('Eigenvalue Estimate'); grid on;
legend(arrayfun(@(x) sprintf('\\lambda_%d',x), 1:n2, 'UniformOutput', false), 'Location', 'best');

% ----- Visualization 2: Convergence Error (Max Off-Diagonal) -----
figure('Name', 'QR Convergence Error');
sgtitle('QR Convergence Error (Max |Off-Diagonal Element|)');

subplot(2, 2, 1);
iters_to_plot = min(results.A1_Full.iter, max_plot_iter);
semilogy(1:iters_to_plot, results.A1_Full.off_diag_hist(1:iters_to_plot));
title(sprintf('A1 Full (First %d of %d Iter)', iters_to_plot, results.A1_Full.iter));
xlabel('Iteration'); ylabel('Error (log scale)'); grid on; ylim([tol*0.1, max(results.A1_Full.off_diag_hist)*10]); % Adjust ylim

subplot(2, 2, 2);
iters_to_plot = min(results.A1_Tri.iter, max_plot_iter);
semilogy(1:iters_to_plot, results.A1_Tri.off_diag_hist(1:iters_to_plot));
title(sprintf('A1 Tridiagonal (First %d of %d Iter)', iters_to_plot, results.A1_Tri.iter));
xlabel('Iteration'); ylabel('Error (log scale)'); grid on; ylim([tol*0.1, max(results.A1_Tri.off_diag_hist)*10]);

subplot(2, 2, 3);
iters_to_plot = min(results.A2_Full.iter, max_plot_iter);
semilogy(1:iters_to_plot, results.A2_Full.off_diag_hist(1:iters_to_plot));
title(sprintf('A2 Full (First %d of %d Iter)', iters_to_plot, results.A2_Full.iter));
xlabel('Iteration'); ylabel('Error (log scale)'); grid on; ylim([tol*0.1, max(results.A2_Full.off_diag_hist)*10]);

subplot(2, 2, 4);
iters_to_plot = min(results.A2_Tri.iter, max_plot_iter);
semilogy(1:iters_to_plot, results.A2_Tri.off_diag_hist(1:iters_to_plot));
title(sprintf('A2 Tridiagonal (First %d of %d Iter)', iters_to_plot, results.A2_Tri.iter));
xlabel('Iteration'); ylabel('Error (log scale)'); grid on; ylim([tol*0.1, max(results.A2_Tri.off_diag_hist)*10]);

% ----- Visualization 3: Final Absolute Error per Eigenvalue -----
figure('Name', 'Final Absolute Error per Eigenvalue');
sgtitle('Final Absolute Error per Eigenvalue (|Ref - Calculated|)');

subplot(2, 2, 1);
final_err_qr1 = abs(lambda_ref1_sorted - results.A1_Full.lambda);
stem(1:n1, final_err_qr1); % Use stem plot
set(gca, 'YScale', 'log'); % Use log scale for y-axis
title(sprintf('A1 Full (Norm: %.2e)', results.A1_Full.error_norm));
xlabel('Eigenvalue Index (Sorted)'); ylabel('Abs Error (log)'); grid on; ylim([min(final_err_qr1(final_err_qr1>0))/10, max(final_err_qr1)*10]);

subplot(2, 2, 2);
final_err_t1 = abs(lambda_ref1_sorted - results.A1_Tri.lambda);
stem(1:n1, final_err_t1);
set(gca, 'YScale', 'log');
title(sprintf('A1 Tridiagonal (Norm: %.2e)', results.A1_Tri.error_norm));
xlabel('Eigenvalue Index (Sorted)'); ylabel('Abs Error (log)'); grid on; ylim([min(final_err_t1(final_err_t1>0))/10, max(final_err_t1)*10]);

subplot(2, 2, 3);
final_err_qr2 = abs(lambda_ref2_sorted - results.A2_Full.lambda);
stem(1:n2, final_err_qr2);
set(gca, 'YScale', 'log');
title(sprintf('A2 Full (Norm: %.2e)', results.A2_Full.error_norm));
xlabel('Eigenvalue Index (Sorted)'); ylabel('Abs Error (log)'); grid on; ylim([min(final_err_qr2(final_err_qr2>0))/10, max(final_err_qr2)*10]);

subplot(2, 2, 4);
final_err_t2 = abs(lambda_ref2_sorted - results.A2_Tri.lambda);
stem(1:n2, final_err_t2);
set(gca, 'YScale', 'log');
title(sprintf('A2 Tridiagonal (Norm: %.2e)', results.A2_Tri.error_norm));
xlabel('Eigenvalue Index (Sorted)'); ylabel('Abs Error (log)'); grid on; ylim([min(final_err_t2(final_err_t2>0))/10, max(final_err_t2)*10]);

disp('Visualizations generated.');
%% LaTeX Table Output

fprintf('\n--- LaTeX Table Results ---\n\n');

% Header
fprintf('\\begin{tabular}{|c|c|c|c|c|}\n');
fprintf('\\hline\n');
fprintf('Matrix & Method & Iterations & Time (s) & Norm Error \\\\ \\hline\n');

% Data Rows
fprintf('A1 & Full QR & %d & %.4f & %.2e \\\\ \\hline\n', ...
    results.A1_Full.iter, results.A1_Full.time, results.A1_Full.error_norm);
fprintf('A1 & Tridiagonal QR & %d & %.4f & %.2e \\\\ \\hline\n', ...
    results.A1_Tri.iter, results.A1_Tri.time, results.A1_Tri.error_norm);
fprintf('A2 & Full QR & %d & %.4f & %.2e \\\\ \\hline\n', ...
    results.A2_Full.iter, results.A2_Full.time, results.A2_Full.error_norm);
fprintf('A2 & Tridiagonal QR & %d & %.4f & %.2e \\\\ \\hline\n', ...
    results.A2_Tri.iter, results.A2_Tri.time, results.A2_Tri.error_norm);

% Footer
fprintf('\\end{tabular}\n');

fprintf('\n--- End of Script ---\n');