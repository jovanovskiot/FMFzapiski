% Parameters for the ellipse and ellipsoid
a = 2;          % Semi-axis along x
b = 1;          % Semi-axis along y
c = 3;          % Semi-axis along z (for ellipsoid)
exact_area = pi * a * b;
exact_vol = (4/3) * pi * a * b * c;

% Monte Carlo parameters
N_values = round(logspace(3, 6, 20));  % N from 1e3 to 1e6
error_area_mc = zeros(size(N_values));
error_vol_mc = zeros(size(N_values));
time_mc_area = zeros(size(N_values));
time_mc_vol = zeros(size(N_values));

% Monte Carlo integration for ellipse area
for i = 1:length(N_values)
    N = N_values(i);
    x = -a + 2*a*rand(N, 1);
    y = -b + 2*b*rand(N, 1);
    tic;
    inside = (x.^2 / a^2 + y.^2 / b^2) <= 1;
    area_est = 4*a*b * sum(inside) / N;
    error_area_mc(i) = abs(area_est - exact_area);
    time_mc_area(i) = toc;
end

% Monte Carlo integration for ellipsoid volume
for i = 1:length(N_values)
    N = N_values(i);
    x = -a + 2*a*rand(N, 1);
    y = -b + 2*b*rand(N, 1);
    z = -c + 2*c*rand(N, 1);
    tic;
    inside = (x.^2 / a^2 + y.^2 / b^2 + z.^2 / c^2) <= 1;
    vol_est = 8*a*b*c * sum(inside) / N;
    error_vol_mc(i) = abs(vol_est - exact_vol);
    time_mc_vol(i) = toc;
end

% Uniform grid method for ellipse (comparison)
n_values = [10, 20, 50, 100, 200, 500, 1000];  % Grid points per axis
N_grid = n_values.^2;
error_area_grid = zeros(size(n_values));
time_grid = zeros(size(n_values));

for i = 1:length(n_values)
    n = n_values(i);
    x = linspace(-a, a, n);
    y = linspace(-b, b, n);
    [X, Y] = meshgrid(x, y);
    tic;
    inside = (X.^2 / a^2 + Y.^2 / b^2) <= 1;
    area_est_grid = 4*a*b * sum(inside(:)) / n^2;
    error_area_grid(i) = abs(area_est_grid - exact_area);
    time_grid(i) = toc;
end

% MC for grid-comparable N values
error_mc_gridN = zeros(size(N_grid));
time_mc_gridN = zeros(size(N_grid));
for i = 1:length(N_grid)
    N = N_grid(i);
    x = -a + 2*a*rand(N, 1);
    y = -b + 2*b*rand(N, 1);
    tic;
    inside = (x.^2 / a^2 + y.^2 / b^2) <= 1;
    area_est = 4*a*b * sum(inside) / N;
    error_mc_gridN(i) = abs(area_est - exact_area);
    time_mc_gridN(i) = toc;
end

% Plotting results
figure;
loglog(N_values, error_area_mc, 'b-o', 'DisplayName', 'MC Area');
hold on;
loglog(N_values, error_vol_mc, 'r-o', 'DisplayName', 'MC Volume');
loglog(N_grid, error_area_grid, 'g-s', 'DisplayName', 'Grid Area');
loglog(N_grid, error_mc_gridN, 'm--d', 'DisplayName', 'MC (Grid N)');
xlabel('Number of Points');
ylabel('Absolute Error');
legend('Location', 'northeast');
title('Error vs Number of Points');
grid on;

figure;
loglog(N_grid, time_grid, 'g-s', 'DisplayName', 'Grid Time');
hold on;
loglog(N_grid, time_mc_gridN, 'm--d', 'DisplayName', 'MC Time (Grid N)');
xlabel('Number of Points');
ylabel('Computation Time (s)');
legend('Location', 'northeast');
title('Computation Time Comparison');
grid on;