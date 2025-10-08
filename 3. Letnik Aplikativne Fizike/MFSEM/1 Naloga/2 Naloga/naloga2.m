format long

% polinomi
p1 = [2, 3, -4, 5];   % P1(x) = 2x³ +3x² -4x +5
p2 = [-1, 0.5, 7, -3]; % P2(x) = -x³ +0.5x² +7x -3
p3 = [3, 0, -5, 2];    % P3(x) = 3x³ -5x +2

% nekaj parametrov
tol = 1e-8;
max_iter = 100;
search_range = [-4, 3];
step = 0.51;


results = cell(0, 6); % Poly, Interval, Method, Root, Iterations, Analytic

polynomials = {p1, p2, p3};
poly_names = {'P1', 'P2', 'P3'};

for p_idx = 1:length(polynomials)
    p = polynomials{p_idx};
    poly_name = poly_names{p_idx};

    % sprememba predznaka
    intervals = find_intervals(p, search_range(1), search_range(2), step);

    % vgrajena
    analytic_roots = roots(p);
    analytic_roots = analytic_roots(abs(imag(analytic_roots)) < 1e-8);
    analytic_roots = real(analytic_roots);

    for i = 1:size(intervals, 1)
        a = intervals(i, 1);
        b = intervals(i, 2);

        % vgrajena
        analytic_in_interval = analytic_roots(analytic_roots >= a & analytic_roots <= b);
        if isempty(analytic_in_interval)
            analytic_root = NaN;
        else
            analytic_root = analytic_in_interval(1); % Take first root in interval
        end

        % bisekcija
        [bis_root, bis_iters] = bisection(p, a, b, tol, max_iter);

        % regula
        [rf_root, rf_iters] = regula_falsi(p, a, b, tol, max_iter);

        results(end+1, :) = {poly_name, [a, b], 'Bisekcija', bis_root, bis_iters, analytic_root};
        results(end+1, :) = {poly_name, [a, b], 'RegulaFalsi', rf_root, rf_iters, analytic_root};
        results(end+1, :) = {poly_name, [a, b], 'Roots', analytic_root, NaN, NaN};
    end
end

% tabela
results_table = cell2table(results, 'VariableNames', ...
    {'Polinom', 'Interval', 'Metoda', 'Nicla', 'Iteracije', 'Vgrajena'});
disp(results_table);

