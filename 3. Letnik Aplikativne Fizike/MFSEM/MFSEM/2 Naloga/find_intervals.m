function intervals = find_intervals(coeff, x_min, x_max, step_size)
    % Finds intervals where f(a)*f(b) < 0 by evaluating the polynomial
    % at intervals of 'step_size' between x_min and x_max.
    %
    % Inputs:
    %   coeff: Polynomial coefficients
    %   x_min: Minimum x-value to check
    %   x_max: Maximum x-value to check
    %   step_size: Step between evaluation points
    %
    % Output:
    %   intervals: List of [a, b] intervals

    x = x_min:step_size:x_max;
    f = polyval(coeff, x);
    intervals = []

    for i = 1:length(x)-1
        a = x(i);
        b = x(i+1);
        if f(i)*f(i+1) < 0
            intervals = [intervals; a, b];
        end
    end
end
