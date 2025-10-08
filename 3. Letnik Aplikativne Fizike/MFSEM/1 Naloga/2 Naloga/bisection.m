function [root, iterations] = bisection(coeff, a, b, tolerance, max_iterations)
    % BISECTION finds a root of a polynomial in the interval [a, b]
    % using the bisection method.
    %
    % Inputs:
    %   coeff: Coefficient vector of the polynomial (e.g., [2, 3, -4, 5])
    %   a: Lower bound of the interval
    %   b: Upper bound of the interval
    %   tolerance: Desired accuracy for the root
    %   max_iterations: Maximum number of iterations allowed
    %
    % Outputs:
    %   root: Estimated root value
    %   iterations: Number of iterations performed

    % Check if the interval brackets a root
    f_a = polyval(coeff, a);
    f_b = polyval(coeff, b);
    if f_a * f_b >= 0
        error('Function does not change sign over [a, b]. Bisection cannot proceed.');
    end

    iterations = 0;
    while (b - a)/2 > tolerance && iterations < max_iterations
        c = (a + b)/2;  % Midpoint
        f_c = polyval(coeff, c);

        if f_c == 0  % Found exact root
            root = c;
            iterations = iterations + 1;
            return;
        elseif f_a * f_c < 0  % Root is in [a, c]
            b = c;
            f_b = f_c;
        else  % Root is in [c, b]
            a = c;
            f_a = f_c;
        end
        iterations = iterations + 1;
    end

    % Final midpoint approximation
    root = (a + b)/2;
end
