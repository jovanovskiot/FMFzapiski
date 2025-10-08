function y = erfc_asymptotic(x)
% Asymptotic expansion of the complementary error function erfc(x) for large x
% y = erfc_asymptotic(x) computes the approximation using an adaptive number of terms

if ~isscalar(x)
    error('x must be a scalar.');
end

if real(x) < 0
    % Use the identity erfc(-x) = 2 - erfc(x)
    y = 2 - erfc_asymptotic(-x);
    return;
end

if x == 0
    y = 1;
    return;
end

prefactor = exp(-x^2) / (x * sqrt(pi));
if isinf(prefactor)
    y = 0;
    return;
end

S = 1;          % Initialize the series sum with the first term (k=0)
term = 1;        % Initial term (k=0)
k = 1;           % Term index
max_terms = 100; % Maximum number of terms to prevent infinite loops
tolerance = 1e-15; % Tolerance for term magnitude
prev_term_mag = abs(term); % Track the magnitude of the previous term

while k <= max_terms
    % Compute the next term in the series
    term = term * (-1) * (2*k - 1) / (2 * x^2);
    term_mag = abs(term);

    % Check if the term magnitude is below the tolerance
    if term_mag < tolerance
        S = S + term;
        break;
    end

    % Check if terms start increasing (divergence)
    if term_mag > prev_term_mag
        break;
    end

    % Add the term to the sum and update previous term magnitude
    S = S + term;
    prev_term_mag = term_mag;
    k = k + 1;
end

y = prefactor * S;

end
