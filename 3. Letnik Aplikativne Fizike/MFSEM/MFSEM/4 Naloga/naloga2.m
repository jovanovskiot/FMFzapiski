% --- Clear Workspace and Command Window ---
clearvars;
close all;
clc;

% --- Define Data ---
% Ensure data are column vectors for fitting functions
x_data = [
    1.00000000; 1.50000000; 2.00000000; 2.50000000; 3.00000000; 3.50000000;
    4.00000000; 4.50000000; 5.00000000; 5.50000000; 6.00000000; 6.50000000;
    7.00000000; 7.50000000; 8.00000000; 8.50000000; 9.00000000; 9.50000000;
    10.00000000; 10.50000000; 11.00000000
];

y_data = [
    0.31700705; 0.43791106; 0.56528271; 0.56102378; 0.63664784; 0.65121353;
    0.63487502; 0.64501481; 0.60942923; 0.62411336; 0.61455575; 0.57226264;
    0.54291294; 0.50329224; 0.50314769; 0.46050043; 0.42461463; 0.40771586;
    0.41605889; 0.36732963; 0.33085992
];

sigma_data = [
    0.01548814; 0.01715189; 0.01602763; 0.01544883; 0.01423655; 0.01645894;
    0.01437587; 0.01891773; 0.01963663; 0.01383442; 0.01791725; 0.01528895;
    0.01568045; 0.01925597; 0.01071036; 0.01087129; 0.01020218; 0.01832620;
    0.01778157; 0.01870012; 0.01978618
];

% --- 1. Define the Model and Fit Options ---

% Define the model function type using fittype
% 'a*x*exp(b*x)' is the formula
% 'independent', 'x' specifies the independent variable
% 'coefficients', {'a', 'b'} specifies the parameters to fit
model_fittype = fittype('a*x*exp(b*x)', 'independent', 'x', 'coefficients', {'a', 'b'});

% Define weights (inverse variance: 1/sigma^2)
weights = 1./(sigma_data.^2);

% Provide initial guesses for parameters a and b (helps convergence)
initial_guesses = [0.5, -0.1]; % [guess_for_a, guess_for_b]

% Set up fit options, including starting points and weights
fit_options = fitoptions(model_fittype);
fit_options.StartPoint = initial_guesses;
fit_options.Weights = weights;

% --- 2. Perform the Curve Fit ---
% Use the 'fit' function
% fit(x, y, model, options)
[fit_result, gof] = fit(x_data, y_data, model_fittype, fit_options);

% fit_result is a 'cfit' object containing the fit parameters
% gof is a structure containing goodness-of-fit statistics

% --- 3. Extract Results ---
a_fit = fit_result.a;
b_fit = fit_result.b;

% Get confidence intervals (95% by default) for parameters
conf_intervals = confint(fit_result); % Returns a 2x2 matrix [lower; upper]
a_lower = conf_intervals(1, 1);
a_upper = conf_intervals(2, 1);
b_lower = conf_intervals(1, 2);
b_upper = conf_intervals(2, 2);

% Estimate standard error from 95% confidence interval (approx.)
% Assuming normal distribution, interval is roughly +/- 1.96 * SE
a_err = (a_upper - a_lower) / (2 * 1.96);
b_err = (b_upper - b_lower) / (2 * 1.96);
% Note: More precise uncertainty calculation might involve the Jacobian
% from the fit, but `confint` is the standard way with `fit`.

fprintf('--- Fit Results ---\n');
fprintf('Optimal parameter a = %.6f +/- %.6f\n', a_fit, a_err);
fprintf('Optimal parameter b = %.6f +/- %.6f\n', b_fit, b_err);
fprintf(' (95%% CI for a: [%.6f, %.6f])\n', a_lower, a_upper);
fprintf(' (95%% CI for b: [%.6f, %.6f])\n', b_lower, b_upper);


% --- 4. Calculate Chi-squared ---
% Calculate fitted y values
y_fit = feval(fit_result, x_data); % Or simply: y_fit = fit_result(x_data);

% Calculate residuals
residuals = y_data - y_fit;

% Calculate Chi-squared statistic
chi_squared = sum((residuals ./ sigma_data).^2);

% Degrees of freedom = number of data points - number of parameters
num_params = numcoeffs(fit_result); % Gets the number of fitted coefficients (2 here)
dof = length(x_data) - num_params;

% Calculate reduced Chi-squared
reduced_chi_squared = chi_squared / dof;

fprintf('\n--- Goodness of Fit ---\n');
fprintf('Chi-squared = %.4f\n', chi_squared);
fprintf('Degrees of Freedom (dof) = %d\n', dof);
fprintf('Reduced Chi-squared = %.4f\n', reduced_chi_squared);
% Also display R^2 from the 'gof' structure
fprintf('R-squared = %.4f\n', gof.rsquare);


% --- 5. Visualize the Fit ---
figure('Name', 'MATLAB Curve Fit'); % Create a new figure window
hold on; % Allow multiple plots on the same axes

% Plot original data with error bars
errorbar(x_data, y_data, sigma_data, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b', 'CapSize', 5, 'DisplayName', 'Data points');

% Plot the fitted curve
% 'fit_result' can be plotted directly
h_fit = plot(fit_result, 'r-'); % Plot the cfit object, returns handle
set(h_fit, 'LineWidth', 1.5); % Make the fit line thicker

% Customize plot appearance
%title('Curve Fit to Data (MATLAB)');
xlabel('x');
ylabel('y');
legend('Location', 'best'); % Automatically uses DisplayName, 'fit_result' gets default name
% Update legend entry for the fit
legend_entries = get(gca,'Legend').String;
legend_entries{end} = sprintf('Fit: %.3f * x * exp(%.3f * x)', a_fit, b_fit); % Update last entry
legend(legend_entries, 'Location', 'best');

grid on;
box on;
hold off; % Release the hold on the axes