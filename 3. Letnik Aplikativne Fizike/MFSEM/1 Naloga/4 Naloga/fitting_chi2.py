#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

###############################################################################
# 1) Define the polynomial model and chi-square
###############################################################################
def poly_model(params, x):
    """
    Polynomial model: y(x) = a + b*x + c*x^2

    params: [a, b, c]
    x: float or array
    """
    a, b, c = params
    return a + b*x + c*x*x

def chi2_poly(params, xvals, yvals, errors):
    """
    chi^2 = Σ( (y_i - model(params,x_i))^2 / err_i^2 )

    params: [a, b, c]
    xvals, yvals, errors: data + uncertainties
    """
    chi2_val = 0.0
    for x, y, err in zip(xvals, yvals, errors):
        y_fit = poly_model(params, x)
        chi2_val += ((y - y_fit)**2) / (err**2)
    return chi2_val

###############################################################################
# 2) Example data
###############################################################################
np.random.seed(1234)
params = [0, np.random.uniform(-10, 10), np.random.uniform(-10, 10)]
print(f"Starting parameters {params}")

xdata = np.linspace(0.0, 10.0, 25)
ydata = poly_model(params, xdata)
print(ydata)

offsets = np.random.normal(size=len(xdata), scale=5)
print(offsets)

ydata += offsets
print(ydata)

sigma = 5 * np.ones(len(xdata))  # uncertainties

# Initial guess for [a, b, c]
init_guess = [0.0, 0.0, 0.0]

###############################################################################
# 3) Fit parameters using downhill simplex
###############################################################################
best_params = scipy.optimize.fmin(chi2_poly, init_guess, args=(xdata, ydata, sigma), ftol=1e-15)
print("Best-fit polynomial parameters [a, b, c] =", best_params)

###############################################################################
# 4) Plot data vs. best-fit polynomial
###############################################################################
plt.figure(figsize=(7,5))
plt.title("Polynomial Fit: y = a + b·x + c·x²")
plt.errorbar(xdata, ydata, yerr=sigma, fmt='ro', label="Data")

# Plot best-fit curve over a dense grid of x
x_fine = np.linspace(np.min(xdata), np.max(xdata), 300)
y_fine = poly_model(best_params, x_fine)
plt.plot(x_fine, y_fine, 'g-', label="Best-fit curve")

plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.show()

###############################################################################
# 5) 2D Chi2 slice: vary (a,b) while holding c fixed
###############################################################################
a_fixed = best_params[0]  # keep 'a' at best-fit
b_vals = np.linspace(best_params[1]-1, best_params[1]+1, 100)
c_vals = np.linspace(best_params[2]-1, best_params[2]+1, 100)

chi2_grid = np.zeros((len(b_vals), len(c_vals)))
for j, b_test in enumerate(b_vals):
    for i, c_test in enumerate(c_vals):
        params_test = [a_fixed, b_test, c_test]
        chi2_grid[j, i] = chi2_poly(params_test, xdata, ydata, sigma)

plt.figure(figsize=(7,5))
plt.title(r"$\chi^2$ Slice at a = {:.3f}".format(a_fixed))
plt.xlabel("b")
plt.ylabel("c")

# Plot 2D heatmap of chi2
extent = [b_vals[0], b_vals[-1], c_vals[0], c_vals[-1]]
plt.imshow(chi2_grid, origin='lower', extent=extent, aspect='auto', norm='log')
cbar = plt.colorbar()
cbar.set_label(r'$\chi^2$')

# Mark the best-fit (a,b) on the heatmap
plt.plot([best_params[1]], [best_params[2]], 'r*', ms=12, label="Best Fit (b, c)", markerfacecolor='none')
plt.legend()
plt.show()

#################################################################################

#################################################################################

# Data for {x_i}
# Using np.arange is efficient for this regular sequence
x_data = np.arange(1.0, 11.5, 0.5)

# Alternatively, by direct transcription:
# x_data = np.array([
#     1.00000000, 1.50000000, 2.00000000, 2.50000000, 3.00000000, 3.50000000,
#     4.00000000, 4.50000000, 5.00000000, 5.50000000, 6.00000000, 6.50000000,
#     7.00000000, 7.50000000, 8.00000000, 8.50000000, 9.00000000, 9.50000000,
#     10.00000000, 10.50000000, 11.00000000
# ])

# Data for {y_i}
y_data = np.array([
    0.31700705, 0.43791106, 0.56528271, 0.56102378, 0.63664784, 0.65121353,
    0.63487502, 0.64501481, 0.60942923, 0.62411336, 0.61455575, 0.57226264,
    0.54291294, 0.50329224, 0.50314769, 0.46050043, 0.42461463, 0.40771586,
    0.41605889, 0.36732963, 0.33085992
])

# Data for {σ_i}
sigma_data = np.array([
    0.01548814, 0.01715189, 0.01602763, 0.01544883, 0.01423655, 0.01645894,
    0.01437587, 0.01891773, 0.01963663, 0.01383442, 0.01791725, 0.01528895,
    0.01568045, 0.01925597, 0.01071036, 0.01087129, 0.01020218, 0.01832620,
    0.01778157, 0.01870012, 0.01978618
])