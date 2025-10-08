#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy

def func(x, a, b):
    return a*(1-np.exp(b*x))

t = np.linspace(0.4, 2.0, 9)
temps = np.array([13.76, 17.47, 19.95, 21.62, 22.73, 23.48, 23.98, 24.32, 24.54])
temps = temps + 0.4*np.random.normal(size=len(temps))
errors = np.array([1.2, 0.9, 1.0, 0.8, 0.7, 0.7, 0.6, 0.5, 0.5])

# Z nedolocenostmi, utezeno
fitParams, fitCovariances = scipy.optimize.curve_fit(func, t, temps, p0 = [25.0, -1.0], sigma=errors)

print(' fit coefficients:\n', fitParams)
print(' Covariance matrix:\n', fitCovariances)

plt.cla()
plt.ylabel('Temperature (C)', fontsize = 16)
plt.xlabel('time (s)', fontsize = 16)
plt.xlim(0,2.2)
plt.ylim(12, 26)
plt.errorbar(t, temps, yerr = errors, xerr = 0.1, fmt='ro')
plt.plot(t, func(t, fitParams[0], fitParams[1]),
     t, func(t, fitParams[0] + np.sqrt(fitCovariances[0,0]), fitParams[1] - np.sqrt(fitCovariances[1,1])),
     t, func(t, fitParams[0] - np.sqrt(fitCovariances[0,0]), fitParams[1] + np.sqrt(fitCovariances[1,1])))

plt.show()
