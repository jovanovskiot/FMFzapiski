# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy


invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2

def zlatirez(f, a, c, tol=1e-5):
    """Golden section search.

    Given a function f with a single local minimum in
    the interval [a,c], gss returns a subset interval
    [b,d] that contains the minimum with d-c <= tol.

    Example:
    >>> f = lambda x: (x-2)**2
    >>> a = 1
    >>> c = 5
    >>> tol = 1e-5
    >>> (b,d) = zlatirez(f, a, c, tol)
    >>> print(b, d)
    1.9999959837979107 2.0000050911830893
    """

    (a, c) = (min(a, c), max(a, c))
    tocke = np.array([a,c])

    h = c - a # zacetna sirina intervala
    if h <= tol:
        return (a, c), tocke

    # ocena potrebnih korakov za doseg natancnosti - interval se vsakic zmanjsa za faktor 1/phi=invphi
    N = int(math.ceil(math.log(tol / h) / math.log(invphi)))

    b = a + invphi * h
    d = a + invphi2 * h
    yb = f(b)
    yd = f(d)

    n = 0
    while(h > tol or n <= N):
      n += 1
      h = invphi * h
      if yd < yb:
          c = b
          b = d
          yb = yd
          d = a + invphi2 * h
          yd = f(d)
          tocke = np.append(tocke, d)
      else:
          a = d
          d = b
          yd = yb
          b = a + invphi * h
          yb = f(b)
          tocke = np.append(tocke, d)
    if yd < yb:
        return (a, d), tocke, n
    else:
        return (d, c), tocke, n

def get_parab_min(f, a, b, c):
    if ((b - a) * (f(b) - f(c)) - (b - c) * (f(b) - f(a))) == 0:
        return None
    d = b - 0.5 * ((b - a)**2 * (f(b) - f(c)) - (b - c)**2 * (f(b) - f(a))) / ((b - a) * (f(b) - f(c)) - (b - c) * (f(b) - f(a)))
    return d

def parabmin(f, a, c, tol):

    (a, c) = (min(a, c), max(a, c))
    tocke = np.array([a,c])

    h = c - a
    if h <= tol:
        return (a, c), tocke

    b = a + invphi * h
    tocke = np.append(tocke, b)
    ya = f(a)
    yc = f(c)
    yb = f(b)

    assert ya > yb and yc > yb, f"Invalid initial interval ({a}, {c})"
    d = get_parab_min(f, a, b, c)
    yd = f(d)

    # Lahko se zgodi, da bo nova točka na intervalu [B, C]...
    if d >= b:
        d, b = b, d
        yd, yb = yb, yd

    n = 0
    while abs(c - a) > tol:
      if yd < yb:
          c = b
          b = d
          yc = yb
          yb = yd
          d = get_parab_min(f, a, b, c)
          if not d:
              break
          yd = f(d)

          if d >= b:
              d, b = b, d
              yd, yb = yb, yd
          tocke = np.append(tocke, d)
      else:
          a = d
          d = b
          ya = yd
          yd = yb
          b = get_parab_min(f, a, d, c)
          if not b:
              break
          yb = f(b)

          if d >= b:
              d, b = b, d
              yd, yb = yb, yd
          tocke = np.append(tocke, d)
      n += 1
    if yd < yb:
        return (a, d), tocke, n
    else:
        return (d, c), tocke, n

if __name__ == '__main__':

  import math

  def func(x):
      return x**3*np.sin(x)
  def odvod(x):
      return 3*x**2*np.sin(x) + x**3*np.cos(x)

  meja_a = 1. #spodnja meja
  meja_b = 10. #zgornja meja

  # "Točna" vrednost
  minimum, result = scipy.optimize.brentq(odvod, a=3, b=7, full_output=True)
  print(f"Found minimum at {minimum}")
  print(f"\n\n{result}\n\n")

  # Uporaba vgrajene funkcije
  xmin, fval, niter, nfunc = scipy.optimize.brent(func, brack=(3, 5, 7), tol=1e-15, full_output=True)
  print("Brent\n----------")
  print(f"  Mininum: {xmin}")
  print(f"  Relative Error: {abs(xmin - minimum) / minimum:.3e}")
  print(f"  Vrednost funkcije: f({xmin}) = {fval}")
  print(f"  Št poskusov: N = {niter}\n\n")

  # Zlati rez
  rezultat, tocke, koraki = zlatirez(func, 2, 8, tol=1e-15)
  mean = sum(rezultat) / len(rezultat)
  print("Zlati rez\n----------")
  print(f"  Mininum: {mean}")
  print(f"  Relative Error: {abs(mean - minimum) / minimum:.3e}")
  print(f"  Vrednost funkcije: f({mean}) = {func(mean)}")
  print(f"  Št poskusov: N = {koraki}\n\n")

plrange = np.linspace(meja_a, meja_b, 100)
fig, axs = plt.subplots(2, 1, figsize=(8, 10))

# Function plot
axs[0].plot(plrange, func(plrange), "b")
axs[0].grid()
axs[0].plot(tocke, func(tocke), "ro")
axs[0].legend(("funkcija", "priblizki"), loc="lower right")
axs[0].set_xlabel("x")
axs[0].set_ylabel("f(x)")
axs[0].set_title("Zlatirez - Funkcija")

# Derivative plot
axs[1].plot(plrange, odvod(plrange), "b")
axs[1].plot(tocke, odvod(tocke), "ro")
axs[1].plot(plrange, np.zeros(len(plrange)), "g-")
axs[1].legend(("odvod", "priblizki"), loc="lower right")
axs[1].set_xlabel("x")
axs[1].set_ylabel("f'(x)")
axs[1].set_title("Zlatirez - Odvod")

plt.tight_layout()
plt.show()

# Parabolicna metoda
rezultat, tocke, koraki = parabmin(func, 2, 8, tol=1e-15)
mean = sum(rezultat) / len(rezultat)
print("Parabolična metoda\n----------")
print(f"  Mininum: {mean}")
print(f"  Relative Error: {abs(mean - minimum) / minimum:.3e}")
print(f"  Vrednost funkcije: f({mean}) = {func(mean)}")
print(f"  Št poskusov: N = {koraki}\n\n")

# Create a 2-row subplot for Parabolic
fig, axs = plt.subplots(2, 1, figsize=(8, 10))  # 2 rows, 1 column

# Plot function and derivative for Parabolic
axs[0].plot(plrange, func(plrange), "b")
axs[0].grid()
axs[0].plot(tocke, func(tocke), "ro")
axs[0].legend(("funkcija", "priblizki"), loc="lower right")
axs[0].set_xlabel("x")
axs[0].set_ylabel("f(x)")
axs[0].set_title("Parabolična - Funkcija")

axs[1].plot(plrange, odvod(plrange), "b")
axs[1].plot(tocke, odvod(tocke), "ro")
axs[1].plot(plrange, np.zeros(len(plrange)), "g-")
axs[1].legend(("odvod", "priblizki"), loc="lower right")
axs[1].set_xlabel("x")
axs[1].set_ylabel("f'(x)")
axs[1].set_title("Parabolična - Odvod")

plt.tight_layout()
plt.show()
