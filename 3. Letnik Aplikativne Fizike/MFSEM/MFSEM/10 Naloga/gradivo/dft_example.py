#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Enable interactive mode for live plotting
plt.ion()


def dft_slow(x: np.ndarray) -> np.ndarray:
    """Compute the Discrete Fourier Transform (DFT) of 1D array x (matrix method)."""
    x = np.asarray(x, dtype=float)
    N = x.size
    n = np.arange(N)
    k = n.reshape((N, 1))
    M = np.exp(-2j * np.pi * k * n / N)
    return M @ x


def dft_direct(x: np.ndarray) -> np.ndarray:
    """Compute the DFT of 1D array x via direct summation (double loop)."""
    x = np.asarray(x, dtype=float)
    N = x.size
    result = np.zeros(N, dtype=complex)
    for k in range(N):
        # sum over n: x[n] * exp(-2πi k n / N)
        angles = np.exp(-2j * np.pi * k * np.arange(N) / N)
        result[k] = np.dot(angles, x)
    return result


def main():
    # --- setup sampling grid and signal ---
    N = 250                # number of points
    T = 250.0              # total period
    dt = T / N             # time step
    fs = 1 / dt            # sampling frequency
    f_nyquist = 0.5 * fs   # Nyquist (critical) frequency

    print(f"Sampling frequency: {fs:.2f} Hz")
    print(f"Nyquist frequency:   {f_nyquist:.2f} Hz")

    t = np.linspace(0, T, N, endpoint=False)

    # signal frequencies (in Hz)
    f1 = 1 / 5.0
    f2 = 1 / 10.0
    print(f"f₁ (sin): {f1:.2f} Hz, f₂ (cos): {f2:.2f} Hz")

    # composite signal: sin(2π f1 t) + 2·cos(2π f2 t)
    h = np.sin(2 * np.pi * f1 * t) + 2 * np.cos(2 * np.pi * f2 * t)

    # --- plot time domain ---
    plt.figure()
    plt.plot(t, h, 'r.-')
    plt.xlabel("t")
    plt.ylabel("h(t)")
    plt.title("h(t) = sin(2πf₁ t) + 2 cos(2πf₂ t)")
    plt.grid(True)
    plt.show()
    input("Press ENTER to compute DFT... ")

    # --- compute DFT and shift zero-frequency to center ---
    H = dft_slow(h)             # or dft_direct(h)
    H_shifted = np.roll(H, int(N/2)) / N

    # Alternatively, use numpy's fft implementation:
    # H_shifted = np.fft.fftshift(H) / N

    freqs = np.linspace(-f_nyquist, f_nyquist, N, endpoint=False)

    # --- plot frequency domain ---
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 6))
    fig.suptitle("Discrete Fourier Transform of h(t)")

    axes[0].plot(freqs, H_shifted.real, 'b')
    axes[0].set_ylabel("Re H(k)")

    axes[1].plot(freqs, H_shifted.imag, 'r')
    axes[1].set_ylabel("Im H(k)")

    axes[2].plot(freqs, np.abs(H_shifted)**2, 'k')
    axes[2].set_ylabel("|H(k)|²")
    axes[2].set_xlabel("Frequency (Hz)")

    for ax in axes:
        ax.grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    input("Press ENTER to exit... ")


if __name__ == "__main__":
    main()