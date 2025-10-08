#!/usr/bin/env python3
import numpy as np
from numpy.fft import fft, fftshift, ifftshift, fftfreq
import matplotlib.pyplot as plt

plt.ion()


def gaussian(x: np.ndarray, center: float, width: float, amplitude: float) -> np.ndarray:
    """Return a Gaussian centered at `center` with given width and amplitude."""
    return amplitude * np.exp(-((x - center) ** 2) / width ** 2)


def compute_fft(signal: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute three variants of the FFT of `signal` sampled with spacing dt:
      1. properly shifted (using ifftshift → fft → fftshift)
      2. unshifted (fft → shift)
      3. explicit phase-correction
    Returns (frequencies, H_shifted, H_unshifted, H_corrected).
    """
    
    # Signal size
    N = signal.size

    # Natural frequency bins
    freqs = fftshift(fftfreq(N, d=dt))

    # (1) correct: undo data shift, then FFT, then center
    H_shifted = fftshift(fft(ifftshift(signal))) / N

    # (2) naive: FFT → center
    H_unshifted = fftshift(fft(signal)) / N

    # (3) explicit correction for a window shift of T
    T = N * dt
    phase = np.exp(-2j * np.pi * freqs * (T / 2))
    H_corrected = fftshift(fft(signal) * phase) / N

    return freqs, H_shifted, H_unshifted, H_corrected


def main():
    # --- parameters ---
    N = 100               # number of samples
    T = 1.0               # total interval length
    dt = T / N            # sampling interval
    fs = 1 / dt           # sampling frequency
    f_nyquist = 0.5 * fs  # Nyquist frequency

    print(f"Samples: {N} over interval T={T}")
    print(f"dt = {dt:.5f}, fs = {fs:.2f}, f_Nyquist = {f_nyquist:.2f}")

    # --- time grid ---
    t = np.linspace(0, T, N, endpoint=False)

    # --- signal: single Gaussian pulse centered at T/2 ---
    amp = 1.0
    sigma = 0.05
    center = 0.5 * T
    h = gaussian(t, center, sigma, amp)

    # --- plot time-domain signal ---
    plt.figure()
    plt.plot(t, h, 'k.-')
    plt.title("Time-domain signal: Gaussian pulse")
    plt.xlabel("t")
    plt.ylabel("h(t)")
    plt.grid(True)
    plt.show()
    input("Press ENTER to compute FFT... ")

    # --- compute FFT variants ---
    freqs, H_shifted, H_unshifted, H_corrected = compute_fft(h, dt)

    # --- plot frequency-domain results ---
    fig, ax = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    fig.suptitle("FFT of the Gaussian Pulse (N points = {})".format(N))

    ax[0].plot(freqs, H_shifted.real, label="shifted", color='C1')
    ax[0].plot(freqs, H_unshifted.real, '--', label="no shift", color='C2')
    ax[0].plot(freqs, H_corrected.real, ':', label="corrected", color='C3')
    ax[0].set_ylabel("Re H(f)")
    ax[0].legend(loc="upper right")
    ax[0].grid(True)

    ax[1].plot(freqs, H_shifted.imag, label="shifted", color='C1')
    ax[1].plot(freqs, H_unshifted.imag, '--', label="no shift", color='C2')
    ax[1].plot(freqs, H_corrected.imag, ':', label="corrected", color='C3')
    ax[1].set_ylabel("Im H(f)")
    ax[1].grid(True)

    ax[2].plot(freqs, np.abs(H_shifted) ** 2, label="|H_shifted|²", color='C4')
    ax[2].set_ylabel("Power |H|²")
    ax[2].set_xlabel("Frequency (Hz)")
    ax[2].grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    input("Press ENTER to exit... ")


if __name__ == "__main__":
    main()