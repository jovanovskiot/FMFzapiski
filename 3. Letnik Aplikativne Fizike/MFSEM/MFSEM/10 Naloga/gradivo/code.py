#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy

# --- DFT and IDFT functions according to problem statement ---

def dft_problem(hk: np.ndarray) -> np.ndarray:
    """
    Computes DFT as defined in problem statement: H_n = sum_{k=0}^{N-1} h_k exp(2*pi*i*k*n/N)
    This is equivalent to N * np.fft.ifft(hk).
    """
    if hk.ndim != 1:
        raise ValueError("Input array hk must be 1-dimensional.")
    N = hk.shape[0]
    # np.fft.ifft(hk) computes (1/N) * sum_{k=0}^{N-1} h_k exp(2*pi*i*k*n/N)
    return N * np.fft.ifft(hk)

def idft_problem(Hn: np.ndarray) -> np.ndarray:
    """
    Computes IDFT as defined in problem statement: h_k = (1/N) sum_{n=0}^{N-1} H_n exp(-2*pi*i*k*n/N)
    This is equivalent to (1/N) * np.fft.fft(Hn).
    """
    if Hn.ndim != 1:
        raise ValueError("Input array Hn must be 1-dimensional.")
    N = Hn.shape[0]
    # np.fft.fft(Hn) computes sum_{n=0}^{N-1} H_n exp(-2*pi*i*k*n/N)
    return (1/N) * np.fft.fft(Hn)

# --- Signal Analysis Helper Function ---

def analyze_signal(hk: np.ndarray, dt: float, signal_name: str, show_plots: bool = True):
    """
    Analyzes a given signal hk: computes DFT, plots time and frequency domains,
    checks IDFT reconstruction, and verifies Parseval's theorem.
    """
    N = hk.shape[0]
    T_total = N * dt
    fs = 1.0 / dt
    
    print(f"\n--- Analyzing: {signal_name} ---")
    print(f"N = {N}, dt = {dt:.4g} s, T_total = {T_total:.3g} s, fs = {fs:.3g} Hz")
    print(f"Nyquist frequency f_c = {fs/2:.3g} Hz")

    # Compute DFT
    Hn = dft_problem(hk)
    
    # Frequencies for plotting (shifted to center zero frequency)
    freqs_shifted = np.fft.fftshift(np.fft.fftfreq(N, d=dt))
    Hn_shifted = np.fft.fftshift(Hn) # Shift Hn to match freqs_shifted
    
    if show_plots:
        fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
        fig.suptitle(f"Signal: {signal_name}\nN={N}, T_total={T_total:.2f}s, fs={fs:.2f}Hz", fontsize=12)
        
        # Time domain plot
        t_axis = np.arange(N) * dt
        axs[0].plot(t_axis, hk.real, '.-', label="Re(h_k)")
        if np.any(np.iscomplex(hk)) and not np.allclose(hk.imag, 0):
             axs[0].plot(t_axis, hk.imag, '.-', label="Im(h_k)")
        axs[0].set_xlabel("Time (s)")
        axs[0].set_ylabel("h_k")
        axs[0].legend()
        axs[0].grid(True)
        axs[0].set_title("Time Domain")
        
        # Frequency domain plot (Real and Imaginary parts of H_n)
        axs[1].plot(freqs_shifted, Hn_shifted.real, '.-', label="Re(H_n)")
        axs[1].plot(freqs_shifted, Hn_shifted.imag, '.-', label="Im(H_n)")
        axs[1].set_ylabel("H_n")
        axs[1].legend()
        axs[1].grid(True)
        axs[1].set_title("Frequency Domain (DFT H_n as per problem statement)")

        # Power spectrum |H_n/N|^2 (Two-sided)
        power_spectrum_val = np.abs(Hn_shifted / N)**2 
        axs[2].plot(freqs_shifted, power_spectrum_val, '.-', label="|H_n/N|^2")
        axs[2].set_xlabel("Frequency (Hz)")
        axs[2].set_ylabel("Power (scaled)")
        axs[2].legend()
        axs[2].grid(True)
        axs[2].set_title("Power Spectrum (Two-sided, |H_n/N|^2)")
        
        # plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Replaced by constrained_layout
        plt.show()

    # Check IDFT
    hk_reconstructed = idft_problem(Hn)
    is_reconstruction_accurate = np.allclose(hk, hk_reconstructed)
    max_diff = np.max(np.abs(hk - hk_reconstructed))
    if is_reconstruction_accurate:
        print(f"IDFT check: PASSED. Max reconstruction error: {max_diff:.2e}")
    else:
        print(f"IDFT check: FAILED. Max reconstruction error: {max_diff:.2e}")
        if show_plots:
            plt.figure(figsize=(8,5), constrained_layout=True)
            t_axis = np.arange(N) * dt
            plt.plot(t_axis, hk.real, label='Original h_k (real)')
            plt.plot(t_axis, hk_reconstructed.real, '--', label='Reconstructed h_k (real)')
            if (np.any(np.iscomplex(hk)) and not np.allclose(hk.imag,0)) or \
               (np.any(np.iscomplex(hk_reconstructed)) and not np.allclose(hk_reconstructed.imag,0)):
                 if not np.allclose(hk.imag,0): plt.plot(t_axis, hk.imag, label='Original h_k (imag)')
                 if not np.allclose(hk_reconstructed.imag,0): plt.plot(t_axis, hk_reconstructed.imag, ':', label='Reconstructed h_k (imag)')
            plt.legend()
            plt.title(f"IDFT Reconstruction Check for {signal_name}")
            plt.xlabel("Time (s)")
            plt.ylabel("h_k")
            plt.grid(True)
            plt.show()

    # Parseval's theorem check
    sum_sq_hk = np.sum(np.abs(hk)**2)
    sum_sq_Hn_div_N = np.sum(np.abs(Hn)**2) / N # Note Hn is from dft_problem
    parseval_ok = np.isclose(sum_sq_hk, sum_sq_Hn_div_N)
    print(f"Parseval's check (sum|h_k|^2 = (1/N)sum|H_n|^2): {'PASSED' if parseval_ok else 'FAILED'}")
    print(f"  sum|h_k|^2 = {sum_sq_hk:.4e}")
    print(f"  (1/N)sum|H_n|^2 = {sum_sq_Hn_div_N:.4e}")


# --- Part 1: Analysis of Simple Signals ---

def run_part1():
    print("====== Starting Part 1: Simple Signals ======")
    # Base parameters for most tests in Part 1
    N_base = 128
    dt_base = 0.01 # s
    
    fs_base = 1.0 / dt_base      # Sampling frequency (100 Hz)
    T_total_base = N_base * dt_base # Total time duration (1.28 s)
    f0_base = fs_base / N_base   # Fundamental frequency resolution (~0.78125 Hz)

    print(f"Base parameters: N={N_base}, dt={dt_base:.2f}s, T_total={T_total_base:.2f}s, fs={fs_base:.2f}Hz, f0={f0_base:.4f}Hz")
    t_base = np.arange(N_base) * dt_base

    # Test 1: Periodic Cosine Signal
    m1 = 5 # Signal frequency is 5 * f0_base
    f_sig1 = m1 * f0_base 
    hk1 = np.cos(2 * np.pi * f_sig1 * t_base)
    analyze_signal(hk1, dt_base, f"Cosine, periodic, f={f_sig1:.2f}Hz (m={m1})")

    # Test 2: Periodic Sine Signal
    m2 = 8 # Signal frequency is 8 * f0_base
    f_sig2 = m2 * f0_base
    hk2 = np.sin(2 * np.pi * f_sig2 * t_base)
    analyze_signal(hk2, dt_base, f"Sine, periodic, f={f_sig2:.2f}Hz (m={m2})")
    
    # Test 3: Sum of Periodic Signals
    f_sig3a = 3 * f0_base # Cosine component frequency
    f_sig3b = 6 * f0_base # Sine component frequency
    hk3 = 0.5 * np.cos(2 * np.pi * f_sig3a * t_base) + 1.5 * np.sin(2 * np.pi * f_sig3b * t_base)
    analyze_signal(hk3, dt_base, f"Sum, periodic, f_cos={f_sig3a:.2f}Hz (m=3), f_sin={f_sig3b:.2f}Hz (m=6)")

    # Test 4: Non-periodic Signal (Spectral Leakage)
    f_sig4 = 4.5 * f0_base # Frequency is not an integer multiple of f0_base
    hk4 = np.sin(2 * np.pi * f_sig4 * t_base)
    analyze_signal(hk4, dt_base, f"Sine, non-periodic, f={f_sig4:.2f}Hz (m=4.5)")

    # Test 5: Aliasing
    f_nyquist_base = fs_base / 2.0 # 50 Hz for fs_base = 100 Hz
    f_sig5_original = 60.0 # Hz (original frequency > Nyquist)
    # Expected aliased frequency: 60 - 100 = -40 Hz
    hk5 = np.sin(2 * np.pi * f_sig5_original * t_base)
    analyze_signal(hk5, dt_base, f"Sine, aliased, f_orig={f_sig5_original:.2f}Hz (expected alias at -40Hz)")
    
    f_sig5b_original = 110.0 # Hz
    # Expected aliased frequency: 110 - 100 = 10 Hz
    hk5b = np.sin(2 * np.pi * f_sig5b_original * t_base)
    analyze_signal(hk5b, dt_base, f"Sine, aliased, f_orig={f_sig5b_original:.2f}Hz (expected alias at 10Hz)")

    # Test 6: Dependence on N (dt fixed, T_total varies, f0 varies)
    # Demonstrates effect on frequency resolution.
    print(f"\n--- Test 6: Vary N, keep dt fixed ({dt_base}s). Signal freq chosen to be periodic for all N. ---")
    f_test6 = 5 * (fs_base / (N_base//2)) # 7.8125 Hz (m=5 for N=64, m=10 for N=128, m=20 for N=256)

    for N_curr in [N_base // 2, N_base, N_base * 2]:
        t_curr = np.arange(N_curr) * dt_base
        hk_curr = np.sin(2 * np.pi * f_test6 * t_curr)
        m_curr = int(round(f_test6 / (fs_base / N_curr))) # Calculate m for current N
        analyze_signal(hk_curr, dt_base, f"Vary N (dt fixed): N={N_curr}, f_sig={f_test6:.3f}Hz (m={m_curr})")

    # Test 7: Dependence on N (T_total fixed, dt varies, fs and f_nyquist vary)
    # Demonstrates effect on Nyquist frequency.
    print(f"\n--- Test 7: Vary N, keep T_total fixed ({T_total_base}s). Signal freq chosen to be periodic. ---")
    f_test7 = 10.0 / T_total_base # 7.8125 Hz (m=10 for all N in this test, as T_total is fixed)

    for N_curr in [N_base // 2, N_base, N_base * 2]:
        dt_curr = T_total_base / N_curr
        t_curr = np.arange(N_curr) * dt_curr
        hk_curr = np.sin(2 * np.pi * f_test7 * t_curr)
        analyze_signal(hk_curr, dt_curr, f"Vary N (T_total fixed): N={N_curr}, f_sig={f_test7:.3f}Hz (m=10)")
    
    print("\n====== Part 1 finished. Press Enter in terminal after closing plots. ======")


# --- Part 2: Analysis of Audio Signal ---

def run_part2():
    print("\n====== Starting Part 2: Audio Signal Analysis ======")
    audio_file_path = r"C:\Users\Teodor\Desktop\Faks\3_Letnik_AF\MFSEM\10 Naloga\gradivo\waveform.txt"
    try:
        hk_audio = np.loadtxt(audio_file_path)
        if hk_audio.ndim > 1: # Handle stereo or multi-column data by taking the first column
             print(f"Audio data has {hk_audio.shape[1]} columns, using the first column.")
             hk_audio = hk_audio[:,0]
    except FileNotFoundError:
        print(f"ERROR: Audio file '{audio_file_path}' not found.")
        print("Please make sure 'waveform.txt' is in the same directory as the script.")
        print("Skipping Part 2.")
        return
    except Exception as e:
        print(f"Error loading '{audio_file_path}': {e}")
        print("Skipping Part 2.")
        return

    N_audio = hk_audio.shape[0]
    if N_audio == 0:
        print("Audio file is empty. Skipping Part 2.")
        return

    fs_audio = 44100.0 # Hz (sampling frequency given in problem)
    dt_audio = 1.0 / fs_audio
    T_total_audio = N_audio * dt_audio
    
    print(f"Loaded audio signal: N = {N_audio} samples")
    print(f"Sampling frequency fs = {fs_audio:.0f} Hz, dt = {dt_audio*1e6:.2f} us")
    print(f"Total duration T = {T_total_audio:.2f} s")
    print(f"Nyquist frequency f_c = {fs_audio/2:.0f} Hz")

    # Plot a segment of the time-domain signal
    plt.figure(figsize=(12, 6), constrained_layout=True)
    plot_N_max = min(N_audio, int(0.05 * fs_audio)) # Plot up to 50ms or whole signal
    t_audio_segment_ms = (np.arange(plot_N_max) * dt_audio) * 1000 # time in ms
    plt.plot(t_audio_segment_ms, hk_audio[:plot_N_max])
    plt.xlabel("Time (ms)")
    plt.ylabel("Amplitude")
    plt.title(f"Audio Signal (first {plot_N_max*dt_audio*1000:.1f} ms)")
    plt.grid(True)
    plt.show()

    # Compute DFT of the audio signal
    print("Computing DFT of audio signal (this may take a moment for large N)...")
    Hn_audio = dft_problem(hk_audio)
    print("DFT computation finished.")

    freqs_audio_shifted = np.fft.fftshift(np.fft.fftfreq(N_audio, d=dt_audio))
    Hn_audio_shifted = np.fft.fftshift(Hn_audio)
    power_spectrum_audio = np.abs(Hn_audio_shifted / N_audio)**2
    
    # Plot Power Spectrum (linear scale)
    plt.figure(figsize=(12, 6), constrained_layout=True)
    plt.plot(freqs_audio_shifted, power_spectrum_audio)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power |H_n/N|^2 (linear scale)")
    plt.title("Power Spectrum of Audio Signal")
    plt.grid(True)
    plt.xlim(0, fs_audio/2) # Focus on positive frequencies up to Nyquist
    plt.show()

    # Plot Power Spectrum (dB scale) for better visualization of peaks
    plt.figure(figsize=(12, 6), constrained_layout=True)
    epsilon = 1e-18 # Small constant to prevent log(0)
    power_db = 10 * np.log10(power_spectrum_audio + epsilon)
    plt.plot(freqs_audio_shifted, power_db)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power (dB scale)")
    plt.title("Power Spectrum of Audio Signal (dB scale)")
    plt.grid(True)
    plt.xlim(0, fs_audio/2)
    if N_audio > 0 and len(power_db[freqs_audio_shifted >= 0]) > 0:
        max_db_in_positive_freqs = np.max(power_db[freqs_audio_shifted >= 0])
        plt.ylim(max_db_in_positive_freqs - 80, max_db_in_positive_freqs + 5) # Dynamic range
    plt.show()

    # Identify prominent frequencies using scipy.signal.find_peaks
    try:
        from scipy.signal import find_peaks
        
        idx_positive_freqs_start = N_audio // 2 
        positive_freqs_hz = freqs_audio_shifted[idx_positive_freqs_start:]
        positive_power_linear = power_spectrum_audio[idx_positive_freqs_start:]
        
        if len(positive_freqs_hz) > 0:
            max_power_val_db = np.max(10 * np.log10(positive_power_linear + epsilon))
            min_height_db = max_power_val_db - 40 # Peaks within 40 dB of the max peak
            min_height_linear = 10**(min_height_db / 10.0)

            df_hz = fs_audio / N_audio # Frequency bin width in Hz
            min_dist_hz = 50.0 # Minimum separation between peaks in Hz
            min_dist_samples = int(min_dist_hz / df_hz) if df_hz > 0 else 1
            min_dist_samples = max(1, min_dist_samples)

            peaks_indices, properties = find_peaks(positive_power_linear, 
                                                 height=min_height_linear, 
                                                 distance=min_dist_samples)
        
            print("\nIdentified prominent frequencies (using scipy.signal.find_peaks):")
            if len(peaks_indices) > 0:
                detected_peak_freqs_hz = positive_freqs_hz[peaks_indices]
                detected_peak_powers_linear = properties['peak_heights']
                
                sorted_indices_by_freq = np.argsort(detected_peak_freqs_hz)
                print(f"{'Frequency (Hz)':<18} {'Power (linear)':<18} {'Power (dB)':<15}")
                print("-" * 55)
                for i in sorted_indices_by_freq:
                    freq_val = detected_peak_freqs_hz[i]
                    pow_lin_val = detected_peak_powers_linear[i]
                    pow_db_val = 10 * np.log10(pow_lin_val + epsilon)
                    print(f"{freq_val:<18.2f} {pow_lin_val:<18.2e} {pow_db_val:<15.1f}")
            else:
                print("  No significant peaks found with current criteria.")
        else:
            print("  No positive frequency data to analyze for peaks.")
            
    except ImportError:
        print("\n`scipy` module not found. Skipping detailed peak analysis with `find_peaks`.")
        print("Please inspect the plots manually to identify prominent frequencies.")
    except Exception as e:
        print(f"\nError during peak finding: {e}")
        print("Please inspect the plots manually to identify prominent frequencies.")

    print("\n====== Part 2 finished. Press Enter in terminal after closing plots. ======")


# --- Main execution ---

if __name__ == "__main__":
    run_part1()
    # input("Press Enter to continue to Part 2...") # Optional pause
    run_part2()
    
    print("\nAll tasks complete. Ensure Matplotlib plots are closed to terminate the script if needed.")