import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import time

# --- Definicija parametrov problema in modela ---
N_STEPS = 100
Z_TARGET = 0.8
Y_LIMIT = 1.5
BETA = 50.0

# --- Definicija funkcij (ostanejo enake) ---

def calculate_physics(y_profile):
    N = len(y_profile)
    if N < 2: return 0, 0
    dx = 1 / (N - 1)
    total_distance = np.trapezoid(y_profile, dx=dx)
    dy = np.diff(y_profile)
    total_effort = np.sum(dy**2) / dx
    return total_distance, total_effort

def lagrangian_with_penalty(y_internal_points, y0, yf, lambda_val, y_limit, beta):
    full_y_profile = np.concatenate(([y0], y_internal_points, [yf]))
    total_distance, total_effort = calculate_physics(full_y_profile)
    speeding_amount = np.maximum(0, full_y_profile - y_limit)
    penalty_per_step = np.exp(beta * speeding_amount) - 1
    total_penalty = np.sum(penalty_per_step)
    return total_effort - lambda_val * total_distance + total_penalty

def find_optimal_profile(lambda_val, y0, yf, y_limit, beta):
    initial_guess = np.linspace(y0, yf, N_STEPS)[1:-1]
    bounds = [(0, None) for _ in range(N_STEPS - 2)]
    result = optimize.minimize(
        fun=lagrangian_with_penalty,
        x0=initial_guess,
        args=(y0, yf, lambda_val, y_limit, beta),
        method='L-BFGS-B',
        bounds=bounds
    )
    return np.concatenate(([y0], result.x, [yf]))

def distance_error(lambda_val, y0, yf, y_limit, beta, z_target):
    y_profile = find_optimal_profile(lambda_val, y0, yf, y_limit, beta)
    distance_travelled, _ = calculate_physics(y_profile)
    return distance_travelled - z_target

def run_simulation(y0, yf):
    try:
        solution = optimize.root_scalar(
            f=distance_error,
            args=(y0, yf, Y_LIMIT, BETA, Z_TARGET),
            bracket=[-200, 200],
            method='brentq'
        )
        lambda_final = solution.root
        return find_optimal_profile(lambda_final, y0, yf, Y_LIMIT, BETA)
    except ValueError:
        print(f"  OPOZORILO: Rešitev za y0={y0:.2f}, yf={yf:.2f} ni bila najdena. Preskakujem.")
        return None

# --- Glavni del skripte: Generiranje grafov (BREZ LEGEND) ---

start_time = time.time()
cmap = plt.get_cmap('viridis')

# --- SET 1: Fiksna začetna hitrost ---
y0_fixed_values = [0.5, 1.0, 1.5]
print("--- Generiram grafe za fiksno začetno hitrost (Set 1) ---")
for y0_val in y0_fixed_values:
    print(f"\nObdelujem primer za fiksni y0 = {y0_val:.2f}...")
    plt.figure(figsize=(10, 7))
    colors = cmap(np.linspace(0, 1, 10))
    yf_varying_values = np.linspace(0.1, 2.0, 10)
    for i, yf_val in enumerate(yf_varying_values):
        y_profile = run_simulation(y0_val, yf_val)
        if y_profile is not None:
            plt.plot(np.linspace(0, 1, N_STEPS), y_profile, color=colors[i])
    plt.axhline(y=Y_LIMIT, color='r', linestyle='--')
    plt.xlabel("Brezdimenzijski čas (t / t₀)")
    plt.ylabel("Brezdimenzijska hitrost (v / v₀)")
    plt.title(f"Družina optimalnih profilov pri fiksni začetni hitrosti $y_0={y0_val:.2f}$")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(f'semafor_fiksni_y0_{y0_val:.1f}.pdf')
    plt.close()

# --- SET 2: Fiksna končna hitrost ---
yf_fixed_values = [0.5, 1.0, 1.5]
print("\n--- Generiram grafe za fiksno končno hitrost (Set 2) ---")
for yf_val in yf_fixed_values:
    print(f"\nObdelujem primer za fiksni yf = {yf_val:.2f}...")
    plt.figure(figsize=(10, 7))
    colors = cmap(np.linspace(0, 1, 10))
    y0_varying_values = np.linspace(0.1, 2.0, 10)
    for i, y0_val in enumerate(y0_varying_values):
        y_profile = run_simulation(y0_val, yf_val) # Uporabimo popravljeno kodo
        if y_profile is not None:
            plt.plot(np.linspace(0, 1, N_STEPS), y_profile, color=colors[i])
    plt.axhline(y=Y_LIMIT, color='r', linestyle='--')
    plt.xlabel("Brezdimenzijski čas (t / t₀)")
    plt.ylabel("Brezdimenzijska hitrost (v / v₀)")
    plt.title(f"Družina optimalnih profilov pri fiksni končni hitrosti $y_f={yf_val:.2f}$")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(f'semafor_fiksni_yf_{yf_val:.1f}.pdf')
    plt.close()

total_time = time.time() - start_time
print(f"\nCelotna analiza je končana v {total_time:.2f} sekundah.")