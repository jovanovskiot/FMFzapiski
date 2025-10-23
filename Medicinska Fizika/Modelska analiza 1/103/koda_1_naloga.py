import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import time

# --- Glavne nastavitve (za originalni Thomsonov problem) ---
MAX_N = 20           # Maksimalno število nabojev za analizo
NUM_REPEATS = 5      # Število ponovitev za vsak N za statistiko

# --- Definicija funkcij (za originalni Thomsonov problem) ---

def calculate_energy(coords, N):
    """Izračuna skupno elektrostatično energijo za N nabojev (potencial 1/r)."""
    positions = np.reshape(coords, (N, 2))
    theta, phi = positions[:, 0], positions[:, 1]
    
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    total_energy = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            r_ij = np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2)
            total_energy += 1 / r_ij
            
    return total_energy

def plot_arrangement(optimal_coords, N, filename, title):
    """Izriše končno 3D razporeditev delcev na sferi."""
    positions = np.reshape(optimal_coords, (N, 2))
    theta, phi = positions[:, 0], positions[:, 1]
    
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    
    points = np.stack([x, y, z], axis=1)
    hull = ConvexHull(points)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('white')

    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    x_sphere = np.cos(u) * np.sin(v)
    y_sphere = np.sin(u) * np.sin(v)
    z_sphere = np.cos(v)
    ax.plot_surface(x_sphere, y_sphere, z_sphere, color='deepskyblue', alpha=0.15)

    ax.plot_trisurf(points[:,0], points[:,1], points[:,2], triangles=hull.simplices, color='cyan', alpha=0.4)
    ax.scatter(x, y, z, color='yellow', s=80, edgecolor='black')

    ax.set_box_aspect([1, 1, 1])
    ax.axis('off')
    ax.set_title(title, fontsize=16, color='black', y=0.95)

    plt.savefig(filename, bbox_inches='tight', pad_inches=0.1, facecolor='white')
    plt.close()

# --- Glavni del skripte (za originalni Thomsonov problem) ---

# ... (VSA OBSTOJEČA KODA ZA THOMSONOV PROBLEM OSTANE TUKAJ NESPREMENJENA) ...
# ... (Zanka, izračuni, izris grafa energija_vs_N.pdf, itd.) ...

# Seznami za shranjevanje rezultatov za končni graf
n_values = []
mean_energies_nm, std_devs_nm = [], []
mean_energies_powell, std_devs_powell = [], []
last_optimal_coords_nm = None

start_time = time.time()

for N in range(2, MAX_N + 1):
    print(f"\n--- Računam (Thomson) za N = {N} ---")
    current_energies_nm, current_energies_powell = [], []
    for i in range(NUM_REPEATS):
        initial_guess = np.random.rand(2 * N)
        initial_guess[0::2] *= np.pi
        initial_guess[1::2] *= 2 * np.pi
        res_nm = optimize.minimize(calculate_energy, initial_guess, args=(N,), method='Nelder-Mead')
        res_powell = optimize.minimize(calculate_energy, initial_guess, args=(N,), method='Powell')
        current_energies_nm.append(res_nm.fun)
        current_energies_powell.append(res_powell.fun)
        if N == MAX_N and res_nm.fun == min(current_energies_nm):
            last_optimal_coords_nm = res_nm.x
    mean_energies_nm.append(np.mean(current_energies_nm))
    std_devs_nm.append(np.std(current_energies_nm))
    mean_energies_powell.append(np.mean(current_energies_powell))
    std_devs_powell.append(np.std(current_energies_powell))
    n_values.append(N)
    print(f"Povprečna energija (Nelder-Mead): {mean_energies_nm[-1]:.4f} ± {std_devs_nm[-1]:.4f}")
    print(f"Povprečna energija (Powell):     {mean_energies_powell[-1]:.4f} ± {std_devs_powell[-1]:.4f}")

total_time = time.time() - start_time
print(f"\nAnaliza (Thomson) končana v {total_time:.2f} sekundah.")

print("Shranjujem graf Energija vs. N (Thomson)...")
plt.style.use('seaborn-v0_8-darkgrid')
fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(n_values, mean_energies_nm, yerr=std_devs_nm, fmt='o-', capsize=5, label='Nelder-Mead', color='cyan', ecolor='lightblue')
ax.errorbar(n_values, mean_energies_powell, yerr=std_devs_powell, fmt='o--', capsize=5, label='Powell', color='magenta', ecolor='lightpink')
ax.set_xlabel("Število nabojev (N)", fontsize=14)
ax.set_ylabel("Minimalna energija E", fontsize=14)
ax.set_title("Minimalna energija sistema v odvisnosti od števila nabojev (Thomsonov problem)", fontsize=16)
ax.legend(fontsize=12)
ax.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig('energija_vs_N.pdf')
plt.close()

print(f"Shranjujem 3D vizualizacijo (Thomson) za N={MAX_N}...")
plot_arrangement(last_optimal_coords_nm, MAX_N, 'koncna_razporeditev_N20.pdf', title=f'Thomsonova razporeditev za N={MAX_N}')

print("\nOriginalna analiza končana. Vsi rezultati so shranjeni v .pdf datotekah.")

# --- Dodatni problem: Lennard-Jonesov potencial ---
# Ta odsek je nov in se izvede po koncu originalne analize.

print("\n" + "="*50)
print("Začenjam dodatni problem: Lennard-Jonesovi delci na sferi")
print("="*50)

# --- Parametri za Lennard-Jonesov problem ---
N_LJ = 40      # Število delcev za prikaz grupiranja
SIGMA = 0.6    # Karakteristična razdalja, pri kateri je potencial 0.
               # To je ključni parameter, ki vpliva na velikost "grudic".

def calculate_LJ_energy(coords, N, sigma):
    """
    Izračuna skupno energijo za delce, med katerimi deluje Lennard-Jonesov potencial.
    E_ij = (sigma/r_ij)^12 - (sigma/r_ij)^6
    """
    positions = np.reshape(coords, (N, 2))
    theta, phi = positions[:, 0], positions[:, 1]
    
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    total_energy = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            r_ij = np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2)
            
            # Izračun Lennard-Jonesovega potenciala
            # Izognemo se deljenju z nič, če sta delca slučajno na istem mestu
            if r_ij > 1e-6:
                inv_r = sigma / r_ij
                term_6 = inv_r**6
                term_12 = term_6**2
                energy_ij = term_12 - term_6
                total_energy += energy_ij
            
    return total_energy

# Priprava in izvedba simulacije
print(f"Iščem optimalno razporeditev za N={N_LJ} delcev s SIGMA={SIGMA}...")
lj_start_time = time.time()

# Generiramo naključno začetno postavitev
initial_guess_lj = np.random.rand(2 * N_LJ)
initial_guess_lj[0::2] *= np.pi
initial_guess_lj[1::2] *= 2 * np.pi

# Kličemo minimizator z novo energijsko funkcijo in njenimi parametri (N, sigma)
result_lj = optimize.minimize(
    calculate_LJ_energy, 
    initial_guess_lj, 
    args=(N_LJ, SIGMA), 
    method='Nelder-Mead'
)
lj_optimal_coords = result_lj.x

lj_total_time = time.time() - lj_start_time
print(f"Minimizacija (Lennard-Jones) končana v {lj_total_time:.2f} sekundah.")
print(f"Končna energija: {result_lj.fun:.4f}")

# Vizualizacija rezultata z obstoječo funkcijo
print("Shranjujem 3D vizualizacijo za Lennard-Jonesov problem...")
plot_arrangement(
    lj_optimal_coords, 
    N_LJ, 
    'lennard_jones_arrangement.pdf', 
    title=f'Razporeditev z Lennard-Jonesovim potencialom (N={N_LJ}, σ={SIGMA})'
)

print("\nDodatni problem je končan. Rezultat je shranjen v 'lennard_jones_arrangement.pdf'.")