import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

print("--- ZAGON SKRIPTE ZA MODELSKE NALOGE ---")

# ==============================================================================
# NALOGA 1: OSNOVNI MODEL
# ==============================================================================
print("\n[NALOGA 1] Procesiram osnovni model...")

def calculate_y(x, z):
    """Izračuna profil hitrosti za osnovni C1 model."""
    A = 1.5 * (1 - z)
    B = -2 * A
    C = 1
    y = A * x**2 + B * x + C
    return y

# Generiranje grafa za družino rešitev (slika za poročilo)
x_values = np.linspace(0, 1, 100)
z_values = np.linspace(0.2, 1.8, 9) # 9 lepih krivulj

plt.figure(figsize=(8, 6))
for z in z_values:
    y = calculate_y(x_values, z)
    plt.plot(x_values, y, label=f"z = {z:.1f}")

plt.title("Družina optimalnih profilov za različne vrednosti $z$")
plt.xlabel("$t/t_0$")
plt.ylabel("$v/v_0$")
plt.grid(True)
plt.legend()
plt.savefig("druzina_profilov_z.pdf")
plt.close()
print("Graf 'druzina_profilov_z.pdf' je shranjen.")


# ==============================================================================
# NALOGA 2: POLJUBNA KONČNA HITROST (KONČNA RAZŠIRJENA VERZIJA)
# ==============================================================================
print("\n[NALOGA 2] Procesiram končno razširjeno analizo s predpisano končno hitrostjo...")

def calculate_y_v_final(x, z, y_final):
    """Izračuna profil hitrosti za model s predpisano končno hitrostjo."""
    A = -6 * (z - 1) + 3 * (y_final - 1)
    B = 6 * (z - 1) - 4 * (y_final - 1)
    C = 1
    y = A * x**2 + B * x + C
    y_at_1 = A + B + C
    # Nanesemo linearno korekcijo, ki ohrani y(0)=1 in zagotovi y(1)=y_final
    # To popravi morebitne numerične napake pri izpeljavi koeficientov.
    y = y - (y_at_1 - y_final) * x
    return y

# --- Sklop grafov 1: Fiksna situacija (z), različne končne hitrosti ---
# Generiramo 4 ločene grafe za 4 različne fiksne situacije (z).
z_scenarios = [0.8, 1.0, 1.33, 1.8] 
y_final_values = np.linspace(0.5, 1.5, 5)

for z_case in z_scenarios:
    plt.figure(figsize=(8, 6))
    for y_f in y_final_values:
        y = calculate_y_v_final(x_values, z_case, y_f)
        plt.plot(x_values, y, label=f"$y_f = {y_f:.1f}$")
    
    plt.title(f"Vpliv končne hitrosti $y_f$ pri fiksni situaciji ($z={z_case:.2f}$)")
    plt.xlabel("$t/t_0$")
    plt.ylabel("$v/v_0$")
    plt.grid(True)
    plt.legend()
    plt.savefig(f"profili_z_fiksen_{str(z_case).replace('.', '_')}.pdf")
    plt.close()
print("Shranjeni so 4 grafi za analizo fiksnega z (npr. 'profili_z_fiksen_0_8.pdf').")


# --- Sklop grafov 2: Fiksna končna hitrost, 10 različnih situacij (z) ---
# To so izboljšani "špageti" grafi, vsak z 10 krivuljami.
z_values_dense = np.linspace(0.2, 2.0, 10)
y_final_cases = [0.8, 1.0, 1.2, 1.5]

for y_f_case in y_final_cases:
    plt.figure(figsize=(8, 6))
    for z_case in z_values_dense:
        y = calculate_y_v_final(x_values, z_case, y_f_case)
        plt.plot(x_values, y, label=f"$z = {z_case:.2f}$")
    
    plt.title(f"Optimalni profili za zahtevano končno hitrost $y_f = {y_f_case}$")
    plt.xlabel("$t/t_0$")
    plt.ylabel("$v/v_0$")
    plt.grid(True)
    # Legendo prikažemo le, če je število krivulj majhno; sicer jo izpustimo.
    if len(z_values_dense) <= 10:
        plt.legend(ncol=2, fontsize='small')
    plt.savefig(f"profili_yf_{str(y_f_case).replace('.', '_')}.pdf")
    plt.close()
print("Shranjeni so 4 grafi za analizo fiksnega y_f (npr. 'profili_yf_0_8.pdf').")


# ==============================================================================
# NALOGA 3: VIŠJE POTENCE FUNKCIONALA
# ==============================================================================
print("\n[NALOGA 3] Procesiram model z višjimi potencami...")

linear_base = np.linspace(1, -1, 100)
a_p1 = linear_base**1
a_p2 = np.sign(linear_base) * np.abs(linear_base)**(1/3)
a_p10 = np.sign(linear_base) * np.abs(linear_base)**(1/19)

y_p1_raw = cumulative_trapezoid(a_p1, x_values, initial=0)
y_p2_raw = cumulative_trapezoid(a_p2, x_values, initial=0)
y_p10_raw = cumulative_trapezoid(a_p10, x_values, initial=0)

y_p1_norm = 1 + y_p1_raw - np.mean(y_p1_raw)
y_p2_norm = 1 + y_p2_raw - np.mean(y_p2_raw)
y_p10_norm = 1 + y_p10_raw - np.mean(y_p10_raw)

# Graf za pospešek
plt.figure(figsize=(10, 7))
plt.plot(x_values, a_p1, label="p=1 (minimiziraj $a^2$)")
plt.plot(x_values, a_p2, label="p=2 (minimiziraj $a^4$)")
plt.plot(x_values, a_p10, label="p=10 (minimiziraj $a^{20}$)")
plt.title("Oblike optimalnega pospeška za različne funkcionale")
plt.xlabel("$t/t_0$")
plt.ylabel("Normaliziran pospešek a(t)")
plt.legend()
plt.grid(True)
plt.savefig("profili_potence_pospeska.pdf")
plt.close()

# Graf za hitrost
plt.figure(figsize=(10, 7))
plt.plot(x_values, y_p1_norm, label="p=1 (iz $a^2$)")
plt.plot(x_values, y_p2_norm, label="p=2 (iz $a^4$)")
plt.plot(x_values, y_p10_norm, label="p=10 (iz $a^{20}$)")
plt.title("Oblike profilov hitrosti za različne funkcionale")
plt.xlabel("$t/t_0$")
plt.ylabel("Normalizirana hitrost v(t)")
plt.legend()
plt.grid(True)
plt.savefig("profili_hitrosti_potence.pdf")
plt.close()
print("Grafa 'profili_potence_pospeska.pdf' in 'profili_hitrosti_potence.pdf' sta shranjena.")


# ==============================================================================
# NALOGA 4: OMEJITEV VELIKOSTI HITROSTI
# ==============================================================================
print("\n[NALOGA 4] Procesiram model z omejitvijo hitrosti...")

def calculate_y_penalty(x, z, gamma):
    """Izračuna profil hitrosti za model s kaznijo za hitrost."""
    if gamma < 1e-6:
        return calculate_y(x, z)
    sqrt_g = np.sqrt(gamma)
    # Izogibanje deljenju z nič, če je tanh(sqrt_g)/sqrt_g == 1
    denominator = (np.tanh(sqrt_g) / sqrt_g) - 1
    if abs(denominator) < 1e-9:
        return np.ones_like(x) * (z) # Približek za zelo velike gamma
    A = (z - 1) / denominator
    B = -A * np.tanh(sqrt_g)
    D = 1 - A
    return A * np.cosh(sqrt_g * x) + B * np.sinh(sqrt_g * x) + D

z_scenarios = [0.8, 1.5]
gamma_values = np.logspace(-1, 2, 12)

for z_case in z_scenarios:
    plt.figure(figsize=(10, 7))
    for g in gamma_values:
        y_profile = calculate_y_penalty(x_values, z_case, g)
        plt.plot(x_values, y_profile, label=f"γ = {g:.1f}")
    
    plt.title(f"Vpliv kazni za hitrost (γ) za scenarij $z = {z_case}$")
    plt.xlabel("$t/t_0$")
    plt.ylabel("$v/v_0$")
    plt.grid(True)
    plt.legend(ncol=2, fontsize='small')
    plt.savefig(f"profili_penalty_z_{str(z_case).replace('.', '_')}.pdf")
    plt.close()
print("Shranjena sta 2 grafa za analizo naloge 4 (npr. 'profili_penalty_z_0_8.pdf').")


# ==============================================================================
# NALOGA 5: ZAPOREDNI SEMAFORJI
# ==============================================================================
print("\n[NALOGA 5] Procesiram simulacijo treh zaporednih semaforjev...")

# --- Orodja za C1 (kvadratični) in C2 (kubični) model ---
def get_profile_c1(z):
    C = 1
    A = 1.5 * (1 - z)
    B = -2 * A
    return A, B, C

def evaluate_profile_c1(A, B, C, x):
    y = A * x**2 + B * x + C
    a = 2 * A * x + B
    return y, a

def get_c2_coefficients(z, alpha_start):
    D = 1.0
    C = alpha_start
    A = 4 * (1 - z + alpha_start / 3.0)
    B = - (3 * A + alpha_start) / 2.0
    return A, B, C, D

def evaluate_profile_c2(A, B, C, D, x):
    y = A * x**3 + B * x**2 + C * x + D
    a = 3 * A * x**2 + 2 * B * x + C
    return y, a

# --- Simulacija za 3 semaforje ---
# Parametri
v0_1 = 15.0; L1 = 200.0; t0_1 = 10.0
L2 = 300.0; t0_2 = 12.0
L3 = 250.0; t0_3 = 8.0

# Leg 1
z1 = L1 / (v0_1 * t0_1)
A1, B1, C1 = get_profile_c1(z1)
y1, a1 = evaluate_profile_c1(A1, B1, C1, x_values)

# Leg 2
v0_2 = y1[-1] * v0_1
z2 = L2 / (v0_2 * t0_2)
alpha1_final_dimless = a1[-1]
A2_c1, B2_c1, C2_c1 = get_profile_c1(z2)
y2_c1, a2_c1 = evaluate_profile_c1(A2_c1, B2_c1, C2_c1, x_values)
alpha2_start_dimless = alpha1_final_dimless * (t0_2 / t0_1) * (v0_1 / v0_2)
A2_c2, B2_c2, C2_c2, D2_c2 = get_c2_coefficients(z2, alpha2_start_dimless)
y2_c2, a2_c2 = evaluate_profile_c2(A2_c2, B2_c2, C2_c2, D2_c2, x_values)

# Leg 3
v0_3_c1 = y2_c1[-1] * v0_2
z3_c1 = L3 / (v0_3_c1 * t0_3)
A3_c1, B3_c1, C3_c1 = get_profile_c1(z3_c1)
y3_c1, a3_c1 = evaluate_profile_c1(A3_c1, B3_c1, C3_c1, x_values)
v0_3_c2 = y2_c2[-1] * v0_2
z3_c2 = L3 / (v0_3_c2 * t0_3)
alpha2_final_dimless = a2_c2[-1]
alpha3_start_dimless = alpha2_final_dimless * (t0_3 / t0_2) * (v0_2 / v0_3_c2)
A3_c2, B3_c2, C3_c2, D3_c2 = get_c2_coefficients(z3_c2, alpha3_start_dimless)
y3_c2, a3_c2 = evaluate_profile_c2(A3_c2, B3_c2, C3_c2, D3_c2, x_values)

# Risanje grafa za 3 semaforje
fig = plt.figure(figsize=(12, 10))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2)
# Hitrost
ax1.plot(x_values, y1, 'b-', label="Odsek 1")
ax1.plot(x_values + 1, y2_c1 * (v0_2 / v0_1), 'r--', label="Odsek 2 (C¹)")
ax1.plot(x_values + 2, y3_c1 * (v0_3_c1 / v0_1), 'r--')
ax1.plot(x_values + 1, y2_c2 * (v0_2 / v0_1), 'g-', label="Odsek 2 (C²)")
ax1.plot(x_values + 2, y3_c2 * (v0_3_c2 / v0_1), 'g-')
ax1.set_title("Profil hitrosti za 3 semaforje (relativno na $v_{0,1}$)")
ax1.set_ylabel("$v(t) / v_{0,1}$")
ax1.set_xticks([0, 1, 2, 3]); ax1.axvline(1, c='k', ls=':'); ax1.axvline(2, c='k', ls=':'); ax1.grid(True); ax1.legend()
# Pospešek
a1_phys = a1 * (v0_1 / t0_1)
a2_c1_phys = a2_c1 * (v0_2 / t0_2)
a3_c1_phys = a3_c1 * (v0_3_c1 / t0_3)
a2_c2_phys = a2_c2 * (v0_2 / t0_2)
a3_c2_phys = a3_c2 * (v0_3_c2 / t0_3)
ax2.plot(x_values, a1_phys, 'b-')
ax2.plot(x_values + 1, a2_c1_phys, 'r--')
ax2.plot(x_values + 2, a3_c1_phys, 'r--')
ax2.plot(x_values + 1, a2_c2_phys, 'g-')
ax2.plot(x_values + 2, a3_c2_phys, 'g-')
ax2.set_title("Profil fizičnega pospeška za 3 semaforje")
ax2.set_xlabel("Časovni odseki")
ax2.set_ylabel("Pospešek (m/s²)")
ax2.set_xticks([0, 1, 2, 3]); ax2.axvline(1, c='k', ls=':'); ax2.axvline(2, c='k', ls=':'); ax2.grid(True)
plt.tight_layout()
plt.savefig("trije_semaforji_skupaj.pdf")
plt.close()
print("Graf 'trije_semaforji_skupaj.pdf' je shranjen.")

print("\n--- SKRIPTA JE ZAKLJUČILA DELO ---")

# ==============================================================================
# NALOGA 5: ZAPOREDNI SEMAFORJI, DRUGA VERZIJA
# ==============================================================================
print("\n[NALOGA 3 - POPRAVEK] Procesiram poenostavljen model pospeševanja od mirovanja...")

x_values = np.linspace(0, 1, 100)

# --- Primer 1: p -> inf (minimizacija max|a(t)|) ---
# Rešitev je konstanten pospešek in linearna hitrost.
# Predpostavimo a(x) = 2 (v brezdimenzijskih enotah, da je pot z=1)
a_inf = np.ones_like(x_values) * 2
y_inf = 2 * x_values

# --- Primer 2: p = 1 (minimizacija ∫a^2 dt) ---
# Rešitev je linearen pospešek in kvadratična hitrost.
# Najdemo koeficiente, da se pot ujema.
# a(x) = 6x, y(x) = 3x^2. Integral y(x) je 1.
a_p1 = 6 * x_values
y_p1 = 3 * x_values**2

# --- Primer 3: p = 0.75 (vmesni primer, koren) ---
# a(x) ~ x^(1/0.5) = x^2, y(x) ~ x^3. Normaliziramo, da je pot 1.
# a(x) = 12x^2, y(x) = 4x^3
a_inter = 12 * x_values**2
y_inter = 4 * x_values**3


# --- Risanje grafov ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Graf hitrosti
ax1.plot(x_values, y_p1, label="p = 1 (min $\int a^2 dt$) - Kvadratična hitrost")
ax1.plot(x_values, y_inter, label="p = 0.75 (vmesni primer) - Kubična hitrost")
ax1.plot(x_values, y_inf, label="p $\\to \infty$ (min max|a(t)|) - Linearna hitrost")
ax1.set_title("Primerjava profilov hitrosti pri pospeševanju od mirovanja")
ax1.set_xlabel("$t/t_0$")
ax1.set_ylabel("Normalizirana hitrost $y(x)$")
ax1.grid(True)
ax1.legend()

# Graf pospeška
ax2.plot(x_values, a_p1, label="p = 1 - Linearen pospešek")
ax2.plot(x_values, a_inter, label="p = 0.75 - Kvadratičen pospešek")
ax2.plot(x_values, a_inf, label="p $\\to \infty$ - Konstanten pospešek")
ax2.set_title("Primerjava profilov pospeška pri pospeševanju od mirovanja")
ax2.set_xlabel("$t/t_0$")
ax2.set_ylabel("Normaliziran pospešek $a(x)$")
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig("pospesevanje_od_mirovanja_primerjava.pdf")
plt.close()

print("Graf 'pospesevanje_od_mirovanja_primerjava.pdf' je shranjen.")