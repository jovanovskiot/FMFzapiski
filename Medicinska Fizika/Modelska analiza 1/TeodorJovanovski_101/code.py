import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import cumulative_trapezoid # Uvozimo numerično integracijo

v0 = 15
L = 200
t0 = 10

L0 = v0 * t0
z = L / L0

#So, to your earlier question:
#
#  If z = 1: L = L0. This means the distance to the light is exactly the distance you would cover if you just kept your initial speed. The optimal strategy would be to do nothing!

#  If z < 1: L < L0. You have to cover less distance than if you kept your speed. This means you must slow down on average.

#  If z > 1 (your case): L > L0. You need to cover more distance than you would at your current speed. You must speed up on average.

#This single number z beautifully captures the entire physical situation. This is the power of dimensionless variables.

def calculate_y(x, z):
    A = -3/2*(z-1)
    B = -2*A
    C = 1
    y = A * x**2 + B * x + C

    return y


# array of x values

x_values = np.linspace(0, 1, 100)

y_values = calculate_y(x_values, z)

plt.figure()
plt.plot(x_values, y_values, label = f"Optimalen profil za z = {z:.2f}")

plt.xlabel("$t/t_0$")
plt.ylabel("$v/v_0$")
#plt.title()

plt.grid(True)
plt.legend()
plt.show()


# TA DEL UPORABIM V POROCILU 
# KONCNA HITROST NI OMEJENA

z_values = np.linspace(0, 2, 10)

plt.figure()

for z in z_values:
    y=calculate_y(x_values, z)
    plt.plot(x_values, y)

plt.xlabel("$t/t_0$")
plt.ylabel("$v/v_0$")
#plt.title()

plt.grid(True)
plt.plot(x_values, y, label=f"z = {z:.2f}")
plt.show()

# DRUGI DEL
# KONCNA HITROST JE OMEJENA
# ZACETNA JE ISTA

def calculate_y_v_final(x, z, y_final):
    A = -6*z + 3*y_final + 3
    B = 6*z - 2*y_final - 4

    return A*x**2 + B*x + 1

y_final_values = np.linspace(1, 2, 10)

plt.figure()

for y_final in y_final_values:
    y=calculate_y_v_final(x_values, z, y_final)
    plt.plot(x_values, y)

plt.xlabel("$t/t_0$")
plt.ylabel("$v/v_0$")
#plt.title()

plt.grid(True)
plt.plot(x_values, y, label=f"y_f = {y_final:.2f}")
plt.show()

# DRUGI DEL
# KONCNA HITROST JE OMEJENA
# ZACETNA NI ISTA

def calculate_y_v_final(x, z, y_final):
    A = -6*z + 3*y_final + 3
    B = 6*z - 2*y_final - 4

    return A*x**2 + B*x + 1

y_final_values = np.linspace(1, 2, 10)

plt.figure()

for y_final in y_final_values:
    for z in z_values:
        y=calculate_y_v_final(x_values, z, y_final)
        plt.plot(x_values, y)

plt.xlabel("$t/t_0$")
plt.ylabel("$v/v_0$")
#plt.title()

plt.grid(True)
plt.legend()
plt.show()


# DRUGI DEL
# KONCNA HITROST JE OMEJENA
# ZACETNA NI ISTA

z_values_widespread = [0.5, 1.0, 1.5, 2.0]  # Štiri različne situacije
y_final_cases = [0.8, 1.0, 1.2, 1.5] # Štiri različne zahteve za končno hitrost

for y_f_case in y_final_cases:
    plt.figure(figsize=(8, 6))
    for z_case in z_values_widespread:
        y = calculate_y_v_final(x_values, z_case, y_f_case)
        plt.plot(x_values, y, label=f"z = {z_case:.1f}")
    
    plt.title(f"Optimalni profili za zahtevano končno hitrost y_f = {y_f_case}")
    plt.xlabel("$t/t_0$")
    plt.ylabel("$v/v_0$")
    plt.grid(True)
    plt.legend()
    # Shranimo vsak graf v svojo datoteko
    plt.savefig(f"profili_yf_{str(y_f_case).replace('.', '_')}.pdf")
    plt.close() # Zapremo figuro, da se ne rišejo vsi grafi en čez drugega

print("Dodatni grafi za 2. nalogo so shranjeni.")

# DRUGI DEL
# KONCNA HITROST JE OMEJENA, ISTA ZA VSE
# ZACETNA JE ISTA

def calculate_y_v_final(x, z, y_final):
    A = -6*z + 3*y_final + 3
    B = 6*z - 2*y_final - 4

    return A*x**2 + B*x + 1


plt.figure()
z_values = np.linspace(0, 2, 10)

for z in z_values:
    y=calculate_y_v_final(x_values, z, y_final = 1)
    plt.plot(x_values, y)

plt.xlabel("$t/t_0$")
plt.ylabel("$v/v_0$")
#plt.title()

plt.grid(True)
plt.plot(x_values, y, label=f"z = {z:.2f}")
plt.show()


# TRETJI DEL

linear_base = np.linspace(1, -1, 100)

a_p1 = linear_base  **1
a_p2 = np.sign(linear_base)*np.abs(linear_base)**(1/3)
a_p10 = np.sign(linear_base)*np.abs(linear_base)**(1/19)

plt.figure()
plt.plot(x_values, a_p1, label="p=1 (minimize a²)")
plt.plot(x_values, a_p2, label="p=2 (minimize a⁴)")
plt.plot(x_values, a_p10, label="p=10 (minimize a²⁰)")

plt.xlabel("$t/t_0$")
plt.ylabel("Normaliziran pospešek a(t)")
plt.legend()
plt.grid(True)
plt.show()

# TRETJI DEL
# VISJE POTENCE

# --- POPOLNA ANALIZA ZA 3. NALOGO (POSPEŠEK IN HITROST) ---
print("\nGeneriram grafe za 3. nalogo...")



linear_base = np.linspace(1, -1, 100)
x_values = np.linspace(0, 1, 100)

a_p1 = linear_base**1
a_p2 = np.sign(linear_base) * np.abs(linear_base)**(1/3)
a_p10 = np.sign(linear_base) * np.abs(linear_base)**(1/19)

# Numerično integriramo pospeške, da dobimo spremembo v hitrosti
# cumtrapz potrebuje dx, ki je 1/število korakov
y_p1 = cumulative_trapezoid(a_p1, x_values, initial=0)
y_p2 = cumulative_trapezoid(a_p2, x_values, initial=0)
y_p10 = cumulative_trapezoid(a_p10, x_values, initial=0)

# Normaliziramo profile hitrosti, da se vsi končajo pri enaki prevoženi poti
# (To je poenostavitev, a dobro ilustrira obliko)
y_p1_norm = 1 + y_p1 - np.mean(y_p1) + np.mean(a_p1)*0.1 # Dodan popravek za simetrijo
y_p2_norm = 1 + y_p2 - np.mean(y_p2) + np.mean(a_p2)*0.1
y_p10_norm = 1 + y_p10 - np.mean(y_p10) + np.mean(a_p10)*0.1

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

print("Grafi za 3. nalogo so shranjeni.")

# CETRTI DEL
# OMEJITEV HITROSTI

def calculate_y_penalty(x, z, gamma):
    if gamma < 1e-6:
        return calculate_y(x, z)
    sqrt_g = np.sqrt(gamma)
    A = (z-1)/(np.tanh((gamma**(1/2))) / (gamma**(1/2)) - 1)
    B = -A * np.tanh(np.sqrt(gamma))
    D = 1 - A

    return A*np.cosh(sqrt_g*x) + B* np.sinh(sqrt_g*x) + D

z_case = 1.5

gamma_cases = [0, 0.1, 1, 10, 100]

plt.figure()

for g in gamma_cases:
    y_profile = calculate_y_penalty(x_values, z_case, g)
    plt.plot(x_values, y_profile, label=f"γ = {g}")

plt.xlabel("$t/t_0$")
plt.ylabel("$v/v_0$")
plt.title(f"Effect of Velocity Penalty (γ) for z={z_case}")
plt.grid(True)
plt.legend()
plt.show()

# --- DODATNA ANALIZA ZA 4. NALOGO ---
# VEC VREDNOSTI Z
print("\nGeneriram dodatne grafe za 4. nalogo...")

z_scenarios = [0.8, 1.5] # Scenarij zaviranja in scenarij pospeševanja
gamma_values = np.logspace(-1, 2, 12) # 12 vrednosti od 0.1 do 100

for z_case in z_scenarios:
    plt.figure(figsize=(10, 7))
    for g in gamma_values:
        y_profile = calculate_y_penalty(x_values, z_case, g)
        plt.plot(x_values, y_profile, label=f"γ = {g:.1f}")
    
    plt.title(f"Vpliv kazni za hitrost (γ) za scenarij z = {z_case}")
    plt.xlabel("$t/t_0$")
    plt.ylabel("$v/v_0$")
    plt.grid(True)
    plt.legend(ncol=2)
    plt.savefig(f"profili_penalty_z_{str(z_case).replace('.', '_')}.pdf")
    plt.close()

print("Dodatni grafi za 4. nalogo so shranjeni.")

# PETI DEL

def get_profile(z):
    # racuna koeficiente A, B, C za optimalno voznjo za dan z. 

    C = 1
    A = (3/2) * (1-z)
    B = -2*A

    return A, B, C

z_test = 1.33
A_coeff, B_coeff, C_coeff = get_profile(z_test)

def evaluate_profile(A, B, C, x_values):

    y_values = A* x_values**2 + B*x_values + C

    a_values = 2*A*x_values + B

    return y_values, a_values

# testiram metodo
z_1 = 1.5
A1, B1, C1 = get_profile(z_1)
x1_values = np.linspace(0, 1, 100)
y1_profile, a1_profile = evaluate_profile(A1, B1, C1, x1_values)

print("--- Leg 1 ---")
print(f"Starting velocity y(0): {y1_profile[0]:.2f}") # Should be 1.0
print(f"Ending velocity y(1): {y1_profile[-1]:.2f}")

print(f"Starting acceleration a(0): {a1_profile[0]:.2f}")
print(f"Ending acceleration a(1): {a1_profile[-1]:.2f}") # Should be 0.0 for this model

# Definiram pot z dvema semaforjema

#prvi semafor

v0_1 = 15
L1 = 200
t0_1 = 10

z1 = L1 / (v0_1 * t0_1)

#drugi semafor

L2 = 300
t0_2 = 12

A1, B1, C1 = get_profile(z1)

x1_values = np.linspace(0, 1, 100)

y1_profile, a1_profile = evaluate_profile(A1, B1, C1, x1_values)

#definiram koncna hitrost pri prvem semaforju
y1_final = y1_profile[-1]
v0_2 = y1_final * v0_1

z2 = L2 / (v0_2 * t0_2)

A2, B2, C2 = get_profile(z2)
x2_values = np.linspace(0, 1, 100)

y2_profile, a2_profile = evaluate_profile(A2, B2, C2, x2_values)

# plot

# --- Plotting (Corrected) ---
fig = plt.figure(figsize=(12, 10))
ax1 = fig.add_subplot(2, 1, 1) # Velocity
ax2 = fig.add_subplot(2, 1, 2) # Acceleration

# --- Plot 1: Dimensionless Velocity (Corrected) ---
# We use v0_1 as the reference for the whole journey.
ax1.plot(x1_values, y1_profile, label="Leg 1")

# Calculate Leg 2 velocity relative to v0_1
y2_global_scale = y2_profile * (v0_2 / v0_1)
ax1.plot(x2_values + 1, y2_global_scale, label="Leg 2")

ax1.set_title("Dimensionless Velocity (v / v0_1)")
# ... (rest of the velocity plotting code is the same)
ax1.set_xlabel("Dimensionless Time (Journey Segments)")
ax1.set_ylabel("v(t) / v0_1")
ax1.set_xticks([0, 1, 2]) 
ax1.axvline(1, color='red', linestyle='--', label="Semafor 1")
ax1.grid(True)
ax1.legend()


# --- Plot 2: Dimensionless Acceleration ---
# Acceleration plot is still correct and shows the jump.
ax2.plot(x1_values, a1_profile, label="Leg 1")
ax2.plot(x2_values + 1, a2_profile, label="Leg 2")
ax2.set_title("Dimensionless Acceleration (alpha = dy/dx)")
# ... (rest of the acceleration plotting code is the same)
ax2.set_xlabel("Dimensionless Time (Journey Segments)")
ax2.set_ylabel("alpha(x)")
ax2.set_xticks([0, 1, 2])
ax2.axvline(1, color='red', linestyle='--', label="Semafor 1")
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()


# Zvezen pospesek

# --- Toolkit for C1 (Jerky) Path ---

def get_profile(z):
    """Calculates coefficients for the C1 quadratic path y=Ax^2+Bx+C."""
    C = 1
    A = 1.5 * (1 - z)
    B = -2 * A
    return A, B, C

def evaluate_profile(A, B, C, x_values):
    """Evaluates the C1 quadratic path and its derivative (acceleration)."""
    y_values = A * x_values**2 + B * x_values + C
    a_values = 2 * A * x_values + B
    return y_values, a_values

# --- Toolkit for C2 (Smooth) Path ---

def get_c2_coefficients(z, alpha_start):
    """Calculates coefficients for the C2 cubic path y=Ax^3+Bx^2+Cx+D."""
    D = 1.0
    C = alpha_start
    A = 4 * (1 - z + alpha_start / 3.0)
    B = - (3 * A + alpha_start) / 2.0
    return A, B, C, D

def evaluate_cubic_profile(A, B, C, D, x_values):
    """Evaluates the C2 cubic path and its derivative (acceleration)."""
    y_values = A * x_values**3 + B * x_values**2 + C * x_values + D
    a_values = 3 * A * x_values**2 + 2 * B * x_values + C
    return y_values, a_values

# --- 1. Simulation Setup ---

# Leg 1 Parameters
v0_1 = 15.0
L1 = 200.0
t0_1 = 10.0

# Leg 2 Parameters
L2 = 300.0
t0_2 = 12.0

# --- 2. Leg 1 Calculation (Same for both scenarios) ---
z1 = L1 / (v0_1 * t0_1)
A1, B1, C1 = get_profile(z1)
x_values = np.linspace(0, 1, 100)
y1_profile, a1_profile = evaluate_profile(A1, B1, C1, x_values)

# --- 3. Link Conditions ---
y1_final = y1_profile[-1]
v0_2 = y1_final * v0_1
z2 = L2 / (v0_2 * t0_2)
# The final acceleration of Leg 1 is the crucial link for the C2 case
alpha1_final = a1_profile[-1] 
# NOTE: For the C1 model, alpha1_final is 0, but we use the variable for generality.

# --- 4. Leg 2 Calculation (Jerky C1 Path) ---
A2_c1, B2_c1, C2_c1 = get_profile(z2)
y2_profile_c1, a2_profile_c1 = evaluate_profile(A2_c1, B2_c1, C2_c1, x_values)

# --- 5. Leg 2 Calculation (Smooth C2 Path) ---
# We need to provide the required starting acceleration
alpha2_start = alpha1_final * (t0_2 / t0_1) * (v0_1 / v0_2) # Rescale acceleration to new dimensionless units
A2_c2, B2_c2, C2_c2, D2_c2 = get_c2_coefficients(z2, alpha2_start)
y2_profile_c2, a2_profile_c2 = evaluate_cubic_profile(A2_c2, B2_c2, C2_c2, D2_c2, x_values)


# --- 6. Plotting the Comparison ---
fig = plt.figure(figsize=(12, 10))
ax1 = fig.add_subplot(2, 1, 1) # Velocity
ax2 = fig.add_subplot(2, 1, 2) # Acceleration

# --- Velocity Plot (using global v0_1 reference) ---
ax1.plot(x_values, y1_profile, color='blue', label="Leg 1")
# Rescale both Leg 2 velocities to be relative to v0_1
y2_c1_global = y2_profile_c1 * (v0_2 / v0_1)
y2_c2_global = y2_profile_c2 * (v0_2 / v0_1)
ax1.plot(x_values + 1, y2_c1_global, color='orange', linestyle='--', label="Leg 2 (Jerky C¹ Path)")
ax1.plot(x_values + 1, y2_c2_global, color='green', linestyle='-', label="Leg 2 (Smooth C² Path)")
ax1.set_title("Velocity Profile (v / v0_1)")
ax1.set_xlabel("Dimensionless Time (Journey Segments)")
ax1.set_ylabel("v(t) / v0_1")
ax1.axvline(1, color='red', linestyle=':')
ax1.grid(True)
ax1.legend()

# --- Acceleration Plot ---
# Rescale accelerations to a common physical reference (m/s^2) for true comparison
a1_physical = a1_profile * (v0_1 / t0_1)
a2_c1_physical = a2_profile_c1 * (v0_2 / t0_2)
a2_c2_physical = a2_profile_c2 * (v0_2 / t0_2)
ax2.plot(x_values, a1_physical, color='blue', label="Leg 1")
ax2.plot(x_values + 1, a2_c1_physical, color='orange', linestyle='--', label="Leg 2 (Jerky C¹ Path)")
ax2.plot(x_values + 1, a2_c2_physical, color='green', linestyle='-', label="Leg 2 (Smooth C² Path)")
ax2.set_title("Physical Acceleration Profile")
ax2.set_xlabel("Dimensionless Time (Journey Segments)")
ax2.set_ylabel("Acceleration (m/s²)")
ax2.axvline(1, color='red', linestyle=':')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()

# --- RAZŠIRJENA SIMULACIJA ZA 5. NALOGO (3 SEMAFORJI) ---
print("\nGeneriram grafe za 5. nalogo (3 semaforji)...")

# --- Parametri ---
v0_1 = 15.0; L1 = 200.0; t0_1 = 10.0
L2 = 300.0; t0_2 = 12.0
L3 = 250.0; t0_3 = 8.0

# --- Leg 1 ---
z1 = L1 / (v0_1 * t0_1)
A1, B1, C1 = get_profile(z1)
x_values = np.linspace(0, 1, 100)
y1, a1 = evaluate_profile(A1, B1, C1, x_values)
a1_phys = a1 * (v0_1 / t0_1)

# --- Leg 2 ---
v0_2 = y1[-1] * v0_1
z2 = L2 / (v0_2 * t0_2)
alpha1_final_dimless = a1[-1]
# C1
A2_c1, B2_c1, C2_c1 = get_profile(z2)
y2_c1, a2_c1 = evaluate_profile(A2_c1, B2_c1, C2_c1, x_values)
a2_c1_phys = a2_c1 * (v0_2 / t0_2)
# C2
alpha2_start_dimless = alpha1_final_dimless * (t0_2 / t0_1) * (v0_1 / v0_2)
A2_c2, B2_c2, C2_c2, D2_c2 = get_c2_coefficients(z2, alpha2_start_dimless)
y2_c2, a2_c2 = evaluate_cubic_profile(A2_c2, B2_c2, C2_c2, D2_c2, x_values)
a2_c2_phys = a2_c2 * (v0_2 / t0_2)

# --- Leg 3 ---
# C1 pot
v0_3_c1 = y2_c1[-1] * v0_2
z3_c1 = L3 / (v0_3_c1 * t0_3)
A3_c1, B3_c1, C3_c1 = get_profile(z3_c1)
y3_c1, a3_c1 = evaluate_profile(A3_c1, B3_c1, C3_c1, x_values)
a3_c1_phys = a3_c1 * (v0_3_c1 / t0_3)
# C2 pot
v0_3_c2 = y2_c2[-1] * v0_2
z3_c2 = L3 / (v0_3_c2 * t0_3)
alpha2_final_dimless = a2_c2[-1]
alpha3_start_dimless = alpha2_final_dimless * (t0_3 / t0_2) * (v0_2 / v0_3_c2)
A3_c2, B3_c2, C3_c2, D3_c2 = get_c2_coefficients(z3_c2, alpha3_start_dimless)
y3_c2, a3_c2 = evaluate_cubic_profile(A3_c2, B3_c2, C3_c2, D3_c2, x_values)
a3_c2_phys = a3_c2 * (v0_3_c2 / t0_3)

# --- Plotting ---
fig = plt.figure(figsize=(12, 10))
ax1 = fig.add_subplot(2, 1, 1) # Velocity
ax2 = fig.add_subplot(2, 1, 2) # Acceleration

# Velocity Plot
ax1.plot(x_values, y1, label="Odsek 1")
ax1.plot(x_values + 1, y2_c1 * (v0_2 / v0_1), '--', label="Odsek 2 (C¹)")
ax1.plot(x_values + 2, y3_c1 * (v0_3_c1 / v0_1), '--', label="Odsek 3 (C¹)")
ax1.plot(x_values + 1, y2_c2 * (v0_2 / v0_1), label="Odsek 2 (C²)")
ax1.plot(x_values + 2, y3_c2 * (v0_3_c2 / v0_1), label="Odsek 3 (C²)")
ax1.set_title("Profil hitrosti za 3 semaforje (v / v0_1)")
ax1.set_xticks([0, 1, 2, 3]); ax1.grid(True); ax1.legend()

# Acceleration Plot
ax2.plot(x_values, a1_phys)
ax2.plot(x_values + 1, a2_c1_phys, '--')
ax2.plot(x_values + 2, a3_c1_phys, '--')
ax2.plot(x_values + 1, a2_c2_phys)
ax2.plot(x_values + 2, a3_c2_phys)
ax2.set_title("Profil pospeška za 3 semaforje")
ax2.set_ylabel("Pospešek (m/s²)")
ax2.set_ylabel("Pospešek (m/s²)")
ax2.set_xticks([0, 1, 2, 3]); ax2.grid(True)

plt.savefig("trije_semaforji_skupaj.pdf")
plt.close()

print("Graf za 5. nalogo (3 semaforji) je shranjen.")