import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# --- 1. Simulation Parameters ---
# Ensure this matches the number of particles in your xyz2.dat file
N_particles = 500 

# --- 2. Carnahan-Starling Equation of State (Theoretical) ---
# We want to plot theoretical P values covering the range of our simulation
P_star_values = np.linspace(0.1, 12.0, 100)

def carnahan_starling_equation(eta, P_star):
    """
    Equation to solve: P_star = (6 * eta / pi) * (1 + eta + eta^2 - eta^3) / (1 - eta)^3
    """
    return P_star - (6.0 * eta / np.pi) * (1 + eta + eta**2 - eta**3) / (1 - eta)**3

eta_values_eos = []
eta_guess = 0.1

# Solve for theoretical eta at each pressure
for P in P_star_values:
    eta_sol = fsolve(carnahan_starling_equation, eta_guess, args=(P,))[0]
    eta_values_eos.append(eta_sol)
    eta_guess = eta_sol # Use previous solution as next guess for stability

# --- 3. Load Simulation Data ---
# Skip the first 2 rows (Headers)
try:
    data = np.genfromtxt("pVdata.csv", delimiter=",", skip_header=2)
except FileNotFoundError:
    print("Error: 'pVdata.csv' not found. Make sure you run the C simulation first.")
    exit()

# Extract liquid columns (Columns index 2 and 3)
pressure_liquid = data[:, 2]
mean_vol_liquid = data[:, 3]

# Calculate simulation packing fraction (eta)
# eta = (N * volume_of_one_particle) / Total_Volume
# For 3D spheres of diameter 1: volume_of_one_particle = pi / 6
pf_liquid = (np.pi * N_particles) / (6.0 * mean_vol_liquid)

# --- 4. Plotting ---
plt.figure(figsize=(8, 6))

# Plot the theoretical EoS curve
plt.plot(P_star_values, eta_values_eos, '-', linewidth=2, color='orange', label="Carnahan-Starling EoS")

# Scatter plot for the Monte Carlo simulation data
plt.plot(pressure_liquid, pf_liquid, 'o', markersize=6, color='blue', label="MC Simulation Data")

# Formatting the plot
plt.xlim(0, max(pressure_liquid) + 1)
plt.ylim(0, max(pf_liquid) + 0.05)
plt.xlabel(r"Reduced Pressure ($\beta P \sigma^3$)", fontsize=12)
plt.ylabel(r"Packing Fraction ($\eta$)", fontsize=12)
plt.title("Equation of State: Hard Sphere Fluid Phase", fontsize=14)
plt.legend(loc="lower right", fontsize=11)
plt.grid(True, linestyle='--', alpha=0.6)

# Display the plot
plt.tight_layout()
plt.show()