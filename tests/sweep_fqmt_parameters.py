#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from itertools import product
import json
from collections import defaultdict

# Add project path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from fqmtmcmc.utils import load_hz_data

# === Load H(z) data ===
z, Hz_obs, sigma = load_hz_data(include_farooq=True, deduplicate=True)
n = len(z)

# === Model definition ===
def Hz_ripple(z, eps, omega, phi, gamma, H0, Om):
    r = eps * np.exp(-gamma * z) * np.cos(omega * z + phi)
    r0 = eps * np.cos(phi)
    return H0 * (1 + r) / (1 + r0) * np.sqrt(Om * (1 + z)**3 + (1 - Om))

def chi2_for_params(eps, omega, phi, gamma, Om):
    def chi2_H0(H0):
        model = Hz_ripple(z, eps, omega, phi, gamma, H0, Om)
        return np.sum(((Hz_obs - model) / sigma) ** 2)
    res = minimize(chi2_H0, x0=[70.0], bounds=[(66.0, 74.0)])
    return res.fun, res.x[0]  # best χ², best H₀

# === Sweep settings ===
eps_vals = np.linspace(0.0, 0.3, 25)
phi_vals = np.linspace(0.0, np.pi, 9)  # NEW: φ now free!
omega_vals = [0.145, 0.150, 0.155, 0.160, 0.165]
gamma_vals = [0.10, 0.12, 0.15]
Om_fixed = 0.4

# === Run full (ε, φ) sweep at each (ω, γ) ===
results = []
best_map = {}

for gamma, omega in product(gamma_vals, omega_vals):
    chi2_list = []
    for eps in eps_vals:
        best_phi = None
        best_chi2 = np.inf
        best_H0 = None
        for phi in phi_vals:
            chi2, H0 = chi2_for_params(eps, omega, phi, gamma, Om_fixed)
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_phi = phi
                best_H0 = H0
        chi2_list.append(best_chi2)
        results.append({
            "gamma": gamma,
            "omega": omega,
            "eps": eps,
            "phi": best_phi,
            "chi2": best_chi2,
            "H0": best_H0
        })
        # Save best overall for heatmap
        key = (gamma, omega)
        if key not in best_map or best_chi2 < best_map[key]["chi2"]:
            best_map[key] = {
                "eps": eps,
                "phi": best_phi,
                "chi2": best_chi2
            }

    # Plot χ² vs ε for this (omega, gamma)
    plt.plot(eps_vals, chi2_list, label=f"ω={omega}, γ={gamma}")

# === Plot χ² vs ε ===
output_dir = os.path.join(os.path.dirname(__file__), "..", "outputs")
os.makedirs(output_dir, exist_ok=True)
plt.xlabel("ε (Ripple Amplitude)")
plt.ylabel("χ²")
plt.title("Ripple Fit χ² vs ε (Optimized φ) across (ω, γ)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "sweep_epsilon_grid_chi2.png"))
plt.close()

# === Build heatmap ===
sorted_gamma = sorted(set(g for g, _ in best_map))
sorted_omega = sorted(set(w for _, w in best_map))
chi2_matrix = np.array([
    [best_map[(g, w)]["chi2"] for w in sorted_omega]
    for g in sorted_gamma
])

fig, ax = plt.subplots(figsize=(8, 5))
c = ax.imshow(
    chi2_matrix,
    cmap="viridis",
    origin="lower",
    extent=[min(sorted_omega), max(sorted_omega), min(sorted_gamma), max(sorted_gamma)],
    aspect="auto"
)
plt.colorbar(c, label="Min χ² at best (ε, φ)")
ax.set_xlabel("ω (Ripple Frequency)")
ax.set_ylabel("γ (Damping Rate)")
ax.set_title("Best χ² across (ω, γ) (Optimized φ, ε)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "sweep_epsilon_grid_heatmap.png"))
plt.close()

# === Save results ===
with open(os.path.join(output_dir, "sweep_epsilon_grid_results.json"), "w") as f:
    json.dump(results, f, indent=2)

print("✅ Saved sweep_epsilon_grid_results.json and updated plots.")
