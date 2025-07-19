import numpy as np
from scipy.integrate import cumulative_trapezoid

c = 299792.458  # km/s

def ripple_Hz(z, Om, eps, omega, phi, gamma, H0):
    r = eps * np.exp(-gamma * z) * np.cos(omega * z + phi)
    r0 = eps * np.cos(phi)
    norm = (1 + r) / (1 + r0)
    if np.any(norm <= 0.5):
        return np.inf
    return H0 * norm * np.sqrt(Om * (1 + z)**3 + (1 - Om))

def dL_ripple(z, Om, eps, omega, phi, gamma, H0):
    zp = np.linspace(0, np.max(z), 5000)
    Hz = ripple_Hz(zp, Om, eps, omega, phi, gamma, H0)
    integral = cumulative_trapezoid(c / Hz, zp, initial=0)
    return (1 + z) * np.interp(z, zp, integral)

def mu_model(z, theta, M_locked):
    Om, eps, omega, phi, gamma, H0 = theta
    dL = dL_ripple(z, Om, eps, omega, phi, gamma, H0)
    return 5 * np.log10(dL) + 25 + M_locked
