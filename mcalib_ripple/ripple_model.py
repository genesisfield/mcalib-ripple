import numpy as np
from scipy.integrate import quad

def ripple_Hz(zp, Om, eps, omega, phi, gamma, H0):
    r    = eps * np.exp(-gamma * zp) * np.cos(omega * zp + phi)
    r0   = eps * np.cos(phi)
    norm = (1 + r) / (1 + r0)
    if np.any(norm <= 0.5):
        return np.inf
    return H0 * norm * np.sqrt(Om * (1 + zp)**3 + (1 - Om))

def dL_ripple(zarr, theta):
    Om, eps, omega, phi, gamma, H0 = theta
    c = 299792.458  # km/s

    def integrand(zp):
        return c / ripple_Hz(zp, Om, eps, omega, phi, gamma, H0)

    return np.array([(1 + zi) * quad(integrand, 0, zi)[0] for zi in zarr])

def mu0_model(zarr, theta):
    return 5 * np.log10(dL_ripple(zarr, theta)) + 25
