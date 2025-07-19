import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
import emcee

c = 299792.458  # speed of light in km/s

def ripple_Hz(z_array, Om, eps, omega, phi, gamma, H0):
    r = eps * np.exp(-gamma * z_array) * np.cos(omega * z_array + phi)
    r0 = eps * np.cos(phi)
    norm = (1 + r) / (1 + r0)
    return H0 * norm * np.sqrt(Om * (1 + z_array)**3 + (1 - Om))

def dL_ripple(z_array, Om, eps, omega, phi, gamma, H0):
    z_grid = np.linspace(0, 2.5, 400)
    Hz_grid = ripple_Hz(z_grid, Om, eps, omega, phi, gamma, H0)
    if np.any(np.isnan(Hz_grid)) or np.any(Hz_grid <= 0):
        return np.full_like(z_array, np.nan)
    D_C = cumulative_trapezoid(c / Hz_grid, z_grid, initial=0.0)
    D_L = (1 + z_grid) * D_C
    interp_fn = interp1d(z_grid, D_L, kind='cubic', fill_value='extrapolate')
    return interp_fn(z_array)

def mu_model(z, Om, eps, omega, phi, gamma, H0, M_locked):
    dL = dL_ripple(z, Om, eps, omega, phi, gamma, H0)
    return 5 * np.log10(dL) + 25 + M_locked

# === Combined log-likelihood for SN + H(z) ===
def log_likelihood_joint(theta, z_sn, mu_obs, sigma_mu, z_hz, Hz_obs, sigma_Hz, M_locked):
    Om, eps, omega, phi, gamma, H0 = theta
    mu_pred = mu_model(z_sn, Om, eps, omega, phi, gamma, H0, M_locked)
    Hz_pred = ripple_Hz(z_hz, Om, eps, omega, phi, gamma, H0)
    if np.any(np.isnan(mu_pred)) or np.any(np.isnan(Hz_pred)) or np.any(Hz_pred <= 0):
        return -np.inf
    chi2_sn = np.sum(((mu_obs - mu_pred) / sigma_mu) ** 2)
    chi2_hz = np.sum(((Hz_obs - Hz_pred) / sigma_Hz) ** 2)
    return -0.5 * (chi2_sn + chi2_hz)

# === Relaxed prior ===
def log_prior(theta):
    Om, eps, omega, phi, gamma, H0 = theta
    if (
        0.1 < Om < 0.5 and
        -0.02 < eps < 0.02 and
        0.01 < omega < 1.00 and
        -np.pi < phi < np.pi and
        0.01 < gamma < 0.5 and
        60.0 < H0 < 77.0
    ):
        return 0.0
    return -np.inf

def log_probability(theta, z_sn, mu_obs, sigma_mu, z_hz, Hz_obs, sigma_Hz, M_locked):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood_joint(theta, z_sn, mu_obs, sigma_mu, z_hz, Hz_obs, sigma_Hz, M_locked)
    return lp + ll if np.isfinite(ll) else -np.inf

# === Unified MCMC runner ===
def run_joint_mcmc(z_sn, mu_obs, sigma_mu, z_hz, Hz_obs, sigma_Hz, M_locked, nwalkers=64, nsteps=30000):
    ndim = 6
    pos = np.column_stack([
        np.random.uniform(0.10, 0.50,  nwalkers),  # Om
        np.random.uniform(-0.01, 0.01, nwalkers),  # eps
        np.random.uniform(0.05, 0.25,  nwalkers),  # omega
        np.random.uniform(-1, 1,       nwalkers),  # phi
        np.random.uniform(0.05, 0.25,  nwalkers),  # gamma
        np.random.uniform(67.0, 71.5,  nwalkers)   # H0
    ])

    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability,
        args=(z_sn, mu_obs, sigma_mu, z_hz, Hz_obs, sigma_Hz, M_locked)
    )
    print(f"Running joint MCMC (SN + H(z)) with M_locked = {M_locked:.5f}")
    sampler.run_mcmc(pos, nsteps, progress=True)
    return sampler
