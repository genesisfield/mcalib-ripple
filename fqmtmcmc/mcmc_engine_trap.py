import numpy as np
import emcee
import multiprocessing
from fqmtmcmc.ripple_model_trap import mu_model

def log_prior(theta):
    Om, eps, omega, phi, gamma, H0 = theta
    if (
        0.05 < Om < 0.5 and
        -0.01 < eps < 0.01 and
        0.01 < omega < 0.3 and
        -np.pi < phi < np.pi and
        0.01 < gamma < 0.3 and
        60.0 < H0 < 75.0
    ):
        return 0.0
    return -np.inf

def log_likelihood(theta, z, mu_obs, sigma_mu, M_locked):
    mu_pred = mu_model(z, theta, M_locked)
    if np.any(np.isnan(mu_pred)) or np.any(mu_pred > 200):
        return -np.inf
    return -0.5 * np.sum(((mu_obs - mu_pred) / sigma_mu) ** 2)

def log_probability(theta, z, mu_obs, sigma_mu, M_locked):
    lp = log_prior(theta)
    return lp + log_likelihood(theta, z, mu_obs, sigma_mu, M_locked) if np.isfinite(lp) else -np.inf

def run_mcmc(z, mu_obs, sigma_mu, M_locked, nwalkers=64, nsteps=30000):
    ndim = 6
    pos = np.column_stack([
        np.random.uniform(0.2, 0.4,   nwalkers),
        np.random.uniform(-0.005, 0.005, nwalkers),
        np.random.uniform(0.05, 0.25,  nwalkers),
        np.random.uniform(-1, 1,       nwalkers),
        np.random.uniform(0.05, 0.2,   nwalkers),
        np.random.uniform(67, 71,      nwalkers)
    ])
    print(f"Running MCMC with M_locked = {M_locked}")
    with multiprocessing.Pool() as pool:
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability,
            args=(z, mu_obs, sigma_mu, M_locked),
            pool=pool
        )
        sampler.run_mcmc(pos, nsteps, progress=True)
    return sampler
