import numpy as np
import emcee
import multiprocessing

from mcalib_ripple.ripple_model import mu0_model

def log_probability(theta, z, mu_obs, invC, M_locked):
    Om, eps, omega, phi, gamma, H0 = theta
    if not (0.05   < Om    < 0.5):    return -np.inf
    if not (-0.01  < eps   < 0.01):   return -np.inf
    if not (0.01   < omega < 0.3):    return -np.inf
    if not (-np.pi < phi   < np.pi):  return -np.inf
    if not (0.01   < gamma < 0.3):    return -np.inf
    if not (60.0   < H0    < 75.0):   return -np.inf

    mu_model = mu0_model(z, theta) + M_locked
    delta = mu_obs - mu_model
    return -0.5 * delta @ invC @ delta

def run_mcmc_theta_fit(log_prob_fn, ndim, nwalkers, nsteps, args, seed=42, nprocs=None):
    """
    Run emcee with optional multiprocessing (Windows-compatible)

    Parameters:
        log_prob_fn : function
        ndim : int
        nwalkers : int
        nsteps : int
        args : tuple, args to log_prob_fn
        seed : int
        nprocs : int or None (auto half of CPU)

    Returns:
        flat MCMC chain (array)
    """
    np.random.seed(seed)

    p0 = np.column_stack([
        np.random.uniform(0.2,   0.4,   nwalkers),
        np.random.uniform(-0.005, 0.005, nwalkers),
        np.random.uniform(0.05,  0.25,  nwalkers),
        np.random.uniform(-1,     1,    nwalkers),
        np.random.uniform(0.05,  0.2,   nwalkers),
        np.random.uniform(67,    70,    nwalkers),
    ])

    if nprocs is None:
        nprocs = max(1, multiprocessing.cpu_count() // 2)

    with multiprocessing.Pool(processes=nprocs) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob_fn, args=args, pool=pool)
        sampler.run_mcmc(p0, nsteps, progress=True)

    return sampler.get_chain(discard=nsteps // 2, flat=True)
