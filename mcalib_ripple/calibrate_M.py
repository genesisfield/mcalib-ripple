import numpy as np
from mcalib_ripple.ripple_model import mu0_model

def compute_M_analytic(z, mu_obs, theta, invC):
    mu_model = mu0_model(z, theta)
    ones = np.ones_like(mu_obs)
    delta = mu_obs - mu_model
    M = float(ones @ (invC @ delta)) / float(ones @ (invC @ ones))
    sigma_M = np.sqrt(1.0 / float(ones @ (invC @ ones)))
    return M, sigma_M

def adjust_M_from_residual(mu_obs, z, theta_med, M_locked):
    mu_model = mu0_model(z, theta_med)
    resid = mu_obs - (mu_model + M_locked)
    return float(M_locked + resid.mean())
