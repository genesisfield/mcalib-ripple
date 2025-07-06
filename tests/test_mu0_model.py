import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
from mcalib_ripple.ripple_model import mu0_model

def test_mu0_model_at_z_0p1():
    """
    Tests the distance modulus μ(z=0.1) under ΛCDM-like parameters.
    Expected result: ~38.3 mag
    """
    theta = [0.3, 0.0, 0.2, 0.0, 0.1, 70.0]  # Standard LCDM-like setup
    z = np.array([0.1])
    mu = mu0_model(z, theta)[0]

    assert 37.5 < mu < 39.0, f"μ(z=0.1) = {mu:.2f} mag is outside expected range"
    print(f"✅ test_mu0_model_at_z_0p1 passed — μ(z=0.1) = {mu:.2f} mag")

if __name__ == "__main__":
    test_mu0_model_at_z_0p1()
