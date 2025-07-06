import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
from mcalib_ripple.ripple_model import dL_ripple

def test_dL_ripple_at_z_0p1():
    """
    Tests the luminosity distance d_L(z=0.1) for a flat ΛCDM-like cosmology.
    Expected result: ~470 Mpc
    """
    theta = [0.3, 0.0, 0.2, 0.0, 0.1, 70.0]  # LCDM-like ripple-off
    z = np.array([0.1])
    dL = dL_ripple(z, theta)[0]

    assert 450 < dL < 490, f"Expected d_L(0.1) in [450, 490] Mpc, got {dL:.2f} Mpc"
    print(f"✅ test_dL_ripple_at_z_0p1 passed — d_L(0.1) = {dL:.2f} Mpc")

if __name__ == "__main__":
    test_dL_ripple_at_z_0p1()
