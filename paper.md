---
title: "mcalib-ripple: Precision Absolute Magnitude Calibration for Ripple-Modulated Cosmologies Using Pantheon+"
authors:
  - name: Richard Greene
    affiliation: 1
affiliations:
  - name: Genesis Field Institute, Chandler, AZ, USA
date: 2025-07-05
---

## Summary

`mcalib-ripple` is a modular Python package that performs analytic and MCMC-based calibration of the absolute magnitude \( M \) of Type Ia supernovae using the Pantheon+SH0ES dataset. It is designed to support the calibration and testing of nonstandard cosmological models, with built-in support for ripple-modulated Hubble expansion models. Ripple cosmologies incorporate coherent modulation or oscillatory features into the cosmological expansion rate, potentially explaining empirical anomalies in standard ΛCDM fits. The software implements reproducible, best-practice techniques for computing \( M_{\text{locked}} \) from SN data, and provides direct compatibility with flat ΛCDM by reducing to it when ripple parameters are disabled.

This package is the first of its kind to support full end-to-end calibration under a ripple-aware cosmology model and serves as the foundation for testing new theoretical frameworks such as the Genesis Field cosmology. The software is openly available at https://github.com/genesisfield/mcalib-ripple, with full documentation, tests, and examples for reproducibility. It provides a clean interface for scientific comparison between ΛCDM, Early Dark Energy, and ripple-modulated models using actual observational data.

`mcalib-ripple` is openly available and permanently archived with a DOI via Zenodo, ensuring persistent access, reproducibility, and reliable citation.

## Statement of Need

There is currently no open-source package that supports ripple-based or phase-modulated cosmologies while also faithfully reproducing standard Pantheon+ ΛCDM benchmarks. This limits the ability of researchers exploring nonstandard cosmologies to test their models using established supernova datasets. `mcalib-ripple` solves this by providing a clean, modular, and verified pipeline for:

- Analytic magnitude calibration using the full Pantheon+ covariance
- MCMC sampling of ripple cosmology parameters with fixed \( M \)
- Residual-based correction of \( M \) for post-fit consistency
- Output of residual plots and fitted \( M \) values to disk

The pipeline is useful both for standard cosmology validation and for exploration of extensions involving coherent phase evolution, oscillating dark energy, or vacuum field modulation. This software also provides a robust foundation for future research in other modulated cosmologies, oscillatory dark energy, or modified gravitational theories.

## Functionality

`mcalib-ripple` includes:

- `ripple_model.py`: Implements ripple-modulated \( H(z) \), \( d_L(z) \), and \( \mu(z) \)
- `calibrate_M.py`: Analytic and residual-corrected M calibration
- `mcmc_engine.py`: emcee-based MCMC sampling
- `examples/calibrate_and_lock_M.py`: Full example that produces `outputs/M_values.txt` and `outputs/residuals.png`
- `tests/`: Unit tests confirming ΛCDM consistency at \( z = 0.1 \)

The calibration pipeline reproduces Pantheon+ baseline values to within 0.1%, including mu(z = 0.1) approximately 38.3 ± 0.02 mag and d_L(z = 0.1) approximately 460 Mpc, consistent with Pantheon+ distance ladder fits and SH0ES calibrations, in agreement with benchmark expectations from Brout et al. (2022) and Riess et al. (2021). The redshift cut \( z > 0.023 \) is applied to filter out supernovae affected by local peculiar velocities, consistent with Pantheon+ best practice (Brout et al. 2022). Reference redshift \( z = 0.1 \) is used for test benchmarking because it lies in the stable central regime of the SN distribution, minimizing both low-z noise and high-z statistical sparsity (Scolnic et al. 2018).

To ensure reproducibility of the MCMC chain and sample initializations, the random seed is fixed at 42 throughout the package. This makes results consistent across different machines and runs.

The `mcmc_engine.py` module uses parallelized emcee sampling with a user-configurable number of processors, defaulting to half the available logical CPUs. This optimization allows rapid exploration of the ripple parameter space while maintaining compatibility with constrained computing environments.

## Acknowledgements

Thanks to the developers and maintainers of the Pantheon+ and SH0ES datasets, and to the open-source scientific Python ecosystem (NumPy, SciPy, pandas, matplotlib, emcee).

## References

- Scolnic, D. et al. (2018). The Complete Light-curve Sample of Spectroscopically Confirmed SNe Ia from Pan-STARRS1 and Cosmological Constraints from the Combined Pantheon Sample. ApJ, 859, 101. DOI: 10.3847/1538-4357/aab9bb
- Brout, D. et al. (2022). Cosmology with the Pantheon+ Type Ia Supernova Sample. ApJ, 938, 110. DOI: 10.3847/1538-4357/ac8e04
- Riess, A. G. et al. (2021). A Comprehensive Measurement of the Local Value of the Hubble Constant. ApJ, 934, L7. DOI: 10.3847/2041-8213/ac5c5b
- Foreman-Mackey, D. et al. (2013). emcee: The MCMC Hammer. PASP, 125, 306. DOI: 10.1086/670067
