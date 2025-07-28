---
title: "GENESISFIELDMCMC: A Full Inference Pipeline for Ripple-Modulated Cosmologies Using Supernova and H(z) Data"
authors:
  - name: Richard Greene
    orcid: https://orcid.org/0009-0002-2430-8184
    affiliation: 1
affiliations:
  - name: Genesis Field Institute, Chandler, AZ, USA
    index: 1
date: 2025-07-17
bibliography: paper.bib
---

## Summary

`GENESISFIELDMCMC` is a modular Python package for complete cosmological inference using ripple-modulated and ΛCDM expansion models. Ripple models introduce oscillatory features into cosmic expansion driven by quantum coherence decay in the vacuum, potentially resolving empirical tensions within standard cosmology. The software supports analytic and MCMC-based calibration of Type Ia supernova absolute magnitude (\$M\$) using Pantheon+SH0ES data and fits both ripple and ΛCDM models to observational datasets, including cosmic chronometer H(z) measurements.

The package provides a reproducible and falsifiable inference pipeline with full statistical comparisons between standard and ripple-based cosmologies. Although initially developed for the Genesis Field cosmology, it is broadly applicable to general phase-modulated and alternative cosmological frameworks. Users run `python run_pipeline.py` to calibrate data, execute model fits, and produce all outputs automatically, typically in under 45 minutes on a standard desktop workstation. The necessary datasets from Pantheon+ are included directly within the repository.

## Statement of Need

There is currently no open, modular pipeline enabling systematic inference and direct comparison of ripple-based cosmologies against standard ΛCDM using real observational datasets (Pantheon+, SH0ES, cosmic chronometers). `GENESISFIELDMCMC` fills this gap, providing a reproducible framework for testing and exploring both standard and nonstandard cosmological models without extensive retooling.

## Functionality

The software includes:

* `mcalib_ripple/`: Absolute magnitude (\$M\$) calibration using analytic and residual-based corrections.
* `fqmtmcmc/`: Modular ripple cosmology models, likelihood functions, and emcee-based MCMC inference engines.
* `examples/`: Scripts executing SN-only, H(z)-only, and joint-model inferences, including ΛCDM baselines.
* `run_pipeline.py`: Automates the full inference process in a single command.
* `tests/`: Unit and validation tests confirming standard cosmological accuracy.
* `outputs/`: Auto-generated corner plots, posteriors, residual diagnostics, and summary JSON files.

Ripple model equation:

$$
H(z) = H_0 \left[1 + \varepsilon \cos(\omega \ln(1+z) + \phi) e^{-\gamma \ln(1+z)} \right]
$$

The pipeline uses fixed random seeds and parallelized emcee sampling, ensuring reproducibility. ΛCDM fits serve as baseline comparisons, employing grid scans over critical parameters such as \$H\_0\$ and \$\Omega\_m\$ to provide robust statistical benchmarks.

The pipeline supports inference across three regimes:

* SN-only inference
* H(z)-only inference (tight and relaxed parameter priors)
* Joint SN + H(z) inference

Post-fit diagnostic tools include:

* Parameter posterior comparisons
* Two-parameter falsifiability maps
* Full parameter-space sweeps (ε, ω, γ, φ)

Parameter ranges and values used in the two-parameter falsifiability tests are guided by prior theoretical and observational studies on oscillatory cosmologies (Zhao et al., 2012; Hu, 2000; Ringermacher & Mead, 2015).

## Archive

The software is permanently archived on Zenodo: [https://doi.org/10.5281/zenodo.15825901](https://doi.org/10.5281/zenodo.15825901)

## Acknowledgements

We acknowledge maintainers of Pantheon+, SH0ES, and cosmic chronometer datasets, and the open-source Python packages NumPy, SciPy, astropy, emcee, and corner.

## References

* Scolnic, D. et al. (2018). The Complete Light-curve Sample of Spectroscopically Confirmed SNe Ia from Pan-STARRS1 and Cosmological Constraints from the Combined Pantheon Sample. *The Astrophysical Journal*, 859, 101. DOI: [https://doi.org/10.3847/1538-4357/aab9bb](https://doi.org/10.3847/1538-4357/aab9bb)
* Brout, D. et al. (2022). Cosmology with the Pantheon+ Type Ia Supernova Sample. *The Astrophysical Journal*, 938, 110. DOI: [https://doi.org/10.3847/1538-4357/ac8e04](https://doi.org/10.3847/1538-4357/ac8e04)
* Riess, A. G. et al. (2021). A Comprehensive Measurement of the Local Value of the Hubble Constant. *The Astrophysical Journal*, 934, L7. DOI: [https://doi.org/10.3847/2041-8213/ac5c5b](https://doi.org/10.3847/2041-8213/ac5c5b)
* Foreman-Mackey, D. et al. (2013). emcee: The MCMC Hammer. *Publications of the Astronomical Society of the Pacific*, 125, 306. DOI: [https://doi.org/10.1086/670067](https://doi.org/10.1086/670067)
* Ringermacher, H. I., & Mead, L. R. (2015). Observation of discrete oscillations in a model-independent plot of cosmological scale factor vs. lookback time and scalar field model. *Astronomical Journal*, 149(4), 137. DOI: [https://doi.org/10.1088/0004-6256/149/4/137](https://doi.org/10.1088/0004-6256/149/4/137)
* Zhao, G.-B., et al. (2012). Probing modifications of General Relativity using current cosmological observations. *Journal of Cosmology and Astroparticle Physics*, 2012(10), 016. DOI: [https://doi.org/10.1088/1475-7516/2012/10/016](https://doi.org/10.1088/1475-7516/2012/10/016)
* Hu, W. (2000). Structure formation with generalized dark matter. *Physical Review D*, 62(4), 043007. DOI: [https://doi.org/10.1103/PhysRevD.62.043007](https://doi.org/10.1103/PhysRevD.62.043007)
