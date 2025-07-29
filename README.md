> 📌 **NOTE TO REVIEWERS**
>
> This repository contains the **official GENESISFIELDMCMC v1.1.4** release, under simultaneous review at:
> - 📝 *Classical and Quantum Gravity* (CQG) — Theory + Results Paper
> - 🧪 *Journal of Open Source Software* (JOSS) — Software Pipeline Submission
>
> 🔗 All figures, tables, and statistical results in the CQG paper are generated from this codebase.  
> 🔁 All JOSS functionality and documentation reflects this version.  
> 📄 Final CQG Submission PDF: [📥 Download `Genesis_Field_CQG_Submission_July23.pdf`](https://github.com/genesisfield/genesisfieldmcmc/raw/main/paper%20CQG/Genesis_Field_CQG_Submission_July23.pdf)
> 📝 Research Notes of the AAS (RNAAS) — Empirical Suppression Result (submitted July 28, 2025)
> 📎 Zenodo DOI: [10.5281/zenodo.16251890](https://doi.org/10.5281/zenodo.16251890)  
> 🏷️ Version: `GENESISFIELDMCMC v1.1.4` (Released July 21, 2025)

# GENESISFIELDMCMC: Ripple-Modulated Cosmological Inference in the Genesis Field Framework
[![GitHub stars](https://img.shields.io/github/stars/genesisfield/genesisfieldmcmc?style=social)](https://github.com/genesisfield/genesisfieldmcmc/stargazers)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](#)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15825900.svg)](https://doi.org/10.5281/zenodo.15825900)

GENESISFIELDMCMC is a full-sequence cosmological inference pipeline for ripple-modulated expansion models derived from the Genesis Field framework. The ripple model arises from vacuum coherence decay described by a generalized Gross–Pitaevskii equation. In its constrained limit, the model reduces exactly to ΛCDM. This repository implements and tests the model using Pantheon+SH0ES supernova data and cosmic chronometer H(z) measurements.

The pipeline is falsifiable, empirical, and fully reproducible. All stages pass outputs (M values, JSON summaries, residuals) forward. Deduplication and fixed seeds ensure reproducibility.

## 🗂 outputs/

- `joint_model_comparison_summary.json`  
- `joint_corner.png`  
- `joint_lcdm_grid_fit_summary.json`  
- `joint_lcdm_grid_fit_chi2.png`

These are the canonical output files referenced in the RNAAS notes and the CQG manuscript for reviwer convenience. All values and figures can be regenerated using the provided pipeline and configuration files.

> **Note:** Raw `.npy` chain files and `.log` outputs are not included in this repository, but are automatically generated when the pipeline is run locally. See the [Usage Instructions](#running-the-pipeline) for details.

Quickstart:
To run the full pipeline in one step:

python run_pipeline.py

This will:
- Calibrate absolute magnitude M
- Run ΛCDM baseline fits for SN, H(z), and joint
- Fit ripple models to each dataset and combined
- Generate posterior plots, residuals, comparisons
- Save everything to the outputs/ folder

Step-by-step usage:

Step 1 — Calibrate and Lock M  
Run: python examples/calibrate_and_lock_M.py  
Computes analytic M from Pantheon+SH0ES and applies post-fit residual correction. Produces M_values.txt and residuals.png.

Step 2 — Sanity Tests  
Run: python tests/test_mu0_model.py and python tests/test_dl_ripple.py  
Confirms μ(z) and d_L(z) accuracy around z ≈ 0.1.

Step 3 — SN ΛCDM Grid Fit Relaxed  
Run: python examples/fit_sn_lcdm_grid.py  
Fits flat ΛCDM to SN data. Outputs include sn_lcdm_grid_fit_chi2.png and sn_lcdm_grid_fit_summary.json.

Step 4 — H(z) ΛCDM Grid Fit Relaxed 
Run: python examples/fit_hz_lcdm_grid_relaxed.py  
Fits flat ΛCDM to H(z) using a grid scan of H₀. Produces hz_lcdm_grid_relaxed_fit_chi2.png and summary.

Step 5 — H(z) ΛCDM Grid Fit Tight
Run: python examples/fit_hz_lcdm_grid_tight.py  
Fits flat ΛCDM to H(z) using a grid scan of H₀. Produces hz_lcdm_grid_tight_fit_chi2.png and summary.

Step 6 — Joint ΛCDM Grid Fit  
Run: python examples/fit_joint_lcdm_grid.py  
Performs ΛCDM fit on SN + H(z) combined. Outputs joint_lcdm_grid_fit_chi2.png and joint_lcdm_grid_fit_summary.json for model comparison.

Step 7 — SN Ripple Fit  
Run: python examples/run_fit_pantheon.py  
Fits ripple model to Pantheon+SH0ES. Produces sn_corner_mcmc_pantheon.png, posteriors, residuals, and summary.

Step 8 — H(z) Ripple Fit (Tight)  
Run: python examples/run_fit_hz_tight.py  
Tests SN-calibrated ripple parameters on H(z) data.

Step 9 — H(z) Ripple Fit (Relaxed)  
Run: python examples/run_fit_hz_relaxed.py  
Fully free ripple fit on H(z). Outputs hz_corner_relaxed.png and model comparison stats.

Step 10 — Joint Ripple Fit  
Run: python examples/run_fit_joint.py  
Final ripple fit across SN + H(z). Produces joint_corner.png and joint_model_comparison_summary.json.

Step 11 — Ripple vs ΛCDM Parameter Comparison  
Run: python examples/plot_ripple_parameter_comparison.py  
Visualizes ripple and ΛCDM posterior differences across SN, H(z), and joint fits.

Step 12 — Two-Parameter Ripple vs ΛCDM  
Run: python examples/ripple_vs_lcdm_2param.py  
Compares ripple vs ΛCDM in ε–ω space. Produces ripple_vs_lcdm_2param_hero.png and ripple_vs_lcdm_2param_summary.json.

Step 13 — Sweep Parameter Space (Post-Fit Exploration)  
Run: python examples/sweep_fqmt_parameters.py  
Sweeps ripple parameters ε, ω, γ, and φ over physical ranges. Produces heatmaps and chi² diagnostics to visualize ripple structure.

Ripple Model:
H(z) = H₀ × [1 + ε cos(ω ln(1+z) + φ) × exp(−γ ln(1+z))]  
where ε = amplitude, ω = frequency, φ = phase offset, γ = damping, and H₀ = Hubble constant

Directory Overview:
examples/ — Pipeline scripts  
mcalib_ripple/ — M-calibration logic  
fqmtmcmc/ — Ripple models and MCMC engines  
data/ — Pantheon+SH0ES SN and H(z) inputs  
tests/ — Sanity tests for distance and μ(z)  
outputs/ — All generated posteriors, plots, and JSONs

Outputs:
All results are saved in outputs/. Key files include:

Calibration:
- M_values.txt  
- residuals.png  

Posteriors & Chains:
- sn_chain_mcmc_pantheon.npy  
- hz_chain_mcmc_tight.npy / relaxed.npy  
- joint_chain_mcmc.npy  
- *_log_prob_mcmc.npy  

Corner Plots:
- sn_corner_mcmc_pantheon.png  
- hz_corner_tight.png / relaxed.png  
- joint_corner.png  

Residuals & Comparisons:
- sn_residuals_pantheon.png  
- ripple_vs_lcdm_2param.png  
- ripple_vs_lcdm_hz_2param_hero.png  
- ripple_param_comparison_barplot.png  
- sweep_epsilon_grid_heatmap.png  

Model Summaries:
- sn_mcmc_pantheon_summary.json  
- hz_model_comparison_tight.json / relaxed.json  
- joint_model_comparison_summary.json  
- ripple_vs_lcdm_2param_summary.json  

ΛCDM Grids:
- sn_lcdm_grid_fit_chi2.png / summary.json  
- hz_lcdm_grid_relaxed_fit_chi2.png / summary.json 
- hz_lcdm_grid_tight_fit_chi2.png / summary.json   
- joint_lcdm_grid_fit_chi2.png / summary.json  

Install:
git clone https://github.com/genesisfield/genesisfieldmcmc.git  
cd genesisfieldmcmc  
pip install -r requirements.txt

Run tests:
python -m unittest discover tests

Citation:
@software{genesisfield2025mcmc,  
author = {Greene, Richard},  
title = {GENESISFIELDMCMC: Ripple-Modulated Cosmological Inference},  
year = {2025},  
doi = {10.5281/zenodo.15825901},  
url = {https://github.com/genesisfield/genesisfieldmcmc}}

@article{genesisfield2025theory,  
author = {Greene, Richard},  
title = {Quantum Coherence and Ripple Structure in the Genesis Field},  
journal = {Classical and Quantum Gravity (submitted)},  
year = {2025}}

References:
[1] Brout et al. (2022), Pantheon+ SN Compilation: https://doi.org/10.3847/1538-4357/ac8e04  
[2] Riess et al. (2021), SH0ES Local H₀ Measurement: https://doi.org/10.3847/2041-8213/ac5c5b  
[3] Ringermacher & Mead (2015), Cosmic Oscillation Evidence: https://doi.org/10.1063/1.4907966  
[4] Zhao et al. (2012), Oscillating Expansion Models: https://doi.org/10.1088/1475-7516/2012/10/016  
[5] Hu (2000), Generalized Dark Energy Framework: https://doi.org/10.1103/PhysRevD.62.043007

Acknowledgements:
Built on the Pantheon+SH0ES SN dataset, cosmic chronometer H(z) observations, and open-source scientific tools: numpy, scipy, matplotlib, pandas, astropy, emcee, corner. The ripple intuition is informed by decades of theoretical and observational work on oscillating cosmology, vacuum coherence, and time-domain structure in expansion history.

This pipeline doesn't assume ripples. It lets the data decide. If they exist, this model reveals them. If not, it returns ΛCDM exactly.
