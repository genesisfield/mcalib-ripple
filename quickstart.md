# Quickstart: GENESISFIELDMCMC

This guide walks you through installing and running the full ripple-modulated cosmology pipeline using Pantheon+SH0ES and H(z) data.

## 1. Clone the Repository

git clone https://github.com/genesisfield/genesisfieldmcmc.git  
cd genesisfieldmcmc

## 2. Create a Clean Environment

Windows (PowerShell):

python -m venv venv  
.\venv\Scripts\activate

macOS/Linux:

python3 -m venv venv  
source venv/bin/activate

## 3. Install Requirements

pip install -r requirements.txt

## 4. Run the Full Pipeline

This will:
- Calibrate SN magnitude M
- Run SN-only, H(z)-only, and joint fits
- Fit ripple and ΛCDM models
- Output chains, plots, and comparison JSONs

python examples/run_pipeline.py

## 5. View Results

All results are saved in the outputs/ directory:

Output Type       → Example File  
----------------- → -----------------------------  
M calibration     → M_values.txt, residuals.png  
Corner plots      → sn_corner_mcmc_pantheon.png  
Residuals         → sn_residuals_pantheon.png  
Posterior JSONs   → *_summary.json  
Comparison plots  → ripple_vs_lcdm_2param.png

## 6. Run Unit Tests (Optional)

python -m unittest discover tests

## 7. Troubleshooting

- Ensure you're using Python 3.9 or higher  
- Run from the repo root: `python examples/run_pipeline.py`  
- All datasets (Pantheon+, SH0ES, H(z)) are included in the repo

## Citation

If you use this software, please cite:

@software{genesisfield2025mcmc,  
  author = {Greene, Richard},  
  title = {GENESISFIELDMCMC: Ripple-Modulated Cosmological Inference},  
  year = {2025},  
  doi = {10.5281/zenodo.15825901},  
  url = {https://github.com/genesisfield/genesisfieldmcmc}  
}
