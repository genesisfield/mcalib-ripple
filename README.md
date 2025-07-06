# mcalib-ripple  
[![GitHub stars](https://img.shields.io/github/stars/genesisfield/mcalib-ripple?style=social)](https://github.com/genesisfield/mcalib-ripple/stargazers)  
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](#)  
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)  
  
---  
  
### Precision Calibration of Supernova Magnitudes Using Ripple-Modulated Cosmologies and the Pantheon+ Dataset  
  
---  
  
## Overview  
  
`mcalib-ripple` is a modular Python package providing an end-to-end pipeline for calibrating the absolute magnitude (M) of Type Ia supernovae using the Pantheon+SH0ES dataset. It supports flat ΛCDM as well as ripple-modulated cosmologies, enabling precise tests of alternative models featuring oscillatory or phase-modulated expansion.  
  
Ideal for researchers exploring coherent cosmologies, Hubble-tension solutions, or modified expansion scenarios, it reproduces Pantheon+ benchmarks to within 0.1% and is optimized for reproducibility, modularity, and clarity.  
  
---  
  
## Features  
  
- **Analytic M Calibration** — Fast, covariance-aware magnitude estimation  
- **Ripple H(z) Models** — Phase and damping parameters for oscillatory fits  
- **Parallel MCMC Sampling** — Multiprocessed `emcee` integration  
- **Diagnostic Outputs** — Residual plots and detailed logs  
- **Test Coverage** — Validates μ(z = 0.1) ≈ 38.3 mag and d<sub>L</sub>(0.1) ≈ 460 Mpc  
  
---  
  
## Installation  
  
```bash  
git clone https://github.com/genesisfield/mcalib-ripple.git  
cd mcalib-ripple  
pip install -r requirements.txt  
```  
  
---  
  
## Example Usage  
  
```bash  
python examples/calibrate_and_lock_M.py  
```  
  
Results are saved in `outputs/`:  
- `M_values.txt` — Locked and adjusted magnitude values  
- `residuals.png` — Analytic vs. MCMC residual comparison  
  
---  
  
## Testing  
  
```bash  
python -m unittest discover tests  
```  
  
---  
  
## Citations  
  
Please cite this software as:  
  
```  
Greene, R. (2025). mcalib-ripple: Precision Absolute Magnitude Calibration for Ripple-Modulated Cosmologies Using Pantheon+. Genesis Field Institute. DOI: <your-zenodo-doi>  
```  
  
---  
  
## License  
  
Distributed under the MIT License. See [LICENSE](LICENSE) for details.  
  
---  
  
## Acknowledgements  
  
Built on the Pantheon+ and SH0ES datasets and the Python ecosystem:  
NumPy, SciPy, pandas, matplotlib, emcee  
  
---  
  
## References  
  
- Brout et al. (2022). [Cosmology with Pantheon+](https://doi.org/10.3847/1538-4357/ac8e04)  
- Riess et al. (2021). [Local Hubble Constant Measurement](https://doi.org/10.3847/2041-8213/ac5c5b)  
- Scolnic et al. (2018). [The Complete Pantheon SN Sample](https://doi.org/10.3847/1538-4357/aab9bb)  
- Foreman-Mackey et al. (2013). [emcee: The MCMC Hammer](https://doi.org/10.1086/670067)
