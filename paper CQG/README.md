# CQG Submission Alignment

📄 [Genesis_Field_CQG_Submission_July23.pdf](./Genesis_Field_CQG_Submission_July23.pdf)

This paper was submitted to *Classical and Quantum Gravity* (Ref: CQG-113295) on July 23, 2025.  
It is fully aligned with this repository and the GENESISFIELDMCMC pipeline currently under review at the *Journal of Open Source Software* (JOSS).

All figures, tables, MCMC diagnostics, and statistical comparisons in the manuscript are reproducible using the code provided here.

To regenerate results:
```bash
python run_pipeline.py
```
This executes the complete inference workflow, including SN fits, H(z) analysis, joint ripple vs. ΛCDM comparisons, and statistical summaries (RMS, AIC, BIC, χ²).  
Outputs will appear in the `/outputs/` directory and are designed to match the figures and tables in the paper exactly.

This structure ensures full transparency between the theoretical model, empirical results, and the published record.
