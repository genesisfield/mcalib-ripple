# run_pipeline.py

import subprocess

steps = [
    "examples/calibrate_and_lock_M.py",
    "tests/test_mu0_model.py",
    "tests/test_dl_ripple.py",
    "examples/fit_sn_lcdm_grid.py",
    "examples/fit_hz_lcdm_grid_relaxed.py",
    "examples/fit_hz_lcdm_grid_tight.py",
    "examples/fit_joint_lcdm_grid.py",
    "examples/run_fit_pantheon.py",
    "examples/run_fit_hz_tight.py",
    "examples/run_fit_hz_relaxed.py",
    "examples/run_fit_joint.py",
    "examples/plot_ripple_parameter_comparison.py",
    "examples/ripple_vs_lcdm_2param.py",
    "tests/sweep_fqmt_parameters.py"
]

print("\n=== Running full GENESISFIELDMCMC pipeline ===\n")

for i, step in enumerate(steps, 1):
    label = f"[{i}/{len(steps)}] {step}".ljust(60, ".")
    try:
        subprocess.run(f"python {step}", shell=True, check=True, stdout=None, stderr=None)
        print(f"{label} done")
    except subprocess.CalledProcessError:
        print(f"{label} FAILED ❌")
        break

else:
    print("\n✅ Pipeline completed successfully.")
    print("📁 All results are in: outputs/")
    print("📄 Key files: joint_corner.png, ripple_vs_lcdm_2param.png, M_values.txt, posteriors/*.json\n")
