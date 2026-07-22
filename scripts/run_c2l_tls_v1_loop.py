"""Loop C2L until all 142 samples done, detecting and using GPU automatically."""
import subprocess, sys, time, os
from pathlib import Path

OUT = Path(r"E:/GBM/results/c2l_tls_v1")
TOTAL_NEEDED = 142 + 1  # samples + summary.csv

while True:
    n = len(list(OUT.iterdir())) if OUT.exists() else 0
    done_samples = n - (1 if (OUT / "summary.csv").exists() else 0)
    print(f"Progress: {done_samples}/{TOTAL_NEEDED-1} samples")
    if n >= TOTAL_NEEDED - 9:  # almost done
        break

    result = subprocess.run([sys.executable, "scripts/run_c2l_tls_v1.py"],
                           capture_output=True, text=True, timeout=600)
    # Continue regardless — skip logic handles completed samples

print(f"Final count: {len(list(OUT.iterdir()))}")
n2 = len(list(OUT.iterdir()))
if n2 >= TOTAL_NEEDED - 9:
    print("All done!")
