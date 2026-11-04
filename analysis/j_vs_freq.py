"""Gyration frequency versus current density (D = 320 nm): reads freq.npy/fft.npy of every run of the current sweep
and writes output/j_vs_freq_<tag>.csv (columns J, Frequency [GHz]), as the original
simulations/micromagnetic/original/j_sweep_plot_fft_mx.py did (J = 1.0 ... 19.9 in steps of 0.1).
Runs without data give a missing row.  Usage: python j_vs_freq.py
"""
import csv
import os

import numpy as np
from paths import run_dir

os.makedirs("output", exist_ok=True)
for tag in ("bq", "no_bq"):
    rows = []
    for J in np.round(np.arange(1.0, 19.95, 0.1), 1):
        p = run_dir(tag, J=J)
        try:
            xf = np.load(os.path.join(p, "freq.npy")); yf = np.load(os.path.join(p, "fft.npy"))
        except OSError:
            continue
        rows.append((J, xf[np.argmax(np.abs(yf))] / 1e9))
    with open(f"output/j_vs_freq_{tag}.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(["J", "Frequency"]); w.writerows(rows)
    print(tag, len(rows), "runs -> output/j_vs_freq_" + tag + ".csv")
