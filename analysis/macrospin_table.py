"""Table S2 (Supplementary Material 2): analytical macrospin prediction versus simulation.

Modelled:   a_J = hbar J P / (2 e mu0 Ms d)          (Lambda = 1)
            f_mac = gamma a_J / (2 pi alpha)           (independent of B12)
            M_z0 = a_J / [alpha (Ms + H_B)],  H_B = 2|B12| / (mu0 Ms d)  (H_B = 0 without BQ)
Simulated:  f      peak of the FFT spectrum (data/derived/j_vs_freq_<tag>.csv)
            <M_z>  free-layer mean of m_z over the last 3 ns (column mz of data/derived/power_balance_J.csv)
Usage: python macrospin_table.py
"""
import csv
import os

import numpy as np
from paths import DERIVED
from sw_theory import gam, mu0, Ms, d, HB

hbar, e, alpha, P = 1.0546e-34, 1.602e-19, 0.035, 0.35


def a_J(J):
    return hbar * J * P / (2 * e * mu0 * Ms * d)


def table(rows=(("bq", 6.0), ("bq", 12.0), ("no_bq", 5.5), ("no_bq", 10.0))):
    jf = {t: {float(r["J"]): float(r["Frequency"]) for r in csv.DictReader(open(os.path.join(DERIVED, f"j_vs_freq_{t}.csv")))}
          for t in ("bq", "no_bq")}
    pb = {(r["tag"], float(r["J"])): r for r in csv.DictReader(open(os.path.join(DERIVED, "power_balance_J.csv")))}
    print(f"{'case':10s} {'J':>5s} {'f_mac':>6s} {'M_z0':>5s} | {'f_sim':>5s} {'<M_z>':>5s}")
    for tag, J in rows:
        aj = a_J(J * 1e11)
        fmac = gam * aj / (2 * np.pi * alpha) / 1e9
        mz0 = aj / (alpha * (Ms + (HB if tag == "bq" else 0.0)))
        print(f"{tag:10s} {J:5.1f} {fmac:6.1f} {mz0:5.2f} | {jf[tag][J]:5.1f} {float(pb[(tag, J)]['mz']):5.2f}")
    jc = {t: alpha * (Ms + (HB if t == "bq" else 0)) / a_J(1e11) for t in ("bq", "no_bq")}
    print("macrospin reversal current (M_z0 = 1): " + ", ".join(f"{t}: {v:.1f}e11 A/m^2" for t, v in jc.items()))


if __name__ == "__main__":
    table()
