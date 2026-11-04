"""Table S4 (Supplementary Material 4): plateaus of the diameter sweep.

The plateau ranges below were chosen by hand from Fig. S1 (a new plateau starts where f jumps up and kR drops);
the script only reads the numbers at the first and last diameter of each range from data/derived/wave_table.csv
(f, k, lambda) and data/derived/power_balance.csv (in-plane wave power <a^2>, column aip1).
Usage: python diameter_plateaus.py
"""
import csv
import os

from paths import DERIVED

PLATEAUS = {"bq": [(2, 180, 240), (3, 250, 320), (4, 330, 390), (5, 400, 460), (6, 470, 520)],
            "no_bq": [(2, 180, 220), (3, 230, 330), (4, 340, 390)]}

wt = {(r["tag"], int(r["D_nm"])): r for r in csv.DictReader(open(os.path.join(DERIVED, "wave_table.csv"))) if r["sweep"] == "D"}
pb = {(r["tag"].replace("no_bq_5.5", "no_bq"), int(float(r["D"]))): r
      for r in csv.DictReader(open(os.path.join(DERIVED, "power_balance.csv")))}

print(f"{'case':6s} {'n':>2s} {'D (nm)':>9s} {'f (GHz)':>12s} {'kR':>12s} {'<a^2>':>12s} {'dD':>4s} {'lambda':>6s}")
for tag, rows in PLATEAUS.items():
    for n, d0, d1 in rows:
        a, b = wt[(tag, d0)], wt[(tag, d1)]
        kR = lambda r, D: abs(float(r["k_rad_per_m"])) * D / 2 * 1e-9
        print(f"{tag:6s} {n:2d} {d0:4d}-{d1:4d} {float(a['f_GHz']):5.1f}->{float(b['f_GHz']):4.1f} "
              f"{kR(a, d0):5.1f}->{kR(b, d1):4.1f} {float(pb[(tag, d0)]['aip1']):5.2f}->{float(pb[(tag, d1)]['aip1']):4.2f} "
              f"{d1 - d0 + 10:4d} {float(a['lambda_nm']):6.0f}")
