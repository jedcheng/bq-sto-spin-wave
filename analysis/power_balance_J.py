"""Exact LLG power balance and wave content for the current sweep (D = 320 nm), same analysis as power_balance.py.
Output: power_balance_J.csv (column J in 1e11 A/m^2)."""
import numpy as np, csv, sys
import os
import power_balance as pb
from paths import run_dir
Js = [round(x, 1) for x in np.arange(3.0, 10.01, 0.1)] + [round(x, 1) for x in np.arange(10.5, 20.01, 0.5)]
tasks = [(tag, J) for tag in ("bq", "no_bq") for J in Js]
def safe(tag, J):
    try:
        r = pb.run(run_dir(tag, J=J), J*1e11, f"{J:.1f}"); r["tag"] = tag; r["J"] = J; return r
    except Exception as ex:
        return dict(tag=tag, J=J, err=str(ex))
if __name__ == "__main__":
    from joblib import Parallel, delayed
    out = Parallel(n_jobs=12)(delayed(safe)(*t) for t in tasks)
    keys = ["tag","J","f","sense","bal","W","mip2","mz","coh","Wz1","Wip1","Wz2","Wip2","Wlow","Whigh","az1","aip1","aip_static","fmac","err"]
    os.makedirs("output", exist_ok=True)
    with open("output/power_balance_J.csv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in out: w.writerow({k: r.get(k, "") for k in keys})
    print("done", len(out))
