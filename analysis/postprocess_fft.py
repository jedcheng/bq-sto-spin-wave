"""FFT post-processing of one or more simulation runs: data.h5 -> freq.npy, fft.npy, fft_avg.npy, heat_map.npy.

Same algorithm as the original scripts (simulations/micromagnetic/original/j_sweep_plot_fft_mx.py and
j_sweep_fft_mz_maps.py), merged into one file and without the diagnostic PNGs:

  freq.npy, fft.npy  rfft of the free-layer m_x averaged over the whole N x N layer (cells outside the disk are zero
                     and only rescale the mean), frames [:-1] of data.h5 (2000 frames = 10 ns, df = 0.1 GHz).
                     The gyration frequency f is the frequency of the largest |fft|.
  fft_avg.npy        |rfft(m_z)| of every free-layer cell, averaged over the layer (spectra of Fig. S1(a,b)).
  heat_map.npy       complex rfft(m_z) map at the gyration frequency f, NaN outside the disk; indexed [x, y].

Usage:  python postprocess_fft.py RUN_DIR [RUN_DIR ...] [--force] [-j 8]
"""
import argparse
import os

import h5py
import numpy as np
from scipy.fft import rfft, rfftfreq, rfftn


def process(p, force=False):
    h5 = os.path.join(p, "data.h5")
    with h5py.File(h5, "r") as f:
        m = f["m"][:-1, :, :, 0, :]          # free layer, all components
        t = f["t"][:-1]
    if force or not os.path.exists(os.path.join(p, "fft.npy")):
        mx = np.mean(m[..., 0], axis=(1, 2))
        xf = rfftfreq(len(mx), t[1] - t[0])
        yf = rfft(mx)
        np.save(os.path.join(p, "fft.npy"), yf)
        np.save(os.path.join(p, "freq.npy"), xf)
    xf = np.load(os.path.join(p, "freq.npy")); yf = np.load(os.path.join(p, "fft.npy"))
    imax = np.argmax(np.abs(yf))
    mz = m[..., 2]
    vacuum = np.where(mz[0] == 0)
    yf_full = rfftn(mz, axes=(0))
    np.save(os.path.join(p, "fft_avg.npy"), np.mean(np.abs(yf_full), axis=(1, 2)))
    heat_map = yf_full[imax]
    heat_map[vacuum] = np.nan
    np.save(os.path.join(p, "heat_map.npy"), heat_map)
    print(f"{p}: f = {xf[imax] / 1e9:.2f} GHz")
    return xf[imax]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--force", action="store_true", help="recompute fft.npy/freq.npy even if present")
    ap.add_argument("-j", type=int, default=1, help="parallel workers")
    a = ap.parse_args()
    if a.j > 1:
        from joblib import Parallel, delayed
        Parallel(n_jobs=a.j)(delayed(process)(p, a.force) for p in a.runs)
    else:
        for p in a.runs:
            process(p, a.force)
