"""Shared loaders, model functions and plotting helpers for the v4 figures.

Conventions
-----------
* magnum.np arrays are indexed [x, y, (z), component].  Every map is shown with `show_map`, which transposes so
  that x is horizontal and y vertical (origin lower).  This fixes the mirrored handedness of the v3 figures.
* Colours: BQ = blue, no BQ = orange.
* The physics functions (dispersion, selection rule) are imported from ../analysis/sw_theory.py so that
  the figures and the analysis use one implementation.
* Data locations come from ../analysis/paths.py (environment variable DATA_DIR, default data/sim); figures are
  written to BQSTO_FIGOUT (default figures/out).
"""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import jn_zeros, jnp_zeros, j1

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                                   # repository root
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from paths import DATA, DERIVED, run_dir                       # noqa: E402
OUT = os.environ.get("BQSTO_FIGOUT", os.path.join(HERE, "out"))
os.makedirs(OUT, exist_ok=True)
from sw_theory import Omega, f_of_k, c_sound, P as Pdip, D as Dex, HB, Ms, gam, mu0, d as dfree   # noqa: E402

DX = 2e-9
ALPHA, POL, HBAR, E = 0.035, 0.35, 1.0546e-34, 1.602e-19

# ----------------------------------------------------------------------------- style
COL = {"bq": "#2a78d6", "no_bq": "#eb6834"}
LAB = {"bq": "with BQ", "no_bq": "without BQ"}
INK, INK2, GRID = "#111111", "#555555", "#d9d8d3"
CM_AMP, CM_PHASE, CM_MZ = "magma", "twilight", "RdBu_r"
plt.rcParams.update({                                          # the author's house style (all work)
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 20,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "mathtext.fontset": "stix",
    "errorbar.capsize": 5,
})
plt.rcParams.update({"figure.dpi": 60, "savefig.dpi": 240, "pdf.fonttype": 42,          # output only, not style
                     "xtick.major.pad": 7, "ytick.major.pad": 7})                 # keeps corner tick labels apart
# Figures are drawn at SCALE x the printed size so that 20 pt text prints at ~8 pt on the APL page; every explicit
# font size, line width and marker size in the scripts is already at this scale.  savefig.dpi = 600 / SCALE.
SCALE = 2.5
COLW, PAGEW = 3.37 * SCALE, 6.69 * SCALE      # APL column / page width (inches, at drawing scale)


def panel_label(ax, s, dx=-0.12, dy=1.04, **kw):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=22.5, fontweight="bold", va="bottom", ha="left", **kw)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    print("saved", name)


# ----------------------------------------------------------------------------- physics helpers
def a_J(J):
    "Slonczewski field amplitude (A/m) for current density J (A/m^2), Lambda = 1 -> eps = P/2"
    return HBAR * J * POL / (2 * E * mu0 * Ms * dfree)


def f_mac(J):
    "macrospin OPP frequency gamma a_J / (2 pi alpha), GHz"
    return gam * a_J(J) / ALPHA / 2 / np.pi / 1e9


def mz_mac(J, bq):
    "macrospin tilt a_J / (alpha (Ms + HB)); reversal where it reaches 1"
    return a_J(J) / (ALPHA * (Ms + (HB if bq else 0.0)))


def ladder(R, bq, nmax=10):
    "plateau frequencies (GHz) for rim antinode (j') and rim node (j) limits"
    return (np.array([f_of_k(x / R, bq) for x in jnp_zeros(1, nmax)]),
            np.array([f_of_k(x / R, bq) for x in jn_zeros(1, nmax)]))


# ----------------------------------------------------------------------------- data loaders
def load_jsweep(tag):
    J, F = [], []
    for r in csv.DictReader(open(f"{DERIVED}/j_vs_freq_{tag}.csv")):
        J.append(float(r["J"])); F.append(float(r["Frequency"]))
    return np.array(J), np.array(F)


def wave_table(sweep=None, tag=None):
    rows = list(csv.DictReader(open(f"{DERIVED}/wave_table.csv")))
    if sweep: rows = [r for r in rows if r["sweep"] == sweep]
    if tag: rows = [r for r in rows if r["tag"] == tag]
    return rows


def power_balance(sweep):
    "sweep = 'J' (power_balance_J.csv) or 'D' (power_balance.csv)"
    fn = "power_balance_J.csv" if sweep == "J" else "power_balance.csv"
    rows = [r for r in csv.DictReader(open(f"{DERIVED}/{fn}")) if not r.get("err")]
    for r in rows:
        if r["tag"] == "no_bq_5.5": r["tag"] = "no_bq"
    return rows


def load_maps(tag, J=None, D=None):
    """Return dict with mz snapshot (final), complex FFT map H at the peak frequency, frequency, mask, R [m].
    Arrays are [x, y]."""
    p = run_dir(tag, J, D)
    fin = np.load(p + "/final.npy")
    mask = np.linalg.norm(fin[:, :, 0, :], axis=-1) > 0.5
    mz = np.where(mask, fin[:, :, 0, 2], np.nan)
    H = np.load(p + "/heat_map.npy")
    H = np.where(mask, np.nan_to_num(H), np.nan + 0j)
    xf = np.load(p + "/freq.npy"); yf = np.load(p + "/fft.npy")
    f = xf[np.argmax(np.abs(yf))] / 1e9
    N = fin.shape[0]
    return dict(mz=mz, H=H, f=f, mask=mask, R=N * DX / 2, N=N)


def azimuthal(H, R, ms=(-3, -2, -1, 0, 1, 2, 3)):
    """Azimuthal decomposition of a complex map H[x, y] (NaN outside the disk).
    Returns radii rb [m], mean amplitude profile, dict m -> complex radial profile, power per m, and k from the
    phase slope of the dominant m over 10-85 % of R."""
    H = np.nan_to_num(H)
    N = H.shape[0]; c = (N - 1) / 2
    x, y = np.mgrid[0:N, 0:N]
    X, Y = (x - c) * DX, (y - c) * DX
    r, phi = np.hypot(X, Y), np.arctan2(Y, X)
    rb = np.arange(0.5 * DX, R, DX)
    idx = np.clip((r / DX).astype(int), 0, len(rb) - 1)
    inside = r < R - DX
    cnt = np.bincount(idx[inside], minlength=len(rb))[:len(rb)]
    Fm = {}
    for m in ms:
        w = H * np.exp(-1j * m * phi)
        re = np.bincount(idx[inside], weights=w[inside].real, minlength=len(rb))[:len(rb)]
        im = np.bincount(idx[inside], weights=w[inside].imag, minlength=len(rb))[:len(rb)]
        Fm[m] = (re + 1j * im) / np.maximum(cnt, 1)
    power = {m: np.sum(np.abs(Fm[m]) ** 2 * rb) for m in ms}
    mdom = max(power, key=power.get)
    amp = np.bincount(idx[inside], weights=np.abs(H[inside]), minlength=len(rb))[:len(rb)] / np.maximum(cnt, 1)
    sel = (rb > 0.10 * R) & (rb < 0.85 * R)
    ph = np.unwrap(np.angle(Fm[mdom][sel]))
    k = np.polyfit(rb[sel], ph, 1)[0]
    return dict(rb=rb, amp=amp, Fm=Fm, power=power, mdom=mdom, k=k, frac=power[mdom] / sum(power.values()))


def show_map(ax, arr, R_nm, cmap, vmin=None, vmax=None, circle=True):
    """imshow of an [x, y] array with x horizontal; returns the image handle."""
    im = ax.imshow(np.asarray(arr).T, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=[-R_nm, R_nm, -R_nm, R_nm], interpolation="nearest")
    if circle:
        ax.add_patch(plt.Circle((0, 0), R_nm, fill=False, lw=1, color=INK2))
    st = 100 if R_nm > 120 else 50                             # ticks inside +-R: labels at +-R collide in the corner
    tk = st * np.arange(-int(0.9 * R_nm // st), int(0.9 * R_nm // st) + 1)
    ax.set_xticks(tk); ax.set_yticks(tk)
    ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
    return im


def cbar(fig, im, ax, label, **kw):
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, **kw)
    cb.set_label(label, fontsize=17.5); cb.ax.tick_params(labelsize=16)
    return cb


def k_of_f(f_GHz, bq):
    "inverse of the self-consistent selection rule: wavenumber (rad/m) whose Omega equals 2 pi f"
    from scipy.optimize import brentq
    return brentq(lambda k: f_of_k(k, bq) - f_GHz, 2e6, 2.3e8)   # f_of_k is defined up to ~45 GHz


def plateau_tops(J, F, jmin=3.5, jmax=8.0, jump=1.0):
    "currents and frequencies just before each upward jump larger than `jump` GHz"
    m = (J >= jmin) & (J <= jmax)
    J, F = J[m], F[m]
    idx = np.where(np.diff(F) > jump)[0]
    idx = np.array([i for n, i in enumerate(idx) if n == 0 or i - idx[n - 1] > 2], dtype=int)
    return J[idx], F[idx]
