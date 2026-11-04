"""Direct observation of the launch: the m=+1 wave leaves the rim at switch-on and travels inward.
For each radial bin, t_arr(r) = first time |A1(r,t)| exceeds a threshold, where A1 is the m=+1 azimuthal
component of m_z. A linear fit of r(t_arr) gives the inward front speed, compared with the sound speed
c = gamma sqrt(D (Ms+HB)) and the group velocity at the final k*. Also saves launch_rt.png (r-t maps)."""
import sys, numpy as np, h5py, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sw_theory import Omega, c_sound
from paths import run_dir
runs = (("bq", 6.0, 13.09), ("no_bq", 5.5, 11.93))
fig, axs = plt.subplots(2, 3, figsize=(14, 7))
for row, (tag, J, f0) in enumerate(runs):
    bq = tag == "bq"
    with h5py.File(os.path.join(run_dir(tag, J=J), "data.h5")) as f:
        m = f["m"][:300, :, :, 0, :]; t = f["t"][:300]
    n = m.shape[1]; x = (np.arange(n) + 0.5 - n / 2) * 2e-9; X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.hypot(X, Y); phi = np.arctan2(Y, X); mask = np.linalg.norm(m[0], axis=-1) > 0.5
    mz = m[..., 2]
    rb = np.arange(0, 162e-9, 4e-9); rc = 0.5 * (rb[1:] + rb[:-1])
    A1 = np.array([np.mean(mz[:, mask & (r >= rb[i]) & (r < rb[i + 1])] * np.exp(-1j * phi[mask & (r >= rb[i]) & (r < rb[i + 1])]), axis=1)
                   for i in range(len(rb) - 1)]).T
    env = np.abs(A1); thr = 0.25 * np.median(env[150:, :])          # quarter of the established level
    tarr = np.array([t[np.argmax(env[:, i] > thr)] for i in range(len(rc))])
    sel = (rc > 20e-9) & (rc < 140e-9)                                  # exclude the centre (regularity) and the driven rim
    p = np.polyfit(tarr[sel], rc[sel], 1); v_front = -p[0]
    # theory: sound speed and group velocity at the final k*
    w0 = 2 * np.pi * f0 * 1e9; kk = np.linspace(1e6, 2e8, 4000); Om = Omega(kk, w0, bq)
    ks = kk[np.argmin(np.abs(Om - w0))]; vg = np.gradient(Om, kk)[np.argmin(np.abs(Om - w0))]; vph = w0 / ks
    print(f"{tag} J={J}: front speed {v_front:.0f} m/s | c = {c_sound(bq):.0f} m/s | at k* = {ks/1e6:.0f} rad/um: v_g = {vg:.0f}, v_ph = {vph:.0f} m/s "
          f"| rim->centre arrival {tarr[rc<20e-9].min()*1e9:.2f} ns, first rim response {tarr[rc>140e-9].min()*1e9:.2f} ns")
    tn = t * 1e9; rn = rc * 1e9
    ax = axs[row, 0]; im = ax.pcolormesh(tn, rn, np.real(A1).T, shading="auto", cmap="RdBu_r", vmin=-0.1, vmax=0.1); ax.set_xlim(0, 0.8)
    ax.plot(tarr[sel] * 1e9, rc[sel] * 1e9, "k.", ms=3); ax.plot(tn, np.polyval(p, t) * 1e9, "k--", lw=1)
    ax.set_ylim(0, 160); ax.set_title(f"{tag} J={J}: Re[m=+1](r,t); front {v_front:.0f} m/s, c = {c_sound(bq):.0f} m/s"); plt.colorbar(im, ax=ax)
    ax = axs[row, 1]; im = ax.pcolormesh(tn, rn, env.T, shading="auto", cmap="magma"); ax.set_xlim(0, 1.5); ax.set_title("|m=+1|(r,t)"); plt.colorbar(im, ax=ax)
    ax = axs[row, 2]
    for ti in (0.1, 0.15, 0.2, 0.3, 0.6, 1.4):
        i = np.argmin(np.abs(tn - ti)); ax.plot(rn, env[i], label=f"{tn[i]:.2f} ns")
    ax.set_xlabel("r (nm)"); ax.set_ylabel("|m=+1|"); ax.legend(fontsize=7); ax.set_title("envelope snapshots")
    for a in axs[row, :2]: a.set_xlabel("t (ns)"); a.set_ylabel("r (nm)")
plt.tight_layout(); os.makedirs("output", exist_ok=True); out = "output/launch_rt.png"; plt.savefig(out, dpi=110); print("saved", out)
