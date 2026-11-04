"""Fig. 2: birth of the wave, the rim launches it (D = 320 nm; J = 6.0 with BQ, 5.5 without).
Rows: with BQ / without BQ.
(a,b) M_z snapshot at t = 0.1 ns: the wave exists only in a ring at the rim
(c,d) M_z snapshot at t = 1.0 ns: the standing pattern is established
(e,f) radius-time map of the m = +1 azimuthal component of M_z for 0 <= t <= 0.8 ns: tilted stripes = inward
      travelling wave; dots = arrival of the front in each radial bin, dashed = linear fit, speed quoted against
      the predicted sound speed c = gamma sqrt(D_ex (Ms + HB))
All maps use the same colour map (CM_MZ); column titles on top, each column shares one colour bar (at the bottom);
the rows are labelled on the left."""
from common import *
import h5py

RUNS = (("bq", 6.0), ("no_bq", 5.5))
T_EARLY, T_LATE, T_MAP = 0.10, 1.00, 0.80          # ns
NFR = 230                                           # frames read (5 ps each) -> 1.15 ns


def radial_m1(mz, mask, r, phi, dr=4e-9, rmax=162e-9):
    """m = +1 azimuthal component of m_z(t, x, y) in radial bins; returns bin centres and A1[t, r]."""
    rb = np.arange(0, rmax, dr); rc = 0.5 * (rb[1:] + rb[:-1])
    A1 = np.empty((mz.shape[0], len(rc)), complex)
    for i in range(len(rc)):
        s = mask & (r >= rb[i]) & (r < rb[i + 1])
        A1[:, i] = np.mean(mz[:, s] * np.exp(-1j * phi[s]), axis=1)
    return rc, A1


# fixed layout in inches: two square snapshots, then a wider radius-time map, one horizontal colour bar per column at the bottom
S, WM = 4.1, 5.4                                    # snapshot side, map width
X0 = (1.6, 1.6 + S + 0.35, 1.6 + 2 * S + 0.35 + 1.1) # left edges of the three columns
CBY = 0.85                                           # colour bars, below the x labels of the bottom row
Y0 = (CBY + 1.05 + S + 0.5, CBY + 1.05)              # bottom edges of the two rows
W, H = PAGEW, Y0[0] + S + 0.75                     # column titles above the panel letters
fig = plt.figure(figsize=(W, H))
box = lambda x, y, w, h: fig.add_axes([x / W, y / H, w / W, h / H])
WID = (S, S, WM)
labels = [["(a)", "(c)", "(e)"], ["(b)", "(d)", "(f)"]]   # lettered down the columns (chronological)
TITLE = (r"$M_z$ at $t$ = %.2f ns" % T_EARLY, r"$M_z$ at $t$ = %.2f ns" % T_LATE, r"Re$[m=+1$ component of $M_z]$")
for col in range(3):
    fig.text((X0[col] + WID[col] / 2) / W, (Y0[0] + S + 0.3) / H, TITLE[col], ha="center", va="baseline", fontsize=20)


def letter(ax, row, col, dx):
    "panel letter; the top row sits on the baseline of the column titles"
    if row == 0:
        ax.text(dx, 1 + 0.3 / S, labels[row][col], transform=ax.transAxes, fontsize=22.5, fontweight="bold", va="baseline", ha="left")
    else:
        panel_label(ax, labels[row][col], dx=dx, dy=1.02)


CB = {0: (r"$M_z$", [0, 0.1, 0.2, 0.3, 0.4, 0.5]),
      1: (r"$M_z$", [-0.2, 0, 0.2, 0.4, 0.6, 0.8]),
      2: (r"$M_z$ ($m=+1$)", [-0.08, -0.04, 0, 0.04, 0.08])}


def colour_bar(im, col):
    cax = box(X0[col], CBY, WID[col], 0.16)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=CB[col][1])
    cax.xaxis.set_ticks_position("bottom"); cax.xaxis.set_label_position("bottom")
    cb.set_label(CB[col][0], fontsize=20, labelpad=8); cax.tick_params(labelsize=16, length=5)


for row, (tag, J) in enumerate(RUNS):
    bq = tag == "bq"
    with h5py.File(f"{run_dir(tag, J=J)}/data.h5") as f:
        m = f["m"][:NFR, :, :, 0, :]; t = f["t"][:NFR]
    n = m.shape[1]; R = n * 1e-9; Rn = R * 1e9
    x = (np.arange(n) + 0.5 - n / 2) * 2e-9; X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.hypot(X, Y); phi = np.arctan2(Y, X); mask = np.linalg.norm(m[0], axis=-1) > 0.5
    mz = np.where(mask, m[..., 2], np.nan); tn = t * 1e9
    bottom = row == 1

    # snapshots
    for tt, col, (vmin, vmax) in ((T_EARLY, 0, (0.0, 0.5)), (T_LATE, 1, (-0.3, 0.9))):
        i = int(np.argmin(np.abs(tn - tt)))
        ax = box(X0[col], Y0[row], S, S)
        im = show_map(ax, mz[i], Rn, CM_MZ, vmin=vmin, vmax=vmax)
        if bottom: ax.set_xlabel("x (nm)")
        else: ax.tick_params(labelbottom=False)
        if col == 0: ax.set_ylabel("y (nm)")
        else: ax.tick_params(labelleft=False)
        if row == 0: colour_bar(im, col)
        letter(ax, row, col, -0.9 / S if col == 0 else 0.0)
    # row label, left of the y axis label of (a) / (d)
    fig.text((X0[0] - 1.45) / W, (Y0[row] + S / 2) / H, f"{LAB[tag]} ($J$ = {J:.1f})", rotation=90, ha="center", va="center",
             color=COL[tag], fontsize=20)

    # radius-time map of the m = +1 component and the arrival front
    rc, A1 = radial_m1(np.nan_to_num(mz), mask, r, phi)
    env = np.abs(A1); thr = 0.25 * np.median(env[150:, :])
    tarr = np.array([t[np.argmax(env[:, i] > thr)] for i in range(len(rc))])
    sel = (rc > 20e-9) & (rc < 140e-9)
    p = np.polyfit(tarr[sel], rc[sel], 1); v_front = -p[0]
    ax = box(X0[2], Y0[row], WM, S)
    im = ax.pcolormesh(tn, rc * 1e9, np.real(A1).T, shading="auto", cmap=CM_MZ, vmin=-0.1, vmax=0.1, rasterized=True)
    ax.plot(tarr[sel] * 1e9, rc[sel] * 1e9, "o", color=INK, ms=3.5, mec="none")
    tf = np.linspace(0, T_MAP * 1e-9, 50); ax.plot(tf * 1e9, np.polyval(p, tf) * 1e9, "--", color=INK, lw=2)
    ax.set_xlim(0, T_MAP); ax.set_ylim(0, Rn); ax.set_ylabel("r (nm)")
    ax.set_yticks([50, 100, 150])                   # no 0: it would sit on the t = 0.0 label in the corner
    if bottom: ax.set_xlabel("t (ns)")
    else: ax.tick_params(labelbottom=False)
    ax.text(0.97, 0.96, f"front {v_front:.0f} m/s\n$c$ = {c_sound(bq):.0f} m/s", transform=ax.transAxes, ha="right", va="top",
            color=INK, fontsize=16, linespacing=1.3, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=INK2, lw=1, alpha=0.9))
    if row == 0: colour_bar(im, 2)
    letter(ax, row, 2, -0.9 / WM)
    print(f"{tag} J={J}: front speed {v_front:.0f} m/s, c = {c_sound(bq):.0f} m/s, rim response {tarr[rc > 140e-9].min()*1e9:.2f} ns, "
          f"centre arrival {tarr[rc < 20e-9].min()*1e9:.2f} ns")
save(fig, "fig2")
