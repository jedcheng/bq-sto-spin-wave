"""Fig. S1: size dependence (diameter sweep D = 100-700 nm; J = 6.0 with BQ, 5.5 without), read with the cavity model.
Columns: with BQ / without BQ.
(a,b) time-averaged FFT power spectrum (log) versus D with the simulated STO frequency and the cavity modes
      Omega(j_{1,n}/R)/2pi, R = D/2 (theory, no fitted parameter)
(c,d) cavity phase kR (k from the measured radial phase slope) with the Bessel zeros j_{1,n} (rim node, solid) and
      j'_{1,n} (rim antinode, dotted); open symbols: pattern not a clean m = +1 spiral
(e,f) in-plane wave power <a^2> versus D
(the f versus <a^2> collapse is in main-text Fig. 4(d), J and D sweeps fitted together)"""
from common import *
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

CASES = (("bq", 6.0), ("no_bq", 5.5))
DMIN_FIT = 180                                        # below this the pattern is not a clean cavity mode
fig = plt.figure(figsize=(PAGEW, 11.7))
gs = fig.add_gridspec(3, 2, height_ratios=[1.45, 1, 0.8], hspace=0.2, wspace=0.08, left=0.08, right=0.95, top=0.94, bottom=0.075)
LET = []                                              # (axes, letter, column); letters placed at the end
fit = {}
for c, (tag, J) in enumerate(CASES):
    bq = tag == "bq"
    rows = sorted(wave_table("D", tag), key=lambda r: float(r["D_nm"]))
    Ds = np.array([float(r["D_nm"]) for r in rows]); f = np.array([float(r["f_GHz"]) for r in rows])
    k = np.abs(np.array([float(r["k_rad_per_m"]) for r in rows])); mf = np.array([float(r["m_frac"]) for r in rows])
    P = sorted([r for r in power_balance("D") if r["tag"] == tag], key=lambda r: float(r["D"]))
    DP = np.array([float(r["D"]) for r in P]); a2 = np.array([float(r["aip1"]) for r in P])
    fP = np.array([float(r["f"]) for r in P]); fm = float(P[0]["fmac"])

    # (a,b) spectra with the cavity fan
    ax = fig.add_subplot(gs[0, c])
    xf = np.load(f"{run_dir(tag, D=Ds[0])}/freq.npy") / 1e9
    spec = np.log10(np.array([np.load(f"{run_dir(tag, D=D)}/fft_avg.npy") for D in Ds]).T)
    lo = xf < 25; mid, top = np.median(spec[lo]), np.percentile(spec[lo], 99.7)      # white = noise floor, red = power
    im = ax.imshow(spec, origin="lower", aspect="auto", cmap=CM_MZ, extent=[Ds[0] - 5, Ds[-1] + 5, xf[0], xf[-1]],
              vmin=2 * mid - top, vmax=top, interpolation="nearest", rasterized=True)
    DD = np.linspace(Ds[0], Ds[-1], 300)
    for n in range(1, 11):
        fn = np.array([f_of_k(jn_zeros(1, n)[-1] / (D / 2 * 1e-9), bq) for D in DD]); s = (fn < 24.5) & (fn > 5)
        ax.plot(DD[s], fn[s], color=INK, lw=1.5, ls=(0, (5, 3)))
        if s.any() and DD[s][0] < 600:                  # mode number just left of where the curve enters
            if fn[s][0] > 23 and DD[s][0] < Ds[0] + 15: ax.text(DD[s][0] + 6, 24.4, f"{n}", fontsize=14, color=INK, va="top", ha="left")
            elif fn[s][0] > 23: ax.text(DD[s][0] - 3, 24.4, f"{n}", fontsize=14, color=INK, va="top", ha="right")
            else: ax.text(DD[s][0] + 4, fn[s][0] + 0.4, f"{n}", fontsize=14, color=INK, va="bottom", ha="left")
    ax.plot(Ds, f, "o", ms=6.5, mfc="white", mec=INK, mew=1.5)
    ax.set_xlim(Ds[0] - 5, Ds[-1] + 5); ax.set_ylim(0, 25); ax.set_xticklabels([])
    ax.set_title(f"{LAB[tag]} ($J$ = {J:.1f})", color=COL[tag], pad=10)
    if c == 0:
        ax.set_ylabel("frequency (GHz)")
        ax.legend([Line2D([], [], ls="none", marker="o", ms=6.5, mfc="white", mec=INK, mew=1.5),
                   Line2D([], [], color=INK, lw=1.5, ls=(0, (5, 3)))],
                  ["STO frequency (simulation)", r"cavity mode $n$: $\Omega(j_{1,n}/R)/2\pi$, no fit"], loc="lower right", fontsize=15,
                  frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, borderpad=0.3, handlelength=1.8, borderaxespad=0.3)
    else:
        ax.set_yticklabels([])
        ax.apply_aspect(); pos = ax.get_position()
        cax = fig.add_axes([pos.x1 + 0.008, pos.y0 + 0.08 * pos.height, 0.01, 0.84 * pos.height])
        cb = fig.colorbar(im, cax=cax, ticks=[]); cb.set_label("FFT power (log)", fontsize=16, labelpad=5)
        cax.text(0.5, 1.02, "high", transform=cax.transAxes, ha="center", va="bottom", fontsize=14, color=INK2)
        cax.text(0.5, -0.02, "low", transform=cax.transAxes, ha="center", va="top", fontsize=14, color=INK2)
    LET.append((ax, ["(a)", "(b)"][c], c))

    # (c,d) cavity phase kR
    ax = fig.add_subplot(gs[1, c])
    for n in range(1, 11):
        ax.axhline(jn_zeros(1, n)[-1], color=GRID, lw=2); ax.axhline(jnp_zeros(1, n)[-1], color=GRID, lw=1, ls=":")
        if c == 1 and jn_zeros(1, n)[-1] < 33: ax.text(Ds[-1] + 12, jn_zeros(1, n)[-1], f"$j_{{1,{n}}}$", fontsize=14, va="center", color=INK2)
    kR = k * Ds / 2 * 1e-9; good = (mf > 0.85) & (f > 7)
    ax.plot(Ds[good], kR[good], "o", ms=6.5, color=COL[tag], mec="none")
    ax.plot(Ds[~good], kR[~good], "o", ms=6.5, mfc="none", mew=1, color=COL[tag], alpha=0.6)
    ax.set_xlim(Ds[0] - 5, Ds[-1] + 5); ax.set_ylim(0, 33); ax.set_xticklabels([])
    if c == 0:
        ax.set_ylabel("cavity phase $kR$ (rad)")
        ax.legend([Line2D([], [], color=INK2, marker="o", ms=6.5, mec="none", ls="none"),
                   Line2D([], [], color=INK2, marker="o", ms=6.5, mfc="none", mew=1, ls="none")],
                  ["clean $m=+1$ pattern", "not a clean $m=+1$ pattern"], loc="upper left", fontsize=15, labelspacing=0.3,
                  handlelength=1.2, borderaxespad=0.3, frameon=True, facecolor="white", edgecolor="none", framealpha=1, borderpad=0.2)
    else:
        ax.set_yticklabels([])
        ax.text(0.98, 0.04, "solid: rim node $j_{1,n}$\ndotted: rim antinode $j'_{1,n}$", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=14, color=INK2, bbox=dict(fc="white", ec="none", pad=0.6))
    LET.append((ax, ["(c)", "(d)"][c], c))

    # (e,f) wave power
    ax = fig.add_subplot(gs[2, c])
    ax.plot(DP, a2, "o-", ms=6, lw=1.5, color=COL[tag])
    ax.axvline(DMIN_FIT, color=GRID, lw=2)
    ax.set_xlim(Ds[0] - 5, Ds[-1] + 5); ax.set_ylim(0, 0.5); ax.set_yticks([0, 0.2, 0.4]); ax.set_xlabel("$D$ (nm)")
    ax.set_xticks(range(200, 701, 100))                 # no 100: it would sit on the 0.0 label in the corner
    if c == 0: ax.set_ylabel(r"wave power $\langle a^2\rangle$")
    else: ax.set_yticklabels([])
    LET.append((ax, ["(e)", "(f)"][c], c))

    s = DP >= DMIN_FIT
    beta = curve_fit(lambda x, b: fm / (1 + b * x), a2[s], fP[s], p0=[2.5])[0][0]
    rms = np.sqrt(np.mean((fP[s] / fm - 1 / (1 + beta * a2[s])) ** 2))     # rms of f/f_mac, unitless as in Fig. 4(d)
    fit[tag] = (beta, fm, a2[s], fP[s], rms)
    print(f"{tag}: f_mac = {fm:.1f} GHz, beta = {beta:.2f}, rms of f/f_mac = {rms:.3f} (N = {s.sum()}), "
          f"corr(f, a2) = {np.corrcoef(a2[s], fP[s])[0, 1]:.3f}")

# panel letters at the left edge of each column (left of the y labels in that column)
fig.canvas.draw(); rend = fig.canvas.get_renderer()
for col in range(2):
    x = min(a.get_tightbbox(rend).x0 for a, _, c in LET if c == col) / fig.bbox.width
    for a, lab, c in LET:
        if c == col:
            a.text(x, 1.04 if lab not in ("(a)", "(b)") else 1.11, lab, transform=blended_transform_factory(fig.transFigure, a.transAxes),
                   fontsize=22.5, fontweight="bold", va="bottom", ha="left")
save(fig, "figS1")
