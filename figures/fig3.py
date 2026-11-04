"""Fig. 3: the selection rule and the anatomy of the cavity mode at D = 320 nm (J = 6.0 with BQ: 13.1 GHz;
J = 5.5 without: 12.0 GHz). Panels are lettered down the columns.
Left column: (a) selection-rule test: simulated f versus measured k for all runs, against Omega(k) with no fitted
parameter (moved here from Fig. 1); (b) FFT amplitude and (c) FFT phase maps (x horizontal, correct handedness).
Right column: (d) azimuthal power spectrum, (e) radial amplitude profile with |J1(k r)| (k from the measured phase
slope; tests the Bessel shape), (f) radial phase of the m = +1 component with the pi staircase of J1."""
from common import *
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

CASES = [("bq", 6.0), ("no_bq", 5.5)]
fig = plt.figure(figsize=(PAGEW, 11.75))
# selection rule top left, maps (2 x 2) below it with colourbars right of the second map column,
# line plots in their own column on the right
gsS = fig.add_gridspec(1, 1, left=0.07, right=0.49, top=0.95, bottom=0.685)
gs = fig.add_gridspec(2, 2, wspace=0.06, hspace=0.12, left=0.07, right=0.415, top=0.575, bottom=0.09)
gsR = fig.add_gridspec(3, 1, hspace=0.36, left=0.585, right=0.985, top=0.95, bottom=0.09)

# (a) selection rule: simulated f versus measured k for every run, against Omega(k) with no fitted parameter
ax = fig.add_subplot(gsS[0])
kk = np.linspace(5e6, 2.2e8, 150)
stats = {}
for tag in ("bq", "no_bq"):
    bq = tag == "bq"
    rows = [r for r in wave_table(tag=tag) if float(r["m_frac"]) > 0.9 and float(r["f_GHz"]) > 7]
    k = np.abs(np.array([float(r["k_rad_per_m"]) for r in rows])); f = np.array([float(r["f_GHz"]) for r in rows])
    sw = np.array([r["sweep"] for r in rows])
    ax.plot(k[sw == "J"] / 1e6, f[sw == "J"], "o", ms=6, color=COL[tag], mec="none", label=LAB[tag])
    ax.plot(k[sw == "D"] / 1e6, f[sw == "D"], "s", ms=6.5, mfc="none", mew=1.5, color=COL[tag])
    ax.plot(kk / 1e6, [f_of_k(x, bq) for x in kk], "-", lw=3, color=COL[tag])
    ax.plot(kk / 1e6, c_sound(bq) * kk / 2 / np.pi / 1e9, ":", lw=2, color=COL[tag])
    rat = np.array([f_of_k(x, bq) for x in k]) / f
    stats[tag] = (len(rows), round(rat.mean(), 3), round(rat.std(), 3))
ax.plot(kk / 1e6, gam * Dex * kk ** 2 / 2 / np.pi / 1e9, "--", lw=2, color=INK2)
ax.text(215, 1.0, "filled: $J$ sweep\nopen: $D$ sweep", fontsize=14.5, color=INK2, ha="right", va="bottom", linespacing=1.3)
ax.set_xlim(0, 220); ax.set_ylim(0, 46); ax.set_yticks([0, 10, 20, 30, 40])   # headroom for the legend
ax.set_xlabel(r"measured $k$ (rad/$\mu$m)"); ax.set_ylabel("STO Freq. $f$ (GHz)")
h, l = ax.get_legend_handles_labels()
h += [Line2D([], [], color=INK2, lw=3), Line2D([], [], color=INK2, lw=2, ls=":"), Line2D([], [], color=INK2, lw=2, ls="--")]
l += [r"$\Omega(k)/2\pi$, no fit", r"$ck/2\pi$", "exchange only"]
ax.legend(h, l, loc="upper left", ncol=2, columnspacing=1.0, frameon=False, handletextpad=0.4, borderaxespad=0.2, fontsize=16, labelspacing=0.3, handlelength=1.6)
panel_label(ax, "(a)", dx=-0.16)
XA = 0.07 - 0.16 * (0.49 - 0.07)                     # figure x of the (a) letter; (b) and (c) share it
print("selection-rule statistics (N, mean ratio, std):", stats)


res, ims = {}, {}
for c, (tag, J) in enumerate(CASES):
    m = load_maps(tag, J=J); Rn = m["R"] * 1e9
    res[tag] = (m, azimuthal(m["H"], m["R"]))
    for r_, (arr, cm, vmin, vmax) in enumerate([(np.abs(m["H"]), CM_MZ, 0, np.nanpercentile(np.abs(m["H"]), 99.5)),
                                                (np.angle(m["H"]), CM_MZ, -np.pi, np.pi)]):
        ax = fig.add_subplot(gs[r_, c])
        ims[(r_, c)] = (ax, show_map(ax, arr, Rn, cm, vmin=vmin, vmax=vmax))
        if r_ == 0: ax.set_title(f"{'BQ' if tag=='bq' else 'no BQ'}, $J$ = {J:.1f}: {m['f']:.1f} GHz", fontsize=17, color=COL[tag])
        if c == 0: ax.set_ylabel("y (nm)")
        else: ax.set_yticklabels([])
        if r_ < 1: ax.set_xticklabels([])
        else: ax.set_xlabel("x (nm)")
        if c == 0: ax.text(XA, 1.04, ["(b)", "(c)"][r_], transform=blended_transform_factory(fig.transFigure, ax.transAxes),
                           fontsize=22.5, fontweight="bold", va="bottom", ha="left")

# colourbars: thin axes to the right of the second map column
for r_, (lab, ticks, tl) in enumerate([("FFT amp. (arb.)", None, None), ("phase", [-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])]):
    ax, im = ims[(r_, 1)]; ax.apply_aspect(); pos = ax.get_position()     # box after the equal-aspect shrink
    cax = fig.add_axes([pos.x1 + 0.008, pos.y0 + 0.08 * pos.height, 0.01, 0.84 * pos.height])
    cb = fig.colorbar(im, cax=cax); cb.set_label(lab, fontsize=16, labelpad=5); cb.ax.tick_params(labelsize=15, length=5)
    if ticks is not None: cb.set_ticks(ticks)
    if tl is not None: cb.set_ticklabels(tl)

# ------------------------------------------------------------------ right column
def right_ax(i):
    return fig.add_subplot(gsR[i])

# (d) azimuthal power spectrum
ax = right_ax(0); ms = [-3, -2, -1, 0, 1, 2, 3]
for i, (tag, J) in enumerate(CASES):
    pw = res[tag][1]["power"]; tot = sum(pw.values())
    ax.bar(np.array(ms) + (-0.18 if i == 0 else 0.18), [pw[m_] / tot for m_ in ms], width=0.34, color=COL[tag], label=LAB[tag])
ax.set_xticks(ms); ax.set_xlabel(r"azimuthal index $m$", labelpad=2.5); ax.set_ylabel("power fraction")
ax.set_ylim(0, 1.05); ax.legend(frameon=False, loc="upper left", fontsize=16)
ax.text(0.98, 0.95, "co-rotating\n$m=+1$", transform=ax.transAxes, ha="right", va="top", fontsize=16, color=INK2)
panel_label(ax, "(d)", dx=-0.21)

# (e) radial amplitude profile with |J1(k r)|
ax = right_ax(1); lam = {}
for tag, J in CASES:
    m, az = res[tag]; rb = az["rb"] * 1e9; amp = az["amp"]; k = abs(az["k"])
    s = rb < 0.85 * m["R"] * 1e9
    prof = np.abs(j1(k * az["rb"])); prof *= amp[s].max() / prof[s].max()
    ax.plot(rb, amp / amp[s].max(), color=COL[tag], lw=3)
    ax.plot(rb, prof / amp[s].max(), color=COL[tag], lw=2, ls="--")
    for n in range(1, 4): ax.plot([jn_zeros(1, n)[-1] / k * 1e9] * 2, [0, 1.1], color=COL[tag], lw=1, alpha=0.6)
    lam[tag] = 2 * np.pi / k * 1e9
ax.set_xlim(0, 160); ax.set_ylim(0, 1.8); ax.set_yticks([0, 0.5, 1]); ax.set_ylabel("amplitude (norm.)"); ax.set_xticklabels([])
ax.legend([Line2D([], [], color=INK2, lw=3), Line2D([], [], color=INK2, lw=2, ls="--"), Line2D([], [], color=INK2, lw=1)],
          ["simulation", r"$|J_1(kr)|$", r"node $j_{1,n}/k$"], frameon=False, loc="upper left", ncol=3, fontsize=15,
          handlelength=1.5, columnspacing=0.9, handletextpad=0.4, borderaxespad=0.2)
for y, tag in ((0.79, "bq"), (0.69, "no_bq")):
    ax.text(0.99, y, rf"$\lambda$ = {lam[tag]:.0f} nm", transform=ax.transAxes, ha="right", va="center", fontsize=15, color=COL[tag])
panel_label(ax, "(e)", dx=-0.21)

# (f) radial phase of the m=+1 component
ax = right_ax(2)
for tag, J in CASES:
    m, az = res[tag]; rb = az["rb"] * 1e9; k = az["k"]
    F1 = az["Fm"][az["mdom"]]; s = (rb > 6) & (rb < 158)
    ph = np.unwrap(np.angle(F1[s])); ph -= ph[0]
    ph *= -np.sign(k); k = abs(k)                       # plot with the phase decreasing outward, same sign for both runs
    ax.plot(rb[s], ph / np.pi, color=COL[tag], lw=3, label=f"{LAB[tag]}")
    stair = -np.pi * np.searchsorted(jn_zeros(1, 6) / k * 1e9, rb[s])
    ax.plot(rb[s], stair / np.pi, color=COL[tag], lw=2, ls="--")
ax.set_xlim(0, 160); ax.set_ylim(-3.6, 0.9); ax.set_yticks([-3, -2, -1, 0]); ax.set_xlabel("r (nm)"); ax.set_ylabel(r"$m=+1$ phase / $\pi$")
ax.legend(frameon=False, loc="lower left", fontsize=16)
ax.text(0.98, 0.95, r"dashed: $\pi$ step at each node of $J_1$", transform=ax.transAxes, ha="right", va="top", fontsize=15, color=INK2)
panel_label(ax, "(f)", dx=-0.21)

for tag, J in CASES:
    m, az = res[tag]
    print(f"{tag}: f={m['f']:.2f} GHz, dominant m={az['mdom']} ({az['frac']:.2f}), k_meas={abs(az['k'])/1e6:.1f} rad/um "
          f"(lambda {2*np.pi/abs(az['k'])*1e9:.0f} nm), k* from Omega(k*)=omega: {k_of_f(m['f'], tag=='bq')/1e6:.1f} rad/um")
save(fig, "fig3")
