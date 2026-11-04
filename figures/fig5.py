"""Fig. 5: high-current regime, J = 12 x 1e11 A/m^2.
Left 2x3 block: m_z snapshot, FFT amplitude, FFT phase for BQ (36 GHz) and no BQ (incoherent).
(d) measured tilt <m_z> versus J with the macrospin tilt a_J/[alpha (Ms+HB)] for both cases
(e) gyration coherence |<m_+>| versus J: the uncoupled disk loses the coherent gyration near J = 10,
    the BQ disk keeps it to J ~ 12-13
(f) radial amplitude profile of the BQ run at 36 GHz with |J1(k r)|, k from the phase slope, and the dispersion check."""
from common import *
from matplotlib.transforms import blended_transform_factory

J_HI = 12.0
fig = plt.figure(figsize=(PAGEW, 11.75))
# maps (2 x 3) on the left, colourbars right of the second map column, line plots in their own column on the right
gs = fig.add_gridspec(3, 2, wspace=0.06, hspace=0.2, left=0.07, right=0.46, top=0.95, bottom=0.09)
gsR = fig.add_gridspec(3, 1, hspace=0.3, left=0.6, right=0.93, top=0.95, bottom=0.09)

res, ims, LET = {}, {}, []                        # LET: (axes, letter, column); letters placed at the end
for c, tag in enumerate(("bq", "no_bq")):
    m = load_maps(tag, J=J_HI); Rn = m["R"] * 1e9
    res[tag] = (m, azimuthal(m["H"], m["R"]))
    for r_, (arr, cm, vmin, vmax) in enumerate([(m["mz"], CM_MZ, -0.2, 1.0),
                                                (np.abs(m["H"]), CM_MZ, 0, np.nanpercentile(np.abs(m["H"]), 99.5)),
                                                (np.angle(m["H"]), CM_MZ, -np.pi, np.pi)]):
        ax = fig.add_subplot(gs[r_, c])
        ims[(r_, c)] = (ax, show_map(ax, arr, Rn, cm, vmin=vmin, vmax=vmax))
        if r_ == 0: ax.set_title(f"{'BQ' if tag=='bq' else 'no BQ'}, $J$ = {J_HI:.0f}: {m['f']:.1f} GHz", fontsize=19, color=COL[tag])
        if c == 0: ax.set_ylabel("y (nm)")
        else: ax.set_yticklabels([])
        if r_ < 2: ax.set_xticklabels([])
        else: ax.set_xlabel("x (nm)")
        if c == 0: LET.append((ax, ["(a)", "(b)", "(c)"][r_], 0))
for r_, (lab, ticks, tl) in enumerate([(r"$m_z$", [0, 0.5, 1], None), ("FFT amp. (arb.)", None, None), ("phase", [-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])]):
    ax, im = ims[(r_, 1)]; ax.apply_aspect(); pos = ax.get_position()     # box after the equal-aspect shrink
    cax = fig.add_axes([pos.x1 + 0.008, pos.y0 + 0.08 * pos.height, 0.01, 0.84 * pos.height])
    cb = fig.colorbar(im, cax=cax); cb.set_label(lab, fontsize=16, labelpad=5); cb.ax.tick_params(labelsize=15, length=5)
    if ticks is not None: cb.set_ticks(ticks)
    if tl is not None: cb.set_ticklabels(tl)

def right_ax(i):
    return fig.add_subplot(gsR[i])

pbJ = power_balance("J")
JJ = np.linspace(3, 20, 200)
# (d) tilt vs J
ax = right_ax(0)
for tag in ("bq", "no_bq"):
    P = [r for r in pbJ if r["tag"] == tag]
    J = np.array([float(r["J"]) for r in P]); mz = np.array([float(r["mz"]) for r in P])
    ax.plot(J, mz, "o", ms=6, color=COL[tag], mec="none")
    ax.plot(JJ, np.minimum(mz_mac(JJ * 1e11, tag == "bq"), 1), "--", lw=2, color=COL[tag])
ax.axhline(1, color=GRID, lw=1.5)
ax.set_xlim(3, 20); ax.set_ylim(0, 1.05); ax.set_ylabel(r"tilt $\langle m_z\rangle$"); ax.set_xticklabels([])
from matplotlib.lines import Line2D
ax.legend([Line2D([], [], color=INK2, marker="o", ms=6, mec="none", ls="none"), Line2D([], [], color=INK2, lw=2, ls="--")],
          ["simulation", r"macrospin $a_J/\alpha(M_s+H_B)$"], frameon=False, loc="lower right", fontsize=14.5,
          handlelength=1.2, handletextpad=0.4, borderaxespad=0.1)
LET.append((ax, "(d)", 1))

# (e) coherence vs J
ax = right_ax(1)
for tag in ("bq", "no_bq"):
    P = [r for r in pbJ if r["tag"] == tag]
    J = np.array([float(r["J"]) for r in P]); coh = np.array([float(r["coh"]) for r in P])
    ax.plot(J, coh, "o-", ms=6, lw=1.5, color=COL[tag], label=LAB[tag])
ax.axvline(J_HI, color=GRID, lw=2); ax.text(J_HI + 0.3, 0.5, rf"$J = {J_HI:.0f}\times10^{{11}}$ A/m$^2$", fontsize=14.5, color=INK2)
ax.set_xlim(3, 20); ax.set_ylim(0, 1); ax.set_xlabel(r"$J$ (10$^{11}$ A/m$^2$)", labelpad=2.5); ax.set_ylabel(r"coherence $|\langle m_+\rangle|$")
ax.legend(frameon=False, loc="upper right", fontsize=16)
LET.append((ax, "(e)", 1))

# (f) BQ 36 GHz: amplitude profile and phase of the m=+1 component (inward travelling wave, linear phase)
ax = right_ax(2); ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
m, az = res["bq"]; rb = az["rb"] * 1e9; k = abs(az["k"]); s = (rb > 6) & (rb < 158)
ax.plot(rb, az["amp"] / az["amp"][s].max(), color=COL["bq"], lw=3, label="amplitude")
F1 = az["Fm"][az["mdom"]]; ph = np.unwrap(np.angle(F1[s])); ph = -np.sign(az["k"]) * (ph - ph[0])
ax2.plot(rb[s], ph / (2 * np.pi), color=INK2, lw=2.5, ls="--", label="phase of $m=+1$ part")
ax.set_xlim(0, 160); ax.set_ylim(0, 1.5); ax.set_xlabel("r (nm)"); ax.set_ylabel("amplitude (norm.)")
ax2.set_ylim(-6, 1.8); ax2.set_ylabel(r"phase / $2\pi$ (dashed)", fontsize=17.5); ax2.tick_params(labelsize=16)
ax.text(0.02, 0.97, rf"BQ, {m['f']:.1f} GHz: $\lambda$ = {2*np.pi/k*1e9:.0f} nm" + "\n" + rf"$\Omega(k)/2\pi$ = {f_of_k(k, True):.1f} GHz",
        transform=ax.transAxes, fontsize=14.5, color=INK2, va="top")
LET.append((ax, "(f)", 1))
kstar = k_of_f(m["f"], True)
# panel letters at the left edge of each column (left of the y labels in that column)
fig.canvas.draw(); rend = fig.canvas.get_renderer()
for col in range(2):
    x = min(a.get_tightbbox(rend).x0 for a, _, c in LET if c == col) / fig.bbox.width
    for a, lab, c in LET:
        if c == col:
            a.text(x, 1.04, lab, transform=blended_transform_factory(fig.transFigure, a.transAxes),
                   fontsize=22.5, fontweight="bold", va="bottom", ha="left")
save(fig, "fig5")
