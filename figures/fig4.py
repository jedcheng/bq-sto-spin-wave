"""Fig. 4: the jump as a change of cavity mode, and the power balance (current sweep, D = 320 nm).
(a) radial amplitude profiles just before / after the jump (BQ 6.0 -> 6.1, no BQ 5.5 -> 5.6): n = 3 -> 4
(b) cavity phase k R versus J with the Bessel zeros j_{1,n}: pinned phase, creep within a plateau, drop at the jump
(c) in-plane wave power <a^2> versus J (rotating-frame mean-square amplitude at 0.5-1.5 omega), with f/f_mac
(d) collapse f/f_mac = 1/(1 + beta <a^2>) for the J sweep and the D (diameter) sweep together; beta is the only fitted number."""
from common import *
from scipy.optimize import curve_fit
from matplotlib.transforms import blended_transform_factory

JUMPS = {"bq": (6.0, 6.1), "no_bq": (5.5, 5.6)}
fig, axs = plt.subplots(2, 2, figsize=(PAGEW, 10.75), gridspec_kw=dict(wspace=0.36, hspace=0.42))

# ------------------------------------------------------------------ (a) profiles across the jump
ax = axs[0, 0]
for tag, (J0, J1_) in JUMPS.items():
    for J, ls, lab in ((J0, "-", "before"), (J1_, "--", "after")):
        m = load_maps(tag, J=J); az = azimuthal(m["H"], m["R"])
        rb = az["rb"] * 1e9; s = rb < 0.85 * m["R"] * 1e9
        ax.plot(rb, az["amp"] / az["amp"][s].max() + (0 if tag == "bq" else 1.3), ls=ls, lw=3, color=COL[tag],
                label=f"{'BQ' if tag == 'bq' else 'no BQ'}, $J$={J:.1f}: {m['f']:.1f} GHz")
ax.set_xlim(0, 160); ax.set_ylim(0, 3.9); ax.set_yticks([]); ax.set_xlabel("r (nm)"); ax.set_ylabel("amplitude at $f$ (offset)")
ax.legend(frameon=False, fontsize=14, loc="upper left", ncol=2, handlelength=1.3, columnspacing=1.0, handletextpad=0.4,
          borderaxespad=0.1)
ax.text(0.5, 0.79, "solid: before the jump, 3 interior antinodes\ndashed: after, 4 (the middle one splits)", transform=ax.transAxes,
        fontsize=14, color=INK2, ha="center", va="top")

# ------------------------------------------------------------------ (b) cavity phase vs J
ax = axs[0, 1]
R = 160e-9
for n in range(1, 9):
    ax.axhline(jn_zeros(1, n)[-1], color=GRID, lw=2); ax.text(10.05, jn_zeros(1, n)[-1], f"$j_{{1,{n}}}$", fontsize=14, va="center", color=INK2)
    ax.axhline(jnp_zeros(1, n)[-1], color=GRID, lw=1, ls=":")
for tag in ("bq", "no_bq"):
    rows = wave_table("J", tag)
    J = np.array([float(r["J"]) for r in rows]); k = np.abs(np.array([float(r["k_rad_per_m"]) for r in rows]))
    mf = np.array([float(r["m_frac"]) for r in rows]); f = np.array([float(r["f_GHz"]) for r in rows])
    good = (mf > 0.85) & (f > 7)
    ax.plot(J[good], k[good] * R, "o", ms=6.5, color=COL[tag], mec="none", label=LAB[tag])
    ax.plot(J[~good], k[~good] * R, "o", ms=6.5, mfc="none", mew=1, color=COL[tag], alpha=0.6)
ax.set_xlim(3, 10); ax.set_ylim(0, 27); ax.set_xlabel(r"$J$ (10$^{11}$ A/m$^2$)"); ax.set_ylabel("cavity phase $kR$ (rad)")
from matplotlib.lines import Line2D
h, l = ax.get_legend_handles_labels()
h += [Line2D([], [], color=INK2, marker="o", ms=6.5, mfc="none", mew=1, ls="none")]
l += ["not a clean $m=1$ pattern"]
ax.legend(h, l, loc="upper left", fontsize=15, labelspacing=0.3, handlelength=1.2, borderaxespad=0.3,
          frameon=True, facecolor="white", edgecolor="none", framealpha=1, borderpad=0.2)
ax.text(0.98, 0.03, "solid: rim node $j_{1,n}$\ndotted: rim antinode $j'_{1,n}$", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=14, color=INK2, bbox=dict(fc="white", ec="none", pad=0.6))

# ------------------------------------------------------------------ (c) wave power and f/f_mac vs J
ax = axs[1, 0]; ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
pbJ = power_balance("J")
fitdata = {}
for tag in ("bq", "no_bq"):
    P = [r for r in pbJ if r["tag"] == tag]
    J = np.array([float(r["J"]) for r in P]); a2 = np.array([float(r["aip1"]) for r in P])
    f = np.array([float(r["f"]) for r in P]); fm = np.array([float(r["fmac"]) for r in P]); coh = np.array([float(r["coh"]) for r in P])
    s = (J >= 3.5) & (J <= 10)
    ax.plot(J[s], a2[s], "o-", ms=6, lw=2, color=COL[tag], label=LAB[tag])
    ax2.plot(J[s], (f / fm)[s], "s--", ms=5, lw=1.5, color=COL[tag], alpha=0.6)
    fitdata[tag] = dict(J=J, a2=a2, f=f, fm=fm, coh=coh)
ax.set_xlim(3, 10); ax.set_ylim(0, 0.56); ax.set_xlabel(r"$J$ (10$^{11}$ A/m$^2$)"); ax.set_ylabel(r"wave power $\langle a^2\rangle$ (circles)")
ax2.set_ylim(0, 1.12); ax2.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])   # 0 omitted: it would touch the x tick "10"; ax2.set_ylabel(r"$f/f_{\rm mac}$ (squares)", color=INK2, labelpad=2.5); ax2.tick_params(labelsize=16)
ax.legend(frameon=False, loc="upper center", ncol=2)

# ------------------------------------------------------------------ (d) collapse
ax = axs[1, 1]
pbD = power_balance("D")
xx = np.linspace(0, 0.45, 100)
for tag in ("bq", "no_bq"):
    dJ = fitdata[tag]; sJ = (dJ["J"] >= 3.5) & (dJ["J"] <= 10) & (dJ["coh"] > 0.3)
    PD = [r for r in pbD if r["tag"] == tag]
    Dn = np.array([float(r["D"]) for r in PD]); a2D = np.array([float(r["aip1"]) for r in PD])
    fD = np.array([float(r["f"]) for r in PD]); fmD = np.array([float(r["fmac"]) for r in PD]); sD = Dn >= 180
    x = np.concatenate([dJ["a2"][sJ], a2D[sD]]); y = np.concatenate([(dJ["f"] / dJ["fm"])[sJ], (fD / fmD)[sD]])
    cfit = curve_fit(lambda a, c: 1 / (1 + c * a), x, y, p0=[2.5])[0][0]
    rms = np.sqrt(np.mean((y - 1 / (1 + cfit * x)) ** 2))
    ax.plot(dJ["a2"][sJ], (dJ["f"] / dJ["fm"])[sJ], "o", ms=6.5, color=COL[tag], mec="none")
    ax.plot(a2D[sD], (fD / fmD)[sD], "s", ms=7, mfc="none", mew=1.5, color=COL[tag])
    ax.plot(xx, 1 / (1 + cfit * xx), "-", lw=3, color=COL[tag])
    print(f"{tag}: beta = {cfit:.2f}, rms of f/f_mac = {rms:.3f}  (N = {len(x)})")
    ax.text(0.98, 0.70 - (0 if tag == "bq" else 0.08), rf"{'BQ' if tag == 'bq' else 'no BQ'}: $\beta$ = {cfit:.1f}, rms {rms:.2f}",
            transform=ax.transAxes, ha="right", fontsize=15, color=COL[tag])
ax.set_xlim(0, 0.45); ax.set_ylim(0.3, 1.0); ax.set_xlabel(r"wave power $\langle a^2\rangle$"); ax.set_ylabel(r"$f/f_{\rm mac}$")
ax.legend([Line2D([], [], color=INK2, marker="o", ms=6.5, mec="none", ls="none"),
           Line2D([], [], color=INK2, marker="s", ms=7, mfc="none", mew=1.5, ls="none")],
          ["$J$ sweep", "$D$ sweep"], frameon=False, loc="lower left", fontsize=16)
ax.text(0.98, 0.98, r"$f=\dfrac{f_{\rm mac}}{1+\beta\langle a^2\rangle}$", transform=ax.transAxes, fontsize=20, ha="right", va="top")
# panel letters at the left edge of each column (left of the widest y label in that column)
fig.canvas.draw(); rend = fig.canvas.get_renderer()
for col in range(2):
    x = min(a.get_tightbbox(rend).x0 for a in axs[:, col]) / fig.bbox.width
    for row in range(2):
        a = axs[row, col]
        a.text(x, 1.04, "(%s)" % "abcd"[2 * row + col], transform=blended_transform_factory(fig.transFigure, a.transAxes),
               fontsize=22.5, fontweight="bold", va="bottom", ha="left")
save(fig, "fig4")
