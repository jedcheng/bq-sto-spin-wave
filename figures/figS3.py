"""Fig. S3: ingredients of the dispersion relation.
(a) the two stiffness brackets N11 (twist) and N22 (tilt) versus k, with and without BQ, at the tilt of the 13 GHz state
(b) Omega(k)/2pi with the linear limit c k and the exchange-only curve; ellipticity |a/b| = sqrt(N22/N11) on a twin axis
(c) ladder f_n = Omega(x_n/R)/2pi versus n for R = 160 nm, both limits, with the simulated plateau tops."""
from common import *
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from sw_theory import Ms, HB, D as Dex, d as dfree

def brackets(k, w, bq):
    hb = HB if bq else 0.0
    mz = np.clip(w / (gam * (Ms + hb)), 0, 0.99); c2, s2 = 1 - mz ** 2, mz ** 2
    p = Pdip(k * dfree)
    return Dex * k ** 2 + 0.5 * Ms * (1 - p), Dex * k ** 2 + (Ms * p + hb) * c2 + 0.5 * Ms * (1 - p) * s2

fig, axs = plt.subplots(1, 3, figsize=(PAGEW, 5.75), gridspec_kw=dict(wspace=0.62))
kk = np.linspace(1e6, 2.2e8, 300)
ax = axs[0]
for tag, f0 in (("bq", 13.09), ("no_bq", 11.93)):
    N11, N22 = brackets(kk, 2 * np.pi * f0 * 1e9, tag == "bq")
    ax.plot(kk / 1e6, N11 / 1e5, color=COL[tag], lw=3)
    ax.plot(kk / 1e6, N22 / 1e5, color=COL[tag], lw=3, ls="--")
ax.plot(kk / 1e6, Dex * kk ** 2 / 1e5, color=INK2, lw=2, ls=":")
ax.set_yscale("log"); ax.set_ylim(0.01, 30); ax.set_xlim(0, 220)
ax.set_xlabel(r"$k$ (rad/$\mu$m)"); ax.set_ylabel(r"stiffness (10$^5$ A/m)")
ax.legend([Line2D([], [], color=INK2, lw=3), Line2D([], [], color=INK2, lw=3, ls="--"), Line2D([], [], color=INK2, lw=2, ls=":")],
          ["$N_{11}$ (twist)", "$N_{22}$ (tilt)", "exchange $Dk^2$"], frameon=False, fontsize=15, loc="lower right",
          handlelength=1.6, borderaxespad=0.2)
panel_label(ax, "(a)", dx=-0.25)

ax = axs[1]; ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
for tag, f0 in (("bq", 13.09), ("no_bq", 11.93)):
    bq = tag == "bq"; w0 = 2 * np.pi * f0 * 1e9
    ax.plot(kk / 1e6, Omega(kk, w0, bq) / 2 / np.pi / 1e9, color=COL[tag], lw=3)
    ax.plot(kk / 1e6, c_sound(bq) * kk / 2 / np.pi / 1e9, color=COL[tag], lw=2, ls=":")
    N11, N22 = brackets(kk, w0, bq); ell = np.sqrt(N22 / N11)
    ax2.plot(kk / 1e6, np.where(ell <= 12, ell, np.nan), color=COL[tag], lw=2, ls="-.")    # clipped below the legend
ax.plot(kk / 1e6, gam * Dex * kk ** 2 / 2 / np.pi / 1e9, color=INK2, lw=2, ls="--")
ax.set_xlim(0, 220); ax.set_ylim(0, 70); ax.set_yticks(range(0, 50, 10)); ax.set_xlabel(r"$k$ (rad/$\mu$m)"); ax.set_ylabel(r"$\Omega/2\pi$ (GHz)")
ax2.set_ylim(0, 70 / 45 * 12); ax2.set_yticks(range(0, 13, 2)); ax2.set_ylabel(r"ellipticity $|a/b|$", fontsize=17.5); ax2.tick_params(labelsize=16)
ax.legend([Line2D([], [], color=INK2, lw=3), Line2D([], [], color=INK2, lw=2, ls=":"), Line2D([], [], color=INK2, lw=2, ls="--"),
           Line2D([], [], color=INK2, lw=2, ls="-.")],
          [r"$\Omega(k)/2\pi$", r"$ck/2\pi$", "exchange only", r"$|a/b|$ (right axis)"], frameon=False, fontsize=15,
          loc="upper right", handlelength=1.6, borderaxespad=0.2, labelspacing=0.3)
panel_label(ax, "(b)", dx=-0.25)

ax = axs[2]; R = 160e-9; n = np.arange(1, 10)
for tag in ("bq", "no_bq"):
    ffree, fpin = ladder(R, tag == "bq", nmax=9)
    ax.fill_between(n, ffree, fpin, color=COL[tag], alpha=0.25, lw=0)
    J, F = load_jsweep(tag); Je, Fe = plateau_tops(J, F)
    ax.plot(np.arange(1, len(Fe) + 1), Fe, "o", ms=7.5, color=COL[tag], mec="none")
ax.set_xlabel("interior antinodes $n$"); ax.set_ylabel("$f_n$ (GHz)"); ax.set_xlim(0.5, 9.5); ax.set_ylim(0, 35)
ax.set_xticks(range(1, 10, 2))
ax.legend([Patch(color=INK2, alpha=0.25, lw=0), Line2D([], [], color=INK2, marker="o", ms=7.5, mec="none", ls="none")],
          [r"$\Omega(x_n/R)/2\pi$," + "\n" + r"$j'_{1,n}\leq x_n\leq j_{1,n}$", "simulated" + "\n" + "plateau tops"], frameon=False, fontsize=15,
          loc="upper left", handlelength=1.4, borderaxespad=0.2)
panel_label(ax, "(c)", dx=-0.25)
fig.legend([Line2D([], [], color=COL["bq"], lw=4), Line2D([], [], color=COL["no_bq"], lw=4)], [LAB["bq"], LAB["no_bq"]],
           loc="lower center", bbox_to_anchor=(0.5, 0.93), ncol=2, frameon=False, fontsize=17.5)
save(fig, "figS3")
