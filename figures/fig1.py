"""Fig. 1 at full page width.
(a) device schematic and (b) torque schematic: the hand-made images fig1a.png / fig1b.png (white margins trimmed).
(c) f vs J for BQ / no BQ with the macrospin line, its reversal threshold and the cavity ladder.
Layout: (a) above (b) in a narrow left column, (c) fills the rest of the width.
The selection-rule test that used to be (d) is now Fig. 3(a).
Output: out/fig1.pdf/png."""
from common import *
from scipy.optimize import brentq
from PIL import Image


def load_trimmed(fn, pad=8):
    "RGB image with the white border cropped off"
    im = np.asarray(Image.open(os.path.join(HERE, "assets", fn)).convert("RGB"))
    ys, xs = np.where((im < 245).any(axis=2))
    return im[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]


R = 160e-9
# compact layout in inches: (a) above (b) in a narrow left column, (c) fills the rest of the width
imgs = [(fn, lab, load_trimmed(fn)) for fn, lab in (("fig1a.png", "(a)"), ("fig1b.png", "(b)"))]
WI, GAP, BOT, TOP = 4.6, 0.35, 0.15, 0.1                      # image column width, gap between images, margins
hs = [WI * im.shape[0] / im.shape[1] for _, _, im in imgs]     # image heights at that width
W, H = PAGEW, BOT + hs[0] + GAP + hs[1] + TOP
fig = plt.figure(figsize=(W, H))
box = lambda x, y, w, h: fig.add_axes([x / W, y / H, w / W, h / H])
axs_img, y = {}, H - TOP
for (fn, lab, im), h in zip(imgs, hs):
    y -= h
    ax = box(0.45, y, WI, h); ax.imshow(im, interpolation="lanczos"); ax.axis("off"); axs_img[lab] = ax
    fig.text(0.05 / W, (y + h) / H, lab, fontsize=22.5, fontweight="bold", va="top", ha="left")
    y -= GAP
XC = 0.45 + WI + 1.45                                          # left edge of (c), after its y tick labels and y label
axb = box(XC, 0.95, W - XC - 0.1, H - 0.95 - 0.45)

# ------------------------------------------------------------------ (c) f vs J
ax = axb
JJ = np.linspace(2, 16, 200)
ax.plot(JJ, f_mac(JJ * 1e11), color=INK2, lw=2, ls="--")
slope = np.degrees(np.arctan((f_mac(7e11) - f_mac(6e11))))     # angle in data coordinates (rotated by transData)
ax.text(6.0, f_mac(6e11) + 1.4, r"macrospin $f_{\rm mac}=\gamma a_J/2\pi\alpha$", fontsize=15, color=INK2, rotation=slope,
        transform_rotates_text=True, rotation_mode="anchor", ha="left", va="bottom")
for tag in ("bq", "no_bq"):
    J, F = load_jsweep(tag)
    s = (J >= 2) & (J <= 14)
    ax.plot(J[s], F[s], "o-", ms=5, lw=1.5, color=COL[tag], label=LAB[tag])
    Jc = brentq(lambda j: mz_mac(j * 1e11, tag == "bq") - 1, 1, 100)       # macrospin reversal (tilt -> 1)
    print(f"{tag}: macrospin reversal current J_c = {Jc:.1f} x1e11 A/m^2")
    if Jc < 15.5:
        ax.axvline(Jc, color=COL[tag], lw=1.5, ls=":")
    ffree, fpin = ladder(R, tag == "bq")
    x0 = 14.3 if tag == "bq" else 14.9
    for a, b in zip(ffree, fpin):
        if b < 40: ax.fill_between([x0, x0 + 0.5], [a, a], [b, b], color=COL[tag], alpha=0.55, lw=0)
ax.text(15.85, 21, "cavity ladder $f_n$", fontsize=14, ha="center", va="center", rotation=90, color=INK2)
ax.set_xlim(2, 16.2); ax.set_ylim(0, 42)
ax.set_xlabel(r"Current Density $J$ (10$^{11}$ A/m$^2$)"); ax.set_ylabel("STO Freq. $f$ (GHz)")
ax.legend(loc="lower right", frameon=False, bbox_to_anchor=(0.83, 0.0), fontsize=16, handlelength=1.6, borderaxespad=0.3)
panel_label(ax, "(c)", dx=-1.2 / (W - XC - 0.1))

save(fig, "fig1")
