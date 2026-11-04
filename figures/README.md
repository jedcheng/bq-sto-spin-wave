# Figure scripts

Each script reproduces one figure of the paper (Figs. 1-5) or of the Supplementary Material (Fig. S1 = SM4 diameter
sweep; figS3 = dispersion ingredients, optional, not included in the manuscript). They read the per-run `.npy`/`.h5`
files under `$DATA_DIR` (default `data/sim`, filled by `data/download.sh`) and the tables in `data/derived`, and write
PDF (vector) and PNG (600 dpi at print size) to `figures/out/` (or `$BQSTO_FIGOUT`).
Run everything with `bash make_all.sh` (set `PY=` to choose the interpreter).
Overlap check: `python check_overlap.py figX.py` runs a script and lists overlapping text before it saves.

| script | output | content |
|---|---|---|
| `fig1.py` | `fig1` | (a) device schematic above (b) torque schematic in a narrow left column, from `assets/fig1a.png` / `assets/fig1b.png` (white margins trimmed); (c) f vs J with macrospin line, reversal thresholds and cavity ladder, filling the rest of the width (selection-rule test moved to fig3(a)) |
| `fig2.py` | `fig2` | birth of the wave: M_z snapshots at 0.1 and 1.0 ns, then the radius-time map of the m=+1 component with the inward front (1116 / 1043 m/s vs c = 1130 / 1006 m/s) |
| `fig3.py` | `fig3` | anatomy of the mode at 13.1 / 12.0 GHz, lettered down the columns: (a) selection-rule test, 177 runs, (b) FFT amplitude and (c) FFT phase maps (correct handedness); (d) azimuthal spectrum, (e) radial amplitude vs |J1(kr)|, (f) radial phase staircase |
| `fig4.py` | `fig4` | profiles before/after the jump; cavity phase kR vs J with j_{1,n}; wave power and f/f_mac vs J; collapse f = f_mac/(1 + beta<a^2>) for J and D (diameter) sweeps |
| `fig5.py` | `fig5` | J = 12: maps for BQ (36 GHz) and no BQ (incoherent); tilt vs J with macrospin tilt; gyration coherence vs J; inward-wave profile and phase at 36 GHz |
| `figS3.py` | `figS3` | dispersion ingredients: stiffness brackets, Omega(k) with limits and ellipticity, ladder vs n |
| `figS1.py` | `figS1` | size dependence (D = 100-700 nm), in the main-figure style: (a,b) FFT spectra with the cavity fan Omega(j_{1,n}/R), (c,d) cavity phase kR with j_{1,n}, (e,f) wave power vs D, (g) f vs <a^2> with f_mac/(1 + beta<a^2>) for D >= 180 nm |

Conventions: magnum.np arrays are `[x, y]`; `show_map` transposes so x is horizontal.
Blue = with BQ, orange = without BQ. Numbers printed by each script (ratios, fitted c, k values) are the ones to quote in the text.

Numbers printed on the last run:
- launch (fig2): rim responds at 0.05 ns, front reaches the centre at 0.18 ns; front speed 1116 (BQ) / 1043 (no BQ) m/s vs c = 1130 / 1006 m/s
- selection rule (fig3(a)): f_theory/f_sim = 1.073 +/- 0.091 (BQ, 126 runs), 1.068 +/- 0.063 (no BQ, 51 runs)
- 13.1 GHz state: m = +1 carries 97 % of the map power, lambda = 94 nm (92 nm without BQ)
- power law: beta = 2.4 (BQ), 2.6 (no BQ), rms of f/f_mac = 0.04, J and D sweeps together (fig4(d)); D sweep alone (figS1(g)): beta = 2.2 / 2.6, f_mac = 20.5 / 18.8 GHz, rms of f/f_mac = 0.02 / 0.02
- macrospin reversal current (tilt -> 1): 13.9 (no BQ), 17.5 (BQ) x 1e11 A/m^2
- J = 12 with BQ: 36.2 GHz, m = +1 fraction 0.81, lambda = 31 nm, Omega(k)/2pi = 37.9 GHz at the measured k
