# Self-generated spin wave in orthogonal biquadratic coupled spin torque oscillators

Code and data pipeline for the paper *"Self-Generated Spin Wave in Orthogonal Biquadratic Coupled Spin Torque
Oscillators"* (DOI / arXiv: `TO BE ADDED`).

**Summary.** A 320 nm CoPt/Cr/CoFe orthogonal spin-torque oscillator (STO) with biquadratic (BQ) interlayer coupling is
simulated micromagnetically (magnum.np). Its frequency lies well below the macrospin prediction and rises in plateaus
with current. The simulations show a spin wave launched at the rim of the free layer by the stray field of the
perpendicular fixed layer. In the spin frame co-rotating with the free layer the wave obeys an analytical dispersion
relation Omega(k), which predicts its wavelength and velocity. The wave forms a standing one-armed Bessel mode
J1(kr)e^{i phi}; the Bessel-zero boundary conditions give bands of allowed cavity frequencies that explain the plateaus
in f(J). BQ coupling stiffens the wave against out-of-plane tilting and keeps the free layer coherent at high current,
which is where the frequency gain comes from.

This repository replaces the "code" supplementary section of the manuscript: it documents how to run the simulations
and how every figure and table was made.

## Repository map

```
simulations/
  micromagnetic/run.py        one magnum.np LLG run (all four sweeps of the paper), parameterised
  micromagnetic/jobs/         PJM job examples: single run, J-sweep and D-sweep bulk jobs
  micromagnetic/original/     verbatim original scripts; CHANGES.md documents the (non-physics) differences
  macrospin/                  two-macrospin LLG without demag (control)
  requirements-sim.txt        simulation environment
analysis/
  paths.py                    data location: one setting, the environment variable DATA_DIR (default data/sim)
  postprocess_fft.py          data.h5 -> freq/fft/fft_avg/heat_map .npy (frequency, spectra, FFT maps)
  j_vs_freq.py                f(J) tables
  wave_extract.py             azimuthal decomposition, dominant m, radial k, lambda  -> wave_table.csv
  power_balance.py, power_balance_J.py   LLG power balance, tilt, coherence, wave power <a^2>  -> power_balance*.csv
  sw_theory.py                analytical model: Omega(k), self-consistent selection rule, cavity ladder
  launch_front.py             inward front speed of the launched wave (Fig. 2 numbers)
  frame_check.py              uniform rotation / co-rotating component check (SM3.4)
  macrospin_table.py          Table S2;   diameter_plateaus.py   Table S4
figures/                      fig1-fig5.py, figS1.py (SM4), figS3.py (extra), make_all.sh, common.py, assets/
data/                         README (layout, figure -> file table), download.sh, MANIFEST (sha256), derived/*.csv
```

## Quick start

Analysis environment (CPU only): Python 3.11 and `pip install -r requirements.txt` (or `conda env create -f environment.yml`).

### A. Download the data and regenerate the figures

The simulation outputs (504 runs, 605 GB) are on a public Nextcloud share:
**https://genkai-nextcloud.hpc.kyushu-u.ac.jp/s/iHs5kSBoLckLR6a** (layout `J_sweep/{bq,no_bq}/<J>/`,
`D_sweep/{bq,no_bq_5.5}/<D>/`, see `data/README.md`). You never need all of it: `data/download.sh` fetches only the
files a given figure needs and checks their sha256 against `data/MANIFEST.tsv`. All scripts read the data from
`$DATA_DIR` (default `data/sim`), which has the same layout as the share.

```bash
bash data/download.sh -n figures   # dry run: what the figures need (270 files, 2.8 GB, nearly all of it the two data.h5 of Fig. 2)
bash data/download.sh figures      # or per figure: fig2 (2.8 GB), fig3, fig4, fig5, figS1 (a few MB each)
bash figures/make_all.sh           # -> figures/out/fig1..fig5, figS1, figS3 (.pdf, .png); ~1 min on a laptop-class CPU
cd analysis
python macrospin_table.py          # Table S2            (no download)
python diameter_plateaus.py        # Table S4            (no download)
python sw_theory.py                # selection-rule statistics, predicted plateau ladder (no download)
python wave_extract.py             # rebuilds output/wave_table.csv   (needs: download.sh npy, 0.83 GB)
python j_vs_freq.py                # rebuilds output/j_vs_freq_*.csv  (needs: download.sh npy)
python launch_front.py             # launch diagnostics, r-t maps     (needs: download.sh ref)
python frame_check.py bq 6.0       # SM3.4                            (needs: download.sh ref)
```
Fig. 1, Fig. S3, Tables S2/S4 and the notebook need no download at all (they read `data/derived/`). The power-balance
tables need every `data.h5` (`download.sh h5-J` 381 GB, `download.sh h5-D` 223 GB). If you have the whole share
mounted or copied somewhere, just `export DATA_DIR=/path/to/share` instead of downloading.

### B. Re-run the simulations

Needs an NVIDIA GPU, PyTorch (CUDA) and the latest magnum.np:

```bash
pip install git+https://gitlab.com/magnum.np/magnum.np@main
```

Details in `simulations/README.md`.

```bash
python simulations/micromagnetic/run.py --D 320 --J 6.0 --bq    --precision float64 --outdir data/sim/J_sweep/bq/6.0
python simulations/micromagnetic/run.py --D 320 --J 5.5 --no-bq                     --outdir data/sim/J_sweep/no_bq/5.5
python analysis/postprocess_fft.py data/sim/J_sweep/bq/6.0 data/sim/J_sweep/no_bq/5.5
```
Full sweeps on a PJM cluster: `pjsub --bulk --sparam 10-200 -x TAG=bq simulations/micromagnetic/jobs/j_sweep.sh` etc.
(see the headers of the job scripts). Cost: one D = 320 nm run takes roughly 15-25 min on a 1/7 H100 (MIG) slice
including the relaxation, about 1 GB of output; the 504 runs of the paper produce 604 GB.

Floating-point note: runs on other GPUs/software versions are not bit-identical; the gyration frequency (FFT bin
0.1 GHz) and the wave pattern are the reproducible quantities (see "Verification").

## Paper item -> script -> input

| paper item | script | input |
|---|---|---|
| Table I | `simulations/micromagnetic/run.py` | - |
| Fig. 1(c) f vs J, macrospin line, cavity ladder | `figures/fig1.py` | `data/derived/j_vs_freq_*.csv` |
| Fig. 2 birth of the wave | `figures/fig2.py` (prints the front speeds; `analysis/launch_front.py` is the diagnostic version) | `data.h5` of `J_sweep/bq/6.0`, `J_sweep/no_bq/5.5` |
| Fig. 3 selection rule, FFT maps, Bessel profile | `figures/fig3.py` | `wave_table.csv`; npy of `J_sweep/bq/6.0`, `J_sweep/no_bq/5.5` |
| Fig. 4 jump, cavity phase, wave power, collapse | `figures/fig4.py` | `wave_table.csv`, `power_balance*.csv`; npy of bq 6.0/6.1, no_bq 5.5/5.6 |
| Fig. 5 J = 12 regime | `figures/fig5.py` | `power_balance_J.csv`; npy of bq 12.0, no_bq 12.0 |
| SM1 effective fields | magnum.np `BiquadraticRKKYField`, `SpinTorqueSlonczewski` | - |
| SM2 macrospin, Table S2 | `analysis/macrospin_table.py` | `j_vs_freq_*.csv`, `power_balance_J.csv` |
| SM3 dispersion / selection rule | `analysis/sw_theory.py`; SM3.4: `analysis/frame_check.py` | `wave_table.csv`, `j_vs_freq_*.csv`; reference `data.h5` for SM3.4 |
| SM4 Fig. S1, Table S4 | `figures/figS1.py`, `analysis/diameter_plateaus.py` | `wave_table.csv`, `power_balance.csv`; `freq.npy`, `fft_avg.npy` of the D sweep (Fig. S1 only) |
| Fig. S3 (extra, not in the manuscript) | `figures/figS3.py` | `j_vs_freq_*.csv` |

Derived tables in `data/derived/` come from: `j_vs_freq.py`, `wave_extract.py` (npy only), `power_balance.py` and
`power_balance_J.py` (need `data.h5`). Full file list per item: `data/README.md`.

## Verification (what was checked for this release)

* **Figures, end to end from the public share.** `download.sh fig3 fig4 fig5 figS1` (268 files, all checksums OK)
  into an empty `DATA_DIR` plus the two reference `data.h5`, then `make_all.sh`: all seven PNGs (Figs. 1-5, S1, S3)
  are pixel-identical to the figures of the manuscript. Tables S2 and S4 and the selection-rule statistics
  (1.073 +/- 0.091, 1.068 +/- 0.063) are reproduced by the scripts above (one rounding difference in Table S4: no BQ,
  n = 2, <a^2> = 0.325 -> 0.215 printed as 0.33 -> 0.22 in the manuscript).
* `wave_extract.py` and `j_vs_freq.py` reproduce `wave_table.csv` and `j_vs_freq_*.csv` exactly from the `.npy` files;
  `power_balance.run` reproduces the rows of `power_balance_J.csv` for the two reference runs.
* **Simulation (`simulations/micromagnetic/run.py`)**, re-run of the two reference cases on one H100 MIG 1g.12gb slice
  (magnum.np 2.0.1 + a local non-physics patch, torch 2.3.1+cu121; the repository now targets magnum.np `main`, which was not re-tested), compared with the published outputs:
  - no BQ, J = 5.5, float32 (15 min): **bit-identical** to the published run (relaxed state, every stored frame of
    `data.h5`, `final.npy`); f = 12.0 GHz, m = +1 (98 %), lambda = 92 nm.
  - BQ, J = 6.0, float64 (25 min): not bit-identical. The relaxed state already differs by 2e-4 (the published run was
    most likely started from a stored relaxed state, `--init-from`; see `simulations/micromagnetic/CHANGES.md`), and the
    difference grows over the 10 ns. The observables agree: f = 13.2 vs 13.1 GHz (one 0.1 GHz FFT bin),
    m = +1 with 97 % of the power in both, k = 67.5 vs 67.0 rad/um (lambda 93 vs 94 nm), 4 amplitude rings in both.
* Not verified: the full sweeps, `power_balance*.py` over all runs, `download.sh` for the large `data.h5` targets
  (the same code path was tested on small files, including resume and checksum-failure handling).

## Citation

See `CITATION.cff`. Until the paper is published:

```
@article{bqsto_selfgenerated_spinwave,
  title   = {Self-Generated Spin Wave in Orthogonal Biquadratic Coupled Spin Torque Oscillators},
  author  = {AUTHORS},
  journal = {Applied Physics Letters},
  year    = {YEAR},
  doi     = {DOI}
}
```
Please also cite magnum.np: F. Bruckner et al., Sci. Rep. 13, 12054 (2023).

## License

Code: MIT (`LICENSE`). Data (Nextcloud share above): CC BY 4.0 suggested (authors' choice, to be stated with the data share).
magnum.np itself is GPL-3.0 and is installed separately. The manuscript source is not included in this repository.
