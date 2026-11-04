# Data

## What is in git

`derived/` holds the small tables made by the analysis scripts from the simulation outputs. They are enough for
Fig. 1, Fig. S3, Tables S2 and S4 and the theory notebook, without downloading anything.

| file | rows | made by |
|---|---|---|
| `derived/j_vs_freq_bq.csv`, `derived/j_vs_freq_no_bq.csv` | gyration frequency f vs J, D = 320 nm (J = 1.0-19.9) | `analysis/j_vs_freq.py` |
| `derived/wave_table.csv` | per run: f, dominant azimuthal index m and its power fraction, radial k from the phase slope, lambda (J sweep 3.0-13.0, D sweep 100-700 nm) | `analysis/wave_extract.py` |
| `derived/power_balance_J.csv` | per run of the J sweep (3.0-10.0 step 0.1, 10.5-20.0 step 0.5): LLG power balance, tilt `mz`, coherence `coh`, in-plane wave power `aip1`, f_mac, ... over the last 3 ns | `analysis/power_balance_J.py` |
| `derived/power_balance.csv` | same for the D sweep (tag `no_bq_5.5` = no BQ at J = 5.5) | `analysis/power_balance.py` |

`MANIFEST.tsv` (path, bytes, sha256) and `MANIFEST.sha256` (same, `sha256sum -c` format) list the 3020 simulation
output files used by the scripts (the five `.npy` files and `data.h5` of each of the 504 runs, 605 GB), with paths
relative to the share / `DATA_DIR`.

## Simulation outputs (public Nextcloud share)

**https://genkai-nextcloud.hpc.kyushu-u.ac.jp/s/iHs5kSBoLckLR6a** (read-only, no account needed). Layout:

```
J_sweep/bq/<J>/           current sweep with BQ, D = 320 nm, J = 1.0 ... 20.0 (1e11 A/m^2, folder name "%.1f")
J_sweep/no_bq/<J>/        current sweep without BQ, J = 1.0 ... 20.0
D_sweep/bq/<D>/           diameter sweep with BQ, J = 6.0, D = 100 ... 700 nm in steps of 10
D_sweep/no_bq_5.5/<D>/    diameter sweep without BQ, J = 5.5
each run folder:
  freq.npy, fft.npy       rfft of the spatially averaged free-layer m_x (gives f)
  fft_avg.npy             cell-averaged |rfft(m_z)| spectrum
  heat_map.npy            complex rfft(m_z) map at f, [x, y], NaN outside the disk
  final.npy, initial.npy  magnetisation at t = 10 ns / after the relaxation, [x, y, z, 3]
  data.h5                 full trajectory: m [2001, x, y, 3, 3] every 5 ps, t [2001]  (0.1-4 GB per run)
  *.png                   diagnostic plots of the original post-processing (not used)
```
Approximate sizes: `J_sweep/bq` 191 GB, `J_sweep/no_bq` 191 GB, `D_sweep/bq` 111 GB, `D_sweep/no_bq_5.5` 112 GB.
`J_sweep/bq/20.1` is empty (run not completed); `J_sweep/*/20.0` has no `heat_map.npy`/`fft_avg.npy`. Neither is used.

The scripts read the data from **`DATA_DIR`** (environment variable, default `data/sim`), which has exactly the share
layout. Either download what you need with `download.sh`, or point `DATA_DIR` at a full copy of the share.

### download.sh: fetch only what a figure needs

```bash
bash data/download.sh -n fig4            # dry run: list files and total size
bash data/download.sh fig3 fig4 fig5     # npy files of the runs used by these figures (~10 MB)
bash data/download.sh figures            # everything any figure needs (2.8 GB, mostly the two data.h5 of Fig. 2)
DATA_DIR=/scratch/bqsto bash data/download.sh npy    # all .npy files of all runs (0.83 GB, ~2500 files)
bash data/download.sh 'J_sweep/bq/6.?/data.h5'       # any glob over the paths in MANIFEST.tsv
```
Every file is checked against `MANIFEST.tsv` (sha256); files already present with the right checksum are skipped, so
an interrupted download can be restarted. Downloads use the share's public WebDAV endpoint
(`https://genkai-nextcloud.hpc.kyushu-u.ac.jp/public.php/webdav/<path>`, user = share token, empty password); a single
file can also be fetched in a browser or with
`curl -L "https://genkai-nextcloud.hpc.kyushu-u.ac.jp/s/iHs5kSBoLckLR6a/download?path=/J_sweep/bq/6.0&files=freq.npy"`.
Throughput is roughly 0.6 s per small file plus transfer time.

### Which files each figure / table needs

| item | script | inputs | `download.sh` target |
|---|---|---|---|
| Table I (parameters) | `simulations/micromagnetic/run.py` | - | - |
| Fig. 1(c) | `figures/fig1.py` | `derived/j_vs_freq_*.csv`, `figures/assets/fig1a.png`, `fig1b.png` (hand-drawn panels a, b) | none |
| Fig. 2 | `figures/fig2.py` | `J_sweep/bq/6.0/data.h5`, `J_sweep/no_bq/5.5/data.h5` (first 230 frames) | `fig2` (2.8 GB) |
| Fig. 3 | `figures/fig3.py` | `derived/wave_table.csv`; `final, heat_map, freq, fft .npy` of `J_sweep/bq/6.0`, `J_sweep/no_bq/5.5` | `fig3` (4 MB) |
| Fig. 4 | `figures/fig4.py` | `derived/wave_table.csv`, `derived/power_balance_J.csv`, `derived/power_balance.csv`; npy of `J_sweep/bq/6.0`, `6.1`, `J_sweep/no_bq/5.5`, `5.6` | `fig4` (6 MB) |
| Fig. 5 | `figures/fig5.py` | `derived/power_balance_J.csv`; npy of `J_sweep/bq/12.0`, `J_sweep/no_bq/12.0` | `fig5` (3 MB) |
| Table S2 | `analysis/macrospin_table.py` | `derived/j_vs_freq_*.csv`, `derived/power_balance_J.csv` | none |
| SM3 selection rule, plateau ladder | `analysis/sw_theory.py` | `derived/wave_table.csv`, `derived/j_vs_freq_*.csv` | none |
| SM3.4 (uniform rotation check) | `analysis/frame_check.py bq 6.0` / `no_bq 5.5` | the two reference `data.h5` | `ref` |
| Fig. S1 (SM4) | `figures/figS1.py` | `derived/wave_table.csv`, `derived/power_balance.csv`; `freq.npy` + `fft_avg.npy` of all 61 + 61 `D_sweep` runs | `figS1` (2 MB) |
| Table S4 (SM4) | `analysis/diameter_plateaus.py` | `derived/wave_table.csv`, `derived/power_balance.csv` | none |
| Fig. S3 (not in the manuscript) | `figures/figS3.py` | `derived/j_vs_freq_*.csv` | none |
| rebuild `wave_table.csv`, `j_vs_freq_*.csv` | `analysis/wave_extract.py`, `analysis/j_vs_freq.py` | npy of every run | `npy` (0.83 GB) |
| rebuild `power_balance*.csv` | `analysis/power_balance.py`, `analysis/power_balance_J.py` | `data.h5` of every run | `h5-D` (223 GB), `h5-J` (381 GB) |

Only Fig. 2 and SM3.4 need `data.h5` files among the figures; everything else runs from a few `.npy` files and the
derived tables.

## For the maintainer

The share holds symbolic links to the original output folders (`cofe_j_sweep/{bq,no_bq}/output`,
`size_CoFe_spectrum_3/{bq,no_bq_5.5}/output`). If a file is replaced, regenerate the manifest from a copy of the
share layout, e.g. `cd $DATA_DIR && find . -name '*.npy' ! -name initial.npy -o -name data.h5 | sed 's#^\./##' | sort |
xargs sha256sum > MANIFEST.sha256` (and the `.tsv` with sizes).

## Licence of the data

Suggested: CC BY 4.0 for the data set (the code is MIT). This is the authors' choice; state it on the Nextcloud share
and in the archive record.
