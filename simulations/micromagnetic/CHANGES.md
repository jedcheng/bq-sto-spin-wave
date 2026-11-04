# run.py versus the original scripts

`original/` holds verbatim copies of the scripts that produced the published data:

| original file | produced | precision | B12 (J/m^2) | command line in the original job script |
|---|---|---|---|---|
| `j_sweep_bq_run.py` | current sweep with BQ, D = 320 nm | float64 | -0.6e-3 | `run.py --d 32 --J <J>` (bulk job) |
| `j_sweep_bq_run_resume.py` | variant of the above that loads `after_relax.npy` before relaxing | float64 | -0.6e-3 | - |
| `j_sweep_no_bq_run.py` | current sweep without BQ, D = 320 nm | float32 | 0 | `run.py --d 32 --J <J>` |
| `d_sweep_bq_run.py` | diameter sweep with BQ, J = 6.0 | float32 | -0.6e-3 | `run.py --d <D/10> --J 6.0` |
| `d_sweep_no_bq_run.py` | diameter sweep without BQ, J = 5.5 | float32 | 0 | `run.py --d <D/10> --J 5.5` |

plus the job scripts (`*_run.sh`) and the FFT post-processing (`*_plot_fft_mx.py`, `*_fft_mz_maps.py`), which are
merged into `analysis/postprocess_fft.py`.

## Physics: unchanged

All five originals share the same mesh, materials, field terms, initial state, relaxation, time step and output;
they differ only in B12, the torch dtype and the output folder name (`diff original/j_sweep_bq_run.py original/d_sweep_bq_run.py`).
`run.py` keeps, line for line:

* mesh `(D/2, D/2, 3)` cells of 2 nm; disk mask `((x-r)/r)^2 + ((y-r)/r)^2 < 1`, free layer z = 0, fixed layer z = 2;
* material dict (A, P, Lambda, epsilon', m_p, d, J) and the per-region alpha, Ku, Ms, Ku_axis (incl. the
  un-normalised free-layer axis (1, 1, 0)); Ku(free) = 2362 J/m^3 in every original, including both diameter sweeps;
* the field terms and their order in both `LLGSolver` calls (relaxation without torque, dynamics with torque),
  the BQ term kept with B12 = 0 for "no BQ";
* `llg.relax(state, maxiter=2000)`, then `llg.step(state, 5e-12)` while `t < 10 ns`, the same arrays and the same
  `data.h5` layout (gzip, datasets `m`, `t`), `initial.npy`, `final.npy`.

## What changed (none of it affects the computed magnetisation)

1. Command line: `--D` in nm (the original `--d` was in units of 10 nm), `--bq/--no-bq`, `--precision`,
   `--outdir`, `--t-end`, `--t-save`, `--relax-maxiter`, `--init-from` (reproduces `run_resume.py`), `--device`.
   Defaults reproduce the original values.
2. Output goes to `--outdir` instead of `output/<J>` or `output/<D>`.
3. Removed: the diagnostic quiver/imshow PNG of the relaxed state and the extra copy `after_relax.npy` in the
   working directory (the same array is saved as `initial.npy`), unused imports, commented-out code.
4. `from magnumnp import *` replaced by explicit imports.

## Precision of the published runs (please read)

The script of the current sweep with BQ was edited over time. The file in `original/` (float64) dates from the re-run
of the reference case J = 6.0. Checked on the stored data (every `data.h5` is written as float64): in 503 of the 504
published runs every stored value of m is exactly representable in float32, i.e. they were computed in float32;
only `J_sweep/bq/6.0` was computed in float64 (also visible in its file size, 1.8 GB instead of ~1.0 GB). To reproduce the published data use
`--precision float64` for BQ J = 6.0 and `--precision float32` otherwise (the default). It is not recorded whether the
J = 6.0 re-run used `run.py` or `run_resume.py`.

Re-run check (H100 MIG slice, see the Verification section of the top README): `run.py --D 320 --J 5.5 --no-bq`
(float32) reproduces the published run bit for bit; `run.py --D 320 --J 6.0 --bq --precision float64` reproduces its
frequency (13.2 vs 13.1 GHz, one FFT bin), mode (m = +1, 97 %) and wavelength (93 vs 94 nm) but not its bits: the
relaxed state differs by 2e-4. The published `J_sweep/bq/6.0/initial.npy` is identical to the `after_relax.npy` left
in the original working directory, so the published J = 6.0 run was most likely made with `run_resume.py`, i.e.
starting the relaxation from a stored state. To try a bit-identical re-run (not tested):
`run.py --D 320 --J 6.0 --bq --precision float64 --init-from $DATA_DIR/J_sweep/bq/6.0/initial.npy`.
