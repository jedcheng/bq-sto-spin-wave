# Simulations

## micromagnetic/ (all data of the paper)

`run.py` is one LLG simulation of the CoPt(2)/Cr(2)/CoFe(2) disk with [magnum.np](https://gitlab.com/magnum.np/magnum.np)
(parameters: Table I of the paper and the docstring of `run.py`). The data set consists of four sweeps:

| sweep | command per run | runs |
|---|---|---|
| current, with BQ | `run.py --D 320 --J <J> --bq` | J = 1.0 ... 20.0 step 0.1 (20.1 was submitted but not completed) |
| current, without BQ | `run.py --D 320 --J <J> --no-bq` | J = 1.0 ... 20.0 step 0.1 |
| diameter, with BQ | `run.py --D <D> --J 6.0 --bq` | D = 100 ... 700 nm step 10 |
| diameter, without BQ | `run.py --D <D> --J 5.5 --no-bq` | D = 100 ... 700 nm step 10 |

All published runs used float32 (the default), except the BQ reference run J = 6.0, which was re-run with
`--precision float64` (checked from the stored values: those of every other run are exactly representable in
float32). See `micromagnetic/CHANGES.md` for the comparison with the original scripts (`micromagnetic/original/`).
After a run, `python analysis/postprocess_fft.py <outdir>` makes the `.npy` files used by the figures.

PJM job examples (Fujitsu PJM scheduler; adapt the resource group and `jobs/env_example.sh` to your cluster):
`jobs/single.sh` (one run), `jobs/j_sweep.sh` and `jobs/d_sweep.sh` (bulk jobs, one sub-job per run).

### Installing magnum.np

Install the latest magnum.np (the runs for the paper used v2.0.1 plus a small local patch that does not change the
physics, see `micromagnetic/CHANGES.md`; the current `main` branch has not been re-tested against the published data):

```bash
pip install git+https://gitlab.com/magnum.np/magnum.np@main
```

PyTorch 2.3.1 (CUDA 12.1 build) and Python 3.11.9 were used; pinned versions of the other packages are in
`requirements-sim.txt`.
`BiquadraticRKKYField` and `SpinTorqueSlonczewski` are part of the upstream v2.0.1 release.

### Cost

One D = 320 nm run (10 ns, 2000 outputs) takes about 8-9 min of dynamics on one H100 slice for float32, the
relaxation adds several minutes; D = 700 nm takes about 40 min. The `data.h5` of a D = 320 nm run is ~1.0 GB (float32
run) or 1.8 GB (float64 run); the four sweeps total 604 GB. A GPU with >= 12 GB memory is sufficient for D = 320 nm.

## macrospin/ (control, not a figure of the paper)

`twomac.py` integrates two coupled macrospins (free + fixed layer) with BQ coupling, anisotropy and the Slonczewski
field but **without** demagnetising fields; `jscan.py` scans J. It was used as a negative control for the launch
mechanism (without the stray/demag fields there is no wave: f = 17.98 GHz, <m_z> = 0.98 at J = 6e11 A/m^2).
The analytical macrospin model of SM2 (f_mac, M_z0, Table S2) is in `analysis/macrospin_table.py`.
