"""Micromagnetic simulation of the orthogonal spin-torque oscillator with biquadratic (BQ) interlayer coupling.

Cleaned, parameterised version of the scripts that produced the data of the paper (verbatim copies are in
`original/`, the differences are listed in `CHANGES.md`).  The physics is identical to the originals:

  * magnum.np finite-difference LLG, cells 2 x 2 x 2 nm, a disk of diameter D and 3 layers along z:
        layer 0  free layer  CoFe   (Ms 1.35e6 A/m, A 1.3e-11 J/m, Ku 2362 J/m^3 along (1,1,0), alpha 0.035)
        layer 1  spacer      empty  (Ms = 0)
        layer 2  fixed layer CoPt   (Ms 0.8e6 A/m, Ku 6e5 J/m^3 along z, alpha 0.1)
    Exchange is evaluated separately in each layer, so the layers interact only through the dipolar (demag)
    field and the BQ coupling.
  * BQ coupling B12 = -0.6e-3 J/m^2 between layers 0 and 2 (0 for the "no BQ" case; the field term is kept in the
    solver with B12 = 0 so that the term list is the same).
  * Slonczewski spin-transfer torque as an effective field, P = 0.35, Lambda = 1, epsilon' = 0, polariser
    m_p = (0, 0.01, 0.99), thickness d = 2 nm, current density J.
  * Initial state: free layer (0, 0.99, 0.01), fixed layer (0, 0.01, 0.99); relaxed without the spin torque
    (llg.relax, maxiter 2000); then the current is switched on at t = 0 and m is stored every 5 ps for 10 ns.

Output in OUTDIR:  initial.npy (relaxed state), data.h5 (datasets m [t, x, y, z, 3] and t), final.npy.
Arrays are indexed [x, y, z, component].  Note: the stored arrays have 2001 slots; with float64 the time loop makes
2000 steps and the last slot stays zero, so analysis scripts use [:-1] (or explicit frame ranges).

Examples
    python run.py --D 320 --J 6.0 --bq    --precision float64 --outdir output/J_sweep/bq/6.0
    python run.py --D 320 --J 5.5 --no-bq --precision float32 --outdir output/J_sweep/no_bq/5.5
"""
import argparse
import os

import h5py
import numpy as np
import torch
from magnumnp import (Mesh, State, Timer, LLGSolver, ExchangeField, BiquadraticRKKYField,
                      UniaxialAnisotropyField, DemagField, SpinTorqueSlonczewski)

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--D", type=int, default=320, help="disk diameter in nm (even; paper: 320, diameter sweep 100-700 in steps of 10)")
p.add_argument("--J", type=float, default=6.0, help="current density in 1e11 A/m^2")
g = p.add_mutually_exclusive_group()
g.add_argument("--bq", dest="bq", action="store_true", default=True, help="B12 = -0.6e-3 J/m^2 (default)")
g.add_argument("--no-bq", dest="bq", action="store_false", help="B12 = 0")
p.add_argument("--precision", choices=("float32", "float64"), default="float32",
               help="torch default dtype; see README (most published runs: float32, BQ J=6.0 reference: float64)")
p.add_argument("--outdir", default=None, help="output directory (default output/<bq|no_bq>/D<D>_J<J>)")
p.add_argument("--t-end", type=float, default=10e-9, help="simulated time after switch-on (s)")
p.add_argument("--t-save", type=float, default=5e-12, help="output interval (s)")
p.add_argument("--relax-maxiter", type=int, default=2000)
p.add_argument("--init-from", default=None,
               help="optional .npy state loaded before the relaxation (as in original/j_sweep_bq_run_resume.py)")
p.add_argument("--device", default="cuda", help="cuda or cpu")
args = p.parse_args()

tag = "bq" if args.bq else "no_bq"
outdir = args.outdir or os.path.join("output", tag, f"D{args.D}_J{args.J:.1f}")
os.makedirs(outdir, exist_ok=True)

Timer.enable()
dtype = getattr(torch, args.precision)
torch.set_default_dtype(dtype)
torch.set_default_device(args.device)

# ------------------------------------------------------------------ parameters (identical to the originals)
B_BQ = -0.6e-3 if args.bq else 0.0          # J/m^2, biquadratic coupling
J = args.J * 1e11                           # A/m^2
diameter = args.D                           # nm
N = diameter // 2                           # cells per side (2 nm cells)

mesh = Mesh((N, N, 3), (2e-9, 2e-9, 2e-9))
state = State(mesh)
state.material = {
    "A": 1.3e-11,
    "P": 0.35,
    "Lambda": 1,
    "epsilon_prime": 0,
    "mp": [0, 0.01, 0.99],
    "d": 2e-9,
    "J": J,
}

r = diameter / 2 * 1e-9
x, y, z = mesh.SpatialCoordinate()
free = ((x - r) / r) ** 2 + ((y - r) / r) ** 2 < 1
fixed = ((x - r) / r) ** 2 + ((y - r) / r) ** 2 < 1
free[:, :, 1:3] = False                     # free layer = z index 0
fixed[:, :, 0:2] = False                    # fixed layer = z index 2 (index 1 is the empty spacer)

state.m = state.Constant([0, 0, 0])
state.m[free] = torch.tensor([0, 0.99, 0.01], dtype=dtype)
state.m[fixed] = torch.tensor([0, 0.01, 0.99], dtype=dtype)

alpha = state.Constant([0.0]); alpha[free] = 0.035; alpha[fixed] = 0.1
Ku = state.Constant([0.0]);    Ku[free] = 2362;     Ku[fixed] = 600000
Ms = state.Constant([0.0]);    Ms[free] = 1350e3;   Ms[fixed] = 800e3
Ku_axis = state.Constant([1, 0, 0])
Ku_axis[free] = torch.tensor([1, 1, 0], dtype=dtype)      # passed un-normalised, as in the originals
Ku_axis[fixed] = torch.tensor([0, 0, 1], dtype=dtype)
state.material["alpha"] = alpha
state.material["Ms"] = Ms
state.material["Ku"] = Ku
state.material["Ku_axis"] = Ku_axis

exchange_free = ExchangeField(free)
exchange_fixed = ExchangeField(fixed)
rkkybq = BiquadraticRKKYField(B_BQ, "z", 0, 2)
aniso = UniaxialAnisotropyField()
demag = DemagField()
torque = SpinTorqueSlonczewski()

nsave = int(args.t_end / args.t_save) + 1
t_save = np.zeros(nsave)
m_save = np.zeros((nsave, N, N, 3, 3))

# ------------------------------------------------------------------ relaxation (no current)
if args.init_from:
    state.m = torch.tensor(np.load(args.init_from), dtype=dtype, device=args.device)
llg = LLGSolver([demag, exchange_free, exchange_fixed, aniso, rkkybq])      # same term order as the originals
llg.relax(state, maxiter=args.relax_maxiter)
np.save(os.path.join(outdir, "initial.npy"), state.m.cpu().numpy())

# ------------------------------------------------------------------ dynamics with spin torque
llg = LLGSolver([demag, aniso, exchange_fixed, exchange_free, torque, rkkybq])
count = 0
while state.t < args.t_end:
    llg.step(state, args.t_save)
    m_save[count] = state.m.cpu().numpy()
    t_save[count] = state.t
    count += 1

with h5py.File(os.path.join(outdir, "data.h5"), "w") as f:
    f.create_dataset("m", data=m_save, compression="gzip")
    f.create_dataset("t", data=t_save)

Timer.print_report()
np.save(os.path.join(outdir, "final.npy"), state.m.cpu().numpy())
