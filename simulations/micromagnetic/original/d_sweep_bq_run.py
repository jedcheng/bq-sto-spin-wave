
from magnumnp import *
import torch
# from matplotlib import pyplot as plt
import numpy as np
from scipy.fft import rfft, rfftfreq
# import pandas as pd
# %config InlineBackend.figure_format = 'retina'
import argparse


from matplotlib import pyplot as plt

import h5py



parser = argparse.ArgumentParser()
parser.add_argument("--d", help="Diameter in 10nm", type=int, default=2)
parser.add_argument("--J", help="Electrical current density in 1e11 A/m^2", type=float, default=14)

args = parser.parse_args()
diameter = int(args.d * 10)

J = args.J * 1e11  # Convert to A/m^2


import os 
try:
        os.mkdir("output/" + str(diameter))
except:
        pass

Timer.enable()

dtype = torch.float32

torch.set_default_dtype(dtype)
torch.set_default_device("cuda")
# torch.set_default_device("cpu")

B_BQ = -0.6e-3
A_BL = -0.0e-3


# initialize mesh
n = (int(diameter/2), int(diameter/2), 3)
dx = (2e-9, 2e-9, 2e-9)
mesh = Mesh(n, dx)

state = State(mesh)

state.material = {
        "A": 1.3e-11,
        "P": 0.35,
        "Lambda": 1,
        "epsilon_prime": 0,
        # "mp": state.Tensor((0, 0.01, 0.99)),
        "mp": [0, 0.01, 0.99],
        "d": 2e-9,
        "J": J
        }


# CoFe
r = diameter/2 * 1e-9
x, y, z = mesh.SpatialCoordinate()
free = ((x-r)/r)**2 + ((y-r)/r)**2 < 1 
fixed = ((x-r)/r)**2 + ((y-r)/r)**2 < 1 

free[:,:,1:3] = False
fixed[:,:,0:2] = False


# Set initial magnetization
state.m = state.Constant([0, 0, 0])
# state.m[free] = state.Tensor([0, 0.99, 0.01])
# state.m[fixed] = state.Tensor([0, 0.01, 0.99])
state.m[free] = torch.tensor([0, 0.99, 0.01], dtype=dtype)
state.m[fixed] = torch.tensor([0, 0.01, 0.99], dtype=dtype)



alpha = state.Constant([0.0])
alpha[free] = 0.035
alpha[fixed] = 0.1


Ku = state.Constant([0.0]) 
Ku[free] = 2362
Ku[fixed] = 600000


Ms = state.Constant([0.0])
Ms[free] = 1350e3
Ms[fixed] = 800e3


Ku_axis = state.Constant([1, 0, 0])
# Ku_axis[free] = state.Tensor([1, 1, 0])
# Ku_axis[fixed] = state.Tensor([0, 0, 1])
Ku_axis[free] = torch.tensor([1, 1, 0], dtype=dtype)
Ku_axis[fixed] = torch.tensor([0, 0, 1], dtype=dtype)


state.material['alpha'] = alpha
state.material['Ms'] = Ms
state.material['Ku'] = Ku
state.material['Ku_axis'] = Ku_axis


exchange_free = ExchangeField(free)
exchange_fixed = ExchangeField(fixed)
# rkky      = RKKYField(A_BL, "z", 0, 2)
rkkybq      = BiquadraticRKKYField(B_BQ, "z", 0, 2)
aniso     = UniaxialAnisotropyField()

demag    = DemagField()

torque   = SpinTorqueSlonczewski()

t_save = np.zeros(int((10e-9) / 5e-12)+1)
m_save = np.zeros((int((10e-9) / 5e-12)+1, int(diameter/2), int(diameter/2), 3, 3))



llg = LLGSolver([demag, exchange_free, exchange_fixed, aniso, rkkybq])

llg.relax(state, maxiter=2000)

m_cpu = state.m.cpu().numpy()

# np.save("after_relax.npy", m_cpu)

plt.quiver(m_cpu[::2,::2,0,0], m_cpu[::2,::2,0,1])
plt.imshow(m_cpu[::2,::2,0,2])
plt.colorbar()
plt.savefig("output/" + str(diameter) + "/initial.png", dpi=300)
plt.close()

np.save("output/" + str(diameter) + "/initial.npy", m_cpu)

llg = LLGSolver([demag, aniso, exchange_fixed, exchange_free, torque, rkkybq])
count = 0


while state.t < 10e-9:
        llg.step(state, 5e-12)
        # logger << state
        # state.material["mp"] = state.m[fixed].average()
        # state.material["mp"] = torch.mean(state.m[fixed], dim=( 0)).cpu().tolist()
        m_cpu = state.m.cpu().numpy()

        m_save[count] = m_cpu
        t_save[count] = state.t

        count += 1


with h5py.File("output/" + str(diameter) + "/data.h5", "w") as f:
        f.create_dataset("m", data=m_save, compression="gzip")
        f.create_dataset("t", data=t_save)

Timer.print_report()

m_cpu = state.m.cpu().numpy()
np.save("output/" + str(diameter) + "/final.npy", m_cpu)