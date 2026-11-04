"""J scan of the two-macrospin model (twomac.py): free-layer tilt, in-plane amplitude, frequency.  Writes jscan.npy."""
import numpy as np
from twomac import run
from scipy.fft import rfft, rfftfreq
res=[]
for J in [0.5,1,1.5,2,2.5,3,3.5,4,5,6,8,10]:
    t,y,pr=run(J*1e11,T=8e-9); h=slice(len(t)//2,None); m=y[h,:3]; m2=y[h,3:]
    xf=rfftfreq(len(t[h]),5e-12); Y=np.abs(rfft(m[:,0]-m[:,0].mean())); ip=np.argmax(Y[1:])+1
    f=xf[ip]/1e9 if Y[ip]>1e-3*len(t[h]) else 0.0
    res.append((J,m[:,2].mean(),np.hypot(m[:,0],m[:,1]).mean(),f,np.hypot(m2[:,0],m2[:,1]).mean()))
    print(f"J={J:4.1f}e11: free <mz>={m[:,2].mean():.3f} in-plane amp={np.hypot(m[:,0],m[:,1]).mean():.3f} f={f:5.2f} GHz | fixed in-plane amp={np.hypot(m2[:,0],m2[:,1]).mean():.3f}", flush=True)
np.save('jscan.npy',np.array(res))
