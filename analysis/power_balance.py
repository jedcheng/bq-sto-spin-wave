"""Exact LLG power balance and wave content over the last 3 ns (frames 1399-1999) of each run of the diameter sweep.
Needs data.h5 of every run.  Main columns: f (mean gyration rate), mz (tilt), coh = |<m_+>| (gyration coherence),
aip1 (in-plane mean-square wave amplitude in the rotating frame at 0.5-1.5 f, <a^2> of the paper), fmac, bal
(damping/torque balance).  Output: output/power_balance.csv (shipped: data/derived/power_balance.csv).
Usage: python power_balance.py   (12 parallel workers)
"""
import numpy as np, h5py, sys, time
import os
from paths import run_dir
mu0=4e-7*np.pi; gam=2.211e5; Ms=1.35e6; d=2e-9; hbar=1.0546e-34; e=1.602e-19; alpha=0.035; Pol=0.35
def run(p,J,Dn):
    """p: run directory containing data.h5; J in A/m^2; Dn: label stored in column D"""
    aJ=hbar*J/(mu0*Ms*e*d)*Pol/2
    with h5py.File(p+"/data.h5") as f:
        t=f["t"][1399:2000]; m=f["m"][1399:2000,:,:,0,:].astype(np.float64)
    dt=t[1]-t[0]
    ins=np.linalg.norm(m[0],axis=-1)>0.5
    m=m[:,ins,:]                                   # (T, Ncell, 3)
    mp=m[...,0]+1j*m[...,1]; mz=m[...,2]
    Phi=np.unwrap(np.angle(mp.mean(1)))            # mean gyration phase
    wf=np.polyfit(t,Phi,1)[0]
    # time derivative (central differences, drop end points)
    md=(m[2:]-m[:-2])/(2*dt); mc=m[1:-1]
    D2=np.mean(np.sum(md**2,-1))                                   # <|mdot|^2>
    Lz=np.mean(mc[...,0]*md[...,1]-mc[...,1]*md[...,0])            # <(m x mdot)_z>
    mip2=np.mean(mc[...,0]**2+mc[...,1]**2)
    # rotating spin frame
    mr=mp*np.exp(-1j*Phi)[:,None]
    T=len(t)-1; mr=mr[:T]; mzz=mz[:T]
    Fr=np.fft.fft(mr-mr.mean(0),axis=0); Fz=np.fft.rfft(mzz-mzz.mean(0),axis=0)
    fr=np.fft.fftfreq(T,dt); frz=np.fft.rfftfreq(T,dt); f0=abs(wf)/2/np.pi
    # <|d/dt|^2> per band (Parseval): sum |F|^2 (2 pi f)^2 / T^2
    def band_c(lo,hi):
        s=(np.abs(fr)>=lo*f0)&(np.abs(fr)<hi*f0); return np.mean(np.sum(np.abs(Fr[s])**2*(2*np.pi*fr[s,None])**2,0))/T**2
    def band_z(lo,hi):
        s=(frz>=lo*f0)&(frz<hi*f0); return 2*np.mean(np.sum(np.abs(Fz[s])**2*(2*np.pi*frz[s,None])**2,0))/T**2
    def amp_c(lo,hi):
        s=(np.abs(fr)>=lo*f0)&(np.abs(fr)<hi*f0); return np.mean(np.sum(np.abs(Fr[s])**2,0))/T**2
    def amp_z(lo,hi):
        s=(frz>=lo*f0)&(frz<hi*f0); return 2*np.mean(np.sum(np.abs(Fz[s])**2,0))/T**2
    w2=wf**2
    res=dict(D=Dn,f=f0/1e9,sense=np.sign(wf),bal=alpha*D2/(gam*aJ*abs(Lz)),
             W=D2/(w2*mip2)-1, mip2=mip2, mz=mz.mean(), coh=np.abs(mp.mean(1)).mean(),
             Wz1=band_z(0.5,1.5)/w2, Wip1=band_c(0.5,1.5)/w2, Wz2=band_z(1.5,2.5)/w2, Wip2=band_c(1.5,2.5)/w2,
             Wlow=(band_z(0,0.5)+band_c(0,0.5))/w2, Whigh=(band_z(2.5,99)+band_c(2.5,99))/w2,
             az1=amp_z(0.5,1.5), aip1=amp_c(0.5,1.5), aip_static=np.mean(np.abs(mr.mean(0)-mr.mean())**2),
             fmac=gam*aJ/alpha/2/np.pi/1e9)
    return res
if __name__=="__main__":
    from joblib import Parallel, delayed
    tasks=[("bq",6.0e11,D) for D in range(100,701,10)]+[("no_bq",5.5e11,D) for D in range(100,701,10)]
    def safe(s,J,D):
        try: r=run(run_dir(s,D=D),J,D); r["tag"]=s; return r
        except Exception as ex: return dict(tag=s,D=D,err=str(ex))
    out=Parallel(n_jobs=12)(delayed(safe)(*tk) for tk in tasks)
    import csv
    keys=["tag","D","f","sense","bal","W","mip2","mz","coh","Wz1","Wip1","Wz2","Wip2","Wlow","Whigh","az1","aip1","aip_static","fmac","err"]
    os.makedirs("output",exist_ok=True)
    with open("output/power_balance.csv","w") as fh:
        w=csv.DictWriter(fh,fieldnames=keys); w.writeheader()
        for r in out: w.writerow({k:r.get(k,"") for k in keys})
    print("done")
