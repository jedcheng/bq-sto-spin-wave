"""Check of the co-rotating frame (SM3.4) on the last 5 ns of a reference run: gyration sense and rate, local rotation
rate and tilt versus r, signed-frequency power of the azimuthal components of M_z, in-plane FFT map.
Needs data.h5.  Usage: python frame_check.py bq 6.0   |   python frame_check.py no_bq 5.5
Fix with respect to the original: the definitions of X, Y, r, phi had been swallowed by a comment (the original
script could not run); they are restored unchanged.
"""
import numpy as np, h5py, sys, time
np.set_printoptions(linewidth=200, precision=3, suppress=True)
import os
from paths import run_dir
tag=sys.argv[1]; J=sys.argv[2]
p=run_dir(tag,J=float(J))
t0=time.time()
with h5py.File(p+"/data.h5") as f:
    t=f["t"][1000:2000]; m=f["m"][1000:2000,:,:,0,:]
print("read",time.time()-t0,"s", m.shape)
N=m.shape[1]; c=(N-1)/2; x,y=np.mgrid[0:N,0:N]  # arrays indexed [x,y]
X=(x-c)*2e-9; Y=(y-c)*2e-9; r=np.hypot(X,Y); phi=np.arctan2(Y,X)
inside=np.linalg.norm(m[0],axis=-1)>0.5
dt=t[1]-t[0]; nt=len(t)
mx=m[...,0]; my=m[...,1]; mz=m[...,2]
# 1. gyration sense & frequency from mean in-plane vector
Mx=mx[:,inside].mean(1); My=my[:,inside].mean(1); Mz=mz[:,inside].mean(1)
ang=np.unwrap(np.arctan2(My,Mx)); wf=np.polyfit(t,ang,1)[0]
print(f"[1] mean in-plane angular velocity/2pi = {wf/2/np.pi/1e9:+.3f} GHz  (positive = counter-clockwise, x->y);  <mz>={Mz.mean():+.4f}, |m_ip|={np.hypot(Mx,My).mean():.3f}")
# complex FFT of mean m+ = mx+i my : sign of peak frequency
fr=np.fft.fftfreq(nt,dt)/1e9
S=np.fft.fft(Mx+1j*My); ipk=np.argmax(np.abs(S)); print(f"    complex FFT of <mx>+i<my> peaks at f={fr[ipk]:+.2f} GHz")
# 2. local rotation rate profile and local mean mz profile
angloc=np.unwrap(np.arctan2(my,mx),axis=0); wloc=(angloc[-1]-angloc[0])/(t[-1]-t[0])/2/np.pi/1e9
mzbar=mz.mean(0)
rb=np.arange(0,N//2)*2e-9; idx=np.minimum((r/2e-9).astype(int),len(rb)-1)
prof=lambda A: np.array([np.nanmean(A[(idx==i)&inside]) if np.any((idx==i)&inside) else np.nan for i in range(len(rb))])
print("[2] r(nm):      ", (rb[::8]*1e9).astype(int))
print("    local f(r): ", prof(wloc)[::8])
print("    <mz>(r):    ", prof(mzbar)[::8])
print("    <|m_ip|>(r):", prof(np.hypot(mx,my).mean(0))[::8])
# 3. azimuthal decomposition of mz(t) -> complex time FFT (signed frequency) for m=-2..2
ms=[-2,-1,0,1,2]
pk=np.argmax(np.abs(np.fft.rfft(Mx-Mx.mean())))  # gyration bin (rfft index)
fpk=np.fft.rfftfreq(nt,dt)[pk]/1e9
print(f"[3] gyration f = {fpk:.2f} GHz. Power of mz azimuthal components at signed frequencies (interior r<0.85R):")
sel=inside&(r<0.85*160e-9)
mzc=mz[:,sel]-mz[:,sel].mean(0)
for mm in ms:
    comp=(mzc*np.exp(-1j*mm*phi[sel])).mean(1)       # e^{i m phi} coefficient vs t
    F=np.fft.fft(comp); Pw=np.abs(F)**2
    out=[]
    for fq in [-2*fpk,-fpk,0,fpk,2*fpk]:
        b=np.argmin(np.abs(fr-fq)); wnd=slice(max(b-3,0),b+4); out.append(Pw[wnd].max())
    print(f"    m={mm:+d}: P(-2f)={out[0]:9.3g} P(-f)={out[1]:9.3g} P(0)={out[2]:9.3g} P(+f)={out[3]:9.3g} P(+2f)={out[4]:9.3g}")
print("    convention: component ~ exp(i m phi) * exp(+i 2pi f t) ; pattern angular velocity = -2pi f/m")
# 4. static (time-averaged) mz map azimuthal content
for mm in ms:
    comp=((mzbar-np.nanmean(mzbar[inside]))[sel]*np.exp(-1j*mm*phi[sel]))
    # radial profile of this static component
    cprof=np.array([comp[(idx[sel]==i)].mean() if np.any(idx[sel]==i) else np.nan for i in range(len(rb))])
    print(f"[4] static <mz> m={mm:+d}: |amp| profile every 16nm:", np.abs(cprof)[::8])
# 5. in-plane FFT map at f: is the in-plane rotation uniform in space (amplitude & phase)?
Mp=mx+1j*my
F=np.fft.fft(Mp-Mp.mean(0),axis=0)
b=np.argmin(np.abs(fr-np.sign(wf)*fpk))
amp=np.abs(F[b]); ph=np.angle(F[b])
print(f"[5] in-plane m+ FFT at signed gyration freq: amplitude profile:", prof(np.where(inside,amp,np.nan))[::8]/nt)
phrel=np.angle(F[b]*np.exp(-1j*np.angle(F[b][inside].mean())))
print(f"    phase (rad) profile rel. to mean:", prof(np.where(inside,phrel,np.nan))[::8])
os.makedirs("output",exist_ok=True); np.savez(f"output/frame_{tag}_{J}.npz",fr=fr,mzbar=mzbar,inside=inside,Fip=F[b],wloc=wloc)
