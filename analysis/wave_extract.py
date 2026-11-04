"""Wave content of every run: reads heat_map.npy / freq.npy / fft.npy, decomposes the complex M_z FFT map into
azimuthal harmonics e^{i m phi} (m = -3..3), and extracts the dominant m, its power fraction, the radial wavenumber k
(slope of the unwrapped phase over 10-85 % of R), lambda, phase turns and amplitude rings.
Output: output/wave_table.csv (shipped copy: data/derived/wave_table.csv).  Usage: python wave_extract.py
"""
import numpy as np, os, sys, csv
np.set_printoptions(linewidth=220, precision=3, suppress=True)
from paths import run_dir
dx=2e-9

def analyse(path, R_nm):
    """Return f_peak [GHz], dominant azimuthal m, k from phase slope [rad/m], n_rings(amp), n_turns(phase)"""
    try:
        H=np.load(path+"/heat_map.npy"); xf=np.load(path+"/freq.npy"); yf=np.load(path+"/fft.npy")
    except Exception as e:
        return None
    f=xf[np.argmax(np.abs(yf))]/1e9
    N=H.shape[0]; c=(N-1)/2
    x,y=np.mgrid[0:N,0:N]; X=(x-c)*dx; Y=(y-c)*dx   # magnum.np arrays are indexed [x, y]
    r=np.hypot(X,Y); phi=np.arctan2(Y,X)
    H=np.where(np.isnan(H),0,H)
    R=R_nm*1e-9
    # azimuthal decomposition in radial bins
    rb=np.arange(0.5*dx, R, dx)  # bin centres
    ms=[-3,-2,-1,0,1,2,3]
    Fm=np.zeros((len(ms),len(rb)),complex); cnt=np.zeros(len(rb))
    idx=np.clip(((r)/dx).astype(int),0,len(rb)-1)
    inside=r<R-dx
    for i,m in enumerate(ms):
        w=H*np.exp(-1j*m*phi)
        Fm[i]=np.bincount(idx[inside],weights=w[inside].real,minlength=len(rb))[:len(rb)]+1j*np.bincount(idx[inside],weights=w[inside].imag,minlength=len(rb))[:len(rb)]
    cnt=np.bincount(idx[inside],minlength=len(rb))[:len(rb)]
    Fm/=np.maximum(cnt,1)
    power=np.sum(np.abs(Fm)**2*rb,axis=1)
    idom=np.argmax(power); mdom=ms[idom]
    frac=power[idom]/power.sum()
    # phase slope of dominant m over interior 10%..85% of radius
    sel=(rb>0.10*R)&(rb<0.85*R)
    ph=np.unwrap(np.angle(Fm[idom][sel]))
    k,b=np.polyfit(rb[sel],ph,1)
    turns=(ph[-1]-ph[0])/(2*np.pi)/(rb[sel][-1]-rb[sel][0])*R  # extrapolated turns over full radius
    # amplitude rings: count local maxima of azimuthal-mean |H| profile
    amp=np.bincount(idx[inside],weights=np.abs(H[inside]),minlength=len(rb))[:len(rb)]/np.maximum(cnt,1)
    a=amp[1:-1]; nmax=int(np.sum((a>amp[:-2])&(a>amp[2:])&(a>0.2*amp.max())))
    return dict(f=f,m=mdom,frac=frac,k=k,lam_nm=2*np.pi/abs(k)*1e9,turns=turns,nrings=nmax,amp=amp,rb=rb,Fm=Fm,ms=ms)

rows=[]
print("### J sweep, D=320 nm")
for tag in ["bq","no_bq"]:
    for J in np.round(np.arange(3.0,13.01,0.1),1):
        res=analyse(run_dir(tag,J=J),160)
        if res is None: continue
        rows.append(("J",tag,J,320,res["f"],res["m"],res["frac"],res["k"],res["lam_nm"],res["turns"],res["nrings"]))
print("### D sweep")
for tag,J in [("bq",6.0),("no_bq",5.5)]:
    for D in range(100,701,10):
        res=analyse(run_dir(tag,D=D),D/2)
        if res is None: continue
        rows.append(("D",tag,J,D,res["f"],res["m"],res["frac"],res["k"],res["lam_nm"],res["turns"],res["nrings"]))
os.makedirs("output",exist_ok=True)
with open("output/wave_table.csv","w") as fh:
    w=csv.writer(fh); w.writerow(["sweep","tag","J","D_nm","f_GHz","m_dom","m_frac","k_rad_per_m","lambda_nm","turns","n_amp_rings"]); w.writerows(rows)
for r in rows:
    print(f"{r[0]} {r[1]:6s} J={r[2]:4.1f} D={r[3]:3d} f={r[4]:5.1f} m={r[5]:+d} ({r[6]:.2f}) k={r[7]/1e6:7.1f} /um  lam={r[8]:5.1f} nm turns={r[9]:5.2f} rings={r[10]}")
