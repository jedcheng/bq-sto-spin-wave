"""
Analytical model of the self-generated standing spin wave in the orthogonal (BQ-coupled) STO.

Physics
-------
1. Free layer gyrates in-plane at angular frequency w with a small tilt mz0 = w / (gamma (Ms + H_B)),
   H_B = 2|B12|/(mu0 Ms d) being the biquadratic field coefficient (acts like extra "demag" along z).
2. In the spin frame co-rotating with the free-layer magnetisation the gyrating state is static and
   the net internal field is zero (fictitious field cancels demag + BQ).  Small excitations obey the
   ordinary thin-film dipole-exchange dispersion of an in-plane film with ZERO internal field:
       Omega(k)^2 = gamma^2 [ D k^2 + Ms (1-P) <sin^2 phi_k> ] [ D k^2 + (Ms P + H_B) cos^2 theta0 + Ms (1-P) <cos^2 phi_k> sin^2 theta0 ]
   with D = 2A/(mu0 Ms), P(kd) = (1 - exp(-kd))/(kd), <sin^2 phi_k> = <cos^2 phi_k> = 1/2 for a radial wave.
   Small k:  Omega ~ c k  (Goldstone / "spin-superfluid" sound),  c = gamma sqrt(D (Ms + H_B)).
3. The static radial stray field of the PMA fixed layer (edge) is, in the co-rotating spin frame, an rf drive
   at frequency w with azimuthal index m = +-1.  It resonantly launches an inward m = 1 spiral spin wave at
   the wavenumber k* fixed by  Omega(k*) = w.
4. Reflection at the disk centre gives a J_1(k r) e^{i phi} standing pattern (nodes every lambda/2).
   Cavity condition k_n R ~ j'_{1,n} (free-spin / Neumann edge) picks the plateau frequencies
       f_n = Omega(j'_{1,n}/R) / 2 pi.
   Step height per extra ring: Delta f ~ c / (2R).
"""
import numpy as np, csv, os
from paths import DERIVED
from scipy.optimize import brentq
from scipy.special import jnp_zeros, jn_zeros
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

mu0=4e-7*np.pi; gam=2.211e5           # m A^-1 s^-1  (omega = gamma * H, H in A/m)
Ms=1.35e6; A=1.3e-11; d=2e-9; B12=0.6e-3
D=2*A/(mu0*Ms); HB=2*B12/(mu0*Ms*d)

def P(kd): return np.where(kd>1e-9,(1-np.exp(-kd))/np.maximum(kd,1e-12),1.0)
def Omega(k,w,bq,sa=0.5,sb=0.5):
    hb=HB if bq else 0.0
    mz=np.clip(w/(gam*(Ms+hb)),0,0.99); c2=1-mz**2; s2=mz**2
    p=P(k*d)
    Haa=D*k**2+Ms*(1-p)*sa
    Hbb=(hb+Ms*p)*c2+D*k**2+Ms*(1-p)*sb*s2
    return gam*np.sqrt(Haa*Hbb)
def f_of_k(k,bq):
    """self-consistent gyration frequency [GHz] at which Omega(k,w)=w"""
    try: w=brentq(lambda w: Omega(k,w,bq)-w,2*np.pi*0.2e9,2*np.pi*80e9)
    except ValueError: return np.nan
    return w/2/np.pi/1e9
def c_sound(bq): return gam*np.sqrt(D*(Ms+(HB if bq else 0)))

if __name__=="__main__":
    print(f"D = {D:.3e} A m ; mu0 H_B = {mu0*HB:.3f} T ; c = {c_sound(False):.0f} m/s (no BQ), {c_sound(True):.0f} m/s (BQ)")
    print(f"step per ring c/(2R), R=160 nm : {c_sound(False)/(2*160e-9)/1e9:.2f} GHz (no BQ), {c_sound(True)/(2*160e-9)/1e9:.2f} GHz (BQ)")
    rows=list(csv.DictReader(open(os.path.join(DERIVED,"wave_table.csv"))))
    fig,axs=plt.subplots(1,3,figsize=(18,5.5))
    kk=np.linspace(5e6,2.2e8,200)
    for bq,col,lab in [(True,"C3","with BQ"),(False,"C0","without BQ")]:
        tag="bq" if bq else "no_bq"
        obs=[(abs(float(r["k_rad_per_m"])),float(r["f_GHz"]),r["sweep"]) for r in rows if r["tag"]==tag and float(r["m_frac"])>0.9 and float(r["f_GHz"])>7]
        kJ=[o[0] for o in obs if o[2]=="J"]; fJ=[o[1] for o in obs if o[2]=="J"]
        kD=[o[0] for o in obs if o[2]=="D"]; fD=[o[1] for o in obs if o[2]=="D"]
        axs[0].plot(np.array(kJ)/1e6,fJ,"o",ms=4,color=col,label=f"{lab}: J sweep (D=320 nm)")
        axs[0].plot(np.array(kD)/1e6,fD,"s",ms=4,mfc="none",color=col,label=f"{lab}: D sweep")
        axs[0].plot(kk/1e6,[f_of_k(k,bq) for k in kk],"-",color=col,lw=2,label=f"{lab}: theory  Ω(k)=ω")
        axs[0].plot(kk/1e6,c_sound(bq)*kk/2/np.pi/1e9,":",color=col,lw=1,label=f"{lab}: c k / 2π")
        ratio=np.array([f_of_k(k,bq)/f for k,f,_ in obs])
        print(f"{lab}: N={len(obs)}  f_theory/f_obs = {ratio.mean():.3f} +- {ratio.std():.3f}")
    axs[0].set_xlabel("radial wavenumber k (rad/µm)"); axs[0].set_ylabel("gyration frequency f (GHz)")
    axs[0].set_xlim(0,220); axs[0].set_ylim(0,40); axs[0].legend(fontsize=8); axs[0].set_title("dispersion: observed (k from phase map) vs theory")
    # plateau ladder
    R=160e-9
    for i,(bq,col,lab) in enumerate([(True,"C3","with BQ"),(False,"C0","without BQ")]):
        tag="bq" if bq else "no_bq"
        J=[];F=[]
        for r in csv.DictReader(open(os.path.join(DERIVED,f"j_vs_freq_{tag}.csv"))):
            J.append(float(r["J"]));F.append(float(r["Frequency"]))
        ax=axs[1+i]; ax.plot(J,F,"o-",ms=3,color=col,lw=0.8,label="simulation")
        for n in range(1,10):
            fn=f_of_k(jnp_zeros(1,n)[-1]/R,bq); fd=f_of_k(jn_zeros(1,n)[-1]/R,bq)
            ax.axhline(fn,color="k",lw=1,ls="--"); ax.text(2.6,fn+0.3,f"n={n}",fontsize=8)
            ax.axhspan(fn,fd,color="gray",alpha=0.12)
        ax.set_xlim(2.5,10); ax.set_ylim(0,32); ax.set_xlabel(r"J (10$^{11}$ A/m$^2$)"); ax.set_ylabel("f (GHz)")
        ax.set_title(f"{lab}: plateaus vs f_n=Ω(j'_1n/R)/2π (dashed), band up to j_1n (grey)")
        ax.legend(loc="upper left",fontsize=8)
    plt.tight_layout(); os.makedirs("output",exist_ok=True); plt.savefig(os.path.join("output","standing_wave_theory.png"),dpi=150)
    print("saved output/standing_wave_theory.png")
    print("\nPredicted plateau (edge-antinode) frequencies f_n [GHz], D=320 nm:")
    for n in range(2,9):
        print(f" n={n}: BQ {f_of_k(jnp_zeros(1,n)[-1]/R,True):5.1f}   noBQ {f_of_k(jnp_zeros(1,n)[-1]/R,False):5.1f}")
