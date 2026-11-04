"""Two-macrospin LLG model (free + fixed layer) without demagnetising fields.

Used as a control (README "Macrospin"): BQ coupling h_i = 2 B12 (m_i.m_j) m_j / (mu0 Ms_i d), Slonczewski field
a_J = hbar J P / (2 e mu0 Ms d) on both layers, uniaxial anisotropy on both layers, NO demag/shape anisotropy.
Parameters as in the micromagnetic runs.  python twomac.py  -> frequency, <m> at J = 6e11 A/m^2.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft, rfftfreq
mu0=4e-7*np.pi; gam=2.211e5; hbar=1.0546e-34; e=1.602e-19
B=-0.6e-3; d=2e-9; P=0.35
Ms1,a1,Ku1=1.35e6,0.035,2362.; Ms2,a2,Ku2=8e5,0.1,6e5
u1=np.array([1,1,0])/np.sqrt(2); u2=np.array([0,0,1.]); p=np.array([0,0.01,0.99]); p/=np.linalg.norm(p)
def run(J, T=10e-9, m10=(0,0.99,0.01), m20=(0,0.01,0.99)):
    HB1=2*abs(B)/(mu0*Ms1*d); HB2=2*abs(B)/(mu0*Ms2*d)
    aJ1=hbar*J*P/(2*e*mu0*Ms1*d); aJ2=hbar*J*P/(2*e*mu0*Ms2*d)
    Hk1=2*Ku1/(mu0*Ms1); Hk2=2*Ku2/(mu0*Ms2)
    def rhs(t,y):
        m1,m2=y[:3],y[3:]; m12=m1@m2
        H1=-HB1*m12*m2 + Hk1*(m1@u1)*u1 + aJ1*np.cross(m1,p)
        H2=-HB2*m12*m1 + Hk2*(m2@u2)*u2 + aJ2*np.cross(m2,p)
        out=[]
        for m,H,a in ((m1,H1,a1),(m2,H2,a2)):
            mxH=np.cross(m,H); out.append(-gam/(1+a*a)*(mxH+a*np.cross(m,mxH)))
        return np.concatenate(out)
    y0=np.concatenate([np.array(m10)/np.linalg.norm(m10),np.array(m20)/np.linalg.norm(m20)])
    t=np.arange(0,T,5e-12)
    s=solve_ivp(rhs,(0,T),y0,t_eval=t,rtol=1e-8,atol=1e-10,method='DOP853')
    return t,s.y.T,dict(HB1=HB1,HB2=HB2,aJ1=aJ1,aJ2=aJ2,Hk1=Hk1,Hk2=Hk2)
if __name__=='__main__':
    for J in (6e11,):
        t,y,pr=run(J); print({k:f'{v:.3e}' for k,v in pr.items()})
        h=slice(len(t)//2,None); xf=rfftfreq(len(t[h]),5e-12)
        for lab,c in (('free',0),('fixed',3)):
            m=y[h,c:c+3]; Y=np.abs(rfft(m[:,0]-m[:,0].mean())); ip=np.argmax(Y[1:])+1
            print(f"J={J:.1e} {lab}: <m>={np.round(m.mean(0),3)}  in-plane amp={np.hypot(m[:,0],m[:,1]).mean():.3f}  f={xf[ip]/1e9:.2f} GHz")
        for t0 in (0.1,0.2,0.5,1,2,5): i=np.argmin(abs(t-t0*1e-9)); print(f"  t={t0} ns free={np.round(y[i,:3],3)} fixed={np.round(y[i,3:],3)}")
