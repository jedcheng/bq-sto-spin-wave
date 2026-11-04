# Section 5 of `derivation.tex`, step by step: linear spin waves of the gyrating state

Starting point: the result of Section 4. In the co-rotating spin frame, and with the macrospin balance of Section 3 in force, the free layer obeys a plain LLG equation

$$
\dot{\mathbf m}'=-\gamma\,\mathbf m'\times\mathbf H'_{\rm tot}+\alpha\,\mathbf m'\times\dot{\mathbf m}',
\qquad
\mathbf H'_{\rm tot}=\mathbf H_{\rm ex}+\mathbf H_{\rm d}+\mathbf H_{\rm bq}+(M_s+H_B)\,m_{z0}\,\hat{\mathbf z},
$$

and for the uniform state $\mathbf m'_0=(\cos\theta_0,0,\sin\theta_0)$ the total field is **zero**. Section 5 asks what small ripples on this state do. Primes are dropped from here on; everything is in the rotating frame.

## 1. Choose coordinates for the ripple

Write $\mathbf m=\mathbf m_0+\delta\mathbf m$. Because $|\mathbf m|=1$,

$$
1=|\mathbf m_0+\delta\mathbf m|^2=1+2\,\mathbf m_0\cdot\delta\mathbf m+|\delta\mathbf m|^2
\quad\Longrightarrow\quad \mathbf m_0\cdot\delta\mathbf m=O(\delta^2),
$$

so to first order $\delta\mathbf m$ lies in the plane perpendicular to $\mathbf m_0$. Two components suffice. The text picks

$$
\hat{\mathbf e}_1=\hat{\mathbf y}\quad(\text{in-plane, transverse to }\mathbf m_0),\qquad
\hat{\mathbf e}_2=(-\sin\theta_0,\,0,\,\cos\theta_0)\quad(\text{polar: tilts }\mathbf m_0\text{ toward }+z),
$$

and $\delta\mathbf m=a\,\hat{\mathbf e}_1+b\,\hat{\mathbf e}_2$. Checks worth doing once by hand:

- $\hat{\mathbf e}_1\cdot\mathbf m_0=0$, $\hat{\mathbf e}_2\cdot\mathbf m_0=-\sin\theta_0\cos\theta_0+\cos\theta_0\sin\theta_0=0$, $\hat{\mathbf e}_1\cdot\hat{\mathbf e}_2=0$, all unit length.
- $\mathbf m_0\times\hat{\mathbf e}_1=(\cos\theta_0,0,\sin\theta_0)\times(0,1,0)=(-\sin\theta_0,\,0,\,\cos\theta_0)=\hat{\mathbf e}_2$.
- $\mathbf m_0\times\hat{\mathbf e}_2=(\cos\theta_0,0,\sin\theta_0)\times(-\sin\theta_0,0,\cos\theta_0)=(0,\,-\sin^2\theta_0-\cos^2\theta_0,\,0)=-\hat{\mathbf e}_1$.

So $(\mathbf m_0,\hat{\mathbf e}_1,\hat{\mathbf e}_2)$ is a right-handed triad and "$\mathbf m_0\times$" acts on the tangent plane as a rotation by $90^\circ$: $\hat{\mathbf e}_1\to\hat{\mathbf e}_2\to-\hat{\mathbf e}_1$. Also note, for later, the projections

$$
\hat{\mathbf e}_1\cdot\hat{\mathbf z}=0,\qquad \hat{\mathbf e}_2\cdot\hat{\mathbf z}=\cos\theta_0
\qquad\Longrightarrow\qquad \delta m_z=b\cos\theta_0 .
$$

The out-of-plane ripple, which is what the $M_z$ maps show, is the polar component $b$ scaled by $\cos\theta_0$.

## 2. Linearise the equation of motion

Insert $\mathbf m=\mathbf m_0+\delta\mathbf m$ and $\mathbf H_{\rm tot}=\mathbf H_0+\delta\mathbf H$, where $\mathbf H_0\equiv\mathbf H_{\rm tot}(\mathbf m_0)$ and $\delta\mathbf H$ is linear in $\delta\mathbf m$. Keep first order:

$$
\delta\dot{\mathbf m}=-\gamma\big[\mathbf m_0\times\delta\mathbf H+\delta\mathbf m\times\mathbf H_0\big]+\alpha\,\mathbf m_0\times\delta\dot{\mathbf m}.
$$

The term $\delta\mathbf m\times\mathbf H_0$ is the one that normally makes the algebra messy (it carries the static field, i.e. the "Zeeman" or "internal field" contribution to the stiffness). Here $\mathbf H_0=0$, so it is simply absent. That is the payoff of Section 4.2. Dropping damping for the dispersion:

$$
\boxed{\;\delta\dot{\mathbf m}=-\gamma\,\mathbf m_0\times\delta\mathbf H\;}
$$

Two remarks:

- Any component of $\delta\mathbf H$ along $\mathbf m_0$ is killed by the cross product, so only the projections of $\delta\mathbf H$ on $\hat{\mathbf e}_1,\hat{\mathbf e}_2$ are needed.
- The fictitious field $(M_s+H_B)m_{z0}\hat{\mathbf z}$ is a constant; it does not depend on $\mathbf m$ and therefore contributes nothing to $\delta\mathbf H$. It has already done its job by cancelling $\mathbf H_0$.

## 3. The stiffness matrix $\mathsf N$

Take a plane wave $\delta\mathbf m\propto e^{\mathrm i(\mathbf k\cdot\mathbf r-\Omega t)}$ with $\mathbf k$ in the film plane. Every field is linear in $\delta\mathbf m$, so one can write

$$
\delta H_i\equiv\hat{\mathbf e}_i\cdot\delta\mathbf H=-\big(N_{i1}\,a+N_{i2}\,b\big),\qquad i=1,2 .
$$

The minus sign is a convention that makes $N_{ij}$ positive for restoring fields. $\mathsf N$ is the Hessian of the energy per $\mu_0 M_s$: $\delta E/V=\tfrac12\mu_0M_s\,(N_{11}a^2+2N_{12}ab+N_{22}b^2)$. Now go through the three fields.

**Exchange.** $\delta\mathbf H_{\rm ex}=D\nabla^2\delta\mathbf m=-Dk^2\,\delta\mathbf m$. Projection on $\hat{\mathbf e}_1$ gives $-Dk^2 a$, on $\hat{\mathbf e}_2$ gives $-Dk^2 b$:

$$
N_{11}=N_{22}=Dk^2,\qquad N_{12}=0 .
$$

**Biquadratic.** $\delta\mathbf H_{\rm bq}=-H_B\,\delta m_z\,\hat{\mathbf z}=-H_B\,b\cos\theta_0\,\hat{\mathbf z}$. Project: $\hat{\mathbf e}_1\cdot\hat{\mathbf z}=0$, $\hat{\mathbf e}_2\cdot\hat{\mathbf z}=\cos\theta_0$, so $\delta H_1=0$ and $\delta H_2=-H_B\cos^2\theta_0\,b$:

$$
N_{22}\mathrel{+}=H_B\cos^2\theta_0 .
$$

The two factors of $\cos\theta_0$ have different origins: one converts $b$ into $\delta m_z$, the other projects the $z$ field back onto the polar direction. For $\theta_0\to90^\circ$ (magnetisation along the hard axis) the easy-plane stiffness would vanish, as it should.

**Dipolar, Eq. (dip) of the derivation.**

$$
\delta\mathbf H_{\rm d}=-M_s\Big[P\,\delta m_z\,\hat{\mathbf z}+(1-P)\,(\delta\mathbf m\cdot\hat{\mathbf k})\,\hat{\mathbf k}\Big],\qquad P=P(kd).
$$

The first piece is the biquadratic term with $H_B\to M_sP$:

$$
N_{22}\mathrel{+}=M_sP\cos^2\theta_0 .
$$

For the second piece we need $\hat{\mathbf k}$. Let $\varphi_k$ be the angle between $\hat{\mathbf k}$ and the in-plane direction of $\mathbf m_0$ (which is $\hat{\mathbf x}$), so $\hat{\mathbf k}=(\cos\varphi_k,\sin\varphi_k,0)$. Then

$$
\hat{\mathbf e}_1\cdot\hat{\mathbf k}=\sin\varphi_k,\qquad
\hat{\mathbf e}_2\cdot\hat{\mathbf k}=-\sin\theta_0\cos\varphi_k,\qquad
\delta\mathbf m\cdot\hat{\mathbf k}=a\sin\varphi_k-b\sin\theta_0\cos\varphi_k .
$$

The field is $-M_s(1-P)\,(a\sin\varphi_k-b\sin\theta_0\cos\varphi_k)\,\hat{\mathbf k}$. Project it:

$$
\begin{aligned}
\delta H_1&=-M_s(1-P)\big(a\sin\varphi_k-b\sin\theta_0\cos\varphi_k\big)\sin\varphi_k
=-\Big[M_s(1-P)\sin^2\varphi_k\;a\;-\;M_s(1-P)\sin\theta_0\sin\varphi_k\cos\varphi_k\;b\Big],\\[2pt]
\delta H_2&=-M_s(1-P)\big(a\sin\varphi_k-b\sin\theta_0\cos\varphi_k\big)(-\sin\theta_0\cos\varphi_k)
=-\Big[-M_s(1-P)\sin\theta_0\sin\varphi_k\cos\varphi_k\;a\;+\;M_s(1-P)\sin^2\theta_0\cos^2\varphi_k\;b\Big].
\end{aligned}
$$

Reading off the coefficients:

$$
N_{11}\mathrel{+}=M_s(1-P)\sin^2\varphi_k,\qquad
N_{22}\mathrel{+}=M_s(1-P)\sin^2\theta_0\cos^2\varphi_k,\qquad
N_{12}=N_{21}=-M_s(1-P)\sin\theta_0\sin\varphi_k\cos\varphi_k .
$$

Collecting everything:

$$
\boxed{\;
\begin{aligned}
N_{11}&=Dk^2+M_s(1-P)\sin^2\varphi_k,\\
N_{22}&=Dk^2+\big(M_sP+H_B\big)\cos^2\theta_0+M_s(1-P)\sin^2\theta_0\cos^2\varphi_k,\\
N_{12}&=N_{21}=-M_s(1-P)\sin\theta_0\sin\varphi_k\cos\varphi_k .
\end{aligned}\;}
$$

Physical reading. $N_{11}$ is the cost of twisting the in-plane angle: exchange, plus the volume charges $\nabla\cdot\mathbf m$ that an in-plane twist creates when $\mathbf k$ has a component along the twist direction ($\sin\varphi_k\ne0$). $N_{22}$ is the cost of tilting out of plane: exchange, plus the surface charges (the thin-film demag $M_sP$ with its biquadratic twin $H_B$). The cross term couples twist and tilt only through the small volume-charge contribution of the tilt, proportional to $\sin\theta_0(1-P)$.

## 4. Equations of motion for $a$ and $b$

Substitute $\delta\mathbf H=-\mathsf N\delta\mathbf m$ (with the component along $\mathbf m_0$ irrelevant):

$$
\delta\dot{\mathbf m}=\gamma\,\mathbf m_0\times(\mathsf N\delta\mathbf m),
\qquad
\mathsf N\delta\mathbf m=(N_{11}a+N_{12}b)\,\hat{\mathbf e}_1+(N_{21}a+N_{22}b)\,\hat{\mathbf e}_2 .
$$

Apply the $90^\circ$ rule $\mathbf m_0\times\hat{\mathbf e}_1=\hat{\mathbf e}_2$, $\mathbf m_0\times\hat{\mathbf e}_2=-\hat{\mathbf e}_1$:

$$
\delta\dot{\mathbf m}=\gamma(N_{11}a+N_{12}b)\,\hat{\mathbf e}_2-\gamma(N_{21}a+N_{22}b)\,\hat{\mathbf e}_1 .
$$

Match the $\hat{\mathbf e}_1$ and $\hat{\mathbf e}_2$ components with $\delta\dot{\mathbf m}=\dot a\,\hat{\mathbf e}_1+\dot b\,\hat{\mathbf e}_2$:

$$
\boxed{\;\dot a=-\gamma\,(N_{21}a+N_{22}b),\qquad \dot b=\gamma\,(N_{11}a+N_{12}b)\;}
$$

This is the standard structure of precession: the restoring field in the polar direction drives the azimuthal velocity and vice versa, with opposite signs.

## 5. Harmonic ansatz and the dispersion relation

Put $a,b\propto e^{-\mathrm i\Omega t}$, so $\partial_t\to-\mathrm i\Omega$:

$$
\begin{pmatrix}-\mathrm i\Omega+\gamma N_{21}&\gamma N_{22}\\ -\gamma N_{11}&-\mathrm i\Omega-\gamma N_{12}\end{pmatrix}
\begin{pmatrix}a\\ b\end{pmatrix}=0 .
$$

A non-trivial solution needs the determinant to vanish:

$$
(-\mathrm i\Omega+\gamma N_{21})(-\mathrm i\Omega-\gamma N_{12})+\gamma^2N_{11}N_{22}=0
\;\Longrightarrow\;
-\Omega^2+\mathrm i\Omega\gamma(N_{12}-N_{21})-\gamma^2N_{12}N_{21}+\gamma^2N_{11}N_{22}=0 .
$$

With $N_{12}=N_{21}$ the imaginary term drops and

$$
\boxed{\;\Omega^2=\gamma^2\big(N_{11}N_{22}-N_{12}^2\big)=\gamma^2\det\mathsf N\;}
\qquad\text{(Eq. Omega2, Smit–Beljers with zero static field)}
$$

Useful by-products of the same $2\times2$ system:

- **Stability.** $\Omega$ is real iff $\det\mathsf N>0$. With $N_{11},N_{22}>0$ and $N_{12}^2\ll N_{11}N_{22}$ the uniform state is a strict energy minimum; this is the statement used in Section 4.3 to exclude a static texture in the rotating frame.
- **Goldstone mode.** At $k\to0$, $P\to1$ and $N_{11}\to0$, so $\Omega\to0$: rotating all in-plane angles together costs nothing. In the lab this is a shift of the gyration phase.
- **Ellipticity.** The second row gives $b/a=\gamma N_{11}/(-\mathrm i\Omega-\gamma N_{12})\simeq\mathrm i\gamma N_{11}/\Omega$; with $\Omega=\gamma\sqrt{N_{11}N_{22}}$ this is $|b/a|=\sqrt{N_{11}/N_{22}}$, i.e. $|a/b|=\sqrt{N_{22}/N_{11}}$. The ripple is mostly an in-plane twist ($a$) with a much smaller out-of-plane part ($b$); $a$ and $b$ are $90^\circ$ out of phase. Numbers below.

## 6. From a plane wave to the radial wave: the azimuthal average

The wave in the disk propagates radially, $\hat{\mathbf k}=\hat{\mathbf r}$, so $\varphi_k=\phi-\phi_0$ runs through all values around the disk while $\mathbf m_0$ (and hence $\phi_0$) is fixed. The entries of $\mathsf N$ therefore depend on $\phi$ through $\sin^2\varphi_k=\tfrac12(1-\cos2\varphi_k)$ and $\sin\varphi_k\cos\varphi_k=\tfrac12\sin2\varphi_k$. In an azimuthal-harmonic expansion these $\cos2\varphi_k$, $\sin2\varphi_k$ factors couple a mode $e^{\mathrm im\phi}$ to $e^{\mathrm i(m\pm2)\phi}$. The one-armed $m=1$ pattern is therefore not an exact eigenmode of the full problem; it is an eigenmode of the $\phi$-averaged problem, which drops the $m\pm2$ admixture. The size of what is dropped is $M_s(1-P)/2$ against $N_{11}$, i.e. tens of per cent of the small bracket; this is one source of the systematic +7 % noted in Section 8.

Averaging: $\langle\sin^2\varphi_k\rangle=\langle\cos^2\varphi_k\rangle=\tfrac12$, $\langle\sin\varphi_k\cos\varphi_k\rangle=0$, and $\langle N_{12}^2\rangle=\tfrac18M_s^2(1-P)^2\sin^2\theta_0$, which is below $10^{-3}$ of $N_{11}N_{22}$ and is dropped. Result:

$$
\boxed{\;
\Omega(k)^2=\gamma^2\Big[Dk^2+\tfrac12M_s\big(1-P(kd)\big)\Big]
\Big[Dk^2+\big(M_sP(kd)+H_B\big)\cos^2\theta_0+\tfrac12M_s\big(1-P(kd)\big)\sin^2\theta_0\Big]\;}
\qquad\text{(Eq. disp)}
$$

with $\sin\theta_0=m_{z0}=\omega/[\gamma(M_s+H_B)]$ from Section 3. The first bracket is the phase (twist) stiffness, the second the tilt stiffness.

## 7. Limits

Small $kd$: $e^{-x}\simeq1-x+x^2/2$ gives $P\simeq1-kd/2$, so $1-P\simeq kd/2$ and the in-plane dipolar term is $\tfrac12M_s(1-P)\simeq M_skd/4$, **linear in $k$**, while exchange is quadratic. Exchange wins when

$$
Dk^2\gg\frac{M_skd}{4}\quad\Longleftrightarrow\quad k\gg\frac{M_sd}{4D}=44\ \mu\mathrm{m}^{-1}\ (\lambda\ll140\ \mathrm{nm}).
$$

In that limit, and with $M_sP+H_B\simeq M_s+H_B$ dominating the second bracket,

$$
\Omega\simeq\gamma\sqrt{Dk^2\,(M_s+H_B)\cos^2\theta_0}=c\,k,\qquad
c=\gamma\sqrt{D(M_s+H_B)}\,\cos\theta_0 ,
$$

a linear ("sound-like") dispersion, Eq. (sound). The actual wavenumbers ($k=40$–$200\ \mu$m$^{-1}$) sit right at the crossover, so the in-plane dipolar term is comparable to exchange and the full Eq. (disp) must be used.

## 8. Numbers at the measured wavelength (D = 320 nm)

| quantity | with BQ, $\lambda=94$ nm, $f_{\rm sim}=13.09$ GHz | without BQ, $\lambda=92$ nm, $f_{\rm sim}=11.93$ GHz |
|---|---|---|
| $k$, $kd$ | $66.8\ \mu$m$^{-1}$, 0.134 | $68.3\ \mu$m$^{-1}$, 0.137 |
| $P$, $1-P$ | 0.936, 0.064 | 0.935, 0.065 |
| $Dk^2$ | $0.68\times10^5$ A/m | $0.71\times10^5$ A/m |
| $\tfrac12M_s(1-P)$ | $0.43\times10^5$ A/m | $0.44\times10^5$ A/m |
| $N_{11}$ (twist stiffness) | $1.12\times10^5$ A/m | $1.16\times10^5$ A/m |
| $N_{22}$ (tilt stiffness) | $16.1\times10^5$ A/m | $12.6\times10^5$ A/m |
| $N_{12}^2/(N_{11}N_{22})$ (rms) | $2.5\times10^{-4}$ | $4.2\times10^{-4}$ |
| $\Omega/2\pi$ from Eq. (disp) | 14.9 GHz (ratio 1.14) | 13.4 GHz (ratio 1.12) |
| $ck/2\pi$, linear limit | 11.7 GHz | 10.6 GHz |
| exchange-only $\gamma Dk^2/2\pi$ | 2.4 GHz | 2.5 GHz |
| $|a/b|=\sqrt{N_{22}/N_{11}}$ | 3.8 | 3.3 |

Constants: $D=2A/(\mu_0M_s)=1.53\times10^{-11}$ A m, $H_B=3.54\times10^5$ A/m, $\gamma=2.211\times10^5$ m A$^{-1}$ s$^{-1}$.

Two things to take from the table. The tilt stiffness is more than ten times the twist stiffness, so the wave is mostly an in-plane twist with an out-of-plane part 3–4 times smaller in amplitude; the $M_z$ maps show the minor axis of the ellipse. And the two entries of the small bracket are of the same size, which is why neither the pure-exchange nor the pure-linear limit is adequate; the exchange-only guess is 5–6 times too low at this $k$. The two runs in the table sit at the upper end of the spread over all 177 runs (mean ratio 1.07, standard deviation 0.08).

## 9. Points that are easy to trip over

- **Where did the static field go?** In an ordinary Smit–Beljers calculation the stiffness matrix contains the internal field $H_0$ added to both diagonal entries. Here that field is exactly zero in the rotating frame, and the diagonal contains only exchange, dipolar and biquadratic curvature terms. If you try to do the linearisation in the lab frame instead, the static field is non-zero and time dependent through the rotating $\mathbf m_0$, and the algebra does not close.

- **Two different $\cos\theta_0$ factors.** $(M_sP+H_B)\cos^2\theta_0$ comes from $\delta m_z=b\cos\theta_0$ together with $\hat{\mathbf e}_2\cdot\hat{\mathbf z}=\cos\theta_0$. It is not an approximation; it is the exact projection for a tilted equilibrium.

- **$P(kd)$ is for an unpinned, thickness-averaged film.** At $kd\simeq0.13$ it is close to 1 and the correction $1-P\simeq kd/2$ is what feeds the in-plane term. Using the uniform-film demag ($P=1$) would set $N_{11}=Dk^2$ and lose the dipolar part of the twist stiffness, which is 40 % of $N_{11}$ at the measured $k$.

- **The $\varphi_k$ average is the main approximation of the section.** The in-plane dipolar field knows the real-space direction of $\mathbf k$ relative to the in-plane magnetisation, and for a radial wave that angle sweeps $2\pi$ around the disk. Averaging makes the one-armed mode an eigenmode and keeps the problem one dimensional. The exact treatment couples $m=1$ to $m=\pm1,3$ and is not needed at the accuracy of the data.

- **Damping was dropped, not forgotten.** Including $\alpha\mathbf m_0\times\delta\dot{\mathbf m}$ gives $\Omega\to\Omega-\mathrm i\Gamma$ with $\Gamma=\alpha\gamma(N_{11}+N_{22})/2$; this decay rate is what Section 9.5 uses for the wave decay length $\ell=v_g/\Gamma$.
