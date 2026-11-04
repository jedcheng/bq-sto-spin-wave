# Section 4.1 of `derivation.tex`, step by step: the co-rotating spin frame

Goal of the section: rewrite the LLG equation in a frame in which the gyrating macrospin of Section 3 is at rest. Everything follows from four properties of one rotation matrix, so those are listed first.

## 0. The rotation matrix and the four facts you need

Rotation about $\hat{\mathbf z}$ by the angle $\Phi$ (counter-clockwise seen from $+z$):

$$
\mathsf R_z(\Phi)=\begin{pmatrix}\cos\Phi&-\sin\Phi&0\\ \sin\Phi&\cos\Phi&0\\ 0&0&1\end{pmatrix},\qquad
\mathsf R_z^{-1}=\mathsf R_z^{\mathsf T}=\mathsf R_z(-\Phi).
$$

Abbreviate $\mathsf R\equiv\mathsf R_z(\Phi)$ with $\Phi=\omega_f t$ (so $\Phi$ depends on time but **not on position**). The facts:

**(F1) $\hat{\mathbf z}$ is left alone.** The third column and third row are $(0,0,1)$, so $\mathsf R\hat{\mathbf z}=\hat{\mathbf z}$ and $\mathsf R^{\mathsf T}\hat{\mathbf z}=\hat{\mathbf z}$. Consequently the $z$ component of any vector is unchanged: $(\mathsf R\mathbf v)_z=v_z$.

**(F2) Dot products are preserved.** $(\mathsf R\mathbf a)\cdot(\mathsf R\mathbf b)=\mathbf a^{\mathsf T}\mathsf R^{\mathsf T}\mathsf R\,\mathbf b=\mathbf a\cdot\mathbf b$.

**(F3) Cross products are preserved** (proper rotation, $\det\mathsf R=+1$):

$$
(\mathsf R\mathbf a)\times(\mathsf R\mathbf b)=\mathsf R\,(\mathbf a\times\mathbf b),
\qquad\text{equivalently}\qquad
\mathsf R^{\mathsf T}(\mathbf a\times\mathbf b)=(\mathsf R^{\mathsf T}\mathbf a)\times(\mathsf R^{\mathsf T}\mathbf b).
$$

This is the one that lets you pull $\mathsf R^{\mathsf T}$ through every $\mathbf m\times(\dots)$ in the LLG equation.

**(F4) Time derivative of the matrix.** Differentiate the entries:

$$
\dot{\mathsf R}=\omega_f\begin{pmatrix}-\sin\Phi&-\cos\Phi&0\\ \cos\Phi&-\sin\Phi&0\\ 0&0&0\end{pmatrix}
=\omega_f\,\mathsf J\,\mathsf R=\omega_f\,\mathsf R\,\mathsf J,
\qquad
\mathsf J=\begin{pmatrix}0&-1&0\\ 1&0&0\\ 0&0&0\end{pmatrix}.
$$

Multiply $\mathsf J\mathsf R$ out to check the middle equality: the first row of $\mathsf J\mathsf R$ is minus the second row of $\mathsf R$, i.e. $(-\sin\Phi,-\cos\Phi,0)$, and the second row is the first row of $\mathsf R$. $\mathsf J$ is the matrix of "$\hat{\mathbf z}\times$": $\mathsf J\mathbf v=(-v_y,v_x,0)=\hat{\mathbf z}\times\mathbf v$. So for any constant vector $\mathbf v$,

$$
\dot{\mathsf R}\,\mathbf v=\omega_f\,\hat{\mathbf z}\times(\mathsf R\mathbf v)=\omega_f\,\mathsf R\,(\hat{\mathbf z}\times\mathbf v).
\tag{F4}
$$

($\mathsf J$ and $\mathsf R$ commute because both are rotations about the same axis; equivalently $\hat{\mathbf z}\times(\mathsf R\mathbf v)=(\mathsf R\hat{\mathbf z})\times(\mathsf R\mathbf v)=\mathsf R(\hat{\mathbf z}\times\mathbf v)$ by F1 and F3.)

## 1. What is being rotated, and what is not

The transformation of the text, Eq. (rot):

$$
\mathbf m(\mathbf r,t)=\mathsf R_z(\omega_f t)\,\mathbf m'(\mathbf r,t)
\qquad\Longleftrightarrow\qquad
\mathbf m'(\mathbf r,t)=\mathsf R_z^{\mathsf T}(\omega_f t)\,\mathbf m(\mathbf r,t).
$$

Read it carefully: the **same** $\mathsf R$ is applied in every cell, and the position argument $\mathbf r=(r,\phi)$ is **not** changed. Only the magnetisation *vectors* are turned back by $-\Phi$. This is the frame in which the gyration (every cell's vector turning in place at the same rate, Section 3) becomes static: for the macrospin, $\mathbf m'=(\cos\theta_0,0,\sin\theta_0)$, a constant.

The alternative, rotating *real space* ($\phi\to\phi-\Phi$), is what one instinctively does for a rotating pattern, and it is the wrong tool here: the spatial pattern of the STO does not move, the vectors do. That wrong choice is what produces the apparent paradox discussed in Section 4.3.

## 2. Consequence 1: $m_z$ is invariant

By F1, $m_z(\mathbf r,t)=m'_z(\mathbf r,t)$ exactly, for any texture. Hence if in the rotating frame there is a wave

$$
m'_z(\mathbf r,t)\propto J_1(kr)\,e^{\mathrm i(m\phi-\Omega t)},
$$

then the laboratory $m_z$ is **the same function**: same azimuthal index $m$, same $k$, same frequency $\Omega$. There is no Doppler shift $\Omega\to\Omega+m\omega_f$ for $m_z$, because $\phi$ was never touched. The FFT maps in the manuscript are maps of $M_z$, so what they show *is* the rotating-frame wave.

## 3. Consequence 2: the in-plane components pick up $e^{\pm\mathrm i\Phi}$

Write the first two rows of $\mathbf m=\mathsf R\mathbf m'$:

$$
m_x=\cos\Phi\,m'_x-\sin\Phi\,m'_y,\qquad m_y=\sin\Phi\,m'_x+\cos\Phi\,m'_y .
$$

Form $m_+=m_x+\mathrm i m_y$:

$$
m_+=(\cos\Phi+\mathrm i\sin\Phi)\,m'_x+(-\sin\Phi+\mathrm i\cos\Phi)\,m'_y
=e^{\mathrm i\Phi}m'_x+\mathrm i\,e^{\mathrm i\Phi}m'_y
=e^{\mathrm i\Phi}\,m'_+ .
$$

The same rule holds for any vector, in particular for fields: $H_+=e^{\mathrm i\Phi}H'_+$, i.e. $H'_+=e^{-\mathrm i\omega_f t}H_+$. A field that is **static** in the laboratory therefore **rotates at $-\omega_f$** in the rotating frame. For the radial edge field of the fixed layer, $H_+=h(r)\,e^{\mathrm i\phi}$ in the lab, so

$$
H'_+=h(r)\,e^{\mathrm i(\phi-\omega_f t)} :
$$

an rf field at the gyration frequency with azimuthal index $m=1$. (This is the source term used in Section 6; in Section 4.1 it is only announced.)

## 4. Time derivative of the transformation

Differentiate $\mathbf m=\mathsf R\mathbf m'$ with the product rule, remembering that both factors depend on $t$:

$$
\dot{\mathbf m}=\dot{\mathsf R}\,\mathbf m'+\mathsf R\,\dot{\mathbf m}'
\overset{\text{F4}}{=}\omega_f\,\mathsf R\,(\hat{\mathbf z}\times\mathbf m')+\mathsf R\,\dot{\mathbf m}'
=\mathsf R\big(\dot{\mathbf m}'+\omega_f\,\hat{\mathbf z}\times\mathbf m'\big).
\tag{B}
$$

Sanity check against Section 3: for the macrospin $\mathbf m'$ is constant, so $\dot{\mathbf m}=\omega_f\mathsf R(\hat{\mathbf z}\times\mathbf m')=\omega_f\,\hat{\mathbf z}\times\mathbf m$ (F4 again), which is relation (A) of the Section 3 notes.

Equivalently, $\mathsf R^{\mathsf T}\dot{\mathbf m}=\dot{\mathbf m}'+\omega_f\hat{\mathbf z}\times\mathbf m'$: the rotating-frame observer sees the true velocity minus the velocity of the rigid rotation. This is the same bookkeeping as $\dot{\mathbf r}_{\rm lab}=\dot{\mathbf r}_{\rm rot}+\boldsymbol\omega\times\mathbf r$ in mechanics, applied to a vector in spin space.

## 5. Transform the LLG equation term by term

Start from the lab equation (Eq. 1 of the derivation), with the spin-torque term already in the form of Eq. (stt):

$$
\dot{\mathbf m}=-\gamma\,\mathbf m\times\mathbf H_{\rm eff}+\alpha\,\mathbf m\times\dot{\mathbf m}+\gamma a_J\big[\hat{\mathbf z}-m_z\mathbf m\big].
$$

Apply $\mathsf R^{\mathsf T}$ to both sides and use $\mathsf R^{\mathsf T}\mathbf m=\mathbf m'$ throughout.

**Left side**, from (B):

$$
\mathsf R^{\mathsf T}\dot{\mathbf m}=\dot{\mathbf m}'+\omega_f\,\hat{\mathbf z}\times\mathbf m'.
$$

**Precession term.** By F3,

$$
\mathsf R^{\mathsf T}\big(-\gamma\,\mathbf m\times\mathbf H_{\rm eff}\big)=-\gamma\,(\mathsf R^{\mathsf T}\mathbf m)\times(\mathsf R^{\mathsf T}\mathbf H_{\rm eff})=-\gamma\,\mathbf m'\times\mathbf H'_{\rm eff},
\qquad \mathbf H'_{\rm eff}\equiv\mathsf R^{\mathsf T}\mathbf H_{\rm eff}.
$$

This is just a definition so far. The useful statement is that $\mathbf H'_{\rm eff}$ is obtained from $\mathbf m'$ by **the same formulas** that give $\mathbf H_{\rm eff}$ from $\mathbf m$ ("covariance"). Check each field:

- Exchange: $\mathsf R$ does not depend on $\mathbf r$, so $\mathsf R^{\mathsf T}D\nabla^2(\mathsf R\mathbf m')=D\nabla^2\mathbf m'$. ✓
- Thin-film demagnetisation: $\mathsf R^{\mathsf T}(-M_s m_z\hat{\mathbf z})=-M_s m'_z\hat{\mathbf z}$ by F1. ✓
- Biquadratic: $\mathsf R^{\mathsf T}(-H_B m_z\hat{\mathbf z})=-H_B m'_z\hat{\mathbf z}$. ✓ (The fixed layer along $\hat{\mathbf z}$ is itself invariant.)
- Fine print: the finite-$k$ in-plane dipolar term $-M_s(1-P)(\delta\mathbf m\cdot\hat{\mathbf k})\hat{\mathbf k}$ of Eq. (dip) is **not** invariant, because $\hat{\mathbf k}$ is a real-space direction that the spin rotation does not turn: in the rotating frame it depends on the angle $\varphi_k$ between $\hat{\mathbf k}$ and the in-plane direction of $\mathbf m'_0$. That is why Section 5.2 averages $\sin^2\varphi_k,\cos^2\varphi_k\to\tfrac12$. The uniform-film $z$ part, which is the dominant one, is exactly invariant.

**Damping term.** By F3 and then (B):

$$
\mathsf R^{\mathsf T}\big(\alpha\,\mathbf m\times\dot{\mathbf m}\big)
=\alpha\,\mathbf m'\times(\mathsf R^{\mathsf T}\dot{\mathbf m})
=\alpha\,\mathbf m'\times\dot{\mathbf m}'+\alpha\omega_f\,\mathbf m'\times(\hat{\mathbf z}\times\mathbf m').
$$

The last piece is the triple product of the Section 3 notes:
$\mathbf m'\times(\hat{\mathbf z}\times\mathbf m')=\hat{\mathbf z}-m'_z\mathbf m'$. So

$$
\mathsf R^{\mathsf T}\big(\alpha\,\mathbf m\times\dot{\mathbf m}\big)=\alpha\,\mathbf m'\times\dot{\mathbf m}'+\alpha\omega_f\big[\hat{\mathbf z}-m'_z\mathbf m'\big].
$$

The first part is ordinary Gilbert damping of whatever moves in the rotating frame. The second part is the damping of the rigid rotation itself: it is there even when $\mathbf m'$ is perfectly static, because the lab velocity never vanishes.

**Spin-torque term.** By F1 and $m_z=m'_z$,

$$
\mathsf R^{\mathsf T}\,\gamma a_J\big[\hat{\mathbf z}-m_z\mathbf m\big]=\gamma a_J\big[\hat{\mathbf z}-m'_z\mathbf m'\big].
$$

$a_J$ itself is unchanged: with $\Lambda=1$ the efficiency does not depend on $\mathbf m\cdot\mathbf p$, and anyway $\mathbf m\cdot\mathbf p=m_z=m'_z$.

**Collect:**

$$
\dot{\mathbf m}'+\omega_f\,\hat{\mathbf z}\times\mathbf m'
=-\gamma\,\mathbf m'\times\mathbf H'_{\rm eff}+\alpha\,\mathbf m'\times\dot{\mathbf m}'
+\big(\alpha\omega_f+\gamma a_J\big)\big[\hat{\mathbf z}-m'_z\mathbf m'\big].
$$

## 6. Turn the leftover rotation term into a fictitious field

Move $\omega_f\hat{\mathbf z}\times\mathbf m'$ to the right-hand side and rewrite it as a precession term:

$$
-\omega_f\,\hat{\mathbf z}\times\mathbf m'=+\omega_f\,\mathbf m'\times\hat{\mathbf z}
=-\gamma\,\mathbf m'\times\Big(-\frac{\omega_f}{\gamma}\hat{\mathbf z}\Big).
$$

So it is exactly what an additional field $-(\omega_f/\gamma)\hat{\mathbf z}$ would do, and the equation becomes Eq. (llgrot) of the text:

$$
\boxed{\;
\dot{\mathbf m}'=-\gamma\,\mathbf m'\times\Big(\mathbf H'_{\rm eff}-\frac{\omega_f}{\gamma}\hat{\mathbf z}\Big)
+\alpha\,\mathbf m'\times\dot{\mathbf m}'
+\big(\alpha\omega_f+\gamma a_J\big)\big[\hat{\mathbf z}-m'_z\mathbf m'\big]\;}
$$

Three things to notice:

1. **Sign of the fictitious field.** In the convention $\dot{\mathbf m}=-\gamma\mathbf m\times\mathbf H$, a field $H\hat{\mathbf z}$ makes $\mathbf m$ turn counter-clockwise at angular velocity $+\gamma H$ about $+z$. Going to a frame that already turns at $\omega_f$ removes that much angular velocity, i.e. subtracts the field $\omega_f/\gamma$. Here $\omega_f<0$ (clockwise), so the fictitious field $-\omega_f/\gamma$ points along $+z$ and is $+3.72\times10^5$ A/m (0.468 T) with BQ, $+3.41\times10^5$ A/m (0.429 T) without. This is Larmor's theorem.

2. **It cancels the easy-plane field.** For the macrospin, Eq. (omegaf) says $\omega_f=-\gamma(M_s+H_B)m_{z0}$, so $-\omega_f/\gamma=+(M_s+H_B)m_{z0}$, which is minus the demag-plus-BQ field $-(M_s+H_B)m_{z0}\hat{\mathbf z}$. Total static field in the rotating frame: zero. That is Section 4.2.

3. **The last bracket vanishes identically.** With the polar balance $\alpha\omega_f+\gamma a_J=0$ of Section 3, the damping of the rotation and the spin torque cancel **for any $\mathbf m'$**, uniform or not, because both have the same vector structure $\hat{\mathbf z}-m'_z\mathbf m'$ pointwise. What survives is a plain conservative-plus-Gilbert LLG equation for a film with no internal field. That is what makes the spin-wave problem of Section 5 tractable.

## 7. The same thing in complex notation (quick cross-check)

Using $m_+=e^{\mathrm i\Phi}m'_+$ from part 3:

$$
\dot m_+=e^{\mathrm i\Phi}\big(\dot m'_++\mathrm i\omega_f m'_+\big).
$$

The precession produced by a $z$ field is $\dot m_+=\mathrm i\gamma H_z m_+$ (Section 3 notes, part 4). Dividing by $e^{\mathrm i\Phi}$:

$$
\dot m'_+=\mathrm i\big(\gamma H_z-\omega_f\big)m'_+=\mathrm i\gamma\Big(H_z-\frac{\omega_f}{\gamma}\Big)m'_+ ,
$$

the same shift $H_z\to H_z-\omega_f/\gamma$, and zero when $\omega_f=\gamma H_z$. One line, but it only covers the precession term; the vector derivation above is needed to see what happens to damping and spin torque.

## 8. Points that are easy to trip over

- **Spin space versus real space.** If you rotate real space instead, $\phi\to\phi-\omega_f t$, an $m=1$ wave at the gyration frequency becomes static in that frame, and one is tempted to call it a "frozen spiral". But a real-space rotation does not make the gyrating vectors static, so it does not remove the time dependence from the LLG equation; nothing is gained. The spin rotation makes the background static *and* leaves the $M_z$ map untouched. Both statements are needed and both hold only for the spin-space choice.

- **$\mathsf R$ is time dependent but position independent.** That is why it commutes with $\nabla^2$ (exchange is covariant) and why F4 produces a single extra term rather than a convective one.

- **Why the damping term is the only one that generates something new.** Precession and spin torque are algebraic in $\mathbf m$ and transform covariantly. The Gilbert term contains $\dot{\mathbf m}$, and the time derivative is the one operation that does not commute with a time-dependent rotation; the extra piece $\alpha\omega_f[\hat{\mathbf z}-m'_z\mathbf m']$ is the price, and it is the damping that the spin torque has to pay for.

- **Direction of $\mathbf m'_0$.** The text takes $\mathbf m'_0=(\cos\theta_0,0,\sin\theta_0)$, i.e. the in-plane part along $x'$. This is just the choice of the origin of $\Phi$; the gyration phase is the Goldstone mode and any in-plane direction is equivalent.

- **The polariser does not move.** $\mathbf p=\hat{\mathbf z}$ satisfies $\mathsf R^{\mathsf T}\mathbf p=\mathbf p$. Had the polariser had an in-plane part, it would have become an rf torque in the rotating frame, exactly like the in-plane stray field in part 3. The actual $0.6^\circ$ tilt of $\mathbf p$ is such a term; it is neglected.
