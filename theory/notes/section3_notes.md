# Section 3 of `derivation.tex`, step by step: the gyrating (macrospin) state

Only one vector identity is needed, used three times:

$$
\mathbf m\times(\mathbf m\times\mathbf v)=\mathbf m\,(\mathbf m\cdot\mathbf v)-\mathbf v\,(\mathbf m\cdot\mathbf m)=\mathbf m\,(\mathbf m\cdot\mathbf v)-\mathbf v
\qquad(|\mathbf m|=1).
$$

## 1. The ansatz and its time derivative

Assume the free layer rotates rigidly in-plane at constant angular velocity $\omega_f$ with a constant tilt $\theta_0$:

$$
\mathbf m(t)=\big(\cos\theta_0\cos\Phi,\ \cos\theta_0\sin\Phi,\ \sin\theta_0\big),\qquad \Phi=\omega_f t,\qquad m_z=\sin\theta_0 .
$$

Only $\Phi$ depends on time, so

$$
\dot{\mathbf m}=\omega_f\cos\theta_0\,(-\sin\Phi,\ \cos\Phi,\ 0).
$$

Compare with $\hat{\mathbf z}\times\mathbf m=(-m_y,\ m_x,\ 0)=\cos\theta_0(-\sin\Phi,\cos\Phi,0)$:

$$
\boxed{\ \dot{\mathbf m}=\omega_f\,\hat{\mathbf z}\times\mathbf m\ }\tag{A}
$$

This is exact for the ansatz. Everything else is substituting (A) into the LLG equation.

## 2. The three terms of the LLG equation

The simulator integrates

$$
\dot{\mathbf m}=-\gamma\,\mathbf m\times\mathbf H+\alpha\,\mathbf m\times\dot{\mathbf m},
$$

with $\mathbf H$ the sum of all effective fields. Three of them matter for the rigid rotation.

**Term 1 – conservative field.** For a thin film the demagnetising field is $-M_s m_z\hat{\mathbf z}$ and the biquadratic field is $-H_B m_z\hat{\mathbf z}$ (Eq. 6 of the derivation). Both are along $z$, so $\mathbf H_c=H_z\hat{\mathbf z}$ with

$$
H_z=-(M_s+H_B)\sin\theta_0 .
$$

The in-plane anisotropy (3.5 mT) is dropped. Then

$$
-\gamma\,\mathbf m\times(H_z\hat{\mathbf z})=-\gamma H_z\,(m_y,-m_x,0)=\gamma H_z\,(-m_y,m_x,0)=\gamma H_z\ \hat{\mathbf z}\times\mathbf m .
$$

The conservative field produces a torque along $\hat{\mathbf z}\times\mathbf m$, the same direction as $\dot{\mathbf m}$ in (A).

**Term 2 – spin torque.** The Slonczewski field is $\mathbf H_{\rm stt}=a_J\,\mathbf m\times\mathbf p$ with $\mathbf p=\hat{\mathbf z}$. Insert and use the identity:

$$
-\gamma\,\mathbf m\times(a_J\,\mathbf m\times\hat{\mathbf z})=-\gamma a_J\big[\mathbf m\,(\mathbf m\cdot\hat{\mathbf z})-\hat{\mathbf z}\big]=\gamma a_J\big[\hat{\mathbf z}-m_z\mathbf m\big].
$$

**Term 3 – Gilbert damping.** Use (A) for $\dot{\mathbf m}$ inside the damping term, then the identity:

$$
\alpha\,\mathbf m\times\dot{\mathbf m}=\alpha\omega_f\,\mathbf m\times(\hat{\mathbf z}\times\mathbf m)=\alpha\omega_f\big[\hat{\mathbf z}\,(\mathbf m\cdot\mathbf m)-\mathbf m\,(\mathbf m\cdot\hat{\mathbf z})\big]=\alpha\omega_f\big[\hat{\mathbf z}-m_z\mathbf m\big].
$$

Terms 2 and 3 have exactly the same vector structure, $\hat{\mathbf z}-m_z\mathbf m$. That is the whole trick of the section.

## 3. Collect and match components

Left side from (A), right side from the three terms:

$$
\omega_f\,(\hat{\mathbf z}\times\mathbf m)=\gamma H_z\,(\hat{\mathbf z}\times\mathbf m)+(\gamma a_J+\alpha\omega_f)\big[\hat{\mathbf z}-m_z\mathbf m\big].
$$

The two vectors on the right are perpendicular to each other:

$$
(\hat{\mathbf z}\times\mathbf m)\cdot(\hat{\mathbf z}-m_z\mathbf m)=(\hat{\mathbf z}\times\mathbf m)\cdot\hat{\mathbf z}-m_z(\hat{\mathbf z}\times\mathbf m)\cdot\mathbf m=0-0=0 .
$$

Geometrically, $\hat{\mathbf z}\times\mathbf m$ is the azimuthal direction (along the rotation) and $\hat{\mathbf z}-m_z\mathbf m$ is the polar direction (tilting $\mathbf m$ toward or away from $z$). An equation between two orthogonal directions can only hold if each coefficient matches separately.

Azimuthal component:

$$
\boxed{\ \omega_f=\gamma H_z=-\gamma(M_s+H_B)\sin\theta_0\ }\qquad\text{(Eq. 9)}
$$

Polar component:

$$
\gamma a_J+\alpha\omega_f=0\quad\Longrightarrow\quad\boxed{\ \omega_f=-\frac{\gamma a_J}{\alpha}\ }\qquad\text{(Eq. 11, first part)}
$$

Combine the two:

$$
\boxed{\ \sin\theta_0=\frac{a_J}{\alpha\,(M_s+H_B)}\ }\qquad\text{(Eq. 11, second part)}
$$

Physical reading: the polar balance says spin-torque input equals damping and fixes the rotation rate; the azimuthal balance says the tilt must be whatever makes the easy-plane fields precess the moment at that rate.

## 4. Points that are easy to trip over

- **The complex shortcut in the text.** With $m_+=m_x+\mathrm i\,m_y$, the vector $\hat{\mathbf z}\times\mathbf m=(-m_y,m_x)$ is just $\mathrm i\,m_+$. The azimuthal equation reads $\dot m_+=\mathrm i\gamma H_z\,m_+$, whose solution is $m_+\propto e^{\mathrm i\gamma H_z t}$, so $\dot\Phi=\gamma H_z$. It contains only the conservative field because in steady state the damping and spin-torque terms cancel and never touch the azimuthal direction.

- **Sign and sense of rotation.** With a positive tilt, $H_z<0$, so $\omega_f<0$: clockwise seen from $+z$. The polar equation says the same from the other side: a positive $a_J$ gives negative $\omega_f$, and then Eq. 9 gives a positive $\sin\theta_0$. The two equations are consistent with each other and with the simulation ($-13.1$ GHz, $m_z\approx+0.26$).

- **Why the biquadratic field raises the frequency.** $H_B$ appears only in the azimuthal equation. For a given tilt the precession is faster by $(M_s+H_B)/M_s=1.26$; for a given frequency the tilt is smaller.

- **Stray field of the fixed layer.** If the fixed layer produces a $z$ field $h_z$ at the free layer, it simply adds to $H_z$: $\omega_f=\gamma[-(M_s+H_B)\sin\theta_0+h_z]$.

## 5. Numbers

| quantity | with BQ (13.1 GHz) | without BQ (12.0 GHz) |
|---|---|---|
| $\omega_f/\gamma=2\pi f/\gamma$ | $3.72\times10^5$ A/m = 0.468 T | $3.41\times10^5$ A/m = 0.429 T |
| $\mu_0(M_s+H_B)$ | 1.696 + 0.444 = 2.14 T | 1.696 T |
| $\sin\theta_0$ from Eq. 9 | 0.22 | 0.25 |
| $\langle m_z\rangle$ in the simulation | 0.26 | 0.28 |

Constants: $\gamma=2.211\times10^5\ \mathrm{m\,A^{-1}s^{-1}}$ (so $\omega=\gamma H$ with $H$ in A/m, and $\gamma/2\pi=28.0$ GHz/T), $\mu_0 M_s=1.696$ T, $\mu_0 H_B=2|B_{12}|/(M_s d)=0.444$ T.

The frequency from Eq. 11 alone, with $a_J=\hbar J P/(2e\mu_0 M_s d)=2.04\times10^4$ A/m at $J=6\times10^{11}$ A/m$^2$, is $f=\gamma a_J/(2\pi\alpha)=20.5$ GHz. The disk runs at 13.1 GHz because the standing wave takes part of the spin-torque work; that is the content of the power-balance section (Sec. 9.2). Section 3 only fixes the relation between tilt and frequency, which is what the rest of the derivation uses.
