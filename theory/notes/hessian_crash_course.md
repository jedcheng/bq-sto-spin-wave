# Crash course: the Hessian, and why spin-wave frequencies are "γ times the square root of its determinant"

You have met the Hessian before, just not under that name. It is the matrix of second derivatives of a function of several variables. In physics it appears every time one expands an energy around an equilibrium and keeps the quadratic term: the spring constant, the second-derivative test from calculus, the normal-mode matrix of coupled oscillators, the Kittel formula. This note goes from the one-variable case to the stiffness matrix $\mathsf N$ of Section 5 in small steps.

## 1. One variable: curvature is the spring constant

Expand an energy $E(x)$ around a point $x_0$:

$$
E(x_0+\xi)=E(x_0)+E'(x_0)\,\xi+\tfrac12E''(x_0)\,\xi^2+O(\xi^3).
$$

At an equilibrium $E'(x_0)=0$, so the first thing that survives is the curvature $E''(x_0)\equiv\kappa$. Then:

- the restoring force is $F=-\partial E/\partial\xi=-\kappa\xi$ (Hooke's law with spring constant $\kappa$);
- the equilibrium is stable if $\kappa>0$, unstable if $\kappa<0$;
- a mass $m$ on it oscillates at $\omega^2=\kappa/m$.

Everything below is this paragraph with $x$ replaced by several variables.

## 2. Several variables: the Hessian

For $E(x_1,\dots,x_n)$ the same expansion around $\mathbf x_0$ reads

$$
E(\mathbf x_0+\boldsymbol\xi)=E(\mathbf x_0)+\sum_i\frac{\partial E}{\partial x_i}\xi_i
+\frac12\sum_{i,j}\frac{\partial^2E}{\partial x_i\partial x_j}\,\xi_i\xi_j+\dots
=E_0+\nabla E\cdot\boldsymbol\xi+\tfrac12\,\boldsymbol\xi^{\mathsf T}\mathsf H\,\boldsymbol\xi+\dots
$$

The **Hessian** is the matrix

$$
\mathsf H_{ij}=\frac{\partial^2E}{\partial x_i\partial x_j}\qquad(\text{symmetric: }\mathsf H_{ij}=\mathsf H_{ji}).
$$

It is the multi-dimensional curvature. Three facts you already use:

**(a) Linearised force.** The force is minus the gradient. Differentiating the quadratic form,

$$
F_i=-\frac{\partial E}{\partial\xi_i}=-\sum_j\mathsf H_{ij}\,\xi_j,\qquad\text{i.e.}\qquad \mathbf F=-\mathsf H\,\boldsymbol\xi .
$$

The Hessian is the matrix of spring constants, including the cross couplings. "Linearising the force" and "taking the Hessian of the energy" are the same operation.

**(b) Stability.** The equilibrium is a minimum iff $\boldsymbol\xi^{\mathsf T}\mathsf H\boldsymbol\xi>0$ for every $\boldsymbol\xi\ne0$ (positive definite), i.e. all eigenvalues positive. For $2\times2$ this is the second-derivative test from calculus:

$$
\det\mathsf H=E_{xx}E_{yy}-E_{xy}^2>0\quad\text{and}\quad E_{xx}>0 .
$$

$\det\mathsf H<0$ is a saddle: one direction curves down.

**(c) Principal axes.** Being symmetric, $\mathsf H$ can be diagonalised by a rotation of the coordinates. In the rotated coordinates the energy is a sum of independent parabolas $\tfrac12\kappa_1\eta_1^2+\tfrac12\kappa_2\eta_2^2+\dots$ with $\kappa_i$ the eigenvalues. The cross term $E_{xy}$ only says that your chosen axes are not the principal ones.

**Coupled masses and springs.** For $n$ masses with kinetic energy $\tfrac12\dot{\boldsymbol\xi}^{\mathsf T}\mathsf M\dot{\boldsymbol\xi}$, Newton's law is $\mathsf M\ddot{\boldsymbol\xi}=-\mathsf H\boldsymbol\xi$, and the normal-mode frequencies are the eigenvalues of $\mathsf M^{-1}\mathsf H$: $\det(\mathsf H-\omega^2\mathsf M)=0$. The "K matrix" of the normal-mode problem is the Hessian of the potential energy. Keep this in mind, because spins do **not** follow this rule, and the difference is the key to Smit–Beljers.

## 3. Energy as a function of the magnetisation direction

The magnetic energy density $E(\mathbf m)$ is a function on the unit sphere. The effective field is defined as minus its gradient,

$$
\mathbf H_{\rm eff}=-\frac{1}{\mu_0M_s}\frac{\partial E}{\partial\mathbf m},
$$

exactly like a force is minus the gradient of a potential, with $\mu_0M_s$ converting J/m$^3$ per unit of $\mathbf m$ into A/m. Examples, all of which you can check by differentiating:

| energy density | field |
|---|---|
| Zeeman $-\mu_0M_s\,\mathbf m\cdot\mathbf H_0$ | $\mathbf H_0$ |
| thin-film demag $\tfrac12\mu_0M_s^2m_z^2$ | $-M_sm_z\hat{\mathbf z}$ |
| biquadratic $-\frac{B_{12}}{d}(\mathbf m\cdot\hat{\mathbf z})^2$ with $B_{12}<0$ | $-H_Bm_z\hat{\mathbf z}$, $H_B=2|B_{12}|/(\mu_0M_sd)$ |
| exchange $A|\nabla\mathbf m|^2$ | $D\nabla^2\mathbf m$, $D=2A/(\mu_0M_s)$ |

Now expand around an equilibrium $\mathbf m_0$ with $\mathbf m=\mathbf m_0+\delta\mathbf m$. The constraint $|\mathbf m|=1$ means that $\delta\mathbf m$ has only two independent components, the ones in the tangent plane: $\delta\mathbf m=a\,\hat{\mathbf e}_1+b\,\hat{\mathbf e}_2$. The energy to second order is

$$
E=E_0+\tfrac12\mu_0M_s\big(N_{11}a^2+2N_{12}ab+N_{22}b^2\big),
\qquad
N_{ij}=\frac{1}{\mu_0M_s}\frac{\partial^2E}{\partial x_i\partial x_j},\quad(x_1,x_2)=(a,b).
$$

So **$\mathsf N$ is the Hessian of the energy density, divided by $\mu_0M_s$ so that it has units of a field**, restricted to the two tangent directions. And by fact (a), the linearised field is

$$
\delta\mathbf H=-\mathsf N\,\delta\mathbf m ,
$$

which is precisely the relation used in Section 5.3. Both routes give the same $\mathsf N$: you can either differentiate the energy twice, or linearise the field expressions directly as the note on Section 5 does. The second route is often shorter when the fields are already known; the first makes the symmetry $N_{12}=N_{21}$ obvious.

One subtlety hides in the word "restricted". When $\mathbf m_0$ is in a field, the energy also has a first-derivative term along $\mathbf m_0$ and the constraint feeds it into the tangent curvature: moving by $a$ along $\hat{\mathbf e}_1$ shortens the projection on $\mathbf m_0$ by $a^2/2$, so a Zeeman energy $-\mu_0M_sH_0m_\parallel$ contributes $+\tfrac12\mu_0M_sH_0a^2$. That is where the "$+H_0$ on the diagonal" of the usual Smit–Beljers matrix comes from. In the rotating frame of Section 4 the total static field is zero, so this term is absent, and $\mathsf N$ contains only the genuine curvatures of exchange, dipolar and biquadratic energy.

## 4. Why the frequency is $\gamma\sqrt{\det\mathsf N}$ and not an eigenvalue

Here is the one point where spins differ from masses on springs. The LLG precession term, linearised, is (Section 5.4)

$$
\dot a=-\gamma\,(N_{21}a+N_{22}b),\qquad \dot b=+\gamma\,(N_{11}a+N_{12}b),
\qquad\text{i.e.}\qquad
\frac{d}{dt}\begin{pmatrix}a\\ b\end{pmatrix}=\gamma\,\mathsf J\,\mathsf N\begin{pmatrix}a\\ b\end{pmatrix},
\quad\mathsf J=\begin{pmatrix}0&-1\\ 1&0\end{pmatrix}.
$$

These are **first-order** equations in which the gradient along one coordinate drives the velocity of the other. That is Hamilton's structure: if you call $\mathcal H=\tfrac12(N_{11}a^2+2N_{12}ab+N_{22}b^2)$, the equations are $\dot a=-\gamma\,\partial\mathcal H/\partial b$, $\dot b=+\gamma\,\partial\mathcal H/\partial a$. The two tangent coordinates are canonically conjugate to each other, like $x$ and $p$. (In spherical angles the conjugate pair is $\phi$ and $\cos\theta$; this is the spin version of "angle and angular momentum".)

Compare with a harmonic oscillator written in the same way: $\mathcal H=\tfrac12\kappa x^2+\tfrac12p^2/m$, $\dot x=\partial\mathcal H/\partial p=p/m$, $\dot p=-\partial\mathcal H/\partial x=-\kappa x$, giving

$$
\omega^2=\kappa\cdot\frac1m=\text{product of the two curvatures}.
$$

Not an eigenvalue of the Hessian, the **product** of its two diagonal entries (its determinant, when there is no cross term). The spin system is identical in structure, so

$$
\boxed{\;\Omega^2=\gamma^2\det\mathsf N=\gamma^2\big(N_{11}N_{22}-N_{12}^2\big)\;}
$$

Derivation in one line: $a,b\propto e^{-\mathrm i\Omega t}$ gives $-\mathrm i\Omega\,\mathbf v=\gamma\mathsf J\mathsf N\mathbf v$, so $-\Omega^2=\gamma^2\det(\mathsf J\mathsf N)=\gamma^2\det\mathsf J\det\mathsf N=-\gamma^2\det\mathsf N$ for a $2\times2$ system (the eigenvalues of $\mathsf J\mathsf N$ come as $\pm\mathrm i\Omega$, their product is $\det$).

Consequences that are easy to remember once you see the product structure:

- **Geometric mean.** $\Omega=\gamma\sqrt{N_{11}N_{22}}$ (ignoring $N_{12}$): the frequency is $\gamma$ times the geometric mean of the stiffness along the two tangent directions. A very soft direction and a very stiff direction give an intermediate frequency and a very **elliptical** orbit, with the large amplitude along the soft axis: $|a/b|=\sqrt{N_{22}/N_{11}}$.
- **Stability test is the same as for any Hessian.** $\det\mathsf N>0$ with positive diagonal means a minimum and a real $\Omega$; $\det\mathsf N<0$ (saddle) gives $\Omega^2<0$, i.e. exponential growth away from the equilibrium. Frequency going to zero is the signature of an instability threshold ("mode softening").
- **Goldstone mode.** If the energy does not depend on one coordinate at all, that curvature is zero, $\det\mathsf N=0$ and $\Omega=0$ at $k=0$. For the easy-plane film this is the free in-plane angle, and it is why the dispersion of Section 5 starts at zero.

## 5. Three worked examples

**(i) Spin in a field.** $E=-\mu_0M_sH_0m_z$, equilibrium $\mathbf m_0=\hat{\mathbf z}$, tangent coordinates $a=m_x$, $b=m_y$, $m_z\simeq1-\tfrac12(a^2+b^2)$. Then $E\simeq E_0+\tfrac12\mu_0M_sH_0(a^2+b^2)$, so $N_{11}=N_{22}=H_0$, $N_{12}=0$, and $\Omega=\gamma\sqrt{H_0\cdot H_0}=\gamma H_0$: Larmor precession, circular orbit.

**(ii) Kittel mode of an in-plane film.** Field $H_0$ along $\hat{\mathbf x}$ in the plane, $E=-\mu_0M_sH_0m_x+\tfrac12\mu_0M_s^2m_z^2$. Take $a=m_y$ (in-plane), $b=m_z$ (out of plane). The Zeeman part gives $\tfrac12\mu_0M_sH_0(a^2+b^2)$, the demag part $\tfrac12\mu_0M_s^2b^2$:

$$
N_{11}=H_0,\qquad N_{22}=H_0+M_s\qquad\Longrightarrow\qquad \Omega=\gamma\sqrt{H_0(H_0+M_s)} ,
$$

the Kittel formula, with an orbit elongated in the plane by $\sqrt{(H_0+M_s)/H_0}$. You have used this result; this is where it comes from.

**(iii) The gyrating state of Section 5.** Equilibrium $\mathbf m_0=(\cos\theta_0,0,\sin\theta_0)$ in zero total field; tangent coordinates $a$ (in-plane twist) and $b$ (polar tilt). The Hessian of exchange is $Dk^2$ on both diagonals (the Fourier transform of $-D\nabla^2$); that of the easy-plane energies $\tfrac12\mu_0M_s(M_sP+H_B)m_z^2$ is $(M_sP+H_B)\cos^2\theta_0$ in the $b$-$b$ slot, because $\delta m_z=b\cos\theta_0$; that of the in-plane dipolar energy $\tfrac12\mu_0M_s^2(1-P)(\delta\mathbf m\cdot\hat{\mathbf k})^2$ gives the $\varphi_k$-dependent entries. No Zeeman-type term appears because $\mathbf H_0=0$. The result is the matrix of Section 5.3 and $\Omega=\gamma\sqrt{\det\mathsf N}$ is Eq. (disp). At the measured wavelength $N_{22}/N_{11}\approx14$: a soft twist direction and a stiff tilt direction, hence the strongly elliptical wave whose small $m_z$ part is what the FFT maps display.

## 6. Smit–Beljers in the form you may have seen in a textbook

In spherical angles $(\theta,\phi)$ of $\mathbf m$, with the energy density $E(\theta,\phi)$, the ferromagnetic-resonance formula reads

$$
\omega=\frac{\gamma}{\mu_0M_s\sin\theta_0}\sqrt{E_{\theta\theta}E_{\phi\phi}-E_{\theta\phi}^2}.
$$

This is the same statement. The orthonormal tangent displacements are $b=\delta\theta$ and $a=\sin\theta_0\,\delta\phi$ (arc length along a parallel), so $E_{\phi\phi}=\mu_0M_s\sin^2\theta_0\,N_{11}$, $E_{\theta\theta}=\mu_0M_sN_{22}$, $E_{\theta\phi}=\mu_0M_s\sin\theta_0N_{12}$, and the square root becomes $\mu_0M_s\sin\theta_0\sqrt{\det\mathsf N}$. The $1/\sin\theta_0$ in front is nothing but this change of variables. The word "Hessian" is rarely used in that literature; "second derivatives of the free energy" or "stiffness" is what you will find, but it is the same object.

## 7. From a macrospin to spin waves

For a continuous texture the energy is a functional, and the Hessian becomes an operator acting on $\delta\mathbf m(\mathbf r)$: $-D\nabla^2$ for exchange, a non-local integral kernel for the dipolar energy. Taking plane waves $e^{\mathrm i\mathbf k\cdot\mathbf r}$ diagonalises both, and the operator reduces to the $2\times2$ matrix $\mathsf N(\mathbf k)$ of Section 5. The dipolar block is the thin-film demagnetising tensor at wavevector $\mathbf k$, which is why the same letter $\mathsf N$ is traditional for both. Everything in Sections 4 to 6 then follows from the macrospin rules: frequency $=\gamma\sqrt{\det\mathsf N(k)}$, stability $=\det\mathsf N>0$ for all $k$, Goldstone mode at $k=0$, ellipticity $\sqrt{N_{22}/N_{11}}$.

## 8. Summary card

| concept | mechanics | magnetism |
|---|---|---|
| energy | $V(\mathbf x)$ | $E(\mathbf m)$ on the unit sphere |
| force / field | $\mathbf F=-\nabla V$ | $\mathbf H_{\rm eff}=-\dfrac{1}{\mu_0M_s}\dfrac{\partial E}{\partial\mathbf m}$ |
| Hessian | $\mathsf K_{ij}=\partial^2V/\partial x_i\partial x_j$ | $N_{ij}=\dfrac{1}{\mu_0M_s}\dfrac{\partial^2E}{\partial x_i\partial x_j}$, tangent coords |
| linearised force / field | $\mathbf F=-\mathsf K\boldsymbol\xi$ | $\delta\mathbf H=-\mathsf N\,\delta\mathbf m$ |
| equations of motion | second order, $\mathsf M\ddot{\boldsymbol\xi}=-\mathsf K\boldsymbol\xi$ | first order, $\dot{\mathbf v}=\gamma\mathsf J\mathsf N\mathbf v$ (canonical pair) |
| frequencies | eigenvalues of $\mathsf M^{-1}\mathsf K$ | $\Omega=\gamma\sqrt{\det\mathsf N}$ (geometric mean) |
| stable equilibrium | $\mathsf K$ positive definite | $\mathsf N$ positive definite |
| zero mode | free translation | free in-plane rotation (Goldstone) |
