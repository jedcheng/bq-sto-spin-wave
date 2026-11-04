# Sections 6 and 7 of `derivation.tex`, step by step: the source, and the standing wave in the disk

What is in hand after Section 5: in the co-rotating spin frame the free layer is a static film with zero internal field, and small ripples on it are spin waves with the dispersion $\Omega(k)$ of Eq. (disp). Section 6 asks what *drives* such a ripple and at what frequency; Section 7 asks what shape the ripple takes in a finite disk and why that produces a staircase of frequencies.

---

# Part A: Section 6, the source

## A1. The stray field, in the laboratory

The fixed layer is magnetised along $\hat{\mathbf z}$ and is a disk of the same radius, so its stray field at the free layer is static and axisymmetric:

$$
\mathbf h(r)=h_r(r)\,\hat{\mathbf r}+h_z(r)\,\hat{\mathbf z},\qquad
\hat{\mathbf r}=(\cos\phi,\sin\phi,0).
$$

$h_r$ is concentrated within about 15 nm of the rim (it is the fringe field of the fixed layer's edge charges); $h_z$ is small and positive inside the disk. In Cartesian components the in-plane part is $h_r(r)\,(\cos\phi,\ \sin\phi,\ 0)$. Note that this is the field of a *rigid* magnet: it does not depend on what the free layer does.

## A2. The same field, in the rotating spin frame

Apply the rule of Section 4.1: a field transforms as $\mathbf h'=\mathsf R_z^{\mathsf T}(\Phi)\,\mathbf h$ with $\Phi=\omega_f t$. The transpose rotates by $-\Phi$:

$$
\mathsf R_z^{\mathsf T}(\Phi)\begin{pmatrix}\cos\phi\\ \sin\phi\\ 0\end{pmatrix}
=\begin{pmatrix}\cos\Phi&\sin\Phi&0\\ -\sin\Phi&\cos\Phi&0\\ 0&0&1\end{pmatrix}\begin{pmatrix}\cos\phi\\ \sin\phi\\ 0\end{pmatrix}
=\begin{pmatrix}\cos\phi\cos\Phi+\sin\phi\sin\Phi\\ \sin\phi\cos\Phi-\cos\phi\sin\Phi\\ 0\end{pmatrix}
=\begin{pmatrix}\cos(\phi-\Phi)\\ \sin(\phi-\Phi)\\ 0\end{pmatrix}.
$$

So

$$
\mathbf h'(r,\phi,t)=h_r(r)\big(\cos(\phi-\omega_f t),\ \sin(\phi-\omega_f t),\ 0\big)+h_z(r)\,\hat{\mathbf z}.
$$

Read it: at a fixed position the in-plane field direction has angle $\phi-\omega_f t$, which changes in time. The static lab field has become a **rotating** field in the frame. Its $z$ part is untouched (F1 of the Section 4.1 note).

## A3. Project on the tangent directions of the equilibrium

Only the components perpendicular to $\mathbf m'_0=(\cos\theta_0,0,\sin\theta_0)$ exert a torque. With $\hat{\mathbf e}_1=\hat{\mathbf y}$ and $\hat{\mathbf e}_2=(-\sin\theta_0,0,\cos\theta_0)$:

$$
h'_1=\mathbf h'\cdot\hat{\mathbf e}_1=h_r\,\sin(\phi-\omega_f t),
\qquad
h'_2=\mathbf h'\cdot\hat{\mathbf e}_2=-\sin\theta_0\,h_r\cos(\phi-\omega_f t)+\cos\theta_0\,h_z .
$$

The $h_z\cos\theta_0$ piece is static and axisymmetric: it only shifts the equilibrium tilt (the $h_z$ correction mentioned at the end of Section 3) and is dropped from the drive. What remains is Eq. (drive) of the text:

$$
\boxed{\;h'_1=h_r(r)\sin(\phi-\omega_f t),\qquad h'_2=-\sin\theta_0\,h_r(r)\cos(\phi-\omega_f t)\;}
$$

Both components depend on $\phi$ and $t$ only through the combination $\phi-\omega_f t$. In complex form, $h'_1=\mathrm{Im}\,[h_re^{\mathrm i(\phi-\omega_ft)}]$ and $h'_2=-\sin\theta_0\,\mathrm{Re}\,[h_re^{\mathrm i(\phi-\omega_ft)}]$: a single azimuthal harmonic $e^{\mathrm i\phi}$ ($m=1$) at the single frequency $\omega_f$.

Three readings of the same formula:

- **Frequency.** The drive oscillates at the gyration frequency $|\omega_f|$, nothing else. There is no second frequency in the problem.
- **Azimuthal index.** Lines of constant phase satisfy $\phi-\omega_ft=\text{const}$, i.e. $\phi=\omega_ft+\text{const}$: the pattern turns at $\omega_f$, the same rate and the same sense as the magnetisation. This is what "co-rotating, $m=+1$" means, and it is what the raw-data check in Section 4 of the README confirms (co-rotating power 100 times the counter-rotating).
- **Physical picture.** The torque of an in-plane field on the magnetisation, $-\gamma\mathbf m\times\mathbf h$, is largest where $\mathbf h\perp\mathbf m$. On the rim, $\mathbf h$ points radially everywhere, so the place where it is perpendicular to the rotating $\mathbf m$ is a point that runs around the rim once per gyration period. Seen from the film, that is a wave source going round at $\omega_f$.

The $\sin\theta_0$ in $h'_2$ says that a *tilted* magnetisation feels the radial field also in the polar direction; for an untilted film ($\theta_0=0$) the radial field would only twist the in-plane angle.

## A4. Linear response: why the ripple has the same $(\phi,t)$ dependence

After the azimuthal averaging of Section 5.2 the medium is linear, time-independent and axisymmetric. For such a medium a drive $\propto e^{\mathrm i(m\phi-\Omega t)}$ produces a response $\propto e^{\mathrm i(m\phi-\Omega t)}$ with the same $m$ and $\Omega$ (there is nothing in the medium that could change either). Hence

$$
\delta\mathbf m(r,\phi,t)=\mathrm{Re}\big[\mathbf u(r)\,e^{\mathrm i(\phi-\omega_ft)}\big],\qquad
m_z-m_{z0}=\cos\theta_0\,\mathrm{Re}\big[u_2(r)\,e^{\mathrm i(\phi-\omega_ft)}\big],
$$

where $\mathbf u=(u_1,u_2)$ are the complex tangent amplitudes and $\delta m_z=b\cos\theta_0$ from the Section 5 note. Only the radial profile $\mathbf u(r)$ remains to be found.

## A5. The radial equation and the selection rule

Away from the rim, where $h_r=0$, $\mathbf u(r)$ obeys the *free* linearised equation of Section 5, but now in polar coordinates instead of for a plane wave. The only place the plane-wave assumption entered was $\nabla^2\to-k^2$ in the exchange term (and $k$ in $P(kd)$). For a function $u(r)e^{\mathrm i\phi}$,

$$
\nabla^2\big[u(r)e^{\mathrm i\phi}\big]=\Big(u''+\frac{u'}{r}-\frac{u}{r^2}\Big)e^{\mathrm i\phi},
$$

the $m=1$ Laplacian ($-m^2/r^2$ with $m=1$). Replace $-k^2$ by this operator in the Section 5 equations and demand a solution at the frequency $\Omega=|\omega_f|\equiv\omega$. For a plane wave the equations were solvable only on the dispersion curve $\Omega=\Omega(k)$; in the disk, the same algebra says $u(r)$ must satisfy the Helmholtz equation

$$
u''+\frac{u'}{r}+\Big(k^{*2}-\frac{1}{r^2}\Big)u=0
\qquad\text{with}\qquad
\boxed{\;\Omega(k^*)=\omega\;}\quad\text{(selection rule, Eq. selection)}
$$

In words: the drive fixes the frequency, the dispersion converts that frequency into a wavenumber, and the geometry (next part) converts that wavenumber into a pattern. There is no free parameter. The direction of the logic is the opposite of a normal cavity problem, where the geometry fixes $k$ and the dispersion then gives the frequency; here the oscillator supplies the frequency and the wavelength adjusts. (Section 7.2 closes the loop: the cavity in turn constrains which frequencies the oscillator can settle on.)

The nonlocal dipolar term, $P(kd)$, is evaluated at $k^*$; since $kd\approx0.13$ and the profile is a Bessel function with that $k$, this is a consistent approximation.

## A6. The two checks against raw data

They are quoted in the text; the point of each:

- **Handedness.** Eq. (drive) predicts a specific sign of $m$ relative to the sense of gyration. A wrong frame choice (real-space rotation, or $\mathsf R$ instead of $\mathsf R^{\mathsf T}$) would give the opposite sign. The raw $M_z(r,\phi,t)$ decomposition gives the co-rotating $m=+1$ at 100 times the counter-rotating power. The counter-rotating partner would need a two-fold anisotropy or an $m=2$ edge texture to be excited and is absent.
- **Uniform gyration rate.** The transformation of Section 4 uses one $\omega_f$ for the whole disk. The measured local gyration frequency is uniform across the disk to 3 %, so the disk is phase-locked and a single rotating frame is legitimate.

---

# Part B: Section 7, the standing wave in the disk

## B1. Radial mode shape: Bessel functions

(For a refresher on Bessel functions themselves, with a figure, see `bessel_crash_course.md`.)

The radial equation of A5 is Bessel's equation of order 1 in the variable $x=k^*r$:

$$
\frac{d^2u}{dx^2}+\frac1x\frac{du}{dx}+\Big(1-\frac1{x^2}\Big)u=0 .
$$

Its two independent solutions are $J_1(x)$ (regular at $x=0$, $J_1\simeq x/2$) and $Y_1(x)$ (singular, $Y_1\simeq-2/(\pi x)$). The centre of the disk is an ordinary point of the film, so $u$ must be finite there and

$$
u_2(r)\propto J_1(k^*r).
$$

That is the whole content of "regularity at $r=0$ selects the Bessel function". The $m_z$ pattern is Eq. (pattern):

$$
m_z(r,\phi,t)-m_{z0}\simeq A\,J_1(k^*r)\cos(\phi-\omega_ft+\delta).
$$

**Why it is a standing wave.** $J_1$ is real. Decompose it into Hankel functions, $J_1=\tfrac12\big(H_1^{(1)}+H_1^{(2)}\big)$, where $H^{(1)}_1(kr)e^{-\mathrm i\omega t}\sim e^{\mathrm i(kr-\omega t)}$ travels outward and $H^{(2)}_1$ inward, with equal amplitudes. The edge launches the inward wave; it converges on the centre, passes through, and comes out as the outward wave. Regularity at $r=0$ *is* total reflection at the centre. The sum of equal inward and outward waves is a standing wave.

**What the FFT maps show.** The FFT at frequency $\omega$ of $\cos(\phi-\omega_ft+\delta)J_1(k^*r)$ gives

- amplitude $\propto|J_1(k^*r)|$: maxima at the extrema of $J_1$, $k^*r=j'_{1,n}=1.84,\,5.33,\,8.54,\dots$, minima (nodes) at its zeros, $k^*r=j_{1,n}=3.83,\,7.02,\,10.17,\dots$. With $k^*=0.070$ nm$^{-1}$: antinodes at 26, 76, 122 nm and nodes at 55, 100, 145 nm, which are the ring positions of Table 1. Consecutive zeros of $J_1$ are $\pi$ apart (asymptotically $J_1(x)\simeq\sqrt{2/\pi x}\cos(x-3\pi/4)$), so **the rings are spaced by $\lambda/2$, not $\lambda$**.
- phase: $J_1$ changes sign at every zero, which in the complex FFT is a jump of $\pi$. Three interior nodes give $3\pi=1.5$ turns of total phase, exactly the measured 1.5 turns. The manuscript's "4 rotations" counted amplitude rings, not phase.

**Why the staircase is rounded into a spiral.** With damping, the inward wave loses amplitude on its way to the centre, so the outward wave that re-emerges is weaker than the inward one. Equal amplitudes give a pure standing wave; unequal ones give a standing wave plus a residual *inward travelling* wave, whose phase advances smoothly with $r$. Superposed on the $\pi$ steps, this smooth drift turns the staircase into the "spiral" look of the phase maps. The size of the drift measures the decay length (Section 9.5).

## B2. Cavity condition and plateau frequencies

$J_1(k^*r)$ with $k^*$ from the selection rule is a solution for *any* $\omega$. What makes some frequencies special is the rim. Whatever the rim does to the wave (the outer 15 nm is a driven, textured region), it imposes a condition on $u$ at $r=R$, and a condition at $r=R$ on a function $J_1(kR)$ restricts $kR$ to discrete values $x_n$. The two textbook limits bracket the real situation:

$$
\begin{aligned}
&\text{free spins at the rim (Neumann, }\partial_ru=0\text{):}\quad && J_1'(k_nR)=0,\quad k_nR=j'_{1,n}=1.84,\ 5.33,\ 8.54,\ 11.71,\ 14.86,\dots\\
&\text{pinned spins at the rim (Dirichlet, }u=0\text{):}\quad && J_1(k_nR)=0,\quad k_nR=j_{1,n}=3.83,\ 7.02,\ 10.17,\ 13.32,\ 16.47,\dots
\end{aligned}
$$

(Neumann is what exchange alone gives at a free surface, $\partial_r\mathbf m=0$; Dirichlet is what a hard edge texture gives. $j'_{1,n}<j_{1,n}<j'_{1,n+1}$, so each standing wave with $n$ interior antinodes lives in the window $j'_{1,n}\lesssim k_nR\lesssim j_{1,n}$.)

Now combine with the selection rule. If the cavity phase is pinned at $x_n$, then $k^*=x_n/R$, and the gyration frequency must be the one the dispersion assigns to that $k$:

$$
\boxed{\;f_n=\frac{1}{2\pi}\,\Omega\!\Big(\frac{x_n}{R}\Big),\qquad j'_{1,n}\lesssim x_n\lesssim j_{1,n}\;}
\qquad\text{(Eq. ladder)}
$$

This is a ladder of allowed frequencies for a given radius, one rung per number of interior antinodes. It is the quantitative version of "the oscillator can only run at frequencies whose wavelength fits the disk".

**Self-consistency.** $\Omega(k)$ contains $\cos^2\theta_0$ with $\sin\theta_0=\omega/[\gamma(M_s+H_B)]$, so the right-hand side depends weakly on the answer. Iterate: start with $\theta_0=0$, compute $\Omega$, update $\theta_0$, repeat. For $x_3=j_{1,3}$ with BQ the sequence is 14.66, 14.24, 14.27, 14.27 GHz: converged in three steps because $\cos^2\theta_0\approx0.95$. `sw_theory.py` does this with a root finder.

**Reading Table 2 honestly.** The simulated end-of-plateau frequencies sit near the $j'_{1,n}$ row. But the cavity phase measured directly from the phase slope of the maps, $kR=6.8,\ 10.7,\ 14.5,\ 17.1,\ 21.4$ for $n=2\dots6$, lies 0 to 1 rad *above* the zeros $j_{1,n}=7.02,\ 10.17,\ 13.32,\ 16.47,\ 19.62$. The fully developed standing wave is therefore close to the pinned (rim-node) condition. The reason the simulated frequencies nevertheless look like the $j'$ row is that the dispersion overestimates $\Omega$ at given $k$ by about 7 % (Section 8), i.e. about 1 GHz, which is roughly the gap between the two rows. The bright ring at the rim in the amplitude maps is the driven edge texture, not an antinode of the free wave.

## B3. Plateaus, jumps and the step height

**The mechanism in words.** Inside a plateau the current rises, the spin torque tries to raise $\omega$, and $k^*$ follows the dispersion upward, so the cavity phase $k^*R$ creeps from near $j'_{1,n}$ toward $j_{1,n}$. As it does so, the standing wave gains amplitude (Section 9 measures this), draws more power from the spin torque, and holds the frequency down: the plateau is flat. When the phase cannot advance further, the pattern re-locks with one more interior antinode, $n\to n+1$: a new ring appears at the centre and the frequency jumps to the next rung.

**The step height.** Two consecutive rungs of the ladder differ by one half wavelength of radial phase, because consecutive Bessel zeros (of either kind) are spaced by $\approx\pi$ (3.18, 3.16, 3.15, … for $j_{1,n}$). Expand $\Omega$ to first order in $k$:

$$
\Delta f_n=f_{n+1}-f_n=\frac{1}{2\pi}\Big[\Omega\Big(\frac{x_{n+1}}R\Big)-\Omega\Big(\frac{x_n}R\Big)\Big]
\simeq\frac{1}{2\pi}\,\frac{\partial\Omega}{\partial k}\,\frac{x_{n+1}-x_n}{R}
\simeq\frac{1}{2\pi}\,v_g\,\frac{\pi}{R}=\frac{v_g}{2R},
$$

with $v_g=\partial\Omega/\partial k$ the group velocity. In the linear-dispersion limit $\Omega\simeq ck$ this is Eq. (step), $\Delta f\simeq c/(2R)$: the free spectral range of a cavity of length $R$ for a wave of speed $c$, exactly as for a string of length $R$ with one fixed end.

| | with BQ | without BQ |
|---|---|---|
| $v_g$ at $k^*=0.067$–$0.070$ nm$^{-1}$ | 1140 m/s | 1015 m/s |
| $c=\gamma\sqrt{D(M_s+H_B)}$ | 1130 m/s | 1006 m/s |
| $v_g/2R$ ($R=160$ nm) | 3.56 GHz | 3.17 GHz |
| rung spacing of Table 2 | 3.3–3.6 GHz | 2.9–3.1 GHz |
| simulated step | 3.7 GHz | 3.4 GHz |

The biquadratic coupling raises the step because it stiffens the tilt bracket of the dispersion and so raises the wave speed by $\approx\sqrt{(M_s+H_B)/M_s}=1.12$.

**The same rule for the diameter sweep.** At fixed $J$ the oscillator sits near a fixed frequency, so $k^*$ is fixed; growing the radius makes the phase $k^*R$ grow, and a new rung is reached every time it advances by $\pi$: $\Delta R=\pi/k^*=\lambda/2$, i.e. $\Delta D=\lambda\approx80$–$95$ nm, and within a plateau $f=\Omega(x_n/R)/2\pi$ *falls* as $R$ grows. That is the sawtooth of Section 9.1.

## B4. Points that are easy to trip over

- **Which is the cause and which the effect.** Section 6: frequency (from the gyration) determines $k$. Section 7.2: $k$ (from the cavity) determines the frequency. Both are true at once; together they say the oscillator and the cavity have to agree, and the allowed agreements are the rungs $f_n$. Where the system sits on a rung is decided by the power balance, Section 9.

- **$\lambda/2$ ring spacing.** Amplitude rings of a standing wave are half a wavelength apart because $|J_1|$ has twice the periodicity of $J_1$. Counting rings and calling each a wavelength doubles $k$.

- **Rim node versus rim antinode.** The two Bessel zero families are each spaced by $\pi$ and interleaved, so the ladder spacing is the same for both; only the offset differs. The data pin the offset near $j_{1,n}$ (rim node). Do not read the bright rim ring as an antinode of the free wave.

- **Why the $m=1$ Laplacian.** The $-1/r^2$ term in the radial equation is the centrifugal term for one unit of angular momentum. It is why $J_1$ (and not $J_0$) appears, why $u(0)=0$, and why the pattern has zero amplitude at the centre with a phase singularity there: the point about which the spiral winds.

- **Self-consistency is weak but not zero.** Forgetting the tilt in $\Omega$ changes $f_n$ by 2–3 %. Not important for the picture; included in the numbers.

- **Sign of $\omega_f$.** Everywhere in Section 7 only $\omega=|\omega_f|$ enters; the sign fixed the handedness in Section 6 and plays no further role.
