# Analytical model of the self-generated standing spin wave

Files
- `wave_extract.py` – reads every `heat_map.npy` of the J sweep (`cofe_j_sweep`) and the diameter sweep
  (`size_CoFe_spectrum_3`), decomposes the complex FFT map into azimuthal harmonics e^{imφ}, and extracts the
  dominant m, the radial wavenumber k (slope of the unwrapped phase, 10–85 % of R) and the number of amplitude rings.
  Output: `wave_table.csv`.
- `sw_theory.py` – the model below, comparison with `wave_table.csv`, figure `standing_wave_theory.png`.

## 1. What the FFT maps actually contain (D = 320 nm, J = 6.0 / 5.5 ×10¹¹ A/m²)

| quantity | with BQ (13.1 GHz) | without BQ (12.0 GHz) |
|---|---|---|
| dominant azimuthal index of the Mz map at f_gyr | m = ±1 (97 % of power) | m = ±1 (98 %) |
| radial wavelength λ (phase slope) | 94 nm | 92 nm |
| amplitude maxima along a radius | r ≈ 28, 78, 123 nm + edge | 25, 78, 123 nm + edge |
| amplitude minima | 55, 100, 143 nm | 55, 100, 143 nm |
| total phase winding along a radius | 1.5 turns | 1.3 turns |

Node positions follow the zeros of J₁(kr) (kr = 3.83, 7.02, 10.17 → k ≈ 0.070 nm⁻¹) and the maxima follow its
extrema (1.84, 5.33, 8.54 → 26, 76, 122 nm).  So the pattern is a **J₁(kr)e^{±iφ} standing wave**: the amplitude
rings are spaced by **λ/2**, and the "spiral" in the phase map is a staircase with a π step at every node plus a
slow drift (net inward propagation), not a 4-turn spiral.  The manuscript's "4 rings" are the 3 interior
antinodes plus the edge antinode.

## 2. Model

### 2.1 Gyrating state
Free layer gyrates in-plane at ω with tilt m_z0.  Balance of the z-fields (demag −M_s m_z, biquadratic
−H_B m_z with H_B = 2|B₁₂|/(μ₀M_s d) = 0.354 MA/m ≙ 0.444 T) against the precession gives

    m_z0 = ω / [γ (M_s + H_B)]        (0.22 with BQ, 0.25 without, at 13 / 12 GHz)

This is the usual macrospin OPP result; BQ raises the frequency for a given tilt.

### 2.2 Co-rotating spin frame
Rotate the *spin space* (not real space) at ω about z.  The gyrating state becomes static and the net internal
field vanishes (the fictitious field −ω/γ cancels demag + BQ).  Because only spin space is rotated, a wave's
azimuthal index m, wavenumber k and frequency are the **same** in both frames for M_z (M_z is invariant under
rotation about z), while lab-static in-plane fields become rf fields.

Small excitations on the static, zero-internal-field in-plane film obey the ordinary thin-film dipole–exchange
dispersion (tilt θ₀, sin θ₀ = m_z0):

    Ω(k)² = γ² [ D k² + M_s(1−P) ⟨sin²φ_k⟩ ] · [ D k² + (M_s P + H_B) cos²θ₀ + M_s(1−P)⟨cos²φ_k⟩ sin²θ₀ ]

    D = 2A/(μ₀M_s) = 1.53×10⁻¹¹ A·m,   P(kd) = (1−e^{−kd})/(kd),   ⟨sin²φ_k⟩ = ⟨cos²φ_k⟩ = ½ (radial wave, azimuthal average)

Small-k limit: **Ω ≈ c k** with the Goldstone ("spin-superfluid sound") velocity

    c = γ √(D (M_s + H_B)) = 1006 m/s (no BQ), 1130 m/s (BQ)

The BQ field enters as extra out-of-plane stiffness; that is why the wave (and the step height) is faster with BQ.

### 2.3 The source and the wavenumber selection
The PMA fixed layer produces a static, radially symmetric in-plane stray field h(r) r̂ at the free-layer edge.
In the co-rotating spin frame this field rotates at −ω, i.e. it is an **rf drive at the gyration frequency with
azimuthal index m = ±1** (sign = handedness of the gyration).  It therefore resonantly launches an inward,
one-armed spiral spin wave at the wavenumber k* obeying

    Ω(k*) = ω                                                           (selection rule)

Test: for all 177 J-sweep and D-sweep runs with a clean m = ±1 pattern (8–35 GHz, with and without BQ), the
frequency predicted from the measured k by Ω(k*) = ω is 1.07 ± 0.08 × the simulated gyration frequency, with no
fitted parameter (left panel of `standing_wave_theory.png`).  The alternatives fail: exchange-only
ω = γDk² is 3–10× too low, and a 2ω (counter-rotating) assignment is 1.8× off.

### 2.4 Standing wave and the frequency staircase
The inward m = 1 wave must be regular at r = 0, so the pattern is J₁(k r) e^{iφ}.  With free spins at the rim
(Neumann, J₁'(kR) = 0; the rim is an antinode in all maps) the cavity admits

    k_n R = j'_{1,n} = 1.84, 5.33, 8.54, 11.71, 14.86, 18.02, 21.16, …

and the plateau frequencies are **f_n = Ω(j'_{1,n}/R)/2π** (n = number of interior antinodes):

| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| f_n theory, BQ (GHz) | 4.5 | 8.9 | 12.5 | 15.9 | 19.2 | 22.5 | 25.8 |
| end-of-plateau, simulation BQ | 5.7 | 9.4 | 13.1 | 16.8 | 19.0 | 22.7 | (25.4) |
| f_n theory, no BQ | 4.0 | 7.9 | 11.0 | 14.0 | 17.0 | 19.9 | 22.8 |
| end-of-plateau, simulation no BQ | 5.3 | 8.7 | 12.0 | 15.4 | 17.9 | 19.4 | 22.1 |

Within a plateau the STT pushes ω up while k* follows the dispersion; when the next half wavelength fits the
cavity the pattern re-locks with one more ring (jump).  The step height is one half-wavelength of phase,

    Δf ≈ (∂f/∂k)(π/R) ≈ c/(2R) = 3.1 GHz (no BQ), 3.5 GHz (BQ)   vs simulated 3.4 and 3.7 GHz.

Diameter sweep at fixed J: the same rule gives jumps every ΔR ≈ λ/2 ≈ 40 nm, i.e. ΔD ≈ 70–80 nm (simulated
jumps at D = 180, 250, 330, 400, 470, 530, 580, 650 nm with BQ), and f decreasing within a plateau because k_n =
j'_{1,n}/R shrinks as R grows.

## 3. Limitations / what is not in the model
- Ω(k) is a thin-film formula with the dipolar orientation factors averaged over the azimuth; the systematic +7 %
  is of the size expected from that averaging, the 2 nm cell, and the large wave amplitude (δm_z ≈ 0.17).
- The absolute f(J) (which plateau the STO sits on at a given J) needs the energy balance including wave damping
  and is not given analytically; the model gives the ladder of allowed plateau frequencies and the k–f relation.
- The exact edge phase (the rim is a driven, textured region ~15 nm wide) shifts k_n R by up to ~10 %; the grey
  bands in the figure span the Neumann (j'_{1,n}) to Dirichlet (j_{1,n}) limits.

## 4. Direct checks on the raw h5 data (last 5 ns, `frame_check.py`, results in `frame_check_results_raw.txt`)

Note: magnum.np arrays are indexed `[x, y]`; the figure scripts and the first version of the extraction used
`[row=y, col=x]`, which mirrors the handedness (m → −m).  Both scripts here now use the correct order.

| | with BQ, J = 6.0 | without BQ, J = 5.5 |
|---|---|---|
| in-plane gyration | −13.09 GHz (clockwise) | −11.93 GHz (clockwise) |
| local gyration rate vs r | uniform to ±3 % (phase-locked disk) | uniform to ±1 % |
| ⟨m_z⟩ (measured / m_z0 = ω/γ(M_s+H_B)) | 0.26 / 0.22 | 0.28 / 0.25 |
| M_z power at f, co-rotating spiral (m = +1) | 218 | 406 |
| M_z power at f, counter-rotating spiral | 2.0 | 2.4 |
| M_z power at f, rings (m = 0) | 0.3 | 0.1 |
| M_z power at 2f (m = 2, second harmonic of the spiral) | 18 | 21 |
| static ⟨M_z⟩ spiral | none (< 0.001) | none |

So the wave is a single one-armed spiral at the gyration frequency that **co-rotates with the magnetisation**, the
handedness the edge stray field drives in the co-rotating spin frame (Sec. 2.3); the counter-rotating partner
(which would need the two-fold anisotropy/edge texture to be excited) is 100× weaker.  The measured tilt is
~15 % larger than the macrospin m_z0 because the fixed layer's stray field inside the disk points along +z
(parallel to its own M) and partly cancels the demag field; including a stray-field term h_z,
ω = γ[(M_s + H_B) m_z − h_z], gives h_z ≈ 0.05–0.09 T.

## 5. Size dependence (supplementary figure 1) — `size_dependence.py`, `power_balance.py`, `size_dependence.png`

`power_balance.py` evaluates the exact LLG power balance (Gilbert dissipation α⟨|ṁ|²⟩ versus spin-torque work
γa_J⟨|(m×ṁ)_z|⟩) and the wave content in the co-rotating spin frame from the last 3 ns of every diameter-sweep
run (`power_balance.csv`).  `size_dependence.py` draws the four-row figure.

1. **Sawtooth = hopping between cavity modes.**  Within a plateau the radial phase kR is pinned
   (6.9, 10.7, 14.3, 17.6, 20.6 for n = 2…6 with BQ), so k = x_n/R and f = Ω(x_n/R)/2π falls as R grows.
   Each jump adds about one half wave; the plateau width is ΔD ≈ λ (60–80 nm against λ = 80–87 nm).
   The pinned phases lie 0–1 rad above the zeros j_{1,n} of J₁, i.e. closer to a rim-node than to the rim-antinode
   condition j'_{1,n} used in Sec. 2.4 (the agreement of Table 2.4 with j'_{1,n} partly reflects the +7 % bias of
   the dispersion).  The dashed curves Ω(j_{1,n}/R)/2π run 4–8 % above the plateaus, with no fit.
2. **Frequency is set by the wave power.**  For D = 180–700 nm, f = f_mac/(1 + c⟨a²⟩) with
   f_mac = γa_J/2πα (20.5 GHz at J = 6.0, 18.8 GHz at J = 5.5, not fitted) and one fitted constant
   c = 2.2 (BQ), 2.6 (no BQ); rms residual 0.4 GHz, correlation −0.88 / −0.84.  ⟨a²⟩ is the mean-square in-plane
   wave amplitude at the gyration frequency.  Inside a plateau ⟨a²⟩ grows (0.16 → 0.23–0.26) and drops by
   30–40 % at each jump: this is the feedback that pins kR and produces the falling branches.
3. **Fading of the steps and slow rise.**  The inward wave decays over 110–190 nm (measured in the
   500–700 nm disks; theory v_g/Γ = 173–193 nm), so beyond D ≈ 400–500 nm almost nothing returns from the
   centre and the cavity resonances disappear.  The area-averaged wave power keeps falling roughly as 1/R
   (0.15 at 500 nm, 0.12 at 700 nm), which is why f keeps creeping up rather than truly saturating.
   Without BQ the large disks are also much less phase-coherent (|⟨m₊⟩| = 0.2–0.5 against 0.7–0.8 with BQ).
4. **Not covered.**  D < 180 nm: the 3–5 GHz states are strongly non-uniform (dissipation 3–6× that of a rigid
   rotation) and are not a coherent gyration; the n = 1 branch appears at D = 120–150 nm (BQ), 130–170 nm (no BQ).
   The power balance closes to 80–86 % only; the missing part is nearly independent of D and not identified.

## 6. Self-contained notebook and study notes

- `standing_wave_model.ipynb` — executed Jupyter notebook that rebuilds the numerical part of the model
  (dispersion → selection rule → Bessel ladder → comparison with the 177 runs → summary figure
  `standing_wave_theory_notebook.png`, plus the diameter-sweep fan) with explanations for a reader or reviewer.
  Inputs are only the three CSV files in `data/` (`wave_table.csv`, `j_vs_freq_bq.csv`, `j_vs_freq_no_bq.csv`);
  needs numpy, scipy, matplotlib. To re-execute headless:
  `python -c "import nbformat; from nbclient import NotebookClient; nb=nbformat.read('standing_wave_model.ipynb',4); NotebookClient(nb,kernel_name='python3').execute(); nbformat.write(nb,'standing_wave_model.ipynb')"`
- Step-by-step notes on the derivation: `section3_notes.md` (macrospin), `section4_1_notes.md` (frame
  transformation), `section5_notes.md` (linearisation), `section6_7_notes.md` (source, cavity),
  `hessian_crash_course.md`, `bessel_crash_course.md` (+ `bessel_primer.png`).
