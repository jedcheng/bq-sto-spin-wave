# Theory

* `derivation.tex` / `derivation.pdf` - full derivation of the analytical model (co-rotating spin frame, dispersion
  Omega(k), stray-field drive, Bessel cavity ladder, power balance). Compiles with `pdflatex derivation.tex` (twice);
  it includes `standing_wave_theory.png` (from `analysis/sw_theory.py`) and `size_dependence.png` (an earlier version
  of Fig. S1). The manuscript's SM3 is a condensed version of it.
* `standing_wave_model.ipynb` - step-by-step numerical walk-through of the model (reads `../data/derived/`; run it from
  this folder). It keeps the outputs of its last execution.
* `notes/` - study notes written while building the model (`model_overview.md` summarises what the FFT maps contain and
  the model; the section notes follow the derivation). They are informal working notes, not part of the paper.
