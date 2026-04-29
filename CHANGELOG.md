# Changelog

All notable changes and findings for the `vfd-b-anomaly` project.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/),
with entries describing what was learned and what artefacts were produced —
not just code changes.

---

## [0.2.0] — 2026-04-29 — Paper preprint

### Added
- **Preprint** at [`paper/main.pdf`](paper/main.pdf) (25 pages). Full LaTeX
  source under `paper/`, three figures (kernel shape, bin pulls,
  cross-dataset amplitudes), three rounds of internal hostile review.
- **Non-linear flavio refit** across all five datasets
  ([`scripts/wo016c_nonlinear_refit.py`](scripts/wo016c_nonlinear_refit.py),
  [`scripts/wo016d_nonlinear_xdataset.py`](scripts/wo016d_nonlinear_xdataset.py)).
- **Akaike-weight stacking** of $\Delta\mathrm{AIC}$ across the five fits
  ([`scripts/wo016a_akaike_stack.py`](scripts/wo016a_akaike_stack.py)).
- **Pure-geometry variant-selection test** for the kernel
  ([`scripts/wo016b_variant_geometry.py`](scripts/wo016b_variant_geometry.py)),
  showing the unweighted Laplacian wins on both pure-geometry
  (correlation with continuum kernel) and LHCb-data criteria — defending
  $k=1$ for the kernel against the Round-1 hostile-review concern.
- **End-to-end reproduction driver** [`repro/run_all.sh`](repro/run_all.sh).
- **Paper figures** F1, F2, F3
  ([`scripts/wo017_paper_figures.py`](scripts/wo017_paper_figures.py)),
  embedded in `paper/figures/` (PDF and PNG).

### Major finding (negative-of-headline)
The earlier linearised Mode-B fit on LHCb 2025 had given
$\Delta\mathrm{AIC} = -1.67$ in favour of the geometry-derived kernel
over a constant Wilson-coefficient shift ($\mathrm{FREE\_C9}$).

The non-linear refit using `flavio.np_prediction` directly gives
$\Delta\mathrm{AIC} = +1.09$ on the same dataset — a $+2.77$ AIC-unit
shift that **flips the sign of the kernel-vs-constant comparison**. The
per-bin linearisation residual reaches $4.3\sigma$ at the linearised
best-fit $\Delta C_{9}=-1.34$, well outside the linear regime.

Stacked Akaike weight across the five non-linear fits:
$w_{\mathrm{VFD}}=0.348$ vs $w_{\mathrm{FREE\_C9}}=0.652$. The kernel
is statistically indistinguishable from a constant shift on AIC.

What survived non-linear evaluation:
- Sign-uniformity: $A>0$ in **5/5** fits, $\Delta C_{9}^{\mathrm{eff}}<0$
  in **5/5** fits.
- The kernel describes the $q^{2}$ shape with one parameter per
  dataset, no shape retuning.

The paper was rewritten around this negative finding, with all
linearised numbers retained as a methodology diagnostic.

### Changed
- Paper headline reframed from "kernel beats $\mathrm{FREE\_C9}$ on AIC"
  (linearised) to "kernel is sign-uniform but AIC-degenerate with
  $\mathrm{FREE\_C9}$" (non-linear).
- Falsification programme item ("non-linear refit gives
  $\Delta\mathrm{AIC} \geq 0$") explicitly noted as fired.

---

## [0.1.0] — 2026-04-28 — Linearised cross-dataset and cross-channel test

### Added
- Five-dataset linearised Mode-B fit
  ([`src/vfd_b_anomaly/wo014_cross_dataset.py`](src/vfd_b_anomaly/wo014_cross_dataset.py)).
- Cross-channel $B_{s}\!\to\!\phi$ fit with basis-effect diagnostic
  ([`src/vfd_b_anomaly/wo015_cross_channel.py`](src/vfd_b_anomaly/wo015_cross_channel.py)).
- Stress-test suite (bin bootstrap, region splits, alternative
  Wilson-coefficient models, BSZ form-factor MC)
  ([`archive/scripts/wo013_stress_test.py`](archive/scripts/wo013_stress_test.py)).
- Multi-observable joint fit (P5', P4', P1, P2)
  ([`src/vfd_b_anomaly/wo010_universality.py`](src/vfd_b_anomaly/wo010_universality.py)).
- 600-cell $V_{600}$ graph construction with binary-icosahedral
  symmetry ([`src/vfd_b_anomaly/wo009_full_lift.py`](src/vfd_b_anomaly/wo009_full_lift.py)).
- Spectral-decomposition diagnostic
  ([`archive/scripts/wo011_spectral.py`](archive/scripts/wo011_spectral.py)).
- flavio + wilson SM and Wilson-coefficient backend
  ([`src/vfd_b_anomaly/flavio_predictor.py`](src/vfd_b_anomaly/flavio_predictor.py)).

### Findings (linearised; later superseded by 0.2.0 non-linear refit)
- All five datasets gave $\Delta\mathrm{AIC} < 0$ vs FREE\_C9, suggesting
  the kernel beat the constant shift in linearised analysis.
- Stacked Akaike weight $w_{\mathrm{VFD}}=0.904$ vs
  $w_{\mathrm{FREE\_C9}}=0.096$ (linearised).
- Cross-channel $B_{s}\!\to\!\phi$: $A=+5.48$ vs $B\!\to\!K^{*}$
  $A\approx +1.6$ — gap explained by Krüger--Matias $P$-basis
  amplification factor $\approx 2.2$.

These linearised numbers are retained in the paper as a methodology
diagnostic and superseded by the 0.2.0 non-linear refit.

---

## [0.0.1] — 2026-04-27 — Initial reproduction

- Reproduced the negative direction of $\Delta C_{9}$ on the LHCb 2025
  $P_{5}'$ data using a one-observable Mode-B linearised fit.
- Hand-tabulated SM baselines (later replaced by flavio in 0.1.0
  because the hand-tabulated tables had wrong-sign
  $\partial P_{5}'/\partial C_{9}$ in several bins).
- 64-test pytest suite covering schema, ingest, fits, and
  parameter-count invariance.
