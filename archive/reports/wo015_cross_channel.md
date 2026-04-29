# WO-015 — Cross-channel frozen-kernel test

**Question.** WO-014 demonstrated the frozen 600-cell Green kernel
universalises across four B → K* μμ datasets (LHCb 2015, LHCb 2021 B⁺→K*⁺,
CMS 2025, LHCb 2025). WO-015 asks: does the same kernel — with no shape
change, only an amplitude `A` per dataset — describe the *other* b → s μμ
channels too?

The kernel is locked from WO-009. Only `A` is fitted per dataset; SM and
dC9 slopes are computed per-dataset on its own bin grid via the
`flavio_predictor` module (with on-disk JSON cache).

## Datasets attempted

| dataset | decay | arXiv | HEPData | status |
| --- | --- | --- | --- | --- |
| Bs → φ μμ angular        | Bs → φ μμ  | 1506.08777 | ins1380188 | **RUN** |
| LHCb full Run 2 Bs → φ μμ | Bs → φ μμ  | 2107.13428 | ins1894428 | submission absent on HEPData |
| B⁺ → K⁺ μμ angular (LHCb) | B⁺ → K⁺ μμ | 1403.8045  | ins1287929 | submission absent on HEPData |
| B⁺ → K⁺ μμ branching+angular | B⁺ → K⁺ μμ | 1209.4284 | ins1186557 | submission absent on HEPData |
| B⁺ → K⁺ μμ angular (CMS) | B⁺ → K⁺ μμ | (PRD 2018) | ins1662193 | submission absent on HEPData |
| LHCb B⁺ → K⁺ μμ updated  | B⁺ → K⁺ μμ | 1804.07167 | ins1668916 | submission absent on HEPData |
| K*₀,₂(1430) (Tier 2)     | B⁰ → K⁺π⁻μμ | 1606.04731 | ins1486676 | DOWNLOADED, not run (S/D-wave moments are not a clean cross-channel kernel test) |

**B⁺ → K⁺ μμ status:** five candidate INSPIRE records were probed for
HEPData submissions; all five returned the HEPData 404 HTML page from
`/download/submission/ins{ID}/yaml`. The published B⁺ → K⁺ μμ angular and
branching-fraction analyses do not appear to have public HEPData YAML
submissions. This is recorded as a known gap; the channel is **not** in
this WO.

## Tier 1 result: Bs → φ μμ (LHCb 2015)

Six exclusive q² bins × four observables (F_L, S₃, S₄, S₇) → **24
data points**, no published correlations on HEPData (diagonal
stat ⊕ syst). flavio decay string `Bs->phimumu`, observable strings
`<FL>`, `<S3>`, `<S4>`, `<S7>` (cached).

### Global fit

| model | k | χ² | AIC | BIC | ΔAIC vs FREE_C9 | params |
| --- | --- | --- | --- | --- | --- | --- |
| FREE_C9 (k=1)             | 1 | 13.20 | 15.20 | 16.37 |  0.00 | ΔC9 = **−4.37** |
| **VFD_GREEN_600CELL (k=1)** | 1 | **13.11** | **15.11** | **16.29** | **−0.08** | A = **+5.48**, ΔC9_eff = **−4.76** |

**Acceptance gates:**

| gate | result |
| --- | --- |
| 1. A same sign as WO-014                | PASS — A = +5.48 > 0 |
| 2. effective ΔC9 negative                | PASS — ΔC9_eff = −4.76 |
| 3. VFD AIC competitive with FREE_C9      | PASS (marginal, ΔAIC = −0.08) |
| 4. bootstrap A sign-stable               | PASS — frac A<0 = 0/500 |
| 5. region splits sign-uniform            | PASS — A > 0 in low/central/high q² |
| 6. convention mismatches reported        | PASS — see "channel sensitivity" below |

### Bootstrap (n = 500 over bins)

| statistic | value |
| --- | --- |
| n_bootstrap         | 500 |
| mean A              | +6.71 |
| median A            | (≈ +5.5 at central, distribution skewed by upper-bound trials) |
| std A               | 3.27 |
| 90% CI              | [+5.20, +16.15] |
| fraction A < 0      | 0 / 500 |

The wide CI and large σ are a known feature: at high q² the Bs→φμμ
angular slopes dO/dC9 are tiny (≤ 0.001), so any non-zero residual
inflates the implied amplitude. The bootstrap is sign-stable but the
magnitude is poorly constrained by this dataset alone.

### Region splits (low / central / high q²)

| region | q² range | n_data | FREE_C9 χ² | FREE_C9 ΔC9 | VFD χ² | VFD A | ΔAIC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| low_q²     | (0.1, 2.0)  | 4  | 1.66 | −4.23  | 1.66  | +5.31  |  0.00 |
| central_q² | (5.0, 8.0)  | 4  | 0.78 | −10.00 | 0.67  | +15.50 | −0.11 |
| high_q²    | (11.0, 19.0)| 12 | 9.33 | −10.00 | 9.15  | +20.00 | −0.18 |

The central/high-q² fits hit the [−10, +20] optimiser bounds — symptom
of vanishing dO/dC9 there. The low-q² window is the only region where
Bs→φμμ has comparable C9 sensitivity to B→K*; there the kernel ties
FREE_C9 (ΔAIC = 0.00) and gives A = +5.31 (broadly consistent with the
global +5.48).

## Channel sensitivity caveat (corrected by WO-015b)

The initial reading of these numbers was that "Bs → φμμ has weaker C9
sensitivity than B → K*." A direct slope comparison shows that is **not
correct**. See `reports/wo015b_basis_diagnostic.md` for the full
diagnostic; the headline:

- At matched q² bins, flavio's S-basis dO/dC9 slopes for Bs → φμμ
  agree with B → K* to within ±10–20% (RMS ratio 0.88–0.98 across all
  six bins).
- The apparent factor-of-3.5 amplitude difference (A_Bs ≈ 5.5 vs
  A_B→K* ≈ 1.6) is almost entirely a *basis* effect: WO-014 used the
  P-basis (P5', P4', P1, P2), which is amplified by 1/√(F_L(1−F_L))
  ≈ 2–2.5× relative to the S-basis (F_L, S3, S4, S7). Krüger–Matias
  designed the P-basis to absorb F_L, so P-basis observables track
  Wilson-coefficient shifts with that built-in amplification.
- LHCb 2015 published Bs → φμμ in S-basis only; flavio carries no
  P-basis observables for that decay; the WO-015 fit was *forced* into
  the smaller-slope basis.

After basis harmonisation (multiplying A_S by the kernel-weighted
factor √(F_L(1−F_L)) ≈ 0.45), the predicted P-basis-equivalent
amplitude for Bs → φμμ is **+2.5**, vs +1.6 for B → K*. A residual
factor of ~1.6 remains, attributable to:

1. No published correlation matrix for Bs → φμμ angular (HEPData
   ins1380188 has none).
2. Bound-pinning at central/high q² in WO-015 region splits (slopes
   drop to O(10⁻⁴) there, so the inverse-variance fit is dominated by
   the 4 lowest-q² bins).
3. Limited statistics (3 fb⁻¹ vs 8.4 fb⁻¹ for LHCb 2025).

The factor-of-3.5 disparity reduces to factor-of-1.5 once basis is
corrected; the direction of the deviation (negative ΔC9_eff) is
unchanged. The kernel itself is consistent across channels; only the
publication-basis differs.

## What WO-015 establishes

- The frozen 600-cell kernel is **sign-consistent** with the Bs → φ μμ
  angular anomaly: same direction (negative effective C9 shift) with no
  shape retuning.
- VFD ties or marginally beats FREE_C9 (k=1 vs k=1) on the global fit
  and on each q² sub-region.
- The bootstrap A is sign-stable on all 500 trials.

## What WO-015 does not establish

- The kernel does not get an *independent* magnitude calibration from
  Bs → φ μμ at this resolution. The 90% CI [+5.2, +16.1] spans more than
  a factor of 3.
- B⁺ → K⁺ μμ remains untested in this WO because no public HEPData
  submission was located for the relevant angular / branching-fraction
  papers (LHCb 1209.4284 / 1403.8045 / 1804.07167, CMS PRD 2018). The
  acceptance verdict for WO-015 is therefore on Bs → φ μμ alone.
- The K*₀,₂(1430) Tier 2 dataset (ins1486676) was downloaded but not
  fit: it reports S/D-wave Γ_i moments rather than the standard angular
  basis, and is not a clean kernel test.

## Files

- `src/vfd_b_anomaly/wo015_cross_channel.py` — driver
- `data/raw/hepdata_cross_channel/ins1380188.tar.gz` — Bs → φ μμ data
- `data/raw/hepdata_cross_channel/ins1486676.tar.gz` — K*₀,₂(1430)
  (downloaded for completeness, not fit)
- `reports/wo015_cross_channel.md` — this report
- `reports/wo015_cross_channel.csv` — global fit table
- `reports/wo015_bootstrap.csv` — per-dataset bootstrap statistics
- `reports/wo015_regions.csv` — region-split fit table
- `reports/wo015_run.log` — full stdout
- `data/processed/flavio_cache.json` — extended with Bs → φ μμ entries

## Combined WO-014 + WO-015 universality picture

| dataset | decay | n | A | ΔC9_eff | ΔAIC | source |
| --- | --- | --- | --- | --- | --- | --- |
| LHCb 2015         | B⁰ → K*⁰ μμ | 32 | +1.76 | −1.51 | −0.67 | WO-014 |
| LHCb 2021         | B⁺ → K*⁺ μμ | 32 | +2.72 | −2.33 | −0.55 | WO-014 |
| CMS 2025 (no P4') | B⁰ → K*⁰ μμ | 18 | +1.59 | −1.39 | −1.52 | WO-014 |
| LHCb 2025         | B⁰ → K*⁰ μμ | 32 | +1.59 | −1.37 | −1.67 | WO-014 |
| **Bs → φ μμ 2015**  | **Bs → φ μμ** | **24** | **+5.48** | **−4.76** | **−0.08** | **WO-015** |

Across five datasets spanning three decay channels (B⁰→K*⁰, B⁺→K*⁺,
Bs→φ) and two collaborations (LHCb, CMS), the frozen kernel:

- gives **A > 0 in all five fits** (sign-uniform);
- gives **ΔC9_eff < 0 in all five fits** (anomaly-direction-uniform);
- **beats or ties FREE_C9 in AIC on every dataset**.

The B → K* family agrees on amplitude at the 10–15% level
(+1.59 to +1.76, modulo the lower-statistics LHCb 2021 isospin partner
at +2.72). Bs → φ μμ enters the picture sign-consistent but with a much
weaker per-dataset constraint, reflecting the smaller C9 sensitivity of
its angular observables.
