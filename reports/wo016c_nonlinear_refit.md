# WO-016c — Non-linear flavio refit on LHCb 2025

Tests whether the linearised Mode-B response is sufficient at the fitted Delta C9 ~ -1.4. Three comparisons:
1. Linear fit (paper headline).
2. Non-linear evaluation at the linearised best-fit point (drift diagnostic — both models pinned at linear best-fit).
3. Non-linear refit (best-fit values found by `flavio.np_prediction`).

## Headline

| model | chi^2 | AIC | Delta AIC vs FREE_C9 | fit param |
|---|---:|---:|---:|---|
| FREE_C9_linear | 39.303 | 41.303 | +0.000 | DC9=-1.340 |
| FREE_C9_nonlinear@linear-best-fit | 66.596 | 68.596 | +0.000 | DC9=-1.340 |
| FREE_C9_nonlinear_refit | 40.891 | 42.891 | +0.000 | DC9=-1.002 |
| VFD_linear | 37.631 | 39.631 | -1.672 | A=+1.594 |
| VFD_nonlinear@linear-best-fit | 82.241 | 84.241 | +15.644 | A=+1.594 |
| VFD_nonlinear_refit | 41.983 | 43.983 | +1.093 | A=+1.135 |

- Linearised Delta AIC (FREE_C9 vs VFD): -1.672
- Non-linear Delta AIC at linear best-fit: +15.644 (diagnostic only; both models held at linear best-fit)
- Non-linear Delta AIC after refit: +1.093 (headline-comparable)
- Drift in headline Delta AIC: +2.765

## Best-fit parameters

- FREE_C9 linear: Delta C9 = -1.3403
- FREE_C9 non-linear refit: Delta C9 = -1.0025
- VFD linear: A = +1.5935
- VFD non-linear refit: A = +1.1350

## Per-bin linearisation residual (|nonlinear - linear| / sigma)

- FREE_C9 at linear best-fit: max = 3.106 sigma, mean = 0.630 sigma
- VFD at linear best-fit: max = 4.269 sigma, mean = 0.714 sigma

**Conclusion.** |Drift in headline Delta AIC| = 2.765 > 0.5 AIC unit. The headline must be updated to the non-linear refit value: Delta AIC_NL = +1.093 (vs linearised -1.672). The non-linear refit best-fit is DC9_FREE = -1.002 (vs linear -1.340) and A = +1.135 (vs linear +1.594).
