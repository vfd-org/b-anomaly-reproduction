# WO-016a — Akaike-weight stack across five fits

Per-dataset AIC deltas and Akaike weights, plus stacked weight.
Stacking assumes independence under the null hypothesis (the five datasets share no observation-level information).

| dataset | FREE_C9 ΔAIC | VFD ΔAIC | w(FREE_C9) | w(VFD) |
|---|---:|---:|---:|---:|
| LHCb-2015 | 0.241 | 0.000 | 0.4699 | 0.5301 |
| LHCb-2021-Kstplus | 0.000 | 0.168 | 0.5210 | 0.4790 |
| CMS-2025-noP4p | 0.000 | 0.473 | 0.5589 | 0.4411 |
| LHCb-2025 | 0.000 | 1.093 | 0.6333 | 0.3667 |
| Bs2phi-LHCb-2015 | 0.240 | 0.000 | 0.4701 | 0.5299 |

## Stacked

- log-evidence(FREE_C9) − log-evidence(VFD) = 0.627
- Total ΔAIC sum (FREE_C9 vs VFD): -1.253

| model | stacked Akaike weight |
|---|---:|
| FREE_C9 | 0.6517 |
| VFD_GREEN_600CELL | 0.3483 |

Auxiliary check: under the null hypothesis P(VFD lower AIC on all 5 fits) = $1/2^{5}$ = 0.0312.
