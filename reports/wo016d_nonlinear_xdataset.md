# WO-016d — Non-linear cross-dataset refit

Re-runs the WO-014 cross-dataset and WO-015 cross-channel fits with `flavio.np_prediction` directly (non-linear) instead of the Mode-B Taylor expansion. The linearised values from reports/wo014_cross_dataset.csv and reports/wo015_cross_channel.csv are quoted in parentheses for comparison.

| dataset | non-linear χ² (FREE) | non-linear χ² (VFD) | ΔAIC (NL) | ΔC9 (NL) | A (NL) |
|---|---:|---:|---:|---:|---:|
| LHCb-2015 | 30.691 | 30.450 | -0.241 | -1.080 | +1.235 |
| LHCb-2021-Kstplus | 22.765 | 22.933 | +0.168 | -1.820 | +2.059 |
| CMS-2025-noP4p | 43.731 | 44.205 | +0.473 | -0.954 | +1.053 |
| LHCb-2025 | 40.891 | 41.983 | +1.093 | -1.003 | +1.135 |
| Bs2phi-LHCb-2015 | 13.201 | 12.962 | -0.240 | -4.122 | +4.984 |

