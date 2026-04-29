# WO-016b — Variant selection by pure-geometry criterion

Pure-geometry criterion: correlation between the discrete shell-mean of the V_600 Green's response and the continuum kernel $\kappa(x) = e^{-|x|/\varphi}$ from Layer 1 of the derivation. This criterion does **not** use LHCb data.

Data criterion: chi^2 against LHCb 2025 P5' on the joint fit.

| variant | corr(κ_continuum) | χ² (LHCb) | geom rank | data rank |
|---|---:|---:|---:|---:|
| FULL_LIFT[UNWEIGHTED]_GREENS | 0.9968 | 13.555 | 1 | 1 |
| FULL_LIFT[PHI_GEOMETRIC]_GREENS | 0.9130 | 14.713 | 2 | 2 |
| FULL_LIFT[PHI_ARITHMETIC]_GREENS | 0.8989 | 14.782 | 3 | 3 |

- Pure-geometry winner: **FULL_LIFT[UNWEIGHTED]_GREENS** (corr with continuum kernel = 0.9968)
- LHCb-data winner: **FULL_LIFT[UNWEIGHTED]_GREENS** (χ² = 13.555)

**Agreement.** The same variant (unweighted Laplacian) wins on both criteria. The variant choice is consistent with pure-geometry selection independent of the data; the LHCb data merely confirms it. Effective parameter count for VFD remains k=1 under the pure-geometry interpretation.
