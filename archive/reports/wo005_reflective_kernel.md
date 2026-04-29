# WO-005 — Reflective / cavity-mode kernel test

**Hypothesis under test.** The B0 -> K*0 mu mu Re(C9) residual is dominated by
*boundary-reflected* / cavity-mode behaviour: the closure window acts as a
bounded resonator, with maximum residual amplitude near the kinematic edges
(low-q^2 and high-q^2) and a node near the J/psi-psi(2S) midpoint.

**Operationalisation.** Three single-amplitude (k = 1) kappa kernels added to
`vfd_closure.py`, all pre-frozen with the same dimensionless coordinate
x = (q^2 - midpoint) / Delta_psi and a kinematic boundary x_max = 2.886
derived from q^2_min = 4 m_mu^2 and q^2_max = (m_B - m_K*)^2:

| mode | shape | peak |
| --- | --- | --- |
| `kappa_reflective_centre` | `cosh((x_max - |x|)/phi) / cosh(x_max/phi)` | midpoint |
| `kappa_reflective_edge` | `sinh(|x|/phi) / sinh(x_max/phi)` | boundaries |
| `kappa_reflective_cavity` | `sin(pi |x| / (2 x_max))` | boundaries |

Compared against the prior champion `kappa_exponential` and the FREE_C9 reference.

## Result on real LHCb config-2 P5' data (8 bins, 8.4 fb^-1)

`FREE_C9` reference: chi^2 = 6.7046, AIC = 8.7046.

| mode | k | chi^2 | dAIC vs FREE_C9 | eff DC9 |
| --- | --- | --- | --- | --- |
| `VFD_KAPPA_EXPONENTIAL`        | 1 | 5.904 | -0.80 | -0.123 |
| `VFD_KAPPA_REFLECTIVE_CENTRE`  | 1 | 5.943 | -0.76 | -0.136 |
| `VFD_KAPPA_YUKAWA_SMOOTH`      | 1 | 6.014 | -0.69 | -0.095 |
| `VFD_KAPPA_GAUSSIAN`           | 1 | 6.142 | -0.56 | -0.107 |
| `VFD_KAPPA_LORENTZIAN`         | 1 | 6.226 | -0.48 | -0.081 |
| `VFD_KAPPA_SECH`               | 1 | 6.281 | -0.42 | -0.135 |
| `VFD_KAPPA_CANONICAL`          | 1 | 6.365 | -0.34 | -0.141 |
| `VFD_KAPPA_REFLECTIVE_CAVITY`  | 1 | 8.636 | **+1.93** | -0.115 |
| `VFD_KAPPA_REFLECTIVE_EDGE`    | 1 | 8.927 | **+2.22** | -0.120 |

## Diagnostic — kappa values at bin centres

```
q^2 centres : [ 0.52  1.80  3.25  5.00  7.00 11.75 16.00 18.00]
x           : [-2.77 -2.45 -2.08 -1.65 -1.15  0.04  1.10  1.60]

kappa_reflective_centre  : [0.33 0.34 0.37 0.43 0.53 0.98 0.55 0.43]
kappa_reflective_edge    : [0.92 0.75 0.58 0.42 0.27 0.01 0.25 0.40]
kappa_reflective_cavity  : [1.00 0.97 0.91 0.78 0.59 0.02 0.57 0.77]
kappa_exponential        : [0.18 0.22 0.28 0.36 0.49 0.98 0.51 0.37]
```

## Conclusion

**The boundary-reflected / cavity-standing-mode hypothesis is REJECTED on
single-observable P5' data.**

- Both edge-peaked kernels (`reflective_edge`, `reflective_cavity`) score
  **WORSE** than the FREE_C9 control (dAIC > +1.9), i.e. worse than a global
  C9 shift with the same k = 1.
- The *centre*-peaked cosh kernel (`reflective_centre`) tracks the prior
  champion `kappa_exponential` to within 0.04 in chi^2 and 0.02 in
  effective DC9 — both peak strongly at the J/psi-psi(2S) midpoint.

The compression centre is **at the resonance midpoint, not at the kinematic
boundaries**. The "kaleidoscope reflecting back at the edges" picture does
not match the P5' data: the residual rises towards the midpoint and decays
towards both boundaries, which is the opposite topology.

This is a meaningful negative result: it rules out a whole family of
boundary-anchored closure kernels. It does not refute the broader VFD
closure picture - it constrains it. The next discriminating test is the
*multi-observable consistency* check (WO-027): if a single midpoint-peaked
kappa with one shared amplitude A explains P5', P4', AFB, FL and S_i
simultaneously, that is the structural signature the model needs.

Run command (for reproduction):

```
python3 -c "
from vfd_b_anomaly import hepdata_ingest, model_compare
archive = hepdata_ingest.hepdata_archive_dir('data/raw/hepdata/extracted')
data = hepdata_ingest.load_config(archive, config_index=2, observables=('P5p',))
print(model_compare.compare_all(data, include_binned=False).to_string())
"
```
