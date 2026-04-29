# Archive

Historical scripts and reports that the paper *cites as supporting
evidence* but does not require for headline reproduction. Kept in the
repo for transparency: a reviewer who wants to read §5's stress-test
language back to the actual outputs can find them here.

Nothing in this directory is on the path of `repro/run_all.sh`. The
headline reproduction works without it.

## Contents

```
archive/
├── scripts/                 # earlier WO modules referenced in §3 / §5
│   ├── wo006_multi_observable.py     # joint P5'/P4'/P1/P2 fit (early)
│   ├── wo007_eigenmodes.py           # Dirichlet eigenmode derivation
│   ├── wo008_discrete_lift.py        # discrete VFD lift of phi-kernel
│   ├── wo011_spectral.py             # spectral decomposition of L_V600
│   └── wo013_stress_test.py          # linearised stress test (bin
│                                       bootstrap, region splits, BSZ)
└── reports/                 # outputs of the above + linearised
                              # cross-dataset (wo014) and cross-channel
                              # (wo015) fits, plus the early reflective-
                              # kernel test (wo005). All of these are
                              # superseded by the non-linear refit
                              # outputs in the top-level reports/
                              # directory.
```

## Why these were moved out of the publish surface

The paper's headline numbers come from the non-linear flavio refit
(`reports/wo016c_*`, `wo016d_*`, `wo016a_*`, `wo016b_*`). The
linearised cross-dataset and stress-test outputs are diagnostic, not
headline; the paper retains them only as a methodology comparison.
Keeping them in `archive/` keeps the top-level `reports/` directory
small and unambiguous about what is current.
