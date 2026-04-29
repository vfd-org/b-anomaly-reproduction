"""Shared fixtures: synthetic LHCb-shaped datasets used across tests.

Synthetic data is built by injecting a known ΔC9_true via the linearised
SM response. Tests use this to verify that fits recover the injected value.
This is deliberate scaffolding - it is NOT a claim of reproducing the paper.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from vfd_b_anomaly.constants import (
    C9_SM,
    PROVENANCE_INFERRED,
    REQUIRED_OBSERVABLE_COLUMNS,
)
from vfd_b_anomaly.paper_tables import DEFAULT_Q2_BIN_EDGES_GEV2
from vfd_b_anomaly.sm_baseline import dO_dC9, sm_value


def _synthetic_rows(delta_c9_true: float, observables: tuple[str, ...], stat_err: float, syst_err: float, rng: np.random.Generator):
    rows = []
    n_bins = len(DEFAULT_Q2_BIN_EDGES_GEV2)
    for obs in observables:
        for b, (q2_lo, q2_hi) in enumerate(DEFAULT_Q2_BIN_EDGES_GEV2):
            sm = sm_value(obs, b)
            slope = dO_dC9(obs, b)
            true_val = sm + slope * delta_c9_true
            noise = rng.normal(0.0, stat_err)
            rows.append(
                {
                    "q2_lo": float(q2_lo),
                    "q2_hi": float(q2_hi),
                    "observable": obs,
                    "value": float(true_val + noise),
                    "stat_err": float(stat_err),
                    "syst_err": float(syst_err),
                    "provenance": PROVENANCE_INFERRED,
                }
            )
    df = pd.DataFrame(rows, columns=list(REQUIRED_OBSERVABLE_COLUMNS))
    return df


@pytest.fixture
def synthetic_data_factory():
    def _make(delta_c9_true: float = -1.0,
              observables: tuple[str, ...] = ("FL", "AFB", "P5p", "BR"),
              stat_err: float = 0.02,
              syst_err: float = 0.01,
              seed: int = 12345) -> dict:
        rng = np.random.default_rng(seed)
        df = _synthetic_rows(delta_c9_true, observables, stat_err, syst_err, rng)
        return {
            "observables": df,
            "metadata": {
                "source": "synthetic",
                "covariance_status": "absent_declared",
                "date_accessed": "synthetic",
                "delta_c9_true": delta_c9_true,
                "notes": "Synthetic data injected via linear MODE_B response.",
            },
        }
    return _make


@pytest.fixture
def synthetic_data(synthetic_data_factory):
    return synthetic_data_factory(delta_c9_true=-1.0)


@pytest.fixture
def c9_sm():
    return C9_SM
