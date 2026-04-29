from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from vfd_b_anomaly import data_ingest
from vfd_b_anomaly.constants import (
    PROVENANCE_INFERRED,
    PROVENANCE_PLACEHOLDER,
    REQUIRED_OBSERVABLE_COLUMNS,
)


def test_placeholder_dataset_has_explicit_missing_flag():
    data = data_ingest.make_placeholder_dataset()
    df = data["observables"]
    assert set(REQUIRED_OBSERVABLE_COLUMNS).issubset(df.columns)
    assert (df["provenance"] == PROVENANCE_PLACEHOLDER).all()
    # Missing values are NaN by design - silently ignoring them is forbidden.
    assert df["value"].isna().all()
    assert data["metadata"]["covariance_status"] == "absent_declared"


def test_placeholder_dataset_passes_schema_validation():
    data = data_ingest.make_placeholder_dataset()
    data_ingest.validate_observable_schema(data)


def test_validate_rejects_missing_top_level_keys():
    with pytest.raises(ValueError, match="missing required top-level keys"):
        data_ingest.validate_observable_schema({"observables": pd.DataFrame()})


def test_validate_rejects_unsupported_observable():
    bad = data_ingest.make_placeholder_dataset()
    bad["observables"].loc[0, "observable"] = "ZZZ"
    with pytest.raises(ValueError, match="unsupported observables"):
        data_ingest.validate_observable_schema(bad)


def test_validate_rejects_nonpositive_stat_err():
    bad = data_ingest.make_placeholder_dataset()
    bad["observables"].loc[0, "stat_err"] = 0.0
    with pytest.raises(ValueError, match="stat_err must be > 0"):
        data_ingest.validate_observable_schema(bad)


def test_validate_rejects_inverted_q2_bin():
    bad = data_ingest.make_placeholder_dataset()
    # Force q2_lo >= q2_hi on a single row.
    bad["observables"].loc[0, "q2_lo"] = bad["observables"].loc[0, "q2_hi"]
    with pytest.raises(ValueError, match="q2_lo < q2_hi"):
        data_ingest.validate_observable_schema(bad)


def test_validate_rejects_silent_missing_covariance():
    bad = data_ingest.make_placeholder_dataset()
    bad["metadata"]["covariance_status"] = "missing"  # not the allowed sentinel
    with pytest.raises(ValueError, match="covariance_status"):
        data_ingest.validate_observable_schema(bad)


def test_validate_rejects_empty_provenance():
    bad = data_ingest.make_placeholder_dataset()
    bad["observables"].loc[0, "provenance"] = ""
    with pytest.raises(ValueError, match="non-empty provenance"):
        data_ingest.validate_observable_schema(bad)


def test_load_covariance_returns_none_for_missing_path(tmp_path):
    assert data_ingest.load_covariance_matrix(None) is None
    assert data_ingest.load_covariance_matrix(tmp_path / "no_such.npy") is None


def test_load_covariance_round_trip(tmp_path):
    cov = np.eye(4) + 0.01
    path = tmp_path / "cov.npy"
    np.save(path, cov)
    loaded = data_ingest.load_covariance_matrix(path)
    assert loaded is not None
    assert loaded.shape == (4, 4)


def test_load_covariance_rejects_non_symmetric(tmp_path):
    cov = np.array([[1.0, 0.0], [0.5, 1.0]])
    path = tmp_path / "cov.npy"
    np.save(path, cov)
    with pytest.raises(ValueError, match="not symmetric"):
        data_ingest.load_covariance_matrix(path)


def test_write_and_load_round_trip(tmp_path, synthetic_data):
    paths = data_ingest.write_processed_dataset(synthetic_data, tmp_path)
    assert paths["csv"].exists() and paths["metadata"].exists()
    loaded = data_ingest.load_lhcb_tables(tmp_path)
    pd.testing.assert_frame_equal(
        loaded["observables"].reset_index(drop=True),
        synthetic_data["observables"].reset_index(drop=True),
    )
    assert loaded["metadata"]["covariance_status"] == synthetic_data["metadata"]["covariance_status"]


def test_write_refuses_inconsistent_covariance(tmp_path, synthetic_data):
    n = len(synthetic_data["observables"])
    cov = np.eye(n)
    # metadata still says 'absent_declared' - the writer must refuse.
    with pytest.raises(ValueError, match="refusing to write inconsistent dataset"):
        data_ingest.write_processed_dataset(synthetic_data, tmp_path, covariance=cov)


def test_provenance_round_trip_synthetic(synthetic_data):
    df = synthetic_data["observables"]
    # Synthetic factory marks every row as 'inferred' - readers must see that.
    assert (df["provenance"] == PROVENANCE_INFERRED).all()
