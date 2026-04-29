from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from vfd_b_anomaly import data_ingest, hepdata_ingest
from vfd_b_anomaly.constants import (
    PROVENANCE_REPRODUCED,
    REQUIRED_OBSERVABLE_COLUMNS,
)
from vfd_b_anomaly.paper_tables import DEFAULT_Q2_BIN_EDGES_GEV2


ARCHIVE = Path(__file__).resolve().parent.parent / "data" / "raw" / "hepdata"

skip_no_archive = pytest.mark.skipif(
    not (ARCHIVE / "extracted" / "submission.yaml").exists(),
    reason="HEPData archive not extracted; run download/extract step first",
)


@skip_no_archive
def test_archive_dir_resolves():
    extracted = hepdata_ingest.hepdata_archive_dir(ARCHIVE)
    assert (extracted / "submission.yaml").exists()


@skip_no_archive
def test_load_config_2_p5p_basic_shape():
    data = hepdata_ingest.load_config(ARCHIVE, config_index=2, observables=("P5p",))
    df = data["observables"]
    # Eight exclusive q² bins, one observable (P5p).
    assert len(df) == len(DEFAULT_Q2_BIN_EDGES_GEV2)
    assert (df["observable"] == "P5p").all()
    assert (df["provenance"] == PROVENANCE_REPRODUCED).all()
    # Public schema columns must be present.
    for c in REQUIRED_OBSERVABLE_COLUMNS:
        assert c in df.columns


@skip_no_archive
def test_load_config_2_uses_lhcb_bin_edges():
    data = hepdata_ingest.load_config(ARCHIVE, config_index=2, observables=("P5p",))
    df = data["observables"].sort_values("q2_lo").reset_index(drop=True)
    expected_edges = list(DEFAULT_Q2_BIN_EDGES_GEV2)
    actual_edges = list(zip(df["q2_lo"].astype(float), df["q2_hi"].astype(float)))
    assert actual_edges == expected_edges


@skip_no_archive
def test_load_config_2_p5p_values_match_known_lhcb_signs():
    """Sanity: known LHCb 2025 P5' is positive at low q² and negative above ~ 4 GeV²."""
    data = hepdata_ingest.load_config(ARCHIVE, config_index=2, observables=("P5p",))
    df = data["observables"].sort_values("q2_lo").reset_index(drop=True)
    # bin 0 [0.06, 0.98]: P5' positive
    assert df.loc[0, "value"] > 0.0
    # bin 1 [1.1, 2.5]: P5' positive
    assert df.loc[1, "value"] > 0.0
    # bin 3 [4.0, 6.0] onward: P5' negative
    assert df.loc[3, "value"] < 0.0
    assert df.loc[4, "value"] < 0.0


@skip_no_archive
def test_load_config_2_covariance_present_and_symmetric():
    data = hepdata_ingest.load_config(
        ARCHIVE, config_index=2, observables=("P5p",), correlation_kind="total"
    )
    cov = data["covariance"]
    n = len(data["observables"])
    assert cov.shape == (n, n)
    assert np.allclose(cov, cov.T, atol=1e-12)
    # PD check (Cholesky must succeed up to the tiny ridge applied internally).
    np.linalg.cholesky(cov)
    # Diagonal entries match per-row stat^2 + syst^2 (within 1%).
    df = data["observables"]
    sigma2_expected = df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2
    np.testing.assert_allclose(np.diag(cov), sigma2_expected, rtol=1e-2)


@skip_no_archive
def test_load_config_2_metadata_declares_covariance_present():
    data = hepdata_ingest.load_config(
        ARCHIVE, config_index=2, observables=("P5p",), correlation_kind="total"
    )
    assert data["metadata"]["covariance_status"] == "present"
    assert data["metadata"]["hepdata_record"] == "ins3094698"


@skip_no_archive
def test_load_config_with_correlation_none_omits_covariance():
    data = hepdata_ingest.load_config(
        ARCHIVE, config_index=2, observables=("P5p",), correlation_kind="none"
    )
    assert "covariance" not in data
    assert data["metadata"]["covariance_status"] == "absent_declared"


@skip_no_archive
def test_public_dataframe_strips_internal_helpers():
    data = hepdata_ingest.load_config(ARCHIVE, config_index=2, observables=("P5p",))
    pub = hepdata_ingest.public_dataframe(data)
    assert not any(c.startswith("_hep_") for c in pub.columns)
    # And the schema validator accepts the public df + metadata together.
    pub_data = {"observables": pub, "metadata": dict(data["metadata"])}
    pub_data["metadata"]["covariance_status"] = "present"
    data_ingest.validate_observable_schema(pub_data)


@skip_no_archive
def test_invalid_config_index_rejected():
    with pytest.raises(ValueError, match="config_index"):
        hepdata_ingest.load_config(ARCHIVE, config_index=99, observables=("P5p",))


@skip_no_archive
def test_unknown_observable_rejected():
    with pytest.raises(KeyError, match="not present"):
        hepdata_ingest.load_config(
            ARCHIVE, config_index=2, observables=("not_a_real_observable",)
        )
