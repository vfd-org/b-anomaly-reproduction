from __future__ import annotations

import numpy as np
import pytest

from vfd_b_anomaly.likelihood import aic, bic, chi2, reduced_chi2


def test_chi2_zero_when_perfect_fit():
    obs = np.array([1.0, 2.0, 3.0])
    pred = obs.copy()
    err = np.array([0.1, 0.1, 0.1])
    assert chi2(obs, pred, errors=err) == 0.0


def test_chi2_diagonal_known_value():
    obs = np.array([1.0, 2.0])
    pred = np.array([1.5, 2.0])
    err = np.array([0.5, 1.0])
    # ((1 - 1.5)/0.5)^2 + ((2 - 2)/1.0)^2 = 1.0
    assert chi2(obs, pred, errors=err) == pytest.approx(1.0)


def test_chi2_full_covariance_matches_diagonal_when_uncorrelated():
    obs = np.array([0.0, 0.0])
    pred = np.array([1.0, 2.0])
    err = np.array([1.0, 0.5])
    cov = np.diag(err**2)
    diag_val = chi2(obs, pred, errors=err)
    cov_val = chi2(obs, pred, covariance=cov)
    assert diag_val == pytest.approx(cov_val)


def test_chi2_requires_exactly_one_of_cov_or_errors():
    obs = np.array([0.0])
    pred = np.array([0.0])
    err = np.array([1.0])
    with pytest.raises(ValueError, match="exactly one"):
        chi2(obs, pred, errors=err, covariance=np.eye(1))
    with pytest.raises(ValueError, match="exactly one"):
        chi2(obs, pred)


def test_chi2_rejects_nonpositive_errors():
    obs = np.array([0.0])
    pred = np.array([0.0])
    with pytest.raises(ValueError, match="positive"):
        chi2(obs, pred, errors=np.array([0.0]))


def test_chi2_singular_covariance_raises():
    obs = np.array([0.0, 0.0])
    pred = np.array([0.0, 0.0])
    cov = np.array([[1.0, 1.0], [1.0, 1.0]])
    with pytest.raises(ValueError, match="singular"):
        chi2(obs, pred, covariance=cov)


def test_aic_bic_penalise_more_params():
    c2 = 10.0
    a0 = aic(c2, 0)
    a3 = aic(c2, 3)
    assert a3 > a0  # AIC penalty grows with k
    b0 = bic(c2, 0, n_data=20)
    b3 = bic(c2, 3, n_data=20)
    assert b3 > b0
    # BIC penalises more harshly than AIC for n>=8
    assert (b3 - b0) > (a3 - a0)


def test_reduced_chi2_dof_guard():
    with pytest.raises(ValueError, match="non-positive dof"):
        reduced_chi2(10.0, n_data=2, n_params=2)
