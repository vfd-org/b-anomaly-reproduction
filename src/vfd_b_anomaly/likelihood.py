"""Gaussian chi-squared likelihood for binned observables.

Supports diagonal-error and full-covariance modes. Missing covariance is
treated as 'diagonal only' but logged via the explicit `covariance` argument
being None - never silently fabricated.
"""

from __future__ import annotations

import numpy as np


def chi2(
    observed: np.ndarray,
    predicted: np.ndarray,
    covariance: np.ndarray | None = None,
    errors: np.ndarray | None = None,
) -> float:
    """Compute chi^2.

    Exactly one of `covariance` or `errors` must be supplied.

    Parameters
    ----------
    observed, predicted : (n,) arrays of observable values.
    covariance : (n, n) symmetric positive-definite matrix or None.
    errors : (n,) array of 1-sigma diagonal errors or None.

    Returns
    -------
    chi^2 value (float).
    """
    obs = np.asarray(observed, dtype=float)
    pred = np.asarray(predicted, dtype=float)
    if obs.shape != pred.shape:
        raise ValueError(f"shape mismatch: observed {obs.shape} vs predicted {pred.shape}")
    if obs.ndim != 1:
        raise ValueError("observed must be 1-D")

    if (covariance is None) == (errors is None):
        raise ValueError("supply exactly one of `covariance` or `errors`")

    residual = obs - pred

    if errors is not None:
        err = np.asarray(errors, dtype=float)
        if err.shape != obs.shape:
            raise ValueError("errors shape mismatch")
        if np.any(err <= 0):
            raise ValueError("all errors must be strictly positive")
        return float(np.sum((residual / err) ** 2))

    cov = np.asarray(covariance, dtype=float)
    if cov.shape != (obs.size, obs.size):
        raise ValueError(f"covariance shape {cov.shape} != ({obs.size}, {obs.size})")
    # Symmetrise to suppress floating drift.
    cov_sym = 0.5 * (cov + cov.T)
    try:
        solved = np.linalg.solve(cov_sym, residual)
    except np.linalg.LinAlgError as e:
        raise ValueError(f"covariance is singular: {e}") from e
    return float(residual @ solved)


def reduced_chi2(chi2_value: float, n_data: int, n_params: int) -> float:
    """Reduced chi^2 with explicit DOF guard."""
    dof = n_data - n_params
    if dof <= 0:
        raise ValueError(f"non-positive dof: n_data={n_data}, n_params={n_params}")
    return chi2_value / dof


def aic(chi2_value: float, n_params: int) -> float:
    """Akaike information criterion under Gaussian likelihood (constant terms dropped)."""
    return chi2_value + 2.0 * n_params


def bic(chi2_value: float, n_params: int, n_data: int) -> float:
    """Bayesian information criterion under Gaussian likelihood (constants dropped)."""
    if n_data <= 0:
        raise ValueError("n_data must be positive")
    return chi2_value + n_params * np.log(n_data)
