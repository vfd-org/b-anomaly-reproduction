"""VFD closure-residual model for ΔC9_VFD(q^2) = -A * R(q^2).

Three constrained residual shapes are exposed:

    R0 'constant'           : R(q^2) = 1                           (1 free parameter: A)
    R1 'log_phase'          : R(q^2) = cos(phi * log(q^2/q^2_ref) + theta)
                              (2 or 3 free parameters: A, theta, optional phi)
    R2 'threshold_weighted' : R(q^2) = sum_k w_k * exp(-(q^2 - q^2_k)^2 / (2 sigma_k^2))
                              (2-3 free parameters: A, sigma, optional centre offset)

Parameter counts are FIXED by mode and do NOT scale with the q^2 grid.
The tests in test_vfd_closure.py enforce this property so per-bin freedom
cannot leak in.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.optimize import minimize

from .constants import (
    C9_SM,
    J_PSI_Q2,
    PHI,
    PROVENANCE_VFD,
    PSI2S_Q2,
    Q2_REF_GEV2,
)
from .likelihood import aic, bic, chi2
from .sm_baseline import predict_vector

# Public mode constants.
MODE_CONSTANT = "constant"
MODE_LOG_PHASE = "log_phase"
MODE_THRESHOLD = "threshold_weighted"
# WO-003 fixed-geometry threshold modes: sigma is FIXED from physics, not fitted.
# Three candidate sigma anchors are committed BEFORE looking at the fit result;
# all three are run and their AIC/BIC are reported so no post-hoc cherry-picking
# is possible.
MODE_THRESHOLD_FIXED_HALFSEP = "threshold_fixed_halfsep"  # sigma = (m^2_psi(2S) - m^2_J/psi)/2
MODE_THRESHOLD_FIXED_PHISEP = "threshold_fixed_phisep"    # sigma = phi * halfsep
MODE_THRESHOLD_FIXED_VETO = "threshold_fixed_veto"        # sigma = avg charmonium-veto half-width
# WO-004A dimensionless closure coordinate kappa(q^2). Single resonance-centred
# coordinate x = (q^2 - midpoint) / Delta_psi where:
#   midpoint = (m^2_J/psi + m^2_psi(2S)) / 2
#   Delta_psi = m^2_psi(2S) - m^2_J/psi
# Three shape families, all FROZEN before fitting. Only amplitude A is fitted.
MODE_KAPPA_GAUSSIAN = "kappa_gaussian"        # exp(-x^2 / phi^2)
MODE_KAPPA_EXPONENTIAL = "kappa_exponential"  # exp(-|x| / phi)
MODE_KAPPA_LORENTZIAN = "kappa_lorentzian"    # 1 / (1 + phi^2 * x^2)
# WO-004A iter-2 (Codex pair-programming derivation 2026-04-29):
# Smooth-core massive closure kernel. Removes the cusp at midpoint of the
# pure exponential while preserving evenness, boundedness, exp tail, and the
# phi-fixed length scale. Derivation: gapped closure operator green-function
# with one regularising derivative. No new fitted scale.
MODE_KAPPA_CANONICAL = "kappa_canonical"      # (1 + |x|/phi) * exp(-|x|/phi)
# Two further single-amplitude variants in the same one-derivative-correction
# family, kept for cross-check; all committed in code BEFORE looking at fit.
MODE_KAPPA_YUKAWA_SMOOTH = "kappa_yukawa_smooth"      # exp(-|x|/phi) / (1 + |x|/phi)
MODE_KAPPA_SECH = "kappa_sech"                        # sech(|x|/phi) = 2 / (e^{|x|/phi} + e^{-|x|/phi})
# WO-005 reflective / cavity-mode kappa kernels. Single amplitude parameter.
# These treat the closure window as a bounded cavity rather than a single
# decaying mode, motivated by the apparent twin tension at low- and high-q^2.
# All shapes are normalised to peak value 1 within the kinematic window.
MODE_KAPPA_REFLECTIVE_CENTRE = "kappa_reflective_centre"  # cosh((x_max - |x|)/phi) / cosh(x_max/phi)
MODE_KAPPA_REFLECTIVE_EDGE = "kappa_reflective_edge"      # sinh(|x|/phi) / sinh(x_max/phi)
MODE_KAPPA_REFLECTIVE_CAVITY = "kappa_reflective_cavity"  # sin(pi |x| / (2 x_max)) standing mode
ALL_MODES = (
    MODE_CONSTANT,
    MODE_LOG_PHASE,
    MODE_THRESHOLD,
    MODE_THRESHOLD_FIXED_HALFSEP,
    MODE_THRESHOLD_FIXED_PHISEP,
    MODE_THRESHOLD_FIXED_VETO,
    MODE_KAPPA_GAUSSIAN,
    MODE_KAPPA_EXPONENTIAL,
    MODE_KAPPA_LORENTZIAN,
    MODE_KAPPA_CANONICAL,
    MODE_KAPPA_YUKAWA_SMOOTH,
    MODE_KAPPA_SECH,
    MODE_KAPPA_REFLECTIVE_CENTRE,
    MODE_KAPPA_REFLECTIVE_EDGE,
    MODE_KAPPA_REFLECTIVE_CAVITY,
)
KAPPA_MODES = (
    MODE_KAPPA_GAUSSIAN,
    MODE_KAPPA_EXPONENTIAL,
    MODE_KAPPA_LORENTZIAN,
    MODE_KAPPA_CANONICAL,
    MODE_KAPPA_YUKAWA_SMOOTH,
    MODE_KAPPA_SECH,
    MODE_KAPPA_REFLECTIVE_CENTRE,
    MODE_KAPPA_REFLECTIVE_EDGE,
    MODE_KAPPA_REFLECTIVE_CAVITY,
)

# Resonance midpoint and separation in q^2 (GeV^2). FROZEN constants, not fit.
KAPPA_MIDPOINT_GEV2 = 0.5 * (J_PSI_Q2 + PSI2S_Q2)        # = 11.59 GeV^2
KAPPA_DELTA_PSI_GEV2 = PSI2S_Q2 - J_PSI_Q2               # = 4.00 GeV^2

# WO-005 cavity boundary in dimensionless x = (q^2 - midpoint)/Delta_psi.
# Bounds derived from B -> K* mu mu kinematic limits, NOT fit:
#   q^2_min = 4 m_mu^2 (di-muon threshold)
#   q^2_max = (m_B - m_K*)^2 (zero-recoil endpoint)
KINEMATIC_Q2_MIN_GEV2 = 4.0 * (0.10566) ** 2             # ~ 0.0447 GeV^2
KINEMATIC_Q2_MAX_GEV2 = (5.27966 - 0.89555) ** 2         # ~ 19.21 GeV^2
KAPPA_X_MAX = max(
    abs(KINEMATIC_Q2_MIN_GEV2 - KAPPA_MIDPOINT_GEV2),
    abs(KINEMATIC_Q2_MAX_GEV2 - KAPPA_MIDPOINT_GEV2),
) / KAPPA_DELTA_PSI_GEV2                                  # ~ 2.886

# Geometry constants for the fixed-sigma variants. All values are derived from
# the standard charmonium-resonance and LHCb-veto topology, NOT from any P5'
# fit, so substituting these does not constitute post-hoc tuning.
_HALF_SEPARATION_GEV2 = 0.5 * (PSI2S_Q2 - J_PSI_Q2)  # = 2.0 GeV^2
_PHI_SCALED_HALFSEP_GEV2 = PHI * _HALF_SEPARATION_GEV2  # = 3.236 GeV^2
# Charmonium veto windows: [8.0, 11.0] (half-width 1.5) and [12.5, 15.0] (half-width 1.25).
_VETO_HALFWIDTH_GEV2 = 0.5 * (1.5 + 1.25)  # = 1.375 GeV^2

FIXED_SIGMA_BY_MODE: dict[str, float] = {
    MODE_THRESHOLD_FIXED_HALFSEP: float(_HALF_SEPARATION_GEV2),
    MODE_THRESHOLD_FIXED_PHISEP: float(_PHI_SCALED_HALFSEP_GEV2),
    MODE_THRESHOLD_FIXED_VETO: float(_VETO_HALFWIDTH_GEV2),
}

# Number of free parameters per mode. The KEY discipline of the VFD model.
PARAM_COUNT: dict[str, int] = {
    MODE_CONSTANT: 1,                   # (A,)
    MODE_LOG_PHASE: 2,                  # (A, theta); phi is fixed at PHI by default
    MODE_THRESHOLD: 2,                  # (A, sigma); centres fixed at J/psi and psi(2S)
    MODE_THRESHOLD_FIXED_HALFSEP: 1,    # (A,); sigma fixed from physics
    MODE_THRESHOLD_FIXED_PHISEP: 1,     # (A,); sigma fixed from physics
    MODE_THRESHOLD_FIXED_VETO: 1,       # (A,); sigma fixed from physics
    MODE_KAPPA_GAUSSIAN: 1,             # (A,); kappa(q^2) frozen
    MODE_KAPPA_EXPONENTIAL: 1,          # (A,); kappa(q^2) frozen
    MODE_KAPPA_LORENTZIAN: 1,           # (A,); kappa(q^2) frozen
    MODE_KAPPA_CANONICAL: 1,            # (A,); smooth-core kernel (1+|x|/phi)*exp(-|x|/phi)
    MODE_KAPPA_YUKAWA_SMOOTH: 1,        # (A,); regularised Yukawa exp(-|x|/phi)/(1+|x|/phi)
    MODE_KAPPA_SECH: 1,                 # (A,); sech(|x|/phi)
    MODE_KAPPA_REFLECTIVE_CENTRE: 1,    # (A,); WO-005 cosh((x_max-|x|)/phi)/cosh(x_max/phi)
    MODE_KAPPA_REFLECTIVE_EDGE: 1,      # (A,); WO-005 sinh(|x|/phi)/sinh(x_max/phi)
    MODE_KAPPA_REFLECTIVE_CAVITY: 1,    # (A,); WO-005 sin(pi |x| / (2 x_max))
}

# Parameter bounds. Bounded amplitude, sign-priored (A >= 0 so that
# ΔC9 = -A * R is negative when R > 0 - this is the physical expectation
# from the literature anomaly direction).
_BOUNDS: dict[str, list[tuple[float, float]]] = {
    MODE_CONSTANT: [(0.0, 3.0)],
    MODE_LOG_PHASE: [(0.0, 3.0), (-np.pi, np.pi)],
    MODE_THRESHOLD: [(0.0, 3.0), (0.5, 6.0)],
    MODE_THRESHOLD_FIXED_HALFSEP: [(0.0, 3.0)],
    MODE_THRESHOLD_FIXED_PHISEP: [(0.0, 3.0)],
    MODE_THRESHOLD_FIXED_VETO: [(0.0, 3.0)],
    MODE_KAPPA_GAUSSIAN: [(0.0, 3.0)],
    MODE_KAPPA_EXPONENTIAL: [(0.0, 3.0)],
    MODE_KAPPA_LORENTZIAN: [(0.0, 3.0)],
    MODE_KAPPA_CANONICAL: [(0.0, 3.0)],
    MODE_KAPPA_YUKAWA_SMOOTH: [(0.0, 3.0)],
    MODE_KAPPA_SECH: [(0.0, 3.0)],
    MODE_KAPPA_REFLECTIVE_CENTRE: [(0.0, 3.0)],
    MODE_KAPPA_REFLECTIVE_EDGE: [(0.0, 3.0)],
    MODE_KAPPA_REFLECTIVE_CAVITY: [(0.0, 3.0)],
}


def kappa_coordinate(q2_mid: np.ndarray) -> np.ndarray:
    """Dimensionless closure coordinate x(q^2) = (q^2 - midpoint) / Delta_psi.

    Frozen constants:
      midpoint = (m^2_J/psi + m^2_psi(2S)) / 2
      Delta_psi = m^2_psi(2S) - m^2_J/psi
    """
    q2 = np.asarray(q2_mid, dtype=float)
    return (q2 - KAPPA_MIDPOINT_GEV2) / KAPPA_DELTA_PSI_GEV2


def kappa_shape(q2_mid: np.ndarray, *, mode: str) -> np.ndarray:
    """Frozen kappa(q^2) shape, peaked at the resonance midpoint and bounded in [0, 1]."""
    if mode not in KAPPA_MODES:
        raise ValueError(f"kappa_shape: mode {mode!r} not in {KAPPA_MODES}")
    x = kappa_coordinate(q2_mid)
    if mode == MODE_KAPPA_GAUSSIAN:
        return np.exp(-(x ** 2) / (PHI ** 2))
    if mode == MODE_KAPPA_EXPONENTIAL:
        return np.exp(-np.abs(x) / PHI)
    if mode == MODE_KAPPA_LORENTZIAN:
        return 1.0 / (1.0 + (PHI ** 2) * (x ** 2))
    if mode == MODE_KAPPA_CANONICAL:
        r = np.abs(x) / PHI
        return (1.0 + r) * np.exp(-r)
    if mode == MODE_KAPPA_YUKAWA_SMOOTH:
        r = np.abs(x) / PHI
        return np.exp(-r) / (1.0 + r)
    if mode == MODE_KAPPA_SECH:
        r = np.abs(x) / PHI
        return 1.0 / np.cosh(r)
    if mode == MODE_KAPPA_REFLECTIVE_CENTRE:
        # cosh((x_max - |x|)/phi) / cosh(x_max/phi). Peaks at midpoint, valleys at edges.
        ax = np.abs(x)
        return np.cosh((KAPPA_X_MAX - ax) / PHI) / np.cosh(KAPPA_X_MAX / PHI)
    if mode == MODE_KAPPA_REFLECTIVE_EDGE:
        # sinh(|x|/phi) / sinh(x_max/phi). Zero at midpoint, peak at edges.
        ax = np.abs(x)
        return np.sinh(ax / PHI) / np.sinh(KAPPA_X_MAX / PHI)
    if mode == MODE_KAPPA_REFLECTIVE_CAVITY:
        # sin(pi |x| / (2 x_max)). Lowest cavity standing mode with antinodes at +/- x_max.
        ax = np.minimum(np.abs(x), KAPPA_X_MAX)  # clip to avoid >1 outside the formal domain
        return np.sin(np.pi * ax / (2.0 * KAPPA_X_MAX))
    raise AssertionError("unreachable")


@dataclass
class ClosureFitResult:
    mode: str
    n_data: int
    n_params: int
    chi2: float
    aic: float
    bic: float
    params: dict[str, float]
    effective_delta_c9_mean: float
    effective_delta_c9_grid: np.ndarray
    success: bool
    message: str
    provenance: str = PROVENANCE_VFD
    extra: dict[str, Any] = field(default_factory=dict)


def closure_residual_q2(
    q2_mid: np.ndarray,
    amplitude: float,
    *,
    mode: str = MODE_LOG_PHASE,
    theta: float = 0.0,
    phase_scale: float = PHI,
    sigma: float = 2.0,
) -> np.ndarray:
    """Compute ΔC9_VFD(q^2) = -A * R(q^2) on a q^2 grid.

    Parameters
    ----------
    q2_mid : (n,) array of q^2 bin centres in GeV^2.
    amplitude : non-negative scalar A.
    mode : 'constant', 'log_phase' or 'threshold_weighted'.
    theta : phase offset (log_phase only).
    phase_scale : multiplier of log(q^2/q^2_ref). Defaults to golden ratio.
    sigma : Gaussian width in GeV^2 (threshold mode only).
    """
    if mode not in ALL_MODES:
        raise ValueError(f"unknown VFD mode {mode!r}; choose from {ALL_MODES}")
    if amplitude < 0:
        raise ValueError("amplitude must be >= 0 (sign prior on ΔC9 = -A R)")

    q2 = np.asarray(q2_mid, dtype=float)
    if np.any(q2 <= 0):
        raise ValueError("q2_mid must be strictly positive")

    if mode == MODE_CONSTANT:
        residual = np.ones_like(q2)

    elif mode == MODE_LOG_PHASE:
        residual = np.cos(phase_scale * np.log(q2 / Q2_REF_GEV2) + theta)

    elif mode == MODE_THRESHOLD:
        if sigma <= 0:
            raise ValueError("sigma must be positive")
        b1 = np.exp(-((q2 - J_PSI_Q2) ** 2) / (2.0 * sigma**2))
        b2 = np.exp(-((q2 - PSI2S_Q2) ** 2) / (2.0 * sigma**2))
        residual = 0.5 * (b1 + b2)

    elif mode in FIXED_SIGMA_BY_MODE:
        sigma_fixed = FIXED_SIGMA_BY_MODE[mode]
        b1 = np.exp(-((q2 - J_PSI_Q2) ** 2) / (2.0 * sigma_fixed**2))
        b2 = np.exp(-((q2 - PSI2S_Q2) ** 2) / (2.0 * sigma_fixed**2))
        residual = 0.5 * (b1 + b2)

    elif mode in KAPPA_MODES:
        residual = kappa_shape(q2, mode=mode)

    else:
        raise ValueError(f"unhandled mode {mode!r}")

    return -amplitude * residual


def _bin_centers(data: dict[str, Any]) -> tuple[list[str], list[int], np.ndarray, np.ndarray]:
    """Return (observables, bin_indices, q2_centres, residual_grid_centres).

    The residual grid centres are per-row (one value per data point), since
    each data point sits in one q^2 bin. We do NOT collapse to unique bins
    here: it keeps the per-row vector aligned with the prediction vector.
    """
    df = data["observables"]
    if df["value"].isna().any():
        raise ValueError("observable values contain NaN; cannot fit VFD closure on placeholder data")
    obs = df["observable"].tolist()
    seen_lo: dict[str, list[float]] = {}
    for _, row in df.iterrows():
        seen_lo.setdefault(row["observable"], [])
        seen_lo[row["observable"]].append(float(row["q2_lo"]))
    for o in seen_lo:
        seen_lo[o] = sorted(set(seen_lo[o]))
    bin_indices: list[int] = []
    centres: list[float] = []
    for _, row in df.iterrows():
        bin_indices.append(seen_lo[row["observable"]].index(float(row["q2_lo"])))
        centres.append(0.5 * (float(row["q2_lo"]) + float(row["q2_hi"])))
    return obs, bin_indices, np.asarray(centres, dtype=float), np.asarray(centres, dtype=float)


def fit_vfd_closure(
    data: dict[str, Any],
    *,
    mode: str = MODE_LOG_PHASE,
    smoothness_lambda: float = 0.0,
) -> ClosureFitResult:
    """Fit the VFD closure-residual model to observed data.

    `smoothness_lambda` penalises mode-specific shape excursion (e.g. large
    amplitudes). It is NOT a per-bin penalty and cannot be used to fake
    flexibility.
    """
    if mode not in ALL_MODES:
        raise ValueError(f"unknown VFD mode {mode!r}")

    df = data["observables"]
    observables, bin_indices, q2_centres, _ = _bin_centers(data)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy(dtype=float) ** 2 + df["syst_err"].to_numpy(dtype=float) ** 2)
    covariance = data.get("covariance")
    n_data = len(values)
    n_params = PARAM_COUNT[mode]

    def _delta_c9_array(theta: np.ndarray) -> np.ndarray:
        if mode == MODE_CONSTANT:
            (a,) = theta
            return closure_residual_q2(q2_centres, a, mode=mode)
        if mode == MODE_LOG_PHASE:
            a, th = theta
            return closure_residual_q2(q2_centres, a, mode=mode, theta=float(th))
        if mode == MODE_THRESHOLD:
            a, sig = theta
            return closure_residual_q2(q2_centres, a, mode=mode, sigma=float(sig))
        if mode in FIXED_SIGMA_BY_MODE:
            (a,) = theta
            return closure_residual_q2(q2_centres, a, mode=mode)
        if mode in KAPPA_MODES:
            (a,) = theta
            return closure_residual_q2(q2_centres, a, mode=mode)
        raise AssertionError("unreachable")

    def neg_log_l(theta: np.ndarray) -> float:
        delta_arr = _delta_c9_array(theta)
        pred = predict_vector(observables, bin_indices, C9_SM + delta_arr)
        c2 = _chi2(values, pred, covariance, errors)
        if smoothness_lambda > 0:
            c2 += smoothness_lambda * float(theta[0]) ** 2  # mild amplitude regulariser only
        return c2

    x0 = _seed(mode)
    bounds = _BOUNDS[mode]
    result = minimize(
        neg_log_l,
        x0=x0,
        method="Powell",
        bounds=bounds,
        options={"xtol": 1e-6, "ftol": 1e-8, "maxiter": 5000},
    )
    theta_hat = np.atleast_1d(result.x).astype(float)
    delta_arr = _delta_c9_array(theta_hat)

    # Recompute pure (un-regularised) chi^2 for AIC/BIC reporting.
    pred = predict_vector(observables, bin_indices, C9_SM + delta_arr)
    c2_pure = _chi2(values, pred, covariance, errors)

    params_named = _name_params(mode, theta_hat)
    return ClosureFitResult(
        mode=mode,
        n_data=n_data,
        n_params=n_params,
        chi2=c2_pure,
        aic=aic(c2_pure, n_params),
        bic=bic(c2_pure, n_params, n_data),
        params=params_named,
        effective_delta_c9_mean=float(np.mean(delta_arr)),
        effective_delta_c9_grid=delta_arr,
        success=bool(result.success),
        message=str(result.message),
    )


def _seed(mode: str) -> np.ndarray:
    if mode == MODE_CONSTANT:
        return np.array([0.5])
    if mode == MODE_LOG_PHASE:
        return np.array([0.5, 0.0])
    if mode == MODE_THRESHOLD:
        return np.array([0.5, 2.0])
    if mode in FIXED_SIGMA_BY_MODE:
        return np.array([0.5])
    if mode in KAPPA_MODES:
        return np.array([0.5])
    raise ValueError(f"unknown mode {mode!r}")


def _name_params(mode: str, theta: np.ndarray) -> dict[str, float]:
    if mode == MODE_CONSTANT:
        return {"amplitude": float(theta[0])}
    if mode == MODE_LOG_PHASE:
        return {"amplitude": float(theta[0]), "theta": float(theta[1]), "phase_scale": PHI}
    if mode == MODE_THRESHOLD:
        return {"amplitude": float(theta[0]), "sigma": float(theta[1])}
    if mode in FIXED_SIGMA_BY_MODE:
        return {"amplitude": float(theta[0]), "sigma_fixed": FIXED_SIGMA_BY_MODE[mode]}
    if mode in KAPPA_MODES:
        return {
            "amplitude": float(theta[0]),
            "kappa_midpoint_gev2": KAPPA_MIDPOINT_GEV2,
            "kappa_delta_psi_gev2": KAPPA_DELTA_PSI_GEV2,
            "phi": PHI,
        }
    raise ValueError(f"unknown mode {mode!r}")


def _chi2(values, pred, covariance, errors):
    if covariance is not None:
        return chi2(values, pred, covariance=covariance)
    return chi2(values, pred, errors=errors)
