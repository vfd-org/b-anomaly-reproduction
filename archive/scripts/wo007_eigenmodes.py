"""WO-007 — Dirichlet eigenmode kernel test.

Derivation: see reports/wo007_eigenvalue_derivation.md.

Models: even-parity Dirichlet eigenfunctions of the segment Laplacian on
x in [-x_max, +x_max], where x is the dimensionless closure coordinate
(q^2 - midpoint)/Delta_psi. The empirical kappa_exponential = exp(-|x|/phi)
is identified as the massive Green's function of the free closure operator;
this script tests its truncated-spectral siblings.

    DIRICHLET_M1  k=1  A1*cos(pi x / (2 L))
    DIRICHLET_M2  k=2  A1*cos(pi x / (2 L)) + A3*cos(3 pi x / (2 L))
    DIRICHLET_M3  k=3  + A5*cos(5 pi x / (2 L))

Compare against:
    SM            k=0
    FREE_C9       k=1
    KAPPA_EXP     k=1  (WO-004A champion, A * exp(-|x|/phi))

Run on real LHCb config-2 P5' data, single-observable. Multi-observable
extension belongs to a later WO once SM slopes for the full angular
P-basis are wired (P4', P1, P2, P8' ...).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from . import hepdata_ingest, vfd_closure
from .constants import C9_SM, PROVENANCE_VFD
from .likelihood import aic, bic, chi2 as chi2_fn
from .sm_baseline import predict_vector


def dirichlet_eigenmode(x: np.ndarray, n: int, *, x_max: float = vfd_closure.KAPPA_X_MAX) -> np.ndarray:
    """Return the n-th even-parity Dirichlet eigenfunction on [-x_max, x_max].

    For Dirichlet BCs psi(+/-L) = 0 with L = x_max, the symmetric (even-parity)
    eigenfunctions are psi_k(x) = cos((2k-1) pi x / (2 L)). These vanish at the
    kinematic boundaries and peak (k=1) or oscillate symmetrically (k>=2)
    inside the closure window.

    n is the principal index (1, 2, 3, ...) so n=1 gives the lowest mode.
    """
    if n < 1:
        raise ValueError(f"eigenmode index must be >=1, got {n}")
    L = float(x_max)
    return np.cos((2 * n - 1) * np.pi * x / (2.0 * L))


def _bin_axis(df: pd.DataFrame) -> tuple[list[str], list[int], np.ndarray]:
    obs = df["observable"].tolist()
    seen: dict[str, list[float]] = {}
    for _, r in df.iterrows():
        seen.setdefault(r["observable"], [])
        if float(r["q2_lo"]) not in seen[r["observable"]]:
            seen[r["observable"]].append(float(r["q2_lo"]))
    bin_indices = [seen[r["observable"]].index(float(r["q2_lo"])) for _, r in df.iterrows()]
    q2_centres = np.array(
        [0.5 * (float(r["q2_lo"]) + float(r["q2_hi"])) for _, r in df.iterrows()],
        dtype=float,
    )
    return obs, bin_indices, q2_centres


def _chi2(values: np.ndarray, pred: np.ndarray, data: dict[str, Any], errors: np.ndarray) -> float:
    cov = data.get("covariance")
    if cov is not None:
        return chi2_fn(values, pred, covariance=cov)
    return chi2_fn(values, pred, errors=errors)


@dataclass
class FitRow:
    model: str
    role: str
    k_params: int
    chi2: float
    aic: float
    bic: float
    delta_aic_vs_FREE_C9: float
    delta_aic_vs_KAPPA_EXP: float
    coefficients: dict[str, float] = field(default_factory=dict)
    notes: str = ""


def _make_fit_row(
    model: str,
    role: str,
    chi2_val: float,
    k: int,
    n: int,
    coefficients: dict[str, float] | None = None,
    notes: str = "",
) -> FitRow:
    return FitRow(
        model=model,
        role=role,
        k_params=k,
        chi2=chi2_val,
        aic=aic(chi2_val, k),
        bic=bic(chi2_val, k, n),
        delta_aic_vs_FREE_C9=np.nan,
        delta_aic_vs_KAPPA_EXP=np.nan,
        coefficients=coefficients or {},
        notes=notes,
    )


def _fit_sm(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    pred = predict_vector(obs, bins, C9_SM)
    c2 = _chi2(values, pred, data, errors)
    return _make_fit_row("SM", "baseline", c2, 0, len(values))


def _fit_free_c9(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)

    def loss(theta):
        return _chi2(values, predict_vector(obs, bins, C9_SM + float(theta[0])), data, errors)

    r = minimize(loss, x0=[-0.5], method="Powell", bounds=[(-3.0, 3.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    return _make_fit_row(
        "FREE_C9", "reference", float(r.fun), 1, len(values),
        coefficients={"delta_c9": float(r.x[0])},
    )


def _fit_kappa_exp(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    kappa = vfd_closure.kappa_shape(q2, mode=vfd_closure.MODE_KAPPA_EXPONENTIAL)

    def loss(theta):
        a = float(theta[0])
        delta = -a * kappa
        return _chi2(values, predict_vector(obs, bins, C9_SM + delta), data, errors)

    r = minimize(loss, x0=[0.5], method="Powell", bounds=[(0.0, 5.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    return _make_fit_row(
        "KAPPA_EXP", "wo004a-champion", float(r.fun), 1, len(values),
        coefficients={"amplitude": float(r.x[0])},
        notes="A * exp(-|x|/phi)",
    )


def _fit_dirichlet_modes(data: dict[str, Any], n_modes: int) -> FitRow:
    """Fit a free linear combination of the lowest n_modes even-parity Dirichlet
    eigenfunctions: residual(x) = sum_{k=1..n} c_k * cos((2k-1) pi x / (2 L)),
    Delta_C9(x) = -residual(x). The amplitudes c_k are unconstrained sign — the
    multi-mode case can have either sign for higher modes. k_params = n_modes.
    """
    df = data["observables"]
    obs, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    x = vfd_closure.kappa_coordinate(q2)
    basis = np.stack([dirichlet_eigenmode(x, n=k) for k in range(1, n_modes + 1)], axis=1)

    # The lowest mode has a ground-state physical sign-prior (positive amplitude
    # so Delta_C9 < 0 to match literature). Higher modes are free in sign.
    bounds = [(0.0, 5.0)] + [(-5.0, 5.0)] * (n_modes - 1)
    x0 = np.concatenate([[0.5], np.zeros(n_modes - 1)])

    def loss(theta):
        c = np.asarray(theta, dtype=float)
        residual = basis @ c
        delta = -residual
        return _chi2(values, predict_vector(obs, bins, C9_SM + delta), data, errors)

    r = minimize(loss, x0=x0, method="Powell", bounds=bounds,
                 options={"xtol": 1e-7, "ftol": 1e-9, "maxiter": 5000})
    coeffs = np.atleast_1d(r.x).astype(float)
    role = "primary" if n_modes <= 2 else "diagnostic"
    return _make_fit_row(
        f"DIRICHLET_M{n_modes}",
        role,
        float(r.fun),
        n_modes,
        len(values),
        coefficients={f"c{2*k - 1}": float(coeffs[k - 1]) for k in range(1, n_modes + 1)},
        notes="even-parity Dirichlet on [-x_max, x_max]; psi_k = cos((2k-1) pi x / (2 L))",
    )


def run(
    *,
    archive_dir: Path | str = "data/raw/hepdata/extracted",
    config_index: int = 2,
    observable: str = "P5p",
    output_dir: Path | str = "reports",
) -> dict[str, Any]:
    archive = hepdata_ingest.hepdata_archive_dir(archive_dir)
    data = hepdata_ingest.load_config(
        archive, config_index=config_index, observables=(observable,)
    )

    rows: list[FitRow] = [
        _fit_sm(data),
        _fit_free_c9(data),
        _fit_kappa_exp(data),
        _fit_dirichlet_modes(data, n_modes=1),
        _fit_dirichlet_modes(data, n_modes=2),
        _fit_dirichlet_modes(data, n_modes=3),
    ]
    free_aic = next(r.aic for r in rows if r.model == "FREE_C9")
    exp_aic = next(r.aic for r in rows if r.model == "KAPPA_EXP")
    for r in rows:
        r.delta_aic_vs_FREE_C9 = r.aic - free_aic
        r.delta_aic_vs_KAPPA_EXP = r.aic - exp_aic

    df_main = pd.DataFrame([r.__dict__ for r in rows])

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df_main.to_csv(out / "wo007_eigenmodes.csv", index=False)
    df_main.to_json(out / "wo007_eigenmodes.json", orient="records", indent=2)

    return {
        "data": data,
        "models": df_main,
        "observable": observable,
        "config_index": config_index,
        "x_max": vfd_closure.KAPPA_X_MAX,
        "provenance": PROVENANCE_VFD,
    }


def main() -> None:
    res = run()
    print(res["models"].to_string(index=False))


if __name__ == "__main__":
    main()
