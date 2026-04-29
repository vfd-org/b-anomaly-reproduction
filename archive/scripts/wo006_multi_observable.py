"""WO-006 — multi-observable shared-kernel test.

Question: does ONE midpoint-peaked kappa shape with ONE shared amplitude A
explain multiple LHCb B0 -> K*0 mu mu observables simultaneously?

Models compared on the joint dataset:
    M0  SM                          k = 0
    M1  FREE_C9 (one global DC9)    k = 1
    M2  VFD_SHARED_KAPPA            k = 1   <- the structural-compression test
    M3  VFD_PER_OBSERVABLE_A        k = N_obs (diagnostic only)

Acceptance gates:
    - M2 AIC and BIC <= M1 AIC and BIC.
    - All major observables prefer the same SIGN of effective DC9 in M3.
    - Leave-one-observable-out remains stable for M2 (sign and magnitude).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from . import hepdata_ingest, vfd_closure, wilson_fit
from .constants import C9_SM, PROVENANCE_VFD
from .likelihood import aic, bic, chi2 as chi2_fn
from .sm_baseline import predict_vector


DEFAULT_KERNEL_MODE = vfd_closure.MODE_KAPPA_EXPONENTIAL  # WO-005 champion
DEFAULT_OBSERVABLES = ("P5p", "dBF/dq2")


@dataclass
class FitRow:
    model: str
    role: str
    k_params: int
    chi2: float
    aic: float
    bic: float
    delta_aic_vs_FREE_C9: float
    delta_bic_vs_FREE_C9: float
    effective_delta_c9_mean: float
    per_observable_delta_c9: dict[str, float] = field(default_factory=dict)
    notes: str = ""


def _bin_axis(df: pd.DataFrame) -> tuple[list[str], list[int], np.ndarray]:
    obs = df["observable"].tolist()
    seen: dict[str, list[float]] = {}
    for _, r in df.iterrows():
        seen.setdefault(r["observable"], [])
        if float(r["q2_lo"]) not in seen[r["observable"]]:
            seen[r["observable"]].append(float(r["q2_lo"]))
    bin_indices = [
        seen[r["observable"]].index(float(r["q2_lo"])) for _, r in df.iterrows()
    ]
    q2_centres = np.array(
        [0.5 * (float(r["q2_lo"]) + float(r["q2_hi"])) for _, r in df.iterrows()],
        dtype=float,
    )
    return obs, bin_indices, q2_centres


def _chi2_with_cov(values: np.ndarray, pred: np.ndarray, data: dict[str, Any], errors: np.ndarray) -> float:
    cov = data.get("covariance")
    if cov is not None:
        return chi2_fn(values, pred, covariance=cov)
    return chi2_fn(values, pred, errors=errors)


def fit_sm(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    pred = predict_vector(obs, bins, C9_SM)
    c2 = _chi2_with_cov(values, pred, data, errors)
    n = len(values)
    return FitRow(
        model="SM",
        role="baseline",
        k_params=0,
        chi2=c2,
        aic=aic(c2, 0),
        bic=bic(c2, 0, n),
        delta_aic_vs_FREE_C9=np.nan,
        delta_bic_vs_FREE_C9=np.nan,
        effective_delta_c9_mean=0.0,
    )


def fit_free_c9(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    n = len(values)

    def loss(theta: np.ndarray) -> float:
        delta = float(theta[0])
        pred = predict_vector(obs, bins, C9_SM + delta)
        return _chi2_with_cov(values, pred, data, errors)

    res = minimize(loss, x0=np.array([-0.5]), method="Powell",
                   bounds=[(-3.0, 3.0)], options={"xtol": 1e-7, "ftol": 1e-9})
    delta_hat = float(res.x[0])
    c2 = float(res.fun)
    return FitRow(
        model="FREE_C9",
        role="reference",
        k_params=1,
        chi2=c2,
        aic=aic(c2, 1),
        bic=bic(c2, 1, n),
        delta_aic_vs_FREE_C9=0.0,
        delta_bic_vs_FREE_C9=0.0,
        effective_delta_c9_mean=delta_hat,
        per_observable_delta_c9={o: delta_hat for o in set(obs)},
    )


def fit_shared_kappa(data: dict[str, Any], *, mode: str = DEFAULT_KERNEL_MODE) -> FitRow:
    df = data["observables"]
    obs, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    n = len(values)
    kappa = vfd_closure.kappa_shape(q2, mode=mode)

    def loss(theta: np.ndarray) -> float:
        a = float(theta[0])
        delta = -a * kappa
        pred = predict_vector(obs, bins, C9_SM + delta)
        return _chi2_with_cov(values, pred, data, errors)

    res = minimize(loss, x0=np.array([0.5]), method="Powell",
                   bounds=[(0.0, 5.0)], options={"xtol": 1e-7, "ftol": 1e-9})
    a_hat = float(res.x[0])
    delta_grid = -a_hat * kappa
    c2 = float(res.fun)
    eff_per_obs: dict[str, float] = {}
    for o in sorted(set(obs)):
        mask = np.array([oo == o for oo in obs])
        eff_per_obs[o] = float(np.mean(delta_grid[mask]))
    return FitRow(
        model=f"VFD_SHARED_KAPPA[{mode}]",
        role="primary",
        k_params=1,
        chi2=c2,
        aic=aic(c2, 1),
        bic=bic(c2, 1, n),
        delta_aic_vs_FREE_C9=np.nan,  # filled in by run()
        delta_bic_vs_FREE_C9=np.nan,
        effective_delta_c9_mean=float(np.mean(delta_grid)),
        per_observable_delta_c9=eff_per_obs,
        notes=f"amplitude={a_hat:.4f}",
    )


def fit_per_observable_a(data: dict[str, Any], *, mode: str = DEFAULT_KERNEL_MODE) -> FitRow:
    """Diagnostic: independent A_o per observable. k = N_obs."""
    df = data["observables"]
    obs, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    n = len(values)
    kappa = vfd_closure.kappa_shape(q2, mode=mode)

    obs_unique = sorted(set(obs))
    obs_idx = {o: i for i, o in enumerate(obs_unique)}
    obs_lookup = np.array([obs_idx[o] for o in obs])
    k = len(obs_unique)

    def loss(theta: np.ndarray) -> float:
        a_per = np.asarray(theta, dtype=float)
        delta = -a_per[obs_lookup] * kappa
        pred = predict_vector(obs, bins, C9_SM + delta)
        return _chi2_with_cov(values, pred, data, errors)

    res = minimize(
        loss,
        x0=np.full(k, 0.5),
        method="Powell",
        bounds=[(0.0, 5.0)] * k,
        options={"xtol": 1e-7, "ftol": 1e-9, "maxiter": 5000},
    )
    a_hat = np.atleast_1d(res.x).astype(float)
    delta_grid = -a_hat[obs_lookup] * kappa
    c2 = float(res.fun)
    eff_per_obs = {o: float(np.mean(delta_grid[obs_lookup == obs_idx[o]])) for o in obs_unique}
    a_per_obs = {o: float(a_hat[obs_idx[o]]) for o in obs_unique}
    return FitRow(
        model=f"VFD_PER_OBSERVABLE_A[{mode}]",
        role="diagnostic",
        k_params=k,
        chi2=c2,
        aic=aic(c2, k),
        bic=bic(c2, k, n),
        delta_aic_vs_FREE_C9=np.nan,
        delta_bic_vs_FREE_C9=np.nan,
        effective_delta_c9_mean=float(np.mean(delta_grid)),
        per_observable_delta_c9=eff_per_obs,
        notes=f"per-observable A={a_per_obs}",
    )


def leave_one_observable_out(
    data: dict[str, Any], *, mode: str = DEFAULT_KERNEL_MODE
) -> pd.DataFrame:
    df = data["observables"]
    obs_unique = sorted(set(df["observable"].tolist()))
    rows = []
    for drop in obs_unique:
        keep_mask = (df["observable"] != drop).to_numpy()
        sub_df = df[keep_mask].reset_index(drop=True)
        sub_data: dict[str, Any] = {
            "observables": sub_df,
            "metadata": data["metadata"],
        }
        if "covariance" in data and data["covariance"] is not None:
            cov = np.asarray(data["covariance"])
            keep = np.where(keep_mask)[0].tolist()
            sub_data["covariance"] = cov[np.ix_(keep, keep)]
        row = fit_shared_kappa(sub_data, mode=mode)
        rows.append(
            {
                "dropped_observable": drop,
                "n_remaining": len(sub_df),
                "amplitude": float(row.notes.split("=")[-1]),
                "effective_delta_c9_mean": row.effective_delta_c9_mean,
                "chi2": row.chi2,
                "aic": row.aic,
            }
        )
    return pd.DataFrame(rows)


def run(
    *,
    archive_dir: Path | str = "data/raw/hepdata/extracted",
    config_index: int = 2,
    observables: Iterable[str] = DEFAULT_OBSERVABLES,
    kernel_mode: str = DEFAULT_KERNEL_MODE,
    output_dir: Path | str = "reports",
) -> dict[str, Any]:
    archive = hepdata_ingest.hepdata_archive_dir(archive_dir)
    data = hepdata_ingest.load_config(
        archive, config_index=config_index, observables=tuple(observables)
    )

    rows = [
        fit_sm(data),
        fit_free_c9(data),
        fit_shared_kappa(data, mode=kernel_mode),
        fit_per_observable_a(data, mode=kernel_mode),
    ]
    free = next(r for r in rows if r.model == "FREE_C9")
    for r in rows:
        if not np.isnan(r.delta_aic_vs_FREE_C9):
            continue
        r.delta_aic_vs_FREE_C9 = r.aic - free.aic
        r.delta_bic_vs_FREE_C9 = r.bic - free.bic

    loo = leave_one_observable_out(data, mode=kernel_mode)

    df_main = pd.DataFrame([
        {
            "model": r.model,
            "role": r.role,
            "k_params": r.k_params,
            "chi2": r.chi2,
            "aic": r.aic,
            "bic": r.bic,
            "delta_aic_vs_FREE_C9": r.delta_aic_vs_FREE_C9,
            "delta_bic_vs_FREE_C9": r.delta_bic_vs_FREE_C9,
            "effective_delta_c9_mean": r.effective_delta_c9_mean,
            "per_observable_delta_c9": r.per_observable_delta_c9,
            "notes": r.notes,
        }
        for r in rows
    ])

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df_main.to_json(out_dir / "wo006_multi_observable.json", orient="records", indent=2)
    loo.to_json(out_dir / "wo006_leave_one_observable_out.json", orient="records", indent=2)
    df_main.to_csv(out_dir / "wo006_multi_observable.csv", index=False)
    loo.to_csv(out_dir / "wo006_leave_one_observable_out.csv", index=False)

    return {
        "data": data,
        "models": df_main,
        "leave_one_observable_out": loo,
        "kernel_mode": kernel_mode,
        "observables": list(observables),
        "config_index": config_index,
        "provenance": PROVENANCE_VFD,
    }


def main() -> None:
    res = run()
    print(res["models"].to_string(index=False))
    print()
    print("Leave-one-observable-out (M2 stability):")
    print(res["leave_one_observable_out"].to_string(index=False))


if __name__ == "__main__":
    main()
