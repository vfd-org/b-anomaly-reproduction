"""WO-010 — Full angular shared-kernel universality test.

The kernel is FROZEN as of WO-009:

    kappa(q^2) = projection-to-bin-centre of
                 ((L_V600 + (1/phi^2) I)^{-1} * source_equatorial)
                 averaged over each inner-product shell of the 600-cell.

Continuum form: exp(-|x|/phi) with x = (q^2 - q^2_mid)/Delta_psi.
No new variants, no fitted shape, no fitted centre, no fitted width.

This script tests whether ONE shared amplitude A across multiple LHCb
angular observables (P5p, P4p, P1, P2) compresses better than:
    - SM (k = 0)
    - FREE_C9 (k = 1, one global Delta_C9)
    - VFD_GREEN_600CELL (k = 1, one A * frozen kappa)
    - VFD_PER_OBSERVABLE_A (k = N_obs, separate A_o per observable; diagnostic only)

BR is intentionally excluded (rate / normalisation observable, dominated by
form-factor / charm-loop physics not captured by a static C_9 kernel; tested
inconclusively in WO-006).

Acceptance:
    - VFD_GREEN_600CELL AIC <= FREE_C9 AIC
    - All major angular observables prefer the same SIGN of effective DC9 in
      the diagnostic per-observable fit.
    - Leave-one-observable-out: shared amplitude A stable in sign and within
      a factor ~2 in magnitude when any one observable is dropped.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from . import hepdata_ingest, vfd_closure, wo009_full_lift
from .constants import C9_SM, PHI, PROVENANCE_VFD
from .likelihood import aic, bic, chi2 as chi2_fn
from .paper_tables import DEFAULT_Q2_BIN_EDGES_GEV2
from .sm_baseline import predict_vector


def _canonical_bin_index(q2_lo: float, q2_hi: float) -> int:
    """Map a (q2_lo, q2_hi) tuple to its absolute index in the project's
    canonical 8-bin grid. Required so that fits on data subsets (region
    splits, leave-one-observable-out) still query sm_baseline at the right
    bin index."""
    for i, (lo, hi) in enumerate(DEFAULT_Q2_BIN_EDGES_GEV2):
        if abs(lo - q2_lo) < 1e-6 and abs(hi - q2_hi) < 1e-6:
            return i
    raise KeyError(f"q2 bin ({q2_lo}, {q2_hi}) not in canonical grid")


DEFAULT_ANGULAR_OBSERVABLES = ("P5p", "P4p", "P1", "P2")


# -----------------------------------------------------------------------------
# Frozen kernel from WO-009: the 600-cell Green's response, projected to
# bin centres via shell-mean + linear interpolation. Cached.
# -----------------------------------------------------------------------------

_FROZEN_KERNEL_CACHE: dict[tuple, np.ndarray] = {}


def frozen_kernel_at_bin_centres(q2_centres: np.ndarray) -> np.ndarray:
    """Return the WO-009 600-cell Green's kernel evaluated at the given q^2 bin
    centres. The kernel is frozen and depends only on x = (q^2 - mid)/Delta_psi.

    Cached on the input vector identity (since callers will typically reuse
    the same bin centres).
    """
    key = tuple(np.round(q2_centres, 9))
    if key in _FROZEN_KERNEL_CACHE:
        return _FROZEN_KERNEL_CACHE[key].copy()
    verts = wo009_full_lift.generate_600_cell_vertices()
    adj = wo009_full_lift.build_adjacency(verts)
    shell = wo009_full_lift.inner_product_shells(verts, base_idx=0)
    kappa_vert = wo009_full_lift.cocycle_kappa(shell)
    A_w = wo009_full_lift.edge_weights(adj, kappa_vert, mode="unweighted")
    L_w = wo009_full_lift.graph_laplacian(A_w)
    n_shells = 9
    centre_shell = (n_shells - 1) // 2
    centre_mask = shell == centre_shell
    source = np.zeros(len(verts))
    source[centre_mask] = 1.0 / centre_mask.sum()
    psi = wo009_full_lift.discrete_greens_response(L_w, source, mass2=1.0 / (PHI ** 2))
    shell_psi = wo009_full_lift.shell_mean_projection(psi, shell)
    x = vfd_closure.kappa_coordinate(q2_centres)
    kernel = wo009_full_lift.project_to_bin_centres(shell_psi, x)
    _FROZEN_KERNEL_CACHE[key] = kernel.copy()
    return kernel


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _bin_axis(df: pd.DataFrame) -> tuple[list[str], list[int], np.ndarray]:
    obs = df["observable"].tolist()
    bin_indices = [
        _canonical_bin_index(float(r["q2_lo"]), float(r["q2_hi"]))
        for _, r in df.iterrows()
    ]
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
    delta_bic_vs_FREE_C9: float
    effective_delta_c9_mean: float
    per_observable_delta_c9: dict[str, float] = field(default_factory=dict)
    notes: str = ""


def _row(name: str, role: str, k: int, c2: float, n: int) -> tuple[float, float]:
    return aic(c2, k), bic(c2, k, n)


def fit_sm(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    pred = predict_vector(obs, bins, C9_SM)
    c2 = _chi2(values, pred, data, errors)
    a_v, b_v = _row("SM", "baseline", 0, c2, len(values))
    return FitRow("SM", "baseline", 0, c2, a_v, b_v, np.nan, np.nan, 0.0)


def fit_free_c9(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)

    def loss(theta):
        return _chi2(values, predict_vector(obs, bins, C9_SM + float(theta[0])), data, errors)

    r = minimize(loss, x0=[-0.5], method="Powell", bounds=[(-3.0, 3.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    a_v, b_v = _row("FREE_C9", "reference", 1, float(r.fun), len(values))
    delta_hat = float(r.x[0])
    return FitRow(
        "FREE_C9", "reference", 1, float(r.fun), a_v, b_v, 0.0, 0.0, delta_hat,
        per_observable_delta_c9={o: delta_hat for o in set(obs)},
    )


def fit_shared_kernel(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    kappa = frozen_kernel_at_bin_centres(q2)

    def loss(theta):
        a = float(theta[0])
        return _chi2(values, predict_vector(obs, bins, C9_SM - a * kappa), data, errors)

    r = minimize(loss, x0=[0.5], method="Powell", bounds=[(0.0, 5.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    a_hat = float(r.x[0])
    delta_grid = -a_hat * kappa
    a_v, b_v = _row("VFD_GREEN_600CELL", "primary", 1, float(r.fun), len(values))
    eff_per_obs = {o: float(np.mean(delta_grid[np.array(obs) == o])) for o in sorted(set(obs))}
    return FitRow(
        "VFD_GREEN_600CELL", "primary", 1, float(r.fun), a_v, b_v, np.nan, np.nan,
        float(np.mean(delta_grid)),
        per_observable_delta_c9=eff_per_obs,
        notes=f"frozen 600-cell kernel; A={a_hat:.4f}",
    )


def fit_per_observable(data: dict[str, Any]) -> FitRow:
    df = data["observables"]
    obs, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    kappa = frozen_kernel_at_bin_centres(q2)

    obs_unique = sorted(set(obs))
    obs_idx = {o: i for i, o in enumerate(obs_unique)}
    obs_lookup = np.array([obs_idx[o] for o in obs])
    k = len(obs_unique)

    def loss(theta):
        a_arr = np.asarray(theta, dtype=float)
        delta = -a_arr[obs_lookup] * kappa
        return _chi2(values, predict_vector(obs, bins, C9_SM + delta), data, errors)

    # Allow per-observable amplitudes to be either sign; structural compression
    # would expect them all positive (matching the WO-009 sign convention),
    # but we let the data choose so we can read off sign-stability.
    r = minimize(
        loss,
        x0=np.full(k, 0.3),
        method="Powell",
        bounds=[(-3.0, 3.0)] * k,
        options={"xtol": 1e-7, "ftol": 1e-9, "maxiter": 5000},
    )
    a_hat = np.atleast_1d(r.x).astype(float)
    delta_grid = -a_hat[obs_lookup] * kappa
    a_v, b_v = _row("VFD_PER_OBSERVABLE_A", "diagnostic", k, float(r.fun), len(values))
    eff_per_obs = {o: float(np.mean(delta_grid[obs_lookup == obs_idx[o]])) for o in obs_unique}
    a_per_obs = {o: float(a_hat[obs_idx[o]]) for o in obs_unique}
    return FitRow(
        "VFD_PER_OBSERVABLE_A", "diagnostic", k, float(r.fun), a_v, b_v, np.nan, np.nan,
        float(np.mean(delta_grid)),
        per_observable_delta_c9=eff_per_obs,
        notes=f"A_per_obs = {a_per_obs}",
    )


def leave_one_observable_out(data: dict[str, Any]) -> pd.DataFrame:
    df = data["observables"]
    obs_unique = sorted(set(df["observable"].tolist()))
    rows = []
    for drop in obs_unique:
        keep_mask = (df["observable"] != drop).to_numpy()
        sub_df = df[keep_mask].reset_index(drop=True)
        sub_data: dict[str, Any] = {"observables": sub_df, "metadata": data["metadata"]}
        if "covariance" in data and data["covariance"] is not None:
            cov = np.asarray(data["covariance"])
            keep = np.where(keep_mask)[0].tolist()
            sub_data["covariance"] = cov[np.ix_(keep, keep)]
        row = fit_shared_kernel(sub_data)
        # Extract amplitude from notes string
        a_str = row.notes.split("A=")[1] if "A=" in row.notes else "nan"
        rows.append(
            {
                "dropped_observable": drop,
                "n_remaining": len(sub_df),
                "amplitude": float(a_str),
                "effective_delta_c9_mean": row.effective_delta_c9_mean,
                "chi2": row.chi2,
                "aic": row.aic,
            }
        )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Top-level run
# -----------------------------------------------------------------------------

def run(
    *,
    archive_dir: Path | str = "data/raw/hepdata/extracted",
    config_index: int = 2,
    observables: Iterable[str] = DEFAULT_ANGULAR_OBSERVABLES,
    output_dir: Path | str = "reports",
) -> dict[str, Any]:
    archive = hepdata_ingest.hepdata_archive_dir(archive_dir)
    data = hepdata_ingest.load_config(
        archive, config_index=config_index, observables=tuple(observables)
    )

    rows = [
        fit_sm(data),
        fit_free_c9(data),
        fit_shared_kernel(data),
        fit_per_observable(data),
    ]
    # Also test the continuum form of the frozen kernel: kappa = exp(-|x|/phi).
    # Same single-amplitude shape, different (continuum) realisation.
    df_main = data["observables"]
    obs, bins, q2 = _bin_axis(df_main)
    values = df_main["value"].to_numpy(dtype=float)
    errors = np.sqrt(df_main["stat_err"].to_numpy() ** 2 + df_main["syst_err"].to_numpy() ** 2)
    kappa_cont = vfd_closure.kappa_shape(q2, mode=vfd_closure.MODE_KAPPA_EXPONENTIAL)
    def _loss_cont(theta):
        a = float(theta[0])
        return _chi2(values, predict_vector(obs, bins, C9_SM - a * kappa_cont), data, errors)
    r_cont = minimize(_loss_cont, x0=[0.5], method="Powell", bounds=[(0.0, 10.0)],
                     options={"xtol": 1e-7, "ftol": 1e-9})
    a_cont = float(r_cont.x[0])
    c2_cont = float(r_cont.fun)
    a_v_cont = aic(c2_cont, 1)
    b_v_cont = bic(c2_cont, 1, len(values))
    delta_grid_cont = -a_cont * kappa_cont
    eff_cont = {o: float(np.mean(delta_grid_cont[np.array(obs) == o])) for o in sorted(set(obs))}
    rows.append(FitRow(
        "VFD_GREEN_CONTINUUM (Layer 1)", "primary", 1, c2_cont, a_v_cont, b_v_cont,
        np.nan, np.nan, float(np.mean(delta_grid_cont)),
        per_observable_delta_c9=eff_cont,
        notes=f"continuum kernel exp(-|x|/phi); A={a_cont:.4f}",
    ))
    free = next(r for r in rows if r.model == "FREE_C9")
    for r in rows:
        if np.isnan(r.delta_aic_vs_FREE_C9):
            r.delta_aic_vs_FREE_C9 = r.aic - free.aic
            r.delta_bic_vs_FREE_C9 = r.bic - free.bic

    loo = leave_one_observable_out(data)

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

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df_main.to_csv(out / "wo010_universality.csv", index=False)
    df_main.to_json(out / "wo010_universality.json", orient="records", indent=2)
    loo.to_csv(out / "wo010_loo.csv", index=False)
    loo.to_json(out / "wo010_loo.json", orient="records", indent=2)

    return {
        "models": df_main,
        "leave_one_observable_out": loo,
        "observables": list(observables),
        "config_index": config_index,
        "n_data": len(data["observables"]),
        "provenance": PROVENANCE_VFD,
    }


def main() -> None:
    res = run()
    pd.set_option("display.max_colwidth", 100)
    print(f"Joint fit on {res['observables']}: n_data = {res['n_data']}")
    print()
    cols = ["model", "role", "k_params", "chi2", "aic",
            "delta_aic_vs_FREE_C9", "effective_delta_c9_mean",
            "per_observable_delta_c9", "notes"]
    print(res["models"][cols].to_string(index=False))
    print()
    print("Leave-one-observable-out (shared-kernel stability):")
    print(res["leave_one_observable_out"].to_string(index=False))


if __name__ == "__main__":
    main()
