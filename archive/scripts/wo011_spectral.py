"""WO-011 — Spectral decomposition of the V_600 Green's response.

The frozen kernel from WO-009 is

    psi(v) = (L_V600 + (1/phi^2) I)^{-1} f,    f = source on equatorial shell.

Spectrally,

    psi(v) = sum_n  <psi_n, f> / (lambda_n + 1/phi^2) * psi_n(v)

where (lambda_n, psi_n) is the spectrum of L_V600. This script:

    1. Computes the full eigenspectrum of L_V600 (120 eigenpairs).
    2. Reconstructs psi using truncated sums sum_{n <= N} for N = 1, 2, ..., 120.
    3. Measures convergence:
         - Pearson r between truncated kernel (shell-mean -> bin centres) and
           the full Green's response.
         - chi^2 of the amplitude-only fit to LHCb P5' data for each truncation.
    4. Checks whether HIGH-mode truncations introduce oscillations that the
       data rejects (i.e. higher dAIC vs FREE_C9 than the rank-1 mode alone).

Independent of the SM tables for the angular observables: P5' alone is
used because that is the only well-anchored observable in the project.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from . import hepdata_ingest, vfd_closure, wo009_full_lift
from .constants import C9_SM, PHI, PROVENANCE_VFD
from .likelihood import aic, bic, chi2 as chi2_fn
from .sm_baseline import predict_vector


def _bin_axis(df):
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


def _chi2(values, pred, data, errors):
    cov = data.get("covariance")
    if cov is not None:
        return chi2_fn(values, pred, covariance=cov)
    return chi2_fn(values, pred, errors=errors)


def amplitude_fit(data, kernel_at_bins):
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)

    def loss(theta):
        a = float(theta[0])
        return _chi2(values, predict_vector(obs, bins, C9_SM - a * kernel_at_bins),
                     data, errors)

    r = minimize(loss, x0=[0.5], method="Powell", bounds=[(-5.0, 5.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    return float(r.x[0]), float(r.fun)


@dataclass
class TruncationResult:
    n_modes: int
    eigenvalues_used: list[float]
    reconstruction_relative_error: float
    correlation_with_full_kernel: float
    correlation_with_continuum_exp: float
    chi2_p5p: float
    aic_p5p: float
    amplitude: float
    shell_psi: list[float]


def run(
    *,
    archive_dir: Path | str = "data/raw/hepdata/extracted",
    config_index: int = 2,
    output_dir: Path | str = "reports",
    n_max_to_record: int = 30,
) -> dict[str, Any]:
    # ---- Build V_600 graph + Laplacian + source ----
    verts = wo009_full_lift.generate_600_cell_vertices()
    adj = wo009_full_lift.build_adjacency(verts)
    shell = wo009_full_lift.inner_product_shells(verts, base_idx=0)
    n_shells = 9
    centre_shell = (n_shells - 1) // 2
    centre_mask = shell == centre_shell
    source = np.zeros(len(verts))
    source[centre_mask] = 1.0 / centre_mask.sum()

    A_w = wo009_full_lift.edge_weights(adj, wo009_full_lift.cocycle_kappa(shell), mode="unweighted")
    L_w = wo009_full_lift.graph_laplacian(A_w)
    mass2 = 1.0 / (PHI ** 2)

    # ---- Full spectrum ----
    eigvals, eigvecs = np.linalg.eigh(L_w)
    # Sort ascending (eigh already does, but be explicit)
    order = np.argsort(eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    # ---- Full Green's response (reference) ----
    full_psi = wo009_full_lift.discrete_greens_response(L_w, source, mass2=mass2)

    # ---- Spectral reconstruction at increasing truncation ----
    coeffs = eigvecs.T @ source            # <psi_n, f>
    weights = coeffs / (eigvals + mass2)   # spectral coefficients

    # ---- LHCb P5' data ----
    archive = hepdata_ingest.hepdata_archive_dir(archive_dir)
    data = hepdata_ingest.load_config(archive, config_index=config_index, observables=("P5p",))
    df = data["observables"]
    obs_list, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    x_centres = vfd_closure.kappa_coordinate(q2)

    # Reference fits
    free_loss = lambda t: _chi2(values, predict_vector(obs_list, bins, C9_SM + float(t[0])),
                                 data, errors)
    free = minimize(free_loss, x0=[-0.5], method="Powell", bounds=[(-3, 3)],
                    options={"xtol": 1e-7, "ftol": 1e-9})
    free_chi2 = float(free.fun)
    free_aic = aic(free_chi2, 1)

    # Continuum benchmark
    cont_exp = np.exp(-np.abs(x_centres) / PHI)

    # ---- Loop over truncation orders ----
    results: list[TruncationResult] = []
    cumulative = np.zeros_like(full_psi)
    for n in range(1, len(eigvals) + 1):
        cumulative = cumulative + weights[n - 1] * eigvecs[:, n - 1]
        # Project to bin centres
        shell_psi = wo009_full_lift.shell_mean_projection(cumulative, shell)
        kernel_at_bins = wo009_full_lift.project_to_bin_centres(shell_psi, x_centres)
        # Skip if kernel is identically zero (would NaN the correlation)
        if np.all(np.abs(kernel_at_bins) < 1e-12):
            continue
        a_hat, c2 = amplitude_fit(data, kernel_at_bins)
        full_kernel = full_psi  # reference
        rel_err = float(np.linalg.norm(cumulative - full_psi) / np.linalg.norm(full_psi))
        # Correlation of bin-centre kernel with FULL bin-centre kernel
        full_shell = wo009_full_lift.shell_mean_projection(full_psi, shell)
        full_at_bins = wo009_full_lift.project_to_bin_centres(full_shell, x_centres)
        # When the truncated kernel is constant, correlation undefined
        if np.std(kernel_at_bins) < 1e-12 or np.std(full_at_bins) < 1e-12:
            r_full = float("nan")
        else:
            r_full = float(np.corrcoef(kernel_at_bins, full_at_bins)[0, 1])
        if np.std(kernel_at_bins) < 1e-12:
            r_cont = float("nan")
        else:
            r_cont = float(np.corrcoef(kernel_at_bins, cont_exp)[0, 1])
        results.append(TruncationResult(
            n_modes=n,
            eigenvalues_used=eigvals[:n].tolist(),
            reconstruction_relative_error=rel_err,
            correlation_with_full_kernel=r_full,
            correlation_with_continuum_exp=r_cont,
            chi2_p5p=c2,
            aic_p5p=aic(c2, 1),
            amplitude=a_hat,
            shell_psi=shell_psi.tolist(),
        ))

    # ---- Summary table (record key truncation milestones) ----
    rows = []
    record_set = set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                      15, 20, 30, 50, 60, 80, 100, len(eigvals)])
    for r in results:
        if r.n_modes in record_set or r.n_modes <= n_max_to_record:
            rows.append({
                "n_modes": r.n_modes,
                "rel_recon_err": r.reconstruction_relative_error,
                "r_vs_full_kernel": r.correlation_with_full_kernel,
                "r_vs_continuum_exp": r.correlation_with_continuum_exp,
                "chi2_p5p": r.chi2_p5p,
                "aic_p5p": r.aic_p5p,
                "delta_aic_vs_FREE_C9": r.aic_p5p - free_aic,
                "amplitude": r.amplitude,
                "lambda_max_used": r.eigenvalues_used[-1] if r.eigenvalues_used else float("nan"),
            })

    df_main = pd.DataFrame(rows)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df_main.to_csv(out / "wo011_spectral.csv", index=False)
    df_main.to_json(out / "wo011_spectral.json", orient="records", indent=2)

    # ---- Spectrum diagnostics ----
    spectrum_df = pd.DataFrame({
        "mode_index": list(range(1, len(eigvals) + 1)),
        "eigenvalue": eigvals.tolist(),
        "coefficient_psi_n_dot_f": coeffs.tolist(),
        "spectral_weight": weights.tolist(),
        "abs_weight": np.abs(weights).tolist(),
    })
    spectrum_df.to_csv(out / "wo011_spectrum.csv", index=False)

    # ---- Find smallest N_min such that r >= 0.99 with full kernel ----
    n_min_99 = next((r.n_modes for r in results
                     if not np.isnan(r.correlation_with_full_kernel)
                     and r.correlation_with_full_kernel >= 0.99), None)
    n_min_999 = next((r.n_modes for r in results
                      if not np.isnan(r.correlation_with_full_kernel)
                      and r.correlation_with_full_kernel >= 0.999), None)

    # The N that minimises chi^2 to LHCb P5' (with k = 1 amplitude)
    chi2_arr = np.array([r.chi2_p5p for r in results])
    best_n = results[int(np.argmin(chi2_arr))].n_modes
    best_chi2 = float(np.min(chi2_arr))

    return {
        "summary": df_main,
        "spectrum": spectrum_df,
        "free_c9_chi2": free_chi2,
        "free_c9_aic": free_aic,
        "n_min_r_99_with_full_kernel": n_min_99,
        "n_min_r_999_with_full_kernel": n_min_999,
        "best_n_modes_for_p5p": best_n,
        "best_chi2_p5p": best_chi2,
        "kernel_continuum_correlation": float(np.corrcoef(
            wo009_full_lift.project_to_bin_centres(
                wo009_full_lift.shell_mean_projection(full_psi, shell), x_centres),
            cont_exp)[0, 1]),
        "provenance": PROVENANCE_VFD,
    }


def main() -> None:
    res = run()
    pd.set_option("display.width", 140)
    pd.set_option("display.max_rows", 50)
    print(f"FREE_C9 reference chi^2 = {res['free_c9_chi2']:.4f}, AIC = {res['free_c9_aic']:.4f}")
    print(f"Full kernel correlation with continuum exp: r = {res['kernel_continuum_correlation']:.4f}")
    print(f"Smallest N for r >= 0.99 with full kernel: {res['n_min_r_99_with_full_kernel']}")
    print(f"Smallest N for r >= 0.999 with full kernel: {res['n_min_r_999_with_full_kernel']}")
    print(f"Best N for P5' chi^2 minimisation: {res['best_n_modes_for_p5p']}, chi^2 = {res['best_chi2_p5p']:.4f}")
    print()
    print("Truncation summary:")
    print(res["summary"].to_string(index=False))


if __name__ == "__main__":
    main()
