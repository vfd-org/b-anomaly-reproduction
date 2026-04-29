"""WO-008 — Discrete VFD lift of the phi-kernel.

Bridges the empirical kernel (Layer 1) and the continuum bounded-mode
derivation (Layer 2) to the discrete VFD substrate. Builds a 9-shell
symmetric path graph (the natural 1-D projection of the 600-cell
9-isotypic decomposition: shell sizes {1, 12, 20, 12, 30, 12, 20, 12, 1},
symmetric around the central shell), computes its lowest even eigenmode
under three operator variants, projects onto the LHCb bin centres, and
fits a single amplitude A to real LHCb config-2 P5' data.

Operator variants (all on the 9-vertex symmetric path with shell index
m in {-4, ..., +4}):

    A. FREE_DIRICHLET
       L = -Laplacian_path,   psi(+/-4) = 0
       (free closure operator with hard kinematic boundary)

    B. PHI_MASS
       L = -Laplacian_path + (1/phi^2) * I,  psi(+/-4) = 0
       (Layer-2 continuum operator, discretised; the mass term shifts
       all eigenvalues uniformly so eigenvectors equal variant A's,
       but their eigenvalues match the analytical lambda_n exactly).

    C. PHI_COCYCLE
       L = -Laplacian_path + V_m,  V_m = (1/phi^2) + alpha * phi^(m^2),
       psi(+/-4) = 0
       (Layer-2 + the framework's pentagonal cocycle weights as a
       diagonal "boundary-suppression" potential; alpha = 1 by default,
       the framework natural value. m^2 in {0, 1, 4, 9, 16} matches the
       VFD pentagonal cocycle exactly.)

For each variant: lowest-eigenvalue eigenvector psi_1 (normalised to peak 1),
linear interpolation onto LHCb bin centres in the dimensionless x coordinate,
single-amplitude fit to P5' data. Compare chi^2 / AIC against:

    - SM (k = 0)
    - FREE_C9 (k = 1)
    - KAPPA_EXP (Layer-1 Green's function, k = 1, the WO-004A champion)
    - DIRICHLET_M1 (Layer-2 ground state, k = 1)

Acceptance per WO spec:
    - Pearson r between projected eigenmode and exp(-|x|/phi) > 0.95
    - Pearson r between projected eigenmode and cos(pi x / (2 L)) > 0.95
    - dAIC vs FREE_C9 <= 0 (competitive with global C_9 shift)
    - No fitted width, no fitted centre, only amplitude A

This is the bridge result. It does NOT prove that the 600-cell V_600
Laplacian projects to L_phi; it shows that the simplest discrete VFD-style
graph (9-shell symmetric path with phi-cocycle weighting) reproduces the
continuum kernel within statistical resolution of the LHCb data.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from . import hepdata_ingest, vfd_closure
from .constants import C9_SM, PHI, PROVENANCE_VFD
from .likelihood import aic, bic, chi2 as chi2_fn
from .sm_baseline import predict_vector


# 9-shell symmetric path graph centred on the J/psi-psi(2S) midpoint.
# Shell index m in {-4, -3, -2, -1, 0, +1, +2, +3, +4}; vertex 0 is the
# central shell, vertices +/-4 are the kinematic-edge shells.
N_SHELLS = 9
SHELL_INDICES = np.arange(N_SHELLS) - (N_SHELLS - 1) // 2  # [-4..+4]


def shell_x_coordinate(*, x_max: float = vfd_closure.KAPPA_X_MAX) -> np.ndarray:
    """Map shell index m -> dimensionless x: x_m = m * (x_max / 4).

    Endpoints sit AT the kinematic limits: shell -4 -> q^2_min, shell +4 -> q^2_max.
    """
    half = (N_SHELLS - 1) // 2
    return SHELL_INDICES * (x_max / half)


def path_laplacian(n: int) -> np.ndarray:
    """Standard 1-D path graph negative Laplacian: -L = (-d^2/dx^2) on uniform grid.

    Tridiagonal, with -1 off-diagonal and +2 on diagonal (for interior nodes,
    boundary diagonals are also +2 here -- the *Dirichlet* boundary is then
    enforced by deletion of the boundary rows/columns when solving).
    """
    L = np.zeros((n, n), dtype=float)
    for i in range(n):
        L[i, i] = 2.0
        if i > 0:
            L[i, i - 1] = -1.0
        if i < n - 1:
            L[i, i + 1] = -1.0
    return L


def cocycle_potential(m: np.ndarray, *, alpha: float = 1.0) -> np.ndarray:
    """V_m = alpha * phi^(m^2). The exponent m^2 in {0,1,4,9,16} matches the
    pentagonal cocycle of the framework's adaptive-closure-transport paper:
    on the 600-cell, kappa(v) = (shell(v) - 4)^2 in {0, 1, 4, 9, 16} and
    omega_+ = phi^kappa is the σ-Galois cocycle weight.
    """
    return alpha * (PHI ** (m ** 2))


@dataclass
class EigenmodeFit:
    variant: str
    chi2: float
    aic: float
    bic: float
    delta_aic_vs_FREE_C9: float
    delta_aic_vs_KAPPA_EXP: float
    amplitude: float
    eigenvalue: float
    correlation_with_exp: float
    correlation_with_cos: float
    coefficients: dict[str, float] = field(default_factory=dict)
    notes: str = ""


def lowest_even_eigenmode(L_full: np.ndarray) -> tuple[np.ndarray, float]:
    """Solve the Dirichlet eigenproblem L psi = lambda psi by removing the
    boundary vertices (m = -4 and m = +4) and finding the lowest even-parity
    eigenvector. Returns (psi_extended, lambda_min) where psi_extended has
    length N_SHELLS with zeros at the boundary.
    """
    interior = slice(1, -1)
    L_int = L_full[interior, interior]
    eigvals, eigvecs = np.linalg.eigh(L_int)
    n_int = L_int.shape[0]
    centre = n_int // 2
    # Among the eigenvectors find the lowest one that is even-parity (psi(-m) = psi(m)).
    # In a symmetric tridiagonal Dirichlet problem the eigenvectors alternate even/odd
    # starting with even at the lowest eigenvalue, but we verify explicitly.
    for idx in range(len(eigvals)):
        v = eigvecs[:, idx]
        # Even-parity check
        if np.allclose(v, v[::-1], atol=1e-8):
            psi = np.zeros(N_SHELLS)
            psi[1:-1] = v
            # Normalise to peak 1 with positive sign at centre
            sign = np.sign(psi[centre + 1]) or 1.0
            psi = psi * sign
            psi = psi / np.max(np.abs(psi))
            return psi, float(eigvals[idx])
    raise RuntimeError("no even-parity eigenvector found")


def project_to_bin_centres(psi_shell: np.ndarray, x_centres: np.ndarray, *,
                           x_max: float = vfd_closure.KAPPA_X_MAX) -> np.ndarray:
    """Linearly interpolate the shell eigenmode to the LHCb bin x-centres."""
    shell_x = shell_x_coordinate(x_max=x_max)
    return np.interp(x_centres, shell_x, psi_shell)


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


def fit_amplitude_only(
    data: dict[str, Any],
    kappa_at_bins: np.ndarray,
) -> tuple[float, float]:
    """Fit only the amplitude A of Delta_C9 = -A * kappa(q^2). Returns (A, chi2)."""
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)

    def loss(theta):
        a = float(theta[0])
        delta = -a * kappa_at_bins
        return _chi2(values, predict_vector(obs, bins, C9_SM + delta), data, errors)

    r = minimize(loss, x0=[0.5], method="Powell", bounds=[(0.0, 5.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    return float(r.x[0]), float(r.fun)


def run_variant(
    name: str,
    L_full: np.ndarray,
    data: dict[str, Any],
    free_aic: float,
    exp_aic: float,
) -> EigenmodeFit:
    psi_shell, lam = lowest_even_eigenmode(L_full)
    df = data["observables"]
    _, _, q2 = _bin_axis(df)
    x_centres = vfd_closure.kappa_coordinate(q2)
    psi_bins = project_to_bin_centres(psi_shell, x_centres)

    a_hat, c2 = fit_amplitude_only(data, psi_bins)
    n_data = len(df)
    aic_v = aic(c2, 1)
    bic_v = bic(c2, 1, n_data)

    # Correlations against the two continuum benchmarks at the same bin centres.
    L = vfd_closure.KAPPA_X_MAX
    exp_kernel = np.exp(-np.abs(x_centres) / PHI)
    cos_kernel = np.cos(np.pi * x_centres / (2.0 * L))
    r_exp = float(np.corrcoef(psi_bins, exp_kernel)[0, 1])
    r_cos = float(np.corrcoef(psi_bins, cos_kernel)[0, 1])

    return EigenmodeFit(
        variant=name,
        chi2=c2,
        aic=aic_v,
        bic=bic_v,
        delta_aic_vs_FREE_C9=aic_v - free_aic,
        delta_aic_vs_KAPPA_EXP=aic_v - exp_aic,
        amplitude=a_hat,
        eigenvalue=lam,
        correlation_with_exp=r_exp,
        correlation_with_cos=r_cos,
        coefficients={"amplitude": a_hat, "eigenvalue_lambda1": lam},
        notes=f"shell eigenmode psi at shells = {psi_shell.tolist()}",
    )


def run(
    *,
    archive_dir: Path | str = "data/raw/hepdata/extracted",
    config_index: int = 2,
    observable: str = "P5p",
    cocycle_alpha: float = 1.0,
    output_dir: Path | str = "reports",
) -> dict[str, Any]:
    archive = hepdata_ingest.hepdata_archive_dir(archive_dir)
    data = hepdata_ingest.load_config(
        archive, config_index=config_index, observables=(observable,)
    )

    # Reference rows. Re-derive FREE_C9 and KAPPA_EXP here so the comparison is local.
    df = data["observables"]
    obs, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)

    def free_loss(theta):
        return _chi2(values, predict_vector(obs, bins, C9_SM + float(theta[0])), data, errors)
    free = minimize(free_loss, x0=[-0.5], method="Powell", bounds=[(-3, 3)],
                    options={"xtol": 1e-7, "ftol": 1e-9})
    free_chi2 = float(free.fun)
    free_aic = aic(free_chi2, 1)

    kappa_exp = vfd_closure.kappa_shape(q2, mode=vfd_closure.MODE_KAPPA_EXPONENTIAL)
    a_exp, exp_chi2 = fit_amplitude_only(data, kappa_exp)
    exp_aic = aic(exp_chi2, 1)

    # Variant A: free Dirichlet
    L_A = path_laplacian(N_SHELLS)
    fit_A = run_variant("FREE_DIRICHLET", L_A, data, free_aic, exp_aic)

    # Variant B: + phi-mass
    L_B = path_laplacian(N_SHELLS) + (1.0 / PHI ** 2) * np.eye(N_SHELLS)
    fit_B = run_variant("PHI_MASS", L_B, data, free_aic, exp_aic)

    # Variant C: + phi-cocycle potential
    V_m = cocycle_potential(SHELL_INDICES, alpha=cocycle_alpha)
    L_C = path_laplacian(N_SHELLS) + np.diag((1.0 / PHI ** 2) + V_m)
    fit_C = run_variant(f"PHI_COCYCLE[alpha={cocycle_alpha}]", L_C, data, free_aic, exp_aic)

    rows = [
        {
            "model": "SM",
            "k_params": 0,
            "chi2": float(_chi2(values, predict_vector(obs, bins, C9_SM), data, errors)),
            "amplitude": np.nan,
            "delta_aic_vs_FREE_C9": np.nan,
            "delta_aic_vs_KAPPA_EXP": np.nan,
            "correlation_with_exp": np.nan,
            "correlation_with_cos": np.nan,
            "eigenvalue": np.nan,
            "notes": "baseline",
        },
        {
            "model": "FREE_C9",
            "k_params": 1,
            "chi2": free_chi2,
            "amplitude": float(free.x[0]),
            "delta_aic_vs_FREE_C9": 0.0,
            "delta_aic_vs_KAPPA_EXP": free_aic - exp_aic,
            "correlation_with_exp": np.nan,
            "correlation_with_cos": np.nan,
            "eigenvalue": np.nan,
            "notes": "global C_9 shift",
        },
        {
            "model": "KAPPA_EXP (continuum Green's function)",
            "k_params": 1,
            "chi2": exp_chi2,
            "amplitude": a_exp,
            "delta_aic_vs_FREE_C9": exp_aic - free_aic,
            "delta_aic_vs_KAPPA_EXP": 0.0,
            "correlation_with_exp": 1.0,
            "correlation_with_cos": float(np.corrcoef(kappa_exp, np.cos(np.pi * vfd_closure.kappa_coordinate(q2) / (2 * vfd_closure.KAPPA_X_MAX)))[0, 1]),
            "eigenvalue": np.nan,
            "notes": "Layer 1 derivation: Green's fn of L_phi on R",
        },
    ]
    for f in (fit_A, fit_B, fit_C):
        rows.append({
            "model": f"DISCRETE_LIFT[{f.variant}]",
            "k_params": 1,
            "chi2": f.chi2,
            "amplitude": f.amplitude,
            "delta_aic_vs_FREE_C9": f.delta_aic_vs_FREE_C9,
            "delta_aic_vs_KAPPA_EXP": f.delta_aic_vs_KAPPA_EXP,
            "correlation_with_exp": f.correlation_with_exp,
            "correlation_with_cos": f.correlation_with_cos,
            "eigenvalue": f.eigenvalue,
            "notes": f.notes,
        })

    df_main = pd.DataFrame(rows)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df_main.to_csv(out / "wo008_discrete_lift.csv", index=False)
    df_main.to_json(out / "wo008_discrete_lift.json", orient="records", indent=2)

    return {
        "data": data,
        "models": df_main,
        "shell_eigenvectors": {
            "FREE_DIRICHLET": lowest_even_eigenmode(path_laplacian(N_SHELLS))[0].tolist(),
            "PHI_MASS": lowest_even_eigenmode(path_laplacian(N_SHELLS) + (1.0 / PHI ** 2) * np.eye(N_SHELLS))[0].tolist(),
            "PHI_COCYCLE": lowest_even_eigenmode(path_laplacian(N_SHELLS) + np.diag((1.0 / PHI ** 2) + V_m))[0].tolist(),
        },
        "shell_x_coordinate": shell_x_coordinate().tolist(),
        "cocycle_potential": V_m.tolist(),
        "cocycle_alpha": cocycle_alpha,
        "provenance": PROVENANCE_VFD,
    }


def main() -> None:
    res = run()
    print(res["models"].to_string(index=False))
    print()
    print("Shell x-coordinates:", [f"{x:+.3f}" for x in res["shell_x_coordinate"]])
    print("FREE_DIRICHLET psi :", [f"{v:+.3f}" for v in res["shell_eigenvectors"]["FREE_DIRICHLET"]])
    print("PHI_MASS       psi :", [f"{v:+.3f}" for v in res["shell_eigenvectors"]["PHI_MASS"]])
    print("PHI_COCYCLE    psi :", [f"{v:+.3f}" for v in res["shell_eigenvectors"]["PHI_COCYCLE"]])
    print(f"cocycle V_m (alpha={res['cocycle_alpha']}):", [f"{v:.3f}" for v in res["cocycle_potential"]])


if __name__ == "__main__":
    main()
