"""WO-009 — Full 2I-equivariant edge-space lift on the 600-cell.

Replaces the 1-D shell path projection of WO-008 with the genuine VFD
graph operator: the 120-vertex 600-cell V_600 (binary icosahedral group 2I),
12-regular, 720 edges, with the pentagonal cocycle applied as edge weights
(not a manually-tuned shell potential).

Pipeline:
    1. Generate the 120 vertices of the 600-cell as unit 4-vectors (icosian
       coordinates: 8 + 16 + 96 = 120).
    2. Build adjacency: two vertices are joined iff <v, w> = phi/2 (the
       canonical 600-cell edge condition; equivalent edge length = 1/phi).
    3. Choose a base vertex v_0 and compute graph distances. The shell
       sizes must come out as the framework values {1, 12, 20, 12, 30, 12,
       20, 12, 1} -- this is sanity check #1.
    4. Pentagonal cocycle on vertices: kappa(v) = (shell(v) - 4)^2 in
       {0, 1, 4, 9, 16} (the framework's pentagonal cocycle weights).
    5. EDGE-SPACE lift: w_{vw} = phi^{(kappa(v) + kappa(w))/2} (geometric
       mean of vertex cocycle weights -- preserves multiplicative phi-rational
       structure). This is the proper edge-space analog of WO-008's
       diagonal vertex potential.
    6. Build the phi-weighted graph Laplacian L_w = D_w - A_w on V_600.
    7. Find the lowest non-constant *even-parity* eigenvector (even = invariant
       under v -> -v, the antipodal map). Even parity is the analog of the
       shell-symmetric "central" mode.
    8. Project to 9-shell coordinate by shell-mean. Map shells to the
       dimensionless x via x_s = (s - 4) * (x_max / 4) where x_max = 2.886
       from kinematics.
    9. Linearly interpolate to LHCb bin centres.
   10. Fit a single amplitude A to P5' (Delta_C9 = -A * eigenmode).

Acceptance:
    - r between projected eigenmode and exp(-|x|/phi) >= 0.95
    - amplitude is the only fitted parameter (no shape, no width, no centre)
    - dAIC vs FREE_C9 <= 0
    - shell sizes verified == {1, 12, 20, 12, 30, 12, 20, 12, 1} (no tuning)
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from . import hepdata_ingest, vfd_closure
from .constants import C9_SM, PHI, PROVENANCE_VFD
from .likelihood import aic, bic, chi2 as chi2_fn
from .sm_baseline import predict_vector


N_VERTICES = 120
EXPECTED_DEGREE = 12
EXPECTED_SHELL_SIZES = (1, 12, 20, 12, 30, 12, 20, 12, 1)
TOL = 1e-9


# -----------------------------------------------------------------------------
# Step 1: 600-cell vertex set
# -----------------------------------------------------------------------------

def generate_600_cell_vertices() -> np.ndarray:
    """Return the 120 vertices of the 600-cell as unit 4-vectors.

    Standard construction (Conway-Sloane, Coxeter):
      Group 1 (8):  unit basis vectors and their negatives, +/-e_i
      Group 2 (16): all (+/-1/2, +/-1/2, +/-1/2, +/-1/2)
      Group 3 (96): even permutations of (0, +/-1/2, +/-1/(2 phi), +/-phi/2)
    """
    verts: list[tuple[float, float, float, float]] = []

    # Group 1
    for i in range(4):
        for sign in (+1.0, -1.0):
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = sign
            verts.append(tuple(v))

    # Group 2
    for s0 in (+1, -1):
        for s1 in (+1, -1):
            for s2 in (+1, -1):
                for s3 in (+1, -1):
                    verts.append((s0 / 2.0, s1 / 2.0, s2 / 2.0, s3 / 2.0))

    # Group 3: even permutations of (0, 1/2, 1/(2 phi), phi/2)
    base_vals = (0.0, 0.5, 1.0 / (2.0 * PHI), PHI / 2.0)
    nonzero_indices = [i for i, v in enumerate(base_vals) if v != 0.0]
    for perm in permutations(range(4)):
        # Even permutation iff inversions even
        invs = sum(
            1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j]
        )
        if invs % 2 != 0:
            continue
        # Permuted base: out[i] = base_vals[perm[i]]
        permuted = [base_vals[perm[i]] for i in range(4)]
        # Find positions of nonzero entries in the permuted vector
        nz_positions = [i for i, x in enumerate(permuted) if x != 0.0]
        # Apply all 8 sign choices to those three positions
        for sign_pattern in range(8):
            v = list(permuted)
            for k, pos in enumerate(nz_positions):
                if (sign_pattern >> k) & 1:
                    v[pos] = -v[pos]
            verts.append(tuple(v))

    # Deduplicate (Group 3 enumeration produces some repeats across permutations).
    unique = sorted(set(verts))
    arr = np.array(unique, dtype=float)
    if arr.shape[0] != N_VERTICES:
        raise RuntimeError(
            f"600-cell construction yielded {arr.shape[0]} vertices, expected {N_VERTICES}"
        )
    # Sanity: all unit-norm
    norms = np.linalg.norm(arr, axis=1)
    if not np.allclose(norms, 1.0, atol=1e-9):
        raise RuntimeError("non-unit norms in 600-cell construction")
    return arr


# -----------------------------------------------------------------------------
# Step 2 & 3: adjacency, shells, antipode map
# -----------------------------------------------------------------------------

def build_adjacency(verts: np.ndarray) -> np.ndarray:
    """Return the symmetric 0/1 adjacency matrix. Edges: <v, w> = phi/2."""
    inner = verts @ verts.T
    target = PHI / 2.0
    A = (np.abs(inner - target) < 1e-9).astype(float)
    np.fill_diagonal(A, 0.0)
    deg = A.sum(axis=1).astype(int)
    if not np.all(deg == EXPECTED_DEGREE):
        raise RuntimeError(f"non-uniform degree: {set(deg.tolist())}")
    return A


def antipode_index(verts: np.ndarray) -> np.ndarray:
    """For each vertex v, return the index of -v."""
    n = len(verts)
    out = -np.ones(n, dtype=int)
    for i in range(n):
        diffs = verts + verts[i]
        idx = np.where(np.linalg.norm(diffs, axis=1) < 1e-9)[0]
        if len(idx) != 1:
            raise RuntimeError(f"antipode lookup failed for vertex {i}")
        out[i] = idx[0]
    return out


def shell_distances(adj: np.ndarray, base_idx: int = 0) -> np.ndarray:
    """BFS from base_idx; return per-vertex graph distance. Note: the 600-cell's
    graph diameter is only 5, while its *Euclidean / inner-product* shell
    decomposition has 9 strata. For the framework's pentagonal cocycle we
    want the inner-product shells (see `inner_product_shells`).
    """
    n = adj.shape[0]
    dist = -np.ones(n, dtype=int)
    dist[base_idx] = 0
    queue = deque([base_idx])
    while queue:
        u = queue.popleft()
        for w in np.where(adj[u] > 0)[0]:
            if dist[w] == -1:
                dist[w] = dist[u] + 1
                queue.append(w)
    if (dist < 0).any():
        raise RuntimeError("graph not connected")
    return dist


def inner_product_shells(verts: np.ndarray, base_idx: int = 0) -> np.ndarray:
    """Return the 9-shell stratification by inner product with the base vertex.

    The 600-cell has 9 distinct inner-product values from any vertex:
        {+1, +phi/2, +1/2, +1/(2 phi), 0, -1/(2 phi), -1/2, -phi/2, -1}
    Shell index s in {0..8} maps in DECREASING inner product so that shell 0
    is the base, shell 4 is orthogonal (the equatorial 30-vertex shell), and
    shell 8 is the antipode. Sizes come out as {1, 12, 20, 12, 30, 12, 20, 12, 1}
    which matches the framework's 2I-isotypic decomposition.
    """
    inner = verts @ verts[base_idx]
    expected_values = np.array([
        +1.0,
        +PHI / 2.0,
        +0.5,
        +1.0 / (2.0 * PHI),
        0.0,
        -1.0 / (2.0 * PHI),
        -0.5,
        -PHI / 2.0,
        -1.0,
    ])
    shells = -np.ones(len(verts), dtype=int)
    for s, val in enumerate(expected_values):
        mask = np.abs(inner - val) < 1e-9
        shells[mask] = s
    if (shells < 0).any():
        bad = inner[shells < 0]
        raise RuntimeError(f"unmatched inner-product values: {sorted(set(bad.tolist()))}")
    return shells


# -----------------------------------------------------------------------------
# Steps 4-5: pentagonal cocycle and edge weighting
# -----------------------------------------------------------------------------

def cocycle_kappa(shell: np.ndarray, n_shells: int = 9) -> np.ndarray:
    """kappa(v) = (shell(v) - midshell)^2. For 9 shells indexed 0..8, midshell = 4.
    Values fall in {0, 1, 4, 9, 16}. This matches the framework's pentagonal
    cocycle exponent.
    """
    midshell = (n_shells - 1) // 2
    return (shell - midshell) ** 2


def edge_weights(
    adj: np.ndarray, kappa: np.ndarray, *, mode: str = "geometric"
) -> np.ndarray:
    """Compute edge-weighted adjacency A_w[i,j] = w_{ij} * adj[i,j].

    mode = 'geometric'  : w_{ij} = phi^{(kappa_i + kappa_j) / 2}
    mode = 'arithmetic' : w_{ij} = (phi^{kappa_i} + phi^{kappa_j}) / 2
    mode = 'unweighted' : w_{ij} = 1
    """
    n = adj.shape[0]
    A_w = np.zeros((n, n))
    if mode == "unweighted":
        return adj.copy()
    for i in range(n):
        for j in range(n):
            if adj[i, j] == 0:
                continue
            if mode == "geometric":
                A_w[i, j] = PHI ** ((kappa[i] + kappa[j]) / 2.0)
            elif mode == "arithmetic":
                A_w[i, j] = 0.5 * (PHI ** kappa[i] + PHI ** kappa[j])
            else:
                raise ValueError(f"unknown edge-weight mode {mode!r}")
    return A_w


def graph_laplacian(A_w: np.ndarray) -> np.ndarray:
    """Weighted graph Laplacian L = D - A_w (D diagonal of row sums)."""
    D = np.diag(A_w.sum(axis=1))
    return D - A_w


def discrete_greens_response(
    L: np.ndarray,
    source: np.ndarray,
    *,
    mass2: float = 1.0 / (PHI ** 2),
) -> np.ndarray:
    """Return phi(v) = (L + mass2 * I)^{-1} source, the discrete analog of
    the Layer-1 Green's function response of L_phi = -d^2/dx^2 + 1/phi^2 to
    a unit source. The mass term regularises the inverse and makes the
    operator strictly positive definite (equivalently, suppresses the
    constant zero-eigenvalue gauge mode).
    """
    n = L.shape[0]
    M = L + mass2 * np.eye(n)
    return np.linalg.solve(M, source)


# -----------------------------------------------------------------------------
# Step 7: lowest non-constant even-parity eigenvector
# -----------------------------------------------------------------------------

def even_subspace_basis(antipode: np.ndarray) -> np.ndarray:
    """Build orthonormal basis B (n x n_even) of the even subspace under v -> -v.

    For the 600-cell (no self-antipodal vertices) the even subspace has
    dimension n / 2 = 60. Each basis vector is (e_i + e_{antipode(i)}) / sqrt(2)
    for one representative per antipodal pair.
    """
    n = len(antipode)
    if np.any(antipode == np.arange(n)):
        raise RuntimeError("self-antipodal vertices not supported")
    representatives = [i for i in range(n) if i < antipode[i]]
    n_even = len(representatives)
    B = np.zeros((n, n_even))
    for j, i in enumerate(representatives):
        B[i, j] = 1.0 / np.sqrt(2.0)
        B[antipode[i], j] = 1.0 / np.sqrt(2.0)
    return B


def lowest_even_eigenvector(
    L: np.ndarray, antipode: np.ndarray, *, shell: np.ndarray | None = None
) -> tuple[np.ndarray, float, dict[str, float]]:
    """Return (eigenvector, eigenvalue, info) for the lowest *non-trivial*
    eigenvector restricted to the even-parity subspace under v -> -v.

    Restriction is exact: we project L onto the even subspace and diagonalise
    the 60-dimensional reduced operator, then lift back. This guarantees the
    returned eigenvector is purely even parity (no contamination from
    eigenvalue-degenerate odd mixtures).
    """
    n = L.shape[0]
    B = even_subspace_basis(antipode)
    L_even = B.T @ L @ B
    eigvals, eigvecs = np.linalg.eigh(L_even)
    info: dict[str, float] = {}
    for idx in range(len(eigvals)):
        v_red = eigvecs[:, idx]
        v = B @ v_red
        # Skip the gauge / constant-on-each-pair eigenvector that yields a
        # spatially uniform mode (eigenvalue ~ 0 for the standard graph
        # Laplacian; for a weighted Laplacian still corresponds to L * 1 = 0
        # because D_w - A_w sums to zero on the constant vector).
        v_normed = v / (np.linalg.norm(v) + 1e-30)
        if np.allclose(np.abs(v_normed - np.mean(v_normed)), 0.0, atol=1e-6):
            continue
        # Normalise to peak 1 with positive sign at the largest-magnitude
        # component. If shell info is provided, force positive sign at the
        # central shell (shell index 4) so the kernel is centre-peaked when
        # interpreted as a closure response.
        scale = np.max(np.abs(v))
        v = v / scale
        if shell is not None:
            centre_mask = shell == ((shell.max() // 2))
            centre_val = float(v[centre_mask].mean()) if centre_mask.any() else float(v[np.argmax(np.abs(v))])
            if centre_val < 0:
                v = -v
        else:
            if v[np.argmax(np.abs(v))] < 0:
                v = -v
        info["selected_idx"] = float(idx)
        info["selected_lambda"] = float(eigvals[idx])
        return v, float(eigvals[idx]), info
    raise RuntimeError("no non-trivial even-parity eigenvector found")


# -----------------------------------------------------------------------------
# Step 8: shell projection
# -----------------------------------------------------------------------------

def shell_mean_projection(psi: np.ndarray, shell: np.ndarray) -> np.ndarray:
    """Average the eigenvector over each shell to get a per-shell value."""
    n_shells = int(shell.max()) + 1
    out = np.zeros(n_shells)
    for s in range(n_shells):
        mask = shell == s
        out[s] = float(psi[mask].mean()) if mask.any() else 0.0
    # Normalise to peak 1 with positive sign at centre
    if out[(n_shells - 1) // 2] < 0:
        out = -out
    out = out / np.max(np.abs(out))
    return out


# -----------------------------------------------------------------------------
# Step 9-10: project to LHCb bin centres and amplitude fit
# -----------------------------------------------------------------------------

def shell_x_positions(n_shells: int = 9, *, x_max: float = vfd_closure.KAPPA_X_MAX) -> np.ndarray:
    half = (n_shells - 1) // 2
    return (np.arange(n_shells) - half) * (x_max / half)


def project_to_bin_centres(shell_psi: np.ndarray, x_centres: np.ndarray) -> np.ndarray:
    shell_x = shell_x_positions(len(shell_psi))
    return np.interp(x_centres, shell_x, shell_psi)


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


def fit_amplitude_only(data: dict[str, Any], kernel_at_bins: np.ndarray) -> tuple[float, float]:
    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)

    def loss(theta):
        a = float(theta[0])
        return _chi2(values, predict_vector(obs, bins, C9_SM - a * kernel_at_bins), data, errors)

    r = minimize(loss, x0=[0.5], method="Powell", bounds=[(0.0, 5.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    return float(r.x[0]), float(r.fun)


# -----------------------------------------------------------------------------
# Top-level run
# -----------------------------------------------------------------------------

@dataclass
class VariantResult:
    name: str
    chi2: float
    aic: float
    amplitude: float
    eigenvalue: float
    correlation_with_exp: float
    correlation_with_cos: float
    delta_aic_vs_FREE_C9: float
    delta_aic_vs_KAPPA_EXP: float
    shell_psi: list[float]


def run(
    *,
    archive_dir: Path | str = "data/raw/hepdata/extracted",
    config_index: int = 2,
    observable: str = "P5p",
    output_dir: Path | str = "reports",
) -> dict[str, Any]:
    # ---- Step 1-3: build the 600-cell ----
    verts = generate_600_cell_vertices()
    adj = build_adjacency(verts)
    n_edges = int(adj.sum() / 2)
    if n_edges != 720:
        raise RuntimeError(f"expected 720 edges, got {n_edges}")

    base = 0  # the (1,0,0,0) vertex (or whichever ends up at index 0)
    # Use inner-product shells (the framework's 9-isotypic decomposition),
    # NOT BFS shells (the 600-cell graph has diameter 5).
    shell = inner_product_shells(verts, base_idx=base)
    n_shells = 9
    shell_sizes = tuple(int((shell == s).sum()) for s in range(n_shells))
    if shell_sizes != EXPECTED_SHELL_SIZES:
        raise RuntimeError(
            f"shell-size mismatch: got {shell_sizes}, expected {EXPECTED_SHELL_SIZES}"
        )

    antipode = antipode_index(verts)
    if not np.array_equal(shell[antipode], 8 - shell):
        raise RuntimeError("antipodal shell symmetry broken")

    kappa = cocycle_kappa(shell)

    # ---- Step 4-7: variants ----
    variants: list[tuple[str, np.ndarray]] = [
        ("UNWEIGHTED",     edge_weights(adj, kappa, mode="unweighted")),
        ("PHI_GEOMETRIC",  edge_weights(adj, kappa, mode="geometric")),
        ("PHI_ARITHMETIC", edge_weights(adj, kappa, mode="arithmetic")),
    ]

    # ---- Load data + reference fits ----
    archive = hepdata_ingest.hepdata_archive_dir(archive_dir)
    data = hepdata_ingest.load_config(archive, config_index=config_index, observables=(observable,))
    df = data["observables"]
    obs_list, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    x_centres = vfd_closure.kappa_coordinate(q2)

    # FREE_C9 reference
    def free_loss(theta):
        return _chi2(values, predict_vector(obs_list, bins, C9_SM + float(theta[0])), data, errors)
    free = minimize(free_loss, x0=[-0.5], method="Powell", bounds=[(-3, 3)],
                    options={"xtol": 1e-7, "ftol": 1e-9})
    free_chi2 = float(free.fun)
    free_aic_v = aic(free_chi2, 1)

    # KAPPA_EXP reference
    kappa_exp = vfd_closure.kappa_shape(q2, mode=vfd_closure.MODE_KAPPA_EXPONENTIAL)
    a_exp, exp_chi2 = fit_amplitude_only(data, kappa_exp)
    exp_aic_v = aic(exp_chi2, 1)

    # Continuum benchmarks at bin centres
    L_kine = vfd_closure.KAPPA_X_MAX
    bench_exp = np.exp(-np.abs(x_centres) / PHI)
    bench_cos = np.cos(np.pi * x_centres / (2.0 * L_kine))

    # Centred source: uniform on shell 4 (the equatorial 30-vertex shell), zero
    # elsewhere. This is the natural discrete analog of a delta source at x = 0
    # in the continuum, normalised so the source vector has unit total mass.
    centre_shell = (n_shells - 1) // 2
    source = np.zeros(N_VERTICES)
    centre_mask = shell == centre_shell
    source[centre_mask] = 1.0 / centre_mask.sum()

    results: list[VariantResult] = []
    for name, A_w in variants:
        L_w = graph_laplacian(A_w)
        # ---- (a) lowest non-trivial even eigenvector ----
        psi_eig, lam_eig, _ = lowest_even_eigenvector(L_w, antipode, shell=shell)
        shell_psi_eig = shell_mean_projection(psi_eig, shell)
        psi_eig_bins = project_to_bin_centres(shell_psi_eig, x_centres)
        a_eig, c2_eig = fit_amplitude_only(data, psi_eig_bins)
        aic_eig = aic(c2_eig, 1)
        results.append(VariantResult(
            name=f"FULL_LIFT[{name}]_EIGENMODE",
            chi2=c2_eig, aic=aic_eig, amplitude=a_eig, eigenvalue=lam_eig,
            correlation_with_exp=float(np.corrcoef(psi_eig_bins, bench_exp)[0, 1]),
            correlation_with_cos=float(np.corrcoef(psi_eig_bins, bench_cos)[0, 1]),
            delta_aic_vs_FREE_C9=aic_eig - free_aic_v,
            delta_aic_vs_KAPPA_EXP=aic_eig - exp_aic_v,
            shell_psi=shell_psi_eig.tolist(),
        ))
        # ---- (b) discrete Green's function response from centred source ----
        psi_g = discrete_greens_response(L_w, source, mass2=1.0 / (PHI ** 2))
        shell_psi_g = shell_mean_projection(psi_g, shell)
        psi_g_bins = project_to_bin_centres(shell_psi_g, x_centres)
        a_g, c2_g = fit_amplitude_only(data, psi_g_bins)
        aic_g = aic(c2_g, 1)
        # Effective "eigenvalue": min eigenvalue of L_w + mass2 (always positive)
        lam_g = float(np.linalg.eigvalsh(L_w).min() + 1.0 / PHI ** 2)
        results.append(VariantResult(
            name=f"FULL_LIFT[{name}]_GREENS",
            chi2=c2_g, aic=aic_g, amplitude=a_g, eigenvalue=lam_g,
            correlation_with_exp=float(np.corrcoef(psi_g_bins, bench_exp)[0, 1]),
            correlation_with_cos=float(np.corrcoef(psi_g_bins, bench_cos)[0, 1]),
            delta_aic_vs_FREE_C9=aic_g - free_aic_v,
            delta_aic_vs_KAPPA_EXP=aic_g - exp_aic_v,
            shell_psi=shell_psi_g.tolist(),
        ))

    # ---- Build summary table ----
    rows = [
        {"model": "FREE_C9", "k_params": 1, "chi2": free_chi2, "amplitude": float(free.x[0]),
         "delta_aic_vs_FREE_C9": 0.0, "delta_aic_vs_KAPPA_EXP": free_aic_v - exp_aic_v,
         "correlation_with_exp": np.nan, "correlation_with_cos": np.nan,
         "eigenvalue": np.nan, "notes": "global C9 shift"},
        {"model": "KAPPA_EXP (Layer 1)", "k_params": 1, "chi2": exp_chi2, "amplitude": a_exp,
         "delta_aic_vs_FREE_C9": exp_aic_v - free_aic_v, "delta_aic_vs_KAPPA_EXP": 0.0,
         "correlation_with_exp": 1.0,
         "correlation_with_cos": float(np.corrcoef(kappa_exp, bench_cos)[0, 1]),
         "eigenvalue": np.nan, "notes": "Green's fn of L_phi on R"},
    ]
    for r in results:
        rows.append({
            "model": r.name, "k_params": 1, "chi2": r.chi2, "amplitude": r.amplitude,
            "delta_aic_vs_FREE_C9": r.delta_aic_vs_FREE_C9,
            "delta_aic_vs_KAPPA_EXP": r.delta_aic_vs_KAPPA_EXP,
            "correlation_with_exp": r.correlation_with_exp,
            "correlation_with_cos": r.correlation_with_cos,
            "eigenvalue": r.eigenvalue,
            "notes": f"shell-mean psi = {[round(x, 3) for x in r.shell_psi]}",
        })

    df_main = pd.DataFrame(rows)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df_main.to_csv(out / "wo009_full_lift.csv", index=False)
    df_main.to_json(out / "wo009_full_lift.json", orient="records", indent=2)

    return {
        "models": df_main,
        "n_vertices": N_VERTICES,
        "n_edges": n_edges,
        "shell_sizes": shell_sizes,
        "expected_shell_sizes": EXPECTED_SHELL_SIZES,
        "antipode_check": "passed",
        "variants": [r.__dict__ for r in results],
        "shell_x_positions": shell_x_positions(n_shells).tolist(),
        "provenance": PROVENANCE_VFD,
    }


def main() -> None:
    res = run()
    print(f"600-cell built: {res['n_vertices']} vertices, {res['n_edges']} edges")
    print(f"Shell sizes: {res['shell_sizes']} (expected {res['expected_shell_sizes']})")
    print(f"Shell x positions: {[round(x, 3) for x in res['shell_x_positions']]}")
    print()
    print(res["models"][["model", "k_params", "chi2", "amplitude", "delta_aic_vs_FREE_C9",
                          "delta_aic_vs_KAPPA_EXP", "correlation_with_exp",
                          "correlation_with_cos", "eigenvalue"]].to_string(index=False))
    print()
    for v in res["variants"]:
        print(f"{v['name']:32s} shell-mean psi = {[round(x, 3) for x in v['shell_psi']]}")


if __name__ == "__main__":
    main()
