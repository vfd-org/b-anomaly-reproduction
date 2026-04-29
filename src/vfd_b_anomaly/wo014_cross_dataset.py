"""WO-014 — Cross-dataset kernel validation.

Goal:
    Test whether the FROZEN VFD_GREEN_600CELL kernel from WO-009 generalises
    across multiple independent datasets without any per-dataset retuning.

Datasets:
    1. LHCb 2015 — arXiv:1512.04442 (3 fb^-1, B0->K*0 mumu) — original P5'
       anomaly. HEPData ins1409497, Table 4 (P-basis CP-averaged).
    2. LHCb 2021 — arXiv:2012.13241 (B+ -> K*+ mumu, isospin partner).
       HEPData ins1838196, data2.yaml + corr_P_bin*.yaml.
    3. CMS 2025 — arXiv:2410.18247 (B0 -> K*0 mumu @ 13 TeV).
       HEPData ins2850101, results_p*.yaml + correlation_matrix_q2_bin_*.yaml.
    4. LHCb 2025 — current dataset (8.4 fb^-1, B0 -> K*0 mumu) — reference.
       Already in `data/raw/hepdata/extracted/config_2_*.yaml`.

Rules:
    - kappa(q^2) is FROZEN (WO-009 600-cell Green's response, projected via
      shell-mean to bin centres). Same kernel function for every dataset.
    - Only the amplitude A is fitted per dataset.
    - SM predictions and dO/dC9 slopes are computed per-dataset on its own
      bin grid via flavio (using the FlavioPredictor cache).
    - For B+->K*+ the flavio decay string is B+->K*mumu (charged isospin partner).

Acceptance:
    - A remains the same SIGN across all datasets (no flips).
    - A remains within an order of magnitude across datasets.
    - VFD_GREEN_600CELL AIC <= FREE_C9 AIC on every dataset (or comparable).
    - Effective DC9 mean is negative and consistent across datasets.
"""

from __future__ import annotations

import warnings
warnings.simplefilter("ignore")

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import minimize

from . import hepdata_ingest, vfd_closure, wo009_full_lift, wo010_universality
from .constants import C9_SM, PHI, PROVENANCE_VFD
from .flavio_predictor import default as flavio_default
from .likelihood import aic, bic, chi2 as chi2_fn


# ---------------------------------------------------------------------------
# Dataset loaders
# ---------------------------------------------------------------------------

LEGACY_DIR = Path("data/raw/hepdata_legacy")
CONFIG_2_DIR = Path("data/raw/hepdata/extracted")


def _err_to_sigma(err_entry: dict) -> tuple[float, float]:
    """Return (stat_sigma, syst_sigma) from a HEPData errors list. Asymmetric
    stat is symmetrised by taking the larger of |minus| and |plus|."""
    stat = 0.0
    syst = 0.0
    for e in err_entry.get("errors", []):
        label = e.get("label", "").lower()
        if "asymerror" in e:
            v = max(abs(e["asymerror"].get("minus", 0.0)),
                    abs(e["asymerror"].get("plus", 0.0)))
        elif "symerror" in e:
            v = abs(e["symerror"])
        else:
            v = 0.0
        if label == "stat":
            stat = v
        elif label in ("sys", "syst", "systematic"):
            syst = v
        else:
            # treat unlabeled as combined
            stat = max(stat, v)
    return stat, syst


def load_lhcb_2015(observables=("P5p", "P4p", "P1", "P2")) -> dict[str, Any]:
    """LHCb 2015 (3 fb^-1) Table 4 P-basis CP-averaged angular observables."""
    yaml_path = LEGACY_DIR / "ins1409497/HEPData-ins1409497-v1-yaml/Table4.yaml"
    d = yaml.safe_load(yaml_path.read_text())
    bins = d["independent_variables"][0]["values"]
    name_map = {"$P_{1}$": "P1", "$P_{2}$": "P2", "$P_{3}$": "P3",
                "$P'_{4}$": "P4p", "$P'_{5}$": "P5p", "$P'_{6}$": "P6p",
                "$P'_{8}$": "P8p"}
    rows = []
    for v in d["dependent_variables"]:
        obs_name = name_map.get(v["header"]["name"], v["header"]["name"])
        if obs_name not in observables:
            continue
        for bi, val in enumerate(v["values"]):
            stat, syst = _err_to_sigma(val)
            rows.append({
                "observable": obs_name,
                "q2_lo": float(bins[bi]["low"]),
                "q2_hi": float(bins[bi]["high"]),
                "value": float(val["value"]),
                "stat_err": stat,
                "syst_err": syst,
            })
    df = pd.DataFrame(rows)
    return {
        "dataset": "LHCb-2015",
        "decay": "B0->K*mumu",
        "n_obs": len(df),
        "covariance": None,
        "observables": df,
        "metadata": {"arxiv": "1512.04442", "luminosity_fb": 3.0,
                     "hepdata": "ins1409497", "table": "Table4"},
    }


def load_lhcb_2021_kstplus(observables=("P5p", "P4p", "P1", "P2")) -> dict[str, Any]:
    """LHCb 2021 (B^+ -> K*+ mumu) data2.yaml + per-bin correlation blocks.

    The 10 q^2 bins include 2 wide combined bins (1.1-6, 15-19) that overlap
    the exclusive ones; we keep only the 8 exclusive bins to avoid double-
    counting. Covariance is block-diagonal: per-bin 7x7 across observables
    (Fl, P1..P8 minus P3 because P3 is dropped for the k=4 angular subset),
    no cross-bin correlations exposed by HEPData.
    """
    base = LEGACY_DIR / "ins1838196/HEPData-ins1838196-v1-yaml"
    d = yaml.safe_load((base / "data2.yaml").read_text())
    bins = d["independent_variables"][0]["values"]
    headers = [v["header"]["name"] for v in d["dependent_variables"]]
    # Map flavio-style "P5p" -> file's "P5"
    name_to_file = {"P1": "P1", "P2": "P2", "P3": "P3",
                    "P4p": "P4", "P5p": "P5", "P6p": "P6", "P8p": "P8",
                    "FL": "Fl"}
    file_to_name = {v: k for k, v in name_to_file.items()}

    # Identify exclusive 8 bins (skip wide combined ones)
    keep_bin_idx = []
    for bi, b in enumerate(bins):
        lo, hi = float(b["low"]), float(b["high"])
        # combined bins are [1.1, 6] and [15, 19] in this file
        if abs(lo - 1.1) < 1e-3 and abs(hi - 6.0) < 1e-3:
            continue
        if abs(lo - 15.0) < 1e-3 and abs(hi - 19.0) < 1e-3:
            continue
        keep_bin_idx.append(bi)

    rows = []
    for var in d["dependent_variables"]:
        file_name = var["header"]["name"]
        obs_name = file_to_name.get(file_name)
        if obs_name is None or obs_name not in observables:
            continue
        for bi in keep_bin_idx:
            val = var["values"][bi]
            stat, syst = _err_to_sigma(val)
            rows.append({
                "observable": obs_name,
                "q2_lo": float(bins[bi]["low"]),
                "q2_hi": float(bins[bi]["high"]),
                "value": float(val["value"]),
                "stat_err": stat,
                "syst_err": syst,
            })
    df = pd.DataFrame(rows).sort_values(["q2_lo", "observable"]).reset_index(drop=True)
    cov = _build_block_diag_corr_lhcb_2021(df, base, file_to_name, observables)
    return {
        "dataset": "LHCb-2021-Kstplus",
        "decay": "B+->K*+mumu",
        "n_obs": len(df),
        "covariance": cov,
        "observables": df,
        "metadata": {"arxiv": "2012.13241", "luminosity_fb": 9.0,
                     "hepdata": "ins1838196", "table": "data2 + corr_P_bin*"},
    }


def _build_block_diag_corr_lhcb_2021(df: pd.DataFrame, base: Path,
                                      file_to_name: dict, observables: tuple) -> np.ndarray:
    """Per-bin 8x8 correlation block from corr_P_bin*.yaml -> block-diagonal
    covariance over the kept rows. Uses (stat^2 + syst^2) for the sigmas.

    Requires `df.index == range(len(df))` because off-diagonal entries are
    written via `r.name`, the original integer position. The caller resets
    the index before invoking us; this assert guards future regressions.
    """
    if not (df.index == np.arange(len(df))).all():
        raise ValueError("LHCb-2021 covariance build requires df with reset_index(drop=True)")
    n = len(df)
    sigma = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    cov = np.diag(sigma ** 2)  # initialise diagonal

    # Map exclusive bin idx (in the kept order) back to corr_P_bin{global}.yaml
    bins_kept = sorted(set(zip(df["q2_lo"], df["q2_hi"])))
    # Need the GLOBAL bin index in the source file
    src_d = yaml.safe_load((base / "data2.yaml").read_text())
    src_bins = src_d["independent_variables"][0]["values"]
    src_idx_for = {(float(b["low"]), float(b["high"])): bi
                   for bi, b in enumerate(src_bins)}
    src_headers = [v["header"]["name"] for v in src_d["dependent_variables"]]
    name_to_src_idx = {v["header"]["name"]: i
                       for i, v in enumerate(src_d["dependent_variables"])}

    for q2_lo, q2_hi in bins_kept:
        gbi = src_idx_for[(q2_lo, q2_hi)]
        corr_file = base / f"corr_P_bin{gbi}.yaml"
        if not corr_file.exists():
            continue
        corr_data = yaml.safe_load(corr_file.read_text())
        corr_vals = corr_data["dependent_variables"][0]["values"]
        # 64 entries = 8x8 over [Fl, P1, P2, P3, P4, P5, P6, P8]
        # But the corr file's row/col labels are in independent_variables.
        idep = corr_data["independent_variables"]
        row_labels = [v["value"] for v in idep[0]["values"]]
        col_labels = [v["value"] for v in idep[1]["values"]]
        # Build 8x8 block
        n_block = int(np.sqrt(len(corr_vals)))
        block_corr = np.zeros((n_block, n_block))
        # Pull labels in order
        block_label_order = []
        seen = set()
        for lab in row_labels:
            if lab not in seen:
                block_label_order.append(lab)
                seen.add(lab)
        for k, val in enumerate(corr_vals):
            i = block_label_order.index(row_labels[k])
            j = block_label_order.index(col_labels[k])
            block_corr[i, j] = float(val["value"])

        # For each pair of obs in our selected `observables` list, look up
        # row/col indices in the `df` and the corresponding corr matrix index.
        sel_rows = df[(df["q2_lo"] == q2_lo) & (df["q2_hi"] == q2_hi)]
        for _, r1 in sel_rows.iterrows():
            for _, r2 in sel_rows.iterrows():
                src1 = next((k for k, v in file_to_name.items() if v == r1["observable"]), None)
                src2 = next((k for k, v in file_to_name.items() if v == r2["observable"]), None)
                if src1 is None or src2 is None:
                    continue
                if src1 not in block_label_order or src2 not in block_label_order:
                    continue
                bi1 = block_label_order.index(src1)
                bi2 = block_label_order.index(src2)
                rho = block_corr[bi1, bi2]
                df_i = r1.name
                df_j = r2.name
                if df_i == df_j:
                    continue
                cov[df_i, df_j] = rho * sigma[df_i] * sigma[df_j]
                cov[df_j, df_i] = cov[df_i, df_j]
    return cov


def load_cms_2025(observables=("P5p", "P4p", "P1", "P2")) -> dict[str, Any]:
    """CMS 2024/2025 results_*.yaml angular analysis (13 TeV, B0->K*0 mumu)."""
    base = LEGACY_DIR / "ins2850101/HEPData-ins2850101-v1-yaml"
    obs_to_file = {"P1": "results_p1.yaml", "P2": "results_p2.yaml",
                   "P3": "results_p3.yaml", "P4p": "results_p4p.yaml",
                   "P5p": "results_p5p.yaml", "P6p": "results_p6p.yaml",
                   "P8p": "results_p8p.yaml", "FL": "results_fl.yaml"}
    rows = []
    for obs in observables:
        if obs not in obs_to_file:
            continue
        d = yaml.safe_load((base / obs_to_file[obs]).read_text())
        bins = d["independent_variables"][0]["values"]
        # First dependent_variable is the measurement; others are theory predictions.
        meas = d["dependent_variables"][0]
        for bi, val in enumerate(meas["values"]):
            stat, syst = _err_to_sigma(val)
            rows.append({
                "observable": obs,
                "q2_lo": float(bins[bi]["low"]),
                "q2_hi": float(bins[bi]["high"]),
                "value": float(val["value"]),
                "stat_err": stat,
                "syst_err": syst,
            })
    df = pd.DataFrame(rows).sort_values(["q2_lo", "observable"]).reset_index(drop=True)
    cov = _build_block_diag_corr_cms(df, base, list(observables))
    return {
        "dataset": "CMS-2025",
        "decay": "B0->K*mumu",
        "n_obs": len(df),
        "covariance": cov,
        "observables": df,
        "metadata": {"hepdata": "ins2850101", "tables": "results_p*.yaml",
                     "note": "13 TeV, full Run 2"},
    }


def _build_block_diag_corr_cms(df: pd.DataFrame, base: Path,
                                observables: list[str]) -> np.ndarray | None:
    """CMS publishes correlation_matrix_q2_bin_{i}.yaml for some bins. Build
    block-diagonal covariance over the requested observables.

    Requires `df.index == range(len(df))` because off-diagonal entries are
    written via `r.name`, the original integer position. The caller resets
    the index before invoking us; this assert guards future regressions.
    """
    if not (df.index == np.arange(len(df))).all():
        raise ValueError("CMS covariance build requires df with reset_index(drop=True)")
    n = len(df)
    sigma = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    cov = np.diag(sigma ** 2)

    # CMS bins index: q2_bin labels are 0,1,2,3,5,7 (skipping 4 and 6 = J/psi
    # and psi(2S) regions). Map each unique (lo,hi) -> file label.
    bins_present = sorted(set(zip(df["q2_lo"], df["q2_hi"])))
    # The bin numbering inside the YAML files is the global CMS bin index
    # (0..7 with 4 and 6 reserved for charm resonances); we infer from order.
    cms_full_bins = [
        (1.1, 2.0), (2.0, 4.3), (4.3, 6.0), (6.0, 8.68),
        (8.68, 10.09), (10.09, 12.86), (12.86, 14.18), (14.18, 16.0),
    ]
    bin_to_global = {}
    for (lo, hi) in bins_present:
        for gi, (glo, ghi) in enumerate(cms_full_bins):
            if abs(lo - glo) < 1e-3 and abs(hi - ghi) < 1e-3:
                bin_to_global[(lo, hi)] = gi
                break

    for (lo, hi), gi in bin_to_global.items():
        corr_path = base / f"correlation_matrix_q2_bin_{gi}.yaml"
        if not corr_path.exists():
            continue
        cd = yaml.safe_load(corr_path.read_text())
        # rows / cols labelled by observable name
        idep = cd["independent_variables"]
        row_labels = [v["value"] for v in idep[0]["values"]]
        col_labels = [v["value"] for v in idep[1]["values"]]
        block_label_order = []
        seen = set()
        for lab in row_labels:
            if lab not in seen:
                block_label_order.append(lab)
                seen.add(lab)
        n_block = len(block_label_order)
        block = np.zeros((n_block, n_block))
        for k, val in enumerate(cd["dependent_variables"][0]["values"]):
            i = block_label_order.index(row_labels[k])
            j = block_label_order.index(col_labels[k])
            block[i, j] = float(val["value"])

        # Map block labels (e.g. "P5p", "P1", "FL") back to our observables
        # CMS uses things like "$P_5'$" or "P5p"; try a tolerant lookup.
        def _normalise_label(s):
            return s.replace("$", "").replace("\\", "").replace("'", "p").replace("_", "").replace("{", "").replace("}", "").lower()

        block_norm = [_normalise_label(s) for s in block_label_order]
        sel_rows = df[(df["q2_lo"] == lo) & (df["q2_hi"] == hi)]
        for _, r1 in sel_rows.iterrows():
            for _, r2 in sel_rows.iterrows():
                if r1.name == r2.name:
                    continue
                k1 = _normalise_label(r1["observable"])
                k2 = _normalise_label(r2["observable"])
                if k1 not in block_norm or k2 not in block_norm:
                    continue
                i = block_norm.index(k1); j = block_norm.index(k2)
                rho = block[i, j]
                cov[r1.name, r2.name] = rho * sigma[r1.name] * sigma[r2.name]
                cov[r2.name, r1.name] = cov[r1.name, r2.name]
    return cov


def load_lhcb_2025(observables=("P5p", "P4p", "P1", "P2")) -> dict[str, Any]:
    """LHCb 2025 8.4 fb^-1 — already in `data/raw/hepdata/extracted` config_2."""
    archive = hepdata_ingest.hepdata_archive_dir(CONFIG_2_DIR)
    data = hepdata_ingest.load_config(archive, config_index=2,
                                      observables=tuple(observables))
    return {
        "dataset": "LHCb-2025",
        "decay": "B0->K*mumu",
        "n_obs": len(data["observables"]),
        "covariance": np.asarray(data["covariance"]) if data.get("covariance") is not None else None,
        "observables": data["observables"],
        "metadata": {"arxiv": "2512.18053", "luminosity_fb": 8.4,
                     "hepdata": "ins3094698", "config": "config_2"},
    }


# ---------------------------------------------------------------------------
# Generic universality fit (any q^2 grid, any decay)
# ---------------------------------------------------------------------------

@dataclass
class FitResult:
    dataset: str
    decay: str
    model: str
    n_data: int
    k_params: int
    chi2: float
    aic: float
    bic: float
    delta_aic_vs_FREE_C9: float
    A: float
    DC9_eff_mean: float
    notes: str = ""


def _chi2_with(values, pred, errors, cov):
    if cov is not None:
        try:
            return chi2_fn(values, pred, covariance=0.5 * (cov + cov.T))
        except np.linalg.LinAlgError:
            eps = 1e-9 * np.trace(cov) / cov.shape[0]
            return chi2_fn(values, pred, covariance=cov + eps * np.eye(cov.shape[0]))
    return chi2_fn(values, pred, errors=errors)


def fit_dataset(data: dict[str, Any]) -> dict[str, FitResult]:
    """Fit FREE_C9 and VFD_GREEN_600CELL on a single dataset. Returns both."""
    df = data["observables"]
    decay = data["decay"]
    obs = df["observable"].tolist()
    q2_lo = df["q2_lo"].tolist()
    q2_hi = df["q2_hi"].tolist()
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    cov = data.get("covariance")
    n = len(values)

    fp = flavio_default()
    # Pre-warm cache
    bins_seen = set(zip(q2_lo, q2_hi))
    obs_seen = set(obs)
    fp.precompute(obs_seen, bins_seen, decay=decay, verbose=False)

    # Pre-compute SM and slope vectors
    sm_vec = np.array([fp.sm_value(o, lo, hi, decay=decay)
                       for o, lo, hi in zip(obs, q2_lo, q2_hi)])
    slope_vec = np.array([fp.dO_dC9(o, lo, hi, decay=decay)
                          for o, lo, hi in zip(obs, q2_lo, q2_hi)])

    # FREE_C9 fit
    def loss_free(theta):
        dC9 = float(theta[0])
        pred = sm_vec + slope_vec * dC9
        return _chi2_with(values, pred, errors, cov)

    r_free = minimize(loss_free, x0=[-1.0], method="Powell", bounds=[(-3.0, 3.0)],
                      options={"xtol": 1e-7, "ftol": 1e-9})
    chi2_free = float(r_free.fun)
    DC9_free = float(r_free.x[0])
    aic_free = aic(chi2_free, 1)
    bic_free = bic(chi2_free, 1, n)

    # VFD shared-kernel fit (frozen kernel)
    q2_centres = 0.5 * (np.array(q2_lo) + np.array(q2_hi))
    kappa = wo010_universality.frozen_kernel_at_bin_centres(q2_centres)

    def loss_vfd(theta):
        a = float(theta[0])
        pred = sm_vec + slope_vec * (-a * kappa)
        return _chi2_with(values, pred, errors, cov)

    r_vfd = minimize(loss_vfd, x0=[1.0], method="Powell", bounds=[(-5.0, 5.0)],
                     options={"xtol": 1e-7, "ftol": 1e-9})
    chi2_vfd = float(r_vfd.fun)
    A_vfd = float(r_vfd.x[0])
    aic_vfd = aic(chi2_vfd, 1)
    bic_vfd = bic(chi2_vfd, 1, n)
    delta_grid = -A_vfd * kappa

    fp.flush()

    return {
        "FREE_C9": FitResult(
            dataset=data["dataset"], decay=decay, model="FREE_C9",
            n_data=n, k_params=1, chi2=chi2_free, aic=aic_free, bic=bic_free,
            delta_aic_vs_FREE_C9=0.0, A=DC9_free,
            DC9_eff_mean=DC9_free, notes=f"DC9={DC9_free:+.3f}",
        ),
        "VFD": FitResult(
            dataset=data["dataset"], decay=decay, model="VFD_GREEN_600CELL",
            n_data=n, k_params=1, chi2=chi2_vfd, aic=aic_vfd, bic=bic_vfd,
            delta_aic_vs_FREE_C9=aic_vfd - aic_free, A=A_vfd,
            DC9_eff_mean=float(np.mean(delta_grid)),
            notes=f"A={A_vfd:+.4f}, DC9_eff={float(np.mean(delta_grid)):+.3f}",
        ),
    }


# ---------------------------------------------------------------------------
# Top-level run
# ---------------------------------------------------------------------------

def load_cms_2025_no_p4p(observables=None) -> dict[str, Any]:
    """CMS 2025 without P4'. CMS publishes P4' values that fall outside the
    [-1, +1] physical range expected by flavio's convention (e.g. P4'(14.18,16)
    = -1.159), suggesting a different normalisation convention. Dropping P4'
    isolates the convention question from the universality test. Always uses
    (P5p, P1, P2); the `observables` argument is intentionally ignored."""
    data = load_cms_2025(observables=("P5p", "P1", "P2"))
    data["dataset"] = "CMS-2025-noP4p"
    data["metadata"]["note"] += "; P4' dropped (convention mismatch)"
    return data


LOADERS = {
    "LHCb-2015":         load_lhcb_2015,
    "LHCb-2021-Kstplus": load_lhcb_2021_kstplus,
    "CMS-2025":          load_cms_2025,
    "CMS-2025-noP4p":    load_cms_2025_no_p4p,
    "LHCb-2025":         load_lhcb_2025,
}


def run(*, observables=("P5p", "P4p", "P1", "P2"),
        output_dir: Path | str = "reports") -> pd.DataFrame:
    rows = []
    for name, loader in LOADERS.items():
        print(f">> loading {name} ...")
        try:
            data = loader(observables=observables)
        except Exception as e:
            print(f"   skip {name}: {type(e).__name__}: {e}")
            continue
        n = data["n_obs"]
        bins = sorted(set(zip(data["observables"]["q2_lo"], data["observables"]["q2_hi"])))
        print(f"   n_data={n}, bins={len(bins)}, decay={data['decay']}")
        print(f">> fitting {name} ...")
        fits = fit_dataset(data)
        for tag, fr in fits.items():
            rows.append({
                "dataset": fr.dataset,
                "decay": fr.decay,
                "model": fr.model,
                "n_data": fr.n_data,
                "k_params": fr.k_params,
                "chi2": fr.chi2,
                "aic": fr.aic,
                "bic": fr.bic,
                "delta_aic_vs_FREE_C9": fr.delta_aic_vs_FREE_C9,
                "A_or_DC9": fr.A,
                "DC9_eff_mean": fr.DC9_eff_mean,
                "notes": fr.notes,
            })
    df = pd.DataFrame(rows)
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "wo014_cross_dataset.csv", index=False)
    return df


def main() -> None:
    pd.set_option("display.width", 160)
    pd.set_option("display.max_colwidth", 60)
    df = run()
    print()
    print("=" * 100)
    print("WO-014 — Cross-dataset frozen-kernel universality")
    print("=" * 100)
    print(df.to_string(index=False))

    print()
    print("Per-dataset summary (VFD vs FREE_C9):")
    pivot = df.pivot_table(index=["dataset", "decay"], columns="model",
                           values=["chi2", "aic", "A_or_DC9"])
    print(pivot.to_string())

    # Universality verdict
    vfd = df[df["model"] == "VFD_GREEN_600CELL"]
    print()
    print(f"VFD A range across {len(vfd)} datasets: "
          f"[{vfd['A_or_DC9'].min():+.3f}, {vfd['A_or_DC9'].max():+.3f}]")
    print(f"All VFD A > 0?  {(vfd['A_or_DC9'] > 0).all()}")
    print(f"All VFD DC9_eff < 0?  {(vfd['DC9_eff_mean'] < 0).all()}")
    print(f"All VFD AIC <= FREE_C9?  "
          f"{(vfd['delta_aic_vs_FREE_C9'] <= 0).all()}")


if __name__ == "__main__":
    main()
