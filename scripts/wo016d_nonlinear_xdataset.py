"""WO-016d — Non-linear cross-dataset and cross-channel refit.

Critical: the linearised Mode-B drift diagnostic (wo016c) showed that
the linearised response gives a materially different result from a
non-linear flavio refit on LHCb 2025 (Delta AIC drift = +2.77 AIC
units). The paper headlines must therefore be re-derived in non-linear
mode across all five fits.

This script re-fits FREE_C9 and VFD_GREEN_600CELL on each of the five
datasets (the four B -> K* exclusive sets plus the Bs -> phi
cross-channel) using `flavio.np_prediction` directly, with no Taylor
expansion.

Optimiser: scipy.optimize.minimize_scalar with brent method, sufficient
for 1-D problems and avoids issues with bound-pinning.

Output: reports/wo016d_nonlinear_xdataset.{csv,md}
"""
from __future__ import annotations

import csv
import pathlib
import time
import warnings

warnings.simplefilter("ignore")

import numpy as np
from scipy.optimize import minimize_scalar

from vfd_b_anomaly import wo010_universality
from vfd_b_anomaly.likelihood import aic, bic, chi2 as chi2_fn
from vfd_b_anomaly.wo014_cross_dataset import (
    load_lhcb_2015, load_lhcb_2021_kstplus, load_lhcb_2025,
    load_cms_2025_no_p4p,
)
from vfd_b_anomaly.wo015_cross_channel import load_bs_phimumu_2015


REPO = pathlib.Path(__file__).resolve().parents[1]
OUT_CSV = REPO / "reports" / "wo016d_nonlinear_xdataset.csv"
OUT_MD = REPO / "reports" / "wo016d_nonlinear_xdataset.md"


OBS_TO_FLAVIO = {
    "P5p": "<P5p>", "P4p": "<P4p>", "P1": "<P1>", "P2": "<P2>",
    "FL": "<FL>", "S3": "<S3>", "S4": "<S4>", "S7": "<S7>",
}


DECAY_TO_FLAVIO = {
    "B0->K*mumu":   "B0->K*mumu",
    "B+->K*+mumu":  "B+->K*mumu",
    "Bs->phimumu":  "Bs->phimumu",
}


def chi2_with(values, pred, errors, cov):
    if cov is not None:
        try:
            return chi2_fn(values, pred, covariance=0.5 * (cov + cov.T))
        except np.linalg.LinAlgError:
            eps = 1e-9 * np.trace(cov) / cov.shape[0]
            return chi2_fn(values, pred, covariance=cov + eps * np.eye(cov.shape[0]))
    return chi2_fn(values, pred, errors=errors)


def fit_dataset_nonlinear(data, *, free_bracket=(-3.0, 0.5),
                          vfd_bracket=(0.0, 8.0), tol=5e-4):
    import flavio
    from wilson import Wilson

    df = data["observables"]
    decay = DECAY_TO_FLAVIO[data["decay"]]
    obs = df["observable"].tolist()
    q2_lo = df["q2_lo"].tolist()
    q2_hi = df["q2_hi"].tolist()
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 +
                      df["syst_err"].to_numpy() ** 2)
    cov = data.get("covariance")
    n = len(values)

    full = lambda o: f"{OBS_TO_FLAVIO[o]}({decay})"

    sm = np.array([
        float(flavio.sm_prediction(full(o), q2min=lo, q2max=hi))
        for o, lo, hi in zip(obs, q2_lo, q2_hi)
    ])
    chi2_sm = chi2_with(values, sm, errors, cov)

    def chi2_free(dC9: float) -> float:
        wc = Wilson({"C9_bsmumu": float(dC9)}, scale=4.8,
                    eft="WET", basis="flavio")
        pred = np.array([
            float(flavio.np_prediction(full(o), wc, q2min=lo, q2max=hi))
            for o, lo, hi in zip(obs, q2_lo, q2_hi)
        ])
        return chi2_with(values, pred, errors, cov)

    rfree = minimize_scalar(chi2_free, bracket=free_bracket,
                             method="brent", options={"xtol": tol})
    dC9_free = float(rfree.x)
    chi2_free_val = float(rfree.fun)

    q2_centres = 0.5 * (np.array(q2_lo) + np.array(q2_hi))
    kappa = wo010_universality.frozen_kernel_at_bin_centres(q2_centres)
    dC9_perbin = lambda A: -A * kappa

    def chi2_vfd(A: float) -> float:
        dperb = dC9_perbin(A)
        pred = np.zeros(n)
        for i, (o, lo, hi) in enumerate(zip(obs, q2_lo, q2_hi)):
            wc_i = Wilson({"C9_bsmumu": float(dperb[i])}, scale=4.8,
                          eft="WET", basis="flavio")
            pred[i] = float(flavio.np_prediction(full(o), wc_i,
                                                  q2min=lo, q2max=hi))
        return chi2_with(values, pred, errors, cov)

    rvfd = minimize_scalar(chi2_vfd, bracket=vfd_bracket,
                            method="brent", options={"xtol": tol})
    A = float(rvfd.x)
    chi2_vfd_val = float(rvfd.fun)
    dC9_eff_mean = float(np.mean(dC9_perbin(A)))

    aic_free = aic(chi2_free_val, 1)
    aic_vfd = aic(chi2_vfd_val, 1)
    bic_free = bic(chi2_free_val, 1, n)
    bic_vfd = bic(chi2_vfd_val, 1, n)

    return {
        "dataset": data["dataset"],
        "decay": data["decay"],
        "n_data": n,
        "chi2_sm": chi2_sm,
        "free": {
            "k_params": 1, "chi2": chi2_free_val,
            "aic": aic_free, "bic": bic_free,
            "dC9": dC9_free,
        },
        "vfd": {
            "k_params": 1, "chi2": chi2_vfd_val,
            "aic": aic_vfd, "bic": bic_vfd,
            "A": A, "dC9_eff_mean": dC9_eff_mean,
        },
        "delta_aic_vfd_minus_free": aic_vfd - aic_free,
    }


LOADERS = [
    ("LHCb-2015",         load_lhcb_2015,         (-3.0, 0.5), (0.0, 8.0)),
    ("LHCb-2021-Kstplus", load_lhcb_2021_kstplus, (-3.0, 0.5), (0.0, 8.0)),
    ("CMS-2025-noP4p",    load_cms_2025_no_p4p,   (-3.0, 0.5), (0.0, 8.0)),
    ("LHCb-2025",         load_lhcb_2025,         (-3.0, 0.5), (0.0, 8.0)),
    ("Bs2phi-LHCb-2015",  load_bs_phimumu_2015,   (-6.0, 0.5), (0.0, 20.0)),
]


def main():
    rows = []
    for ds, loader, fb, vb in LOADERS:
        t0 = time.time()
        data = loader()
        # Some loaders set their own dataset name; force consistency
        data["dataset"] = ds
        try:
            result = fit_dataset_nonlinear(data, free_bracket=fb, vfd_bracket=vb)
            rows.append(result)
            print(f"[{time.time()-t0:6.1f}s] {ds}: "
                  f"FREE chi2={result['free']['chi2']:.3f} "
                  f"DC9={result['free']['dC9']:+.3f} | "
                  f"VFD chi2={result['vfd']['chi2']:.3f} "
                  f"A={result['vfd']['A']:+.3f} | "
                  f"DAIC={result['delta_aic_vfd_minus_free']:+.3f}")
        except Exception as exc:
            print(f"[FAIL] {ds}: {exc}")
            rows.append({"dataset": ds, "error": str(exc)})

    with OUT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "dataset", "decay", "model", "n_data", "k_params",
            "chi2", "aic", "bic", "delta_aic_vs_FREE_C9",
            "A_or_DC9", "DC9_eff_mean", "evaluation",
        ])
        for r in rows:
            if "error" in r:
                continue
            w.writerow([
                r["dataset"], r["decay"], "FREE_C9", r["n_data"],
                r["free"]["k_params"], r["free"]["chi2"],
                r["free"]["aic"], r["free"]["bic"], 0.0,
                r["free"]["dC9"], r["free"]["dC9"], "nonlinear",
            ])
            w.writerow([
                r["dataset"], r["decay"], "VFD_GREEN_600CELL", r["n_data"],
                r["vfd"]["k_params"], r["vfd"]["chi2"],
                r["vfd"]["aic"], r["vfd"]["bic"],
                r["delta_aic_vfd_minus_free"],
                r["vfd"]["A"], r["vfd"]["dC9_eff_mean"], "nonlinear",
            ])

    lines = [
        "# WO-016d — Non-linear cross-dataset refit",
        "",
        "Re-runs the WO-014 cross-dataset and WO-015 cross-channel fits "
        "with `flavio.np_prediction` directly (non-linear) instead of the "
        "Mode-B Taylor expansion. The linearised values from "
        "reports/wo014_cross_dataset.csv and reports/wo015_cross_channel.csv "
        "are quoted in parentheses for comparison.",
        "",
        "| dataset | non-linear χ² (FREE) | non-linear χ² (VFD) | "
        "ΔAIC (NL) | ΔC9 (NL) | A (NL) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        if "error" in r:
            lines.append(f"| {r['dataset']} | ERROR: {r['error']} | | | | |")
            continue
        lines.append(
            f"| {r['dataset']} "
            f"| {r['free']['chi2']:.3f} "
            f"| {r['vfd']['chi2']:.3f} "
            f"| {r['delta_aic_vfd_minus_free']:+.3f} "
            f"| {r['free']['dC9']:+.3f} "
            f"| {r['vfd']['A']:+.3f} |"
        )
    lines.append("")

    OUT_MD.write_text("\n".join(lines) + "\n")
    print()
    print(OUT_MD.read_text())


if __name__ == "__main__":
    main()
