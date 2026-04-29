"""WO-016c — Non-linear flavio refit on LHCb 2025.

Reviewer concern: the linearised response model (Mode B) extrapolates a
central-difference slope at delta C9 = +/- 0.5 to the fitted shift
~ -1.4. At |Delta C9| ~ 1.4 the linearisation is approximate; the
question is whether the headline Delta AIC (-1.67) survives a fully
non-linear flavio refit.

Method:
    1. Load LHCb 2025 4-observable joint dataset (32 points).
    2. Evaluate both models in non-linear mode at the linearised
       best-fit point (drift diagnostic).
    3. Re-fit both FREE_C9 (one global Delta C9) and VFD (one A;
       per-bin Delta C9_eff = -A * kappa) using `flavio.np_prediction`
       directly. Compare non-linear best-fit Delta AIC to the
       linearised value (-1.67).
    4. Report whether the headline Delta AIC survives non-linear
       evaluation.

Output: reports/wo016c_nonlinear_refit.{csv,md}
"""
from __future__ import annotations

import csv
import pathlib
import warnings

warnings.simplefilter("ignore")

import numpy as np
from scipy.optimize import minimize_scalar

from vfd_b_anomaly import wo010_universality
from vfd_b_anomaly.flavio_predictor import default as flavio_default
from vfd_b_anomaly.likelihood import aic, chi2 as chi2_fn
from vfd_b_anomaly.wo014_cross_dataset import load_lhcb_2025


REPO = pathlib.Path(__file__).resolve().parents[1]
OUT_CSV = REPO / "reports" / "wo016c_nonlinear_refit.csv"
OUT_MD = REPO / "reports" / "wo016c_nonlinear_refit.md"

# Linearised best-fit values for the LHCb-2025 row (the linearised
# cross-dataset run that produced these is now archived; the values are
# pinned here so this drift diagnostic remains reproducible from a clean
# repo without rerunning the linearised fit).
DC9_FREE_BEST = -1.340331410060973
A_VFD_BEST = 1.5935042407680777


def chi2_with(values, pred, errors, cov):
    if cov is not None:
        try:
            return chi2_fn(values, pred, covariance=0.5 * (cov + cov.T))
        except np.linalg.LinAlgError:
            eps = 1e-9 * np.trace(cov) / cov.shape[0]
            return chi2_fn(values, pred, covariance=cov + eps * np.eye(cov.shape[0]))
    return chi2_fn(values, pred, errors=errors)


def main():
    import flavio
    from wilson import Wilson

    data = load_lhcb_2025()
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

    OBS_TO_FLAVIO = {
        "P5p": "<P5p>",
        "P4p": "<P4p>",
        "P1":  "<P1>",
        "P2":  "<P2>",
    }
    FULL = lambda o: f"{OBS_TO_FLAVIO[o]}(B0->K*mumu)"

    # ---- FREE_C9 non-linear ----
    wc_free = Wilson({"C9_bsmumu": DC9_FREE_BEST}, scale=4.8,
                     eft="WET", basis="flavio")
    pred_free_nl = np.array([
        float(flavio.np_prediction(FULL(o), wc_free, q2min=lo, q2max=hi))
        for o, lo, hi in zip(obs, q2_lo, q2_hi)
    ])
    pred_free_lin = np.array([
        fp.predict(o, lo, hi, c9_value=4.27 + DC9_FREE_BEST,
                   c9_sm=4.27, decay=decay)
        for o, lo, hi in zip(obs, q2_lo, q2_hi)
    ])

    chi2_free_nl = chi2_with(values, pred_free_nl, errors, cov)
    chi2_free_lin = chi2_with(values, pred_free_lin, errors, cov)

    # ---- VFD non-linear ----
    q2_centres = 0.5 * (np.array(q2_lo) + np.array(q2_hi))
    kappa = wo010_universality.frozen_kernel_at_bin_centres(q2_centres)
    dC9_eff_perbin = -A_VFD_BEST * kappa

    pred_vfd_nl = np.zeros(n)
    pred_vfd_lin = np.zeros(n)
    for i, (o, lo, hi) in enumerate(zip(obs, q2_lo, q2_hi)):
        wc_i = Wilson({"C9_bsmumu": float(dC9_eff_perbin[i])}, scale=4.8,
                      eft="WET", basis="flavio")
        pred_vfd_nl[i] = float(flavio.np_prediction(FULL(o), wc_i,
                                                     q2min=lo, q2max=hi))
        pred_vfd_lin[i] = fp.predict(o, lo, hi,
                                     c9_value=4.27 + float(dC9_eff_perbin[i]),
                                     c9_sm=4.27, decay=decay)

    chi2_vfd_nl = chi2_with(values, pred_vfd_nl, errors, cov)
    chi2_vfd_lin = chi2_with(values, pred_vfd_lin, errors, cov)

    # ---- Non-linear refits (find best-fit DC9 / A in non-linear model) ----
    def chi2_free_nonlinear(dC9: float) -> float:
        wc = Wilson({"C9_bsmumu": float(dC9)}, scale=4.8,
                    eft="WET", basis="flavio")
        pred = np.array([
            float(flavio.np_prediction(FULL(o), wc, q2min=lo, q2max=hi))
            for o, lo, hi in zip(obs, q2_lo, q2_hi)
        ])
        return chi2_with(values, pred, errors, cov)

    def chi2_vfd_nonlinear(A: float) -> float:
        dC9_perbin = -A * kappa
        pred = np.zeros(n)
        for i, (o, lo, hi) in enumerate(zip(obs, q2_lo, q2_hi)):
            wc_i = Wilson({"C9_bsmumu": float(dC9_perbin[i])}, scale=4.8,
                          eft="WET", basis="flavio")
            pred[i] = float(flavio.np_prediction(FULL(o), wc_i,
                                                  q2min=lo, q2max=hi))
        return chi2_with(values, pred, errors, cov)

    rfree = minimize_scalar(chi2_free_nonlinear, bracket=(-3.0, 0.0),
                             method="brent", options={"xtol": 1e-4})
    DC9_FREE_NL = float(rfree.x)
    chi2_free_nl_refit = float(rfree.fun)

    rvfd = minimize_scalar(chi2_vfd_nonlinear, bracket=(0.0, 5.0),
                            method="brent", options={"xtol": 1e-4})
    A_VFD_NL = float(rvfd.x)
    chi2_vfd_nl_refit = float(rvfd.fun)

    aic_free_nl = aic(chi2_free_nl, 1)
    aic_vfd_nl = aic(chi2_vfd_nl, 1)
    delta_aic_nl = aic_vfd_nl - aic_free_nl
    aic_free_lin = aic(chi2_free_lin, 1)
    aic_vfd_lin = aic(chi2_vfd_lin, 1)
    delta_aic_lin = aic_vfd_lin - aic_free_lin
    aic_free_nl_refit = aic(chi2_free_nl_refit, 1)
    aic_vfd_nl_refit = aic(chi2_vfd_nl_refit, 1)
    delta_aic_nl_refit = aic_vfd_nl_refit - aic_free_nl_refit

    rows = [
        ("FREE_C9_linear", chi2_free_lin, aic_free_lin, 0.0,
         f"DC9={DC9_FREE_BEST:+.3f}"),
        ("FREE_C9_nonlinear@linear-best-fit", chi2_free_nl, aic_free_nl, 0.0,
         f"DC9={DC9_FREE_BEST:+.3f}"),
        ("FREE_C9_nonlinear_refit", chi2_free_nl_refit, aic_free_nl_refit, 0.0,
         f"DC9={DC9_FREE_NL:+.3f}"),
        ("VFD_linear", chi2_vfd_lin, aic_vfd_lin,
         aic_vfd_lin - aic_free_lin, f"A={A_VFD_BEST:+.3f}"),
        ("VFD_nonlinear@linear-best-fit", chi2_vfd_nl, aic_vfd_nl,
         aic_vfd_nl - aic_free_nl, f"A={A_VFD_BEST:+.3f}"),
        ("VFD_nonlinear_refit", chi2_vfd_nl_refit, aic_vfd_nl_refit,
         aic_vfd_nl_refit - aic_free_nl_refit, f"A={A_VFD_NL:+.3f}"),
    ]

    with OUT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "chi2", "aic", "delta_aic_vs_FREE_C9", "fit_param"])
        for r in rows:
            w.writerow([r[0], f"{r[1]:.4f}", f"{r[2]:.4f}", f"{r[3]:.4f}", r[4]])

    pred_diff_free = np.abs(pred_free_nl - pred_free_lin) / np.abs(errors)
    pred_diff_vfd = np.abs(pred_vfd_nl - pred_vfd_lin) / np.abs(errors)

    lines = [
        "# WO-016c — Non-linear flavio refit on LHCb 2025",
        "",
        "Tests whether the linearised Mode-B response is sufficient at the "
        "fitted Delta C9 ~ -1.4. Three comparisons:",
        "1. Linear fit (paper headline).",
        "2. Non-linear evaluation at the linearised best-fit point "
        "(drift diagnostic — both models pinned at linear best-fit).",
        "3. Non-linear refit (best-fit values found by `flavio.np_prediction`).",
        "",
        "## Headline",
        "",
        "| model | chi^2 | AIC | Delta AIC vs FREE_C9 | fit param |",
        "|---|---:|---:|---:|---|",
    ]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:.3f} | {r[2]:.3f} | "
                     f"{r[3]:+.3f} | {r[4]} |")
    lines.append("")
    lines.append(
        f"- Linearised Delta AIC (FREE_C9 vs VFD): {delta_aic_lin:+.3f}"
    )
    lines.append(
        f"- Non-linear Delta AIC at linear best-fit: {delta_aic_nl:+.3f} "
        f"(diagnostic only; both models held at linear best-fit)"
    )
    lines.append(
        f"- Non-linear Delta AIC after refit: {delta_aic_nl_refit:+.3f} "
        f"(headline-comparable)"
    )
    lines.append(
        f"- Drift in headline Delta AIC: "
        f"{delta_aic_nl_refit - delta_aic_lin:+.3f}"
    )
    lines.append("")
    lines.append("## Best-fit parameters")
    lines.append("")
    lines.append(f"- FREE_C9 linear: Delta C9 = {DC9_FREE_BEST:+.4f}")
    lines.append(f"- FREE_C9 non-linear refit: Delta C9 = {DC9_FREE_NL:+.4f}")
    lines.append(f"- VFD linear: A = {A_VFD_BEST:+.4f}")
    lines.append(f"- VFD non-linear refit: A = {A_VFD_NL:+.4f}")
    lines.append("")
    lines.append("## Per-bin linearisation residual (|nonlinear - linear| / sigma)")
    lines.append("")
    lines.append(
        f"- FREE_C9 at linear best-fit: max = {pred_diff_free.max():.3f} sigma, "
        f"mean = {pred_diff_free.mean():.3f} sigma"
    )
    lines.append(
        f"- VFD at linear best-fit: max = {pred_diff_vfd.max():.3f} sigma, "
        f"mean = {pred_diff_vfd.mean():.3f} sigma"
    )
    lines.append("")
    threshold = 0.5
    drift = delta_aic_nl_refit - delta_aic_lin
    if abs(drift) <= threshold:
        lines.append(
            f"**Conclusion.** |Drift in headline Delta AIC| = "
            f"{abs(drift):.3f} <= {threshold} AIC unit. The linearised "
            f"and non-linear fits agree within tolerance; the headline "
            f"Delta AIC is robust to non-linear flavio evaluation."
        )
    else:
        lines.append(
            f"**Conclusion.** |Drift in headline Delta AIC| = "
            f"{abs(drift):.3f} > {threshold} AIC unit. The headline must be "
            f"updated to the non-linear refit value: "
            f"Delta AIC_NL = {delta_aic_nl_refit:+.3f} (vs linearised "
            f"{delta_aic_lin:+.3f}). The non-linear refit best-fit is "
            f"DC9_FREE = {DC9_FREE_NL:+.3f} (vs linear {DC9_FREE_BEST:+.3f}) "
            f"and A = {A_VFD_NL:+.3f} (vs linear {A_VFD_BEST:+.3f})."
        )

    OUT_MD.write_text("\n".join(lines) + "\n")
    print(OUT_MD.read_text())


if __name__ == "__main__":
    main()
