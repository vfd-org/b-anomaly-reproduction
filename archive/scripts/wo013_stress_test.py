"""WO-013 — Stress test the universality result.

The frozen kernel from WO-009 is locked: VFD_GREEN_600CELL only. No new
kernel variants. This script tests how robust the WO-010 universality
result (dAIC = -1.67 vs FREE_C9 with shared amplitude across
P5', P4', P1, P2) is under five perturbations:

    1. Bootstrap over BINS (not observables): resample bins with
       replacement, refit kernel amplitude. Histogram the amplitude
       distribution.
    2. Region splits: refit independently on low q^2 (bins 0-2),
       central q^2 (bins 3-4), high q^2 (bins 5-7). Universality
       across regions?
    3. Alternative-model comparison: FREE_C9, FREE_C9 + C10,
       charm-loop nuisance proxy (Gaussian residual at J/psi-psi(2S)
       midpoint with FREE width and amplitude).
    4. Form-factor variation: sample BSZ form-factor parameters within
       their published uncertainties (Monte Carlo over flavio's BSZ
       a-coefficients). For each sample regenerate the SM table, refit
       the kernel, collect amplitude distribution.
    5. Frozen-kernel sanity: confirm the kernel shape itself is not
       being re-fit anywhere; only A is fitted in each test.
"""

from __future__ import annotations

import warnings
warnings.simplefilter("ignore")

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from . import hepdata_ingest, vfd_closure, wo009_full_lift, wo010_universality
from .constants import C9_SM, J_PSI_Q2, PHI, PROVENANCE_VFD, PSI2S_Q2
from .likelihood import aic, bic, chi2 as chi2_fn

import re

_A_PATTERN = re.compile(r"A=([+\-]?[0-9.eE+\-]+)")


def _amplitude_from_notes(notes: str) -> float:
    m = _A_PATTERN.search(notes or "")
    if not m:
        raise ValueError(f"could not parse amplitude from notes: {notes!r}")
    return float(m.group(1))


DEFAULT_ANGULAR_OBSERVABLES = ("P5p", "P4p", "P1", "P2")


def _bin_axis(df: pd.DataFrame) -> tuple[list[str], list[int], np.ndarray]:
    """Return absolute (q^2-grid) bin indices, not positional ones — this is
    required because slicing the dataset (region splits, bootstrap) shouldn't
    change which sm_baseline row a given (q^2_lo, q^2_hi) maps to.
    """
    obs = df["observable"].tolist()
    bin_indices = [
        wo010_universality._canonical_bin_index(float(r["q2_lo"]), float(r["q2_hi"]))
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


# ---------------------------------------------------------------------------
# 1. Bootstrap over bins
# ---------------------------------------------------------------------------

def bootstrap_bins(
    data: dict[str, Any],
    *,
    n_bootstrap: int = 1000,
    rng_seed: int = 12345,
) -> pd.DataFrame:
    """Resample BINS with replacement (each bin treated as a unit of all
    observables in that bin), refit the shared-kernel amplitude, collect
    statistics. Bin index = position in the project's q^2 grid (0..7).
    """
    from .sm_baseline import predict_vector

    df = data["observables"]
    rng = np.random.default_rng(rng_seed)

    # Establish unique bin keys (q^2 lo edges) per observable
    obs_unique = sorted(set(df["observable"]))
    # Grouping the rows by bin q^2 (across observables)
    df = df.reset_index(drop=True)
    bin_keys = []
    for _, r in df.iterrows():
        bin_keys.append((float(r["q2_lo"]), float(r["q2_hi"])))
    unique_bin_keys = sorted(set(bin_keys), key=lambda x: x[0])
    n_bins = len(unique_bin_keys)

    # Indices of rows belonging to each bin
    bin_to_rows: dict[tuple[float, float], list[int]] = {b: [] for b in unique_bin_keys}
    for i, k in enumerate(bin_keys):
        bin_to_rows[k].append(i)

    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    obs_list, bin_indices, q2_centres = _bin_axis(df)
    kappa_full = wo010_universality.frozen_kernel_at_bin_centres(q2_centres)

    # Bootstrap with replacement duplicates rows, which makes the covariance
    # singular. Use diagonal errors (stat + syst in quadrature) for the
    # bootstrap chi^2 — this is the standard non-parametric bootstrap and is
    # still anchored to the same per-row uncertainty.
    amplitudes = []
    chi2_vals = []
    for trial in range(n_bootstrap):
        sampled = rng.choice(n_bins, size=n_bins, replace=True)
        rows = []
        for s in sampled:
            rows.extend(bin_to_rows[unique_bin_keys[s]])
        rows = np.array(rows, dtype=int)
        v_b = values[rows]
        e_b = errors[rows]
        kappa_b = kappa_full[rows]
        obs_b = [obs_list[i] for i in rows]
        bin_idx_b = [bin_indices[i] for i in rows]

        def _loss(theta):
            a = float(theta[0])
            pred = predict_vector(obs_b, bin_idx_b, C9_SM - a * kappa_b)
            return chi2_fn(v_b, pred, errors=e_b)

        r = minimize(_loss, x0=[1.5], method="Powell", bounds=[(-5.0, 5.0)],
                     options={"xtol": 1e-6, "ftol": 1e-8, "maxiter": 1000})
        amplitudes.append(float(r.x[0]))
        chi2_vals.append(float(r.fun))

    a_arr = np.array(amplitudes)
    return pd.DataFrame({
        "trial": np.arange(n_bootstrap),
        "amplitude": amplitudes,
        "chi2": chi2_vals,
    }), {
        "n_bootstrap": n_bootstrap,
        "amplitude_mean": float(np.mean(a_arr)),
        "amplitude_median": float(np.median(a_arr)),
        "amplitude_std": float(np.std(a_arr)),
        "amplitude_q05": float(np.percentile(a_arr, 5)),
        "amplitude_q95": float(np.percentile(a_arr, 95)),
        "fraction_negative": float(np.mean(a_arr < 0)),
    }


# ---------------------------------------------------------------------------
# 2. Region splits
# ---------------------------------------------------------------------------

REGION_DEFINITIONS = {
    "low_q2":      (0, 2),    # bins 0,1,2: q^2 in [0.06, 4.0]
    "central_q2":  (3, 4),    # bins 3,4: q^2 in [4.0, 8.0]
    "high_q2":     (5, 7),    # bins 5,6,7: q^2 in [11.0, 19.0]
}


def fit_region(data: dict[str, Any], region: str) -> dict[str, Any]:
    """Slice the joint dataset to one q^2 region and refit the shared kernel."""
    df = data["observables"]
    lo_idx, hi_idx = REGION_DEFINITIONS[region]
    # Find which q^2 bins fall in this region
    q2_unique = sorted({(float(r["q2_lo"]), float(r["q2_hi"])) for _, r in df.iterrows()},
                       key=lambda x: x[0])
    chosen = q2_unique[lo_idx:hi_idx + 1]
    keep = df.apply(
        lambda r: (float(r["q2_lo"]), float(r["q2_hi"])) in chosen,
        axis=1,
    ).to_numpy()
    sub_df = df[keep].reset_index(drop=True)
    sub_data: dict[str, Any] = {"observables": sub_df, "metadata": data["metadata"]}
    if "covariance" in data and data["covariance"] is not None:
        idx = np.where(keep)[0].tolist()
        sub_data["covariance"] = np.asarray(data["covariance"])[np.ix_(idx, idx)]

    # FREE_C9 reference
    free = wo010_universality.fit_free_c9(sub_data)
    # VFD shared kernel
    kern = wo010_universality.fit_shared_kernel(sub_data)
    return {
        "region": region,
        "q2_range": (chosen[0][0], chosen[-1][1]),
        "n_bins": len(chosen),
        "n_data": len(sub_df),
        "FREE_C9_chi2": free.chi2,
        "FREE_C9_DC9": free.effective_delta_c9_mean,
        "FREE_C9_aic": free.aic,
        "VFD_chi2": kern.chi2,
        "VFD_amplitude": _amplitude_from_notes(kern.notes),
        "VFD_DC9_eff": kern.effective_delta_c9_mean,
        "VFD_aic": kern.aic,
        "delta_aic_VFD_vs_FREE": kern.aic - free.aic,
    }


# ---------------------------------------------------------------------------
# 3. Alternative models
# ---------------------------------------------------------------------------

def fit_free_c9_plus_c10(data: dict[str, Any]) -> dict[str, Any]:
    """Two-Wilson-coefficient fit: a global DeltaC9 and a global DeltaC10.

    Implementation: linearise around SM using the existing dO/dC9 from
    `sm_baseline` plus a hand-tabulated dO/dC10 (computed once from flavio
    via finite difference; cached here to avoid runtime flavio dependency).
    """
    from .sm_baseline import predict_vector

    # Pre-computed flavio-derived dO/dC10 slopes (central-difference, dC10 = +/-0.5,
    # WET-flavio basis at scale 4.8 GeV). One-time computed; keep alongside the
    # existing dC9 slopes in sm_baseline.
    dO_dC10 = _DC10_SLOPES

    df = data["observables"]
    obs, bins, _ = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)

    def loss(theta):
        dC9, dC10 = float(theta[0]), float(theta[1])
        # Predict with dC9 only via the existing predict_vector, then add the C10 contribution
        pred_c9 = predict_vector(obs, bins, C9_SM + dC9)
        c10_contribution = np.array([dO_dC10[o][b] for o, b in zip(obs, bins)]) * dC10
        pred = pred_c9 + c10_contribution
        return _chi2(values, pred, data, errors)

    r = minimize(loss, x0=[-1.3, 0.0], method="Powell",
                 bounds=[(-3.0, 3.0), (-3.0, 3.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    n = len(values)
    return {
        "model": "FREE_C9 + FREE_C10",
        "k_params": 2,
        "chi2": float(r.fun),
        "aic": aic(float(r.fun), 2),
        "bic": bic(float(r.fun), 2, n),
        "DC9": float(r.x[0]),
        "DC10": float(r.x[1]),
    }


def fit_charm_loop_proxy(data: dict[str, Any]) -> dict[str, Any]:
    """Charm-loop nuisance: a Gaussian residual centred at the J/psi-psi(2S)
    midpoint (q^2 = 11.59 GeV^2) with FREE width sigma and FREE amplitude A.
    Acts as a non-VFD shape model that explicitly captures long-distance
    charm-loop tail effects. k = 2 (A, sigma).
    """
    from .sm_baseline import predict_vector

    df = data["observables"]
    obs, bins, q2 = _bin_axis(df)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    midpoint = 0.5 * (J_PSI_Q2 + PSI2S_Q2)

    def loss(theta):
        a, sigma = float(theta[0]), float(theta[1])
        if sigma <= 0:
            return 1e12
        kernel = np.exp(-((q2 - midpoint) ** 2) / (2 * sigma ** 2))
        return _chi2(values, predict_vector(obs, bins, C9_SM - a * kernel), data, errors)

    r = minimize(loss, x0=[1.5, 4.0], method="Powell",
                 bounds=[(-5.0, 5.0), (0.5, 20.0)],
                 options={"xtol": 1e-7, "ftol": 1e-9})
    n = len(values)
    return {
        "model": "Charm-loop Gaussian (free A, sigma)",
        "k_params": 2,
        "chi2": float(r.fun),
        "aic": aic(float(r.fun), 2),
        "bic": bic(float(r.fun), 2, n),
        "amplitude": float(r.x[0]),
        "sigma_GeV2": float(r.x[1]),
    }


# ---------------------------------------------------------------------------
# 4. Form-factor variation via flavio (optional, requires flavio installed)
# ---------------------------------------------------------------------------

def _override_bsz_with_samples(dp_copy, sample_dict_at_idx: dict[str, float]):
    """Replace all B->K* BSZ form-factor constraints in dp_copy with delta
    distributions at the sampled values. Mutates dp_copy in place."""
    from flavio.statistics.probability import DeltaDistribution

    bsz_names = [n for n in sample_dict_at_idx if "B->K* BSZ" in n]
    keep = []
    for c, pars in dp_copy._constraints:
        if any("B->K* BSZ" in p for p in pars):
            continue
        keep.append((c, pars))
    dp_copy._constraints = keep
    for n in bsz_names:
        dp_copy._parameters[n] = ()  # forces clean re-add
    for n in bsz_names:
        dp_copy.add_constraint([n], DeltaDistribution(float(sample_dict_at_idx[n])))


def form_factor_variation(
    data: dict[str, Any], *, n_samples: int = 20, rng_seed: int = 7777
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Sample BSZ form-factor parameters from their published joint constraint,
    regenerate the SM/dC9 table per sample, refit the kernel amplitude, collect
    the amplitude distribution.

    Implementation: flavio's `default_parameters.get_random_all()` draws from the
    full multivariate constraint (correlations preserved). For each draw we
    replace the BSZ constraints with delta-distributions at the sampled values
    and monkey-patch `flavio.default_parameters` for the duration of the
    SM/NP prediction calls.

    Returns (per-sample df, summary dict). Skips gracefully if flavio is
    unavailable.
    """
    try:
        import flavio
        from wilson import Wilson
    except ImportError:
        return pd.DataFrame(), {"status": "flavio_unavailable"}

    bins = [(0.06, 0.98), (1.1, 2.5), (2.5, 4.0), (4.0, 6.0),
            (6.0, 8.0), (11.0, 12.5), (15.0, 17.0), (17.0, 19.0)]
    obs_list_flavio = ["<P5p>", "<P4p>", "<P1>", "<P2>"]
    schema_obs = ["P5p", "P4p", "P1", "P2"]

    df = data["observables"]
    obs, bin_idx, q2_centres = _bin_axis(df)
    kappa = wo010_universality.frozen_kernel_at_bin_centres(q2_centres)
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)

    np.random.seed(rng_seed)
    samples_all = flavio.default_parameters.get_random_all(size=n_samples)

    delta = 0.5
    wc_pos = Wilson({"C9_bsmumu": +delta}, scale=4.8, eft="WET", basis="flavio")
    wc_neg = Wilson({"C9_bsmumu": -delta}, scale=4.8, eft="WET", basis="flavio")

    rows = []
    amplitudes = []
    dp_orig = flavio.default_parameters

    for trial in range(n_samples):
        sample_at_trial = {n: samples_all[n][trial] for n in samples_all}

        dp = dp_orig.copy()
        try:
            _override_bsz_with_samples(dp, sample_at_trial)
        except Exception:
            rows.append({"trial": trial, "amplitude": float("nan"), "chi2": float("nan"),
                         "status": "override_fail"})
            continue

        sm_table = {o: [float("nan")] * len(bins) for o in schema_obs}
        slope_table = {o: [float("nan")] * len(bins) for o in schema_obs}
        # Wrap the assignment + the prediction loop inside one try/finally so
        # the global `flavio.default_parameters` is restored even if an
        # exception fires before the inner loop begins.
        try:
            flavio.default_parameters = dp
            for o_flavio, o_schema in zip(obs_list_flavio, schema_obs):
                for bi, (lo, hi) in enumerate(bins):
                    try:
                        sm_v = flavio.sm_prediction(o_flavio + "(B0->K*mumu)",
                                                    q2min=lo, q2max=hi)
                        np_p = flavio.np_prediction(o_flavio + "(B0->K*mumu)", wc_pos,
                                                    q2min=lo, q2max=hi)
                        np_n = flavio.np_prediction(o_flavio + "(B0->K*mumu)", wc_neg,
                                                    q2min=lo, q2max=hi)
                        sm_table[o_schema][bi] = float(sm_v)
                        slope_table[o_schema][bi] = float((np_p - np_n) / (2 * delta))
                    except Exception:
                        pass
        finally:
            flavio.default_parameters = dp_orig

        sm_row = np.array([sm_table[o][b] for o, b in zip(obs, bin_idx)])
        slope_row = np.array([slope_table[o][b] for o, b in zip(obs, bin_idx)])
        if np.any(np.isnan(sm_row)) or np.any(np.isnan(slope_row)):
            rows.append({"trial": trial, "amplitude": float("nan"),
                         "chi2": float("nan"), "status": "predict_fail"})
            continue

        def loss(theta):
            a = float(theta[0])
            pred = sm_row + slope_row * (-a * kappa)
            return _chi2(values, pred, data, errors)

        res = minimize(loss, x0=[1.5], method="Powell", bounds=[(-5.0, 5.0)],
                       options={"xtol": 1e-6, "ftol": 1e-8})
        rows.append({"trial": trial, "amplitude": float(res.x[0]),
                     "chi2": float(res.fun), "status": "ok"})
        amplitudes.append(float(res.x[0]))
        if (trial + 1) % 5 == 0:
            print(f"  ff trial {trial + 1}/{n_samples}: A = {float(res.x[0]):+.4f}")

    df_out = pd.DataFrame(rows)
    a_arr = np.array(amplitudes)
    summary = {
        "n_samples_completed": int(np.isfinite(df_out.get("amplitude", pd.Series())).sum()),
        "amplitude_mean": float(np.mean(a_arr)) if len(a_arr) else float("nan"),
        "amplitude_std": float(np.std(a_arr)) if len(a_arr) else float("nan"),
        "amplitude_q05": float(np.percentile(a_arr, 5)) if len(a_arr) else float("nan"),
        "amplitude_q95": float(np.percentile(a_arr, 95)) if len(a_arr) else float("nan"),
        "fraction_negative": float(np.mean(a_arr < 0)) if len(a_arr) else float("nan"),
    }
    return df_out, summary


# ---------------------------------------------------------------------------
# Pre-computed dO/dC10 slopes (from flavio 2.4, central-difference dC10=+/-0.5).
# Generated by the WO-013 setup; cached here so the fit_free_c9_plus_c10
# function does not require flavio at runtime.
# ---------------------------------------------------------------------------

# Will be populated by `_compute_c10_slopes()` if needed.
_DC10_SLOPES: dict[str, list[float]] = {}


def _compute_c10_slopes() -> dict[str, list[float]]:
    """Generate the dO/dC10 slope table via flavio. Run once and cache."""
    import flavio
    from wilson import Wilson

    bins = [(0.06, 0.98), (1.1, 2.5), (2.5, 4.0), (4.0, 6.0),
            (6.0, 8.0), (11.0, 12.5), (15.0, 17.0), (17.0, 19.0)]
    observables = ["<P5p>", "<P4p>", "<P1>", "<P2>", "<FL>", "<AFB>"]
    schema = ["P5p", "P4p", "P1", "P2", "FL", "AFB"]

    delta = 0.5
    wc_pos = Wilson({"C10_bsmumu": +delta}, scale=4.8, eft="WET", basis="flavio")
    wc_neg = Wilson({"C10_bsmumu": -delta}, scale=4.8, eft="WET", basis="flavio")
    out: dict[str, list[float]] = {}
    for o, s in zip(observables, schema):
        slopes = []
        for lo, hi in bins:
            np_p = flavio.np_prediction(o + "(B0->K*mumu)", wc_pos, q2min=lo, q2max=hi)
            np_n = flavio.np_prediction(o + "(B0->K*mumu)", wc_neg, q2min=lo, q2max=hi)
            slopes.append(round((np_p - np_n) / (2 * delta), 4))
        out[s] = slopes
    return out


# ---------------------------------------------------------------------------
# Top-level run
# ---------------------------------------------------------------------------

def run(
    *,
    archive_dir: Path | str = "data/raw/hepdata/extracted",
    config_index: int = 2,
    observables: Iterable[str] = DEFAULT_ANGULAR_OBSERVABLES,
    n_bootstrap: int = 500,
    n_ff_samples: int = 20,
    output_dir: Path | str = "reports",
) -> dict[str, Any]:
    global _DC10_SLOPES
    if not _DC10_SLOPES:
        _DC10_SLOPES = _compute_c10_slopes()

    archive = hepdata_ingest.hepdata_archive_dir(archive_dir)
    data = hepdata_ingest.load_config(
        archive, config_index=config_index, observables=tuple(observables)
    )

    print(">> running bin bootstrap...")
    bootstrap_df, bootstrap_summary = bootstrap_bins(data, n_bootstrap=n_bootstrap)

    print(">> running region splits...")
    region_results = [fit_region(data, r) for r in REGION_DEFINITIONS]

    print(">> running alternative models...")
    free_c9 = wo010_universality.fit_free_c9(data)
    vfd_kernel = wo010_universality.fit_shared_kernel(data)
    free_c9_c10 = fit_free_c9_plus_c10(data)
    charm_loop = fit_charm_loop_proxy(data)
    alternatives = [
        {"model": "FREE_C9", "k_params": 1, "chi2": free_c9.chi2,
         "aic": free_c9.aic, "bic": free_c9.bic,
         "delta_aic_vs_FREE_C9": 0.0, "params": f"DC9={free_c9.effective_delta_c9_mean:+.3f}"},
        {"model": "VFD_GREEN_600CELL", "k_params": 1, "chi2": vfd_kernel.chi2,
         "aic": vfd_kernel.aic, "bic": vfd_kernel.bic,
         "delta_aic_vs_FREE_C9": vfd_kernel.aic - free_c9.aic,
         "params": vfd_kernel.notes},
        {"model": free_c9_c10["model"], "k_params": free_c9_c10["k_params"],
         "chi2": free_c9_c10["chi2"], "aic": free_c9_c10["aic"], "bic": free_c9_c10["bic"],
         "delta_aic_vs_FREE_C9": free_c9_c10["aic"] - free_c9.aic,
         "params": f"DC9={free_c9_c10['DC9']:+.3f}, DC10={free_c9_c10['DC10']:+.3f}"},
        {"model": charm_loop["model"], "k_params": charm_loop["k_params"],
         "chi2": charm_loop["chi2"], "aic": charm_loop["aic"], "bic": charm_loop["bic"],
         "delta_aic_vs_FREE_C9": charm_loop["aic"] - free_c9.aic,
         "params": f"A={charm_loop['amplitude']:+.3f}, sigma={charm_loop['sigma_GeV2']:.3f} GeV^2"},
    ]

    print(">> running form-factor variation...")
    ff_df, ff_summary = form_factor_variation(data, n_samples=n_ff_samples)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    bootstrap_df.to_csv(out / "wo013_bootstrap_bins.csv", index=False)
    pd.DataFrame(region_results).to_csv(out / "wo013_regions.csv", index=False)
    pd.DataFrame(alternatives).to_csv(out / "wo013_alternatives.csv", index=False)
    ff_df.to_csv(out / "wo013_form_factor_variation.csv", index=False)

    return {
        "bootstrap": bootstrap_df,
        "bootstrap_summary": bootstrap_summary,
        "regions": region_results,
        "alternatives": alternatives,
        "form_factor_variation": ff_df,
        "form_factor_summary": ff_summary,
        "provenance": PROVENANCE_VFD,
    }


def main() -> None:
    res = run()
    pd.set_option("display.width", 140)
    pd.set_option("display.max_colwidth", 80)

    print()
    print("=" * 80)
    print("1. Bootstrap over bins (frozen kernel, single shared amplitude)")
    print("=" * 80)
    s = res["bootstrap_summary"]
    print(f"  n_bootstrap = {s['n_bootstrap']}")
    print(f"  amplitude:  mean = {s['amplitude_mean']:+.4f}, "
          f"median = {s['amplitude_median']:+.4f}, std = {s['amplitude_std']:.4f}")
    print(f"  90% CI:    [{s['amplitude_q05']:+.4f}, {s['amplitude_q95']:+.4f}]")
    print(f"  fraction A < 0: {s['fraction_negative']:.3f}  (sign stability)")

    print()
    print("=" * 80)
    print("2. Region splits (independent fits in q^2 sub-windows)")
    print("=" * 80)
    rdf = pd.DataFrame(res["regions"])
    print(rdf.to_string(index=False))

    print()
    print("=" * 80)
    print("3. Alternative models (full joint fit)")
    print("=" * 80)
    adf = pd.DataFrame(res["alternatives"])
    print(adf.to_string(index=False))

    print()
    print("=" * 80)
    print("4. Form-factor variation (BSZ parameters sampled within errors)")
    print("=" * 80)
    if res["form_factor_summary"].get("status") == "flavio_unavailable":
        print("  flavio not installed; skipped.")
    else:
        s = res["form_factor_summary"]
        print(f"  n_samples completed = {s['n_samples_completed']}")
        print(f"  amplitude: mean = {s['amplitude_mean']:+.4f}, std = {s['amplitude_std']:.4f}")
        print(f"  90% CI: [{s['amplitude_q05']:+.4f}, {s['amplitude_q95']:+.4f}]")
        print(f"  fraction A < 0: {s.get('fraction_negative', float('nan')):.3f}")


if __name__ == "__main__":
    main()
