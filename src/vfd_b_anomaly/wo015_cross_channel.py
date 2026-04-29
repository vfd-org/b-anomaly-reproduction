"""WO-015 — Cross-channel kernel validation.

Goal:
    Test the FROZEN VFD_GREEN_600CELL kernel on b -> s mu+mu- channels other
    than B -> K* mu+mu-. The kernel function is *not* changed; only one
    amplitude A per dataset is fitted.

Tier 1 datasets:
    1. B_s^0 -> phi mu+mu- angular (LHCb 2015, arXiv:1506.08777,
       HEPData ins1380188). Observables: F_L, S3, S4, S7 on 6 exclusive
       q^2 bins (excluding the 2 wide combined bins).

Tier 1 not run:
    2. B+ -> K+ mu+mu-: HEPData submissions for the relevant LHCb (1209.4284,
       1403.8045) and CMS (1662193) papers are NOT available — querying
       /download/submission/ins{ID}/yaml returns the HEPData 404 HTML page.
       Recorded as a known gap. Belle / Belle II angular analyses also lack
       public correlation matrices on HEPData and were skipped in WO-014.

Tier 2 (deferred):
    3. B^0 -> K+ pi- mu+mu- in the K*_{0,2}(1430) region (ins1486676).
       The final-state hadronic system is more complicated (S-wave + D-wave
       moments); not a clean first cross-channel test.

Rules (locked from WO-013/014):
    - kappa(q^2) frozen (600-cell discrete Green's response).
    - one amplitude A per dataset, no kernel reshape.
    - SM and dC9 slopes regenerated per dataset on its own bin grid via flavio
      (cached).
    - Per-dataset diagnostic: bootstrap over bins (sign stability, magnitude
      CI) and q^2 region split (low / central / high) where applicable.

Acceptance gates (from the WO-015 spec):
    1. A has the same sign as WO-014 (positive).
    2. Effective DC9 remains negative.
    3. VFD_GREEN_600CELL is AIC/BIC competitive with FREE_C9.
    4. Bootstrap A is sign-stable.
    5. Region splits are sign-uniform.
    6. Convention mismatches reported transparently.
"""

from __future__ import annotations

import warnings
warnings.simplefilter("ignore")

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import minimize

from . import wo009_full_lift, wo010_universality
from .constants import C9_SM, PHI, PROVENANCE_VFD
from .flavio_predictor import default as flavio_default
from .likelihood import aic, bic, chi2 as chi2_fn


CROSS_CHANNEL_DIR = Path("data/raw/hepdata_cross_channel")


# ---------------------------------------------------------------------------
# Dataset loaders
# ---------------------------------------------------------------------------

def _err_to_sigma(err_entry: dict) -> tuple[float, float]:
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
            stat = max(stat, v)
    return stat, syst


def load_bs_phimumu_2015(observables=("FL", "S3", "S4", "S7")) -> dict[str, Any]:
    """LHCb 2015 (3 fb^-1) B_s -> phi mu+mu- angular: F_L, S3, S4, S7.

    Drops the 2 wide combined q^2 bins ([1, 6] and [15, 19]) — they overlap
    the exclusive ones. Keeps the 6 exclusive bins.
    """
    base = CROSS_CHANNEL_DIR / "ins1380188/HEPData-ins1380188-v1-yaml"
    d = yaml.safe_load((base / "Table2.yaml").read_text())
    bin_pairs = [(float(b["low"]), float(b["high"]))
                 for b in d["independent_variables"][0]["values"]]
    # Exclusive bins (drop wide combined ones at indices 6, 7 by their q^2 ranges)
    keep_bin = [i for i, (lo, hi) in enumerate(bin_pairs)
                if not (abs(lo - 1.0) < 1e-3 and abs(hi - 6.0) < 1e-3)
                and not (abs(lo - 15.0) < 1e-3 and abs(hi - 19.0) < 1e-3)]
    name_map = {r"$F_{\rm L}$": "FL", r"$S_3$": "S3",
                r"$S_4$": "S4", r"$S_7$": "S7"}
    rows = []
    for v in d["dependent_variables"]:
        obs_name = name_map.get(v["header"]["name"], v["header"]["name"])
        if obs_name not in observables:
            continue
        for bi in keep_bin:
            val = v["values"][bi]
            stat, syst = _err_to_sigma(val)
            rows.append({
                "observable": obs_name,
                "q2_lo": bin_pairs[bi][0],
                "q2_hi": bin_pairs[bi][1],
                "value": float(val["value"]),
                "stat_err": stat,
                "syst_err": syst,
            })
    df = pd.DataFrame(rows).sort_values(["q2_lo", "observable"]).reset_index(drop=True)
    return {
        "dataset": "Bs2phi-LHCb-2015",
        "decay": "Bs->phimumu",
        "n_obs": len(df),
        "covariance": None,  # not published in HEPData submission
        "observables": df,
        "metadata": {"arxiv": "1506.08777", "hepdata": "ins1380188",
                     "luminosity_fb": 3.0, "table": "Table2",
                     "covariance": "diagonal (HEPData has none for angular)"},
    }


# ---------------------------------------------------------------------------
# Generic universality fit + bootstrap + region split
# ---------------------------------------------------------------------------

@dataclass
class FitOutcome:
    dataset: str; decay: str; model: str
    n_data: int; k_params: int; chi2: float; aic: float; bic: float
    delta_aic_vs_FREE_C9: float; A: float; DC9_eff_mean: float; notes: str = ""


def _chi2_with(values, pred, errors, cov):
    if cov is not None:
        try:
            return chi2_fn(values, pred, covariance=0.5 * (cov + cov.T))
        except np.linalg.LinAlgError:
            eps = 1e-9 * np.trace(cov) / cov.shape[0]
            return chi2_fn(values, pred, covariance=cov + eps * np.eye(cov.shape[0]))
    return chi2_fn(values, pred, errors=errors)


def _prepare(data):
    df = data["observables"]
    decay = data["decay"]
    obs = df["observable"].tolist()
    q2_lo = df["q2_lo"].tolist()
    q2_hi = df["q2_hi"].tolist()
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 + df["syst_err"].to_numpy() ** 2)
    cov = data.get("covariance")
    fp = flavio_default()
    fp.precompute(set(obs), set(zip(q2_lo, q2_hi)), decay=decay)
    sm_vec = np.array([fp.sm_value(o, lo, hi, decay=decay)
                       for o, lo, hi in zip(obs, q2_lo, q2_hi)])
    slope_vec = np.array([fp.dO_dC9(o, lo, hi, decay=decay)
                          for o, lo, hi in zip(obs, q2_lo, q2_hi)])
    q2_centres = 0.5 * (np.array(q2_lo) + np.array(q2_hi))
    kappa = wo010_universality.frozen_kernel_at_bin_centres(q2_centres)
    fp.flush()
    return {
        "df": df, "decay": decay, "obs": obs, "q2_lo": q2_lo, "q2_hi": q2_hi,
        "values": values, "errors": errors, "cov": cov,
        "sm_vec": sm_vec, "slope_vec": slope_vec, "kappa": kappa,
    }


def _at_bound(value: float, lo: float, hi: float, tol: float = 1e-3) -> str | None:
    if abs(value - lo) < tol:
        return f"pinned at lower bound {lo}"
    if abs(value - hi) < tol:
        return f"pinned at upper bound {hi}"
    return None


def _warn_bound(label: str, msg: str) -> None:
    """Print directly (rather than warnings.warn) because the module-level
    `warnings.simplefilter('ignore')` would swallow it. We want bound-pinning
    to be visible in the run log."""
    print(f"  WARN [{label}]: {msg}")


def fit_dataset(data: dict[str, Any]) -> dict[str, FitOutcome]:
    p = _prepare(data)
    n = len(p["values"])
    cov = p["cov"]

    free_bounds = (-10.0, 10.0)
    vfd_bounds = (-20.0, 20.0)

    def loss_free(theta):
        pred = p["sm_vec"] + p["slope_vec"] * float(theta[0])
        return _chi2_with(p["values"], pred, p["errors"], cov)

    r_free = minimize(loss_free, x0=[-1.0], method="Powell", bounds=[free_bounds],
                      options={"xtol": 1e-7, "ftol": 1e-9})
    chi2_free = float(r_free.fun); DC9_free = float(r_free.x[0])
    a_free = aic(chi2_free, 1); b_free = bic(chi2_free, 1, n)
    free_warn = _at_bound(DC9_free, *free_bounds)
    if free_warn:
        _warn_bound(data["dataset"], f"FREE_C9 ΔC9 {free_warn}")

    def loss_vfd(theta):
        a = float(theta[0])
        pred = p["sm_vec"] + p["slope_vec"] * (-a * p["kappa"])
        return _chi2_with(p["values"], pred, p["errors"], cov)

    r_vfd = minimize(loss_vfd, x0=[1.0], method="Powell", bounds=[vfd_bounds],
                     options={"xtol": 1e-7, "ftol": 1e-9})
    chi2_vfd = float(r_vfd.fun); A = float(r_vfd.x[0])
    a_vfd = aic(chi2_vfd, 1); b_vfd = bic(chi2_vfd, 1, n)
    vfd_warn = _at_bound(A, *vfd_bounds)
    if vfd_warn:
        _warn_bound(data["dataset"], f"VFD A {vfd_warn}")
    DC9_eff = float(np.mean(-A * p["kappa"]))

    return {
        "FREE_C9": FitOutcome(
            data["dataset"], p["decay"], "FREE_C9", n, 1,
            chi2_free, a_free, b_free, 0.0, DC9_free, DC9_free,
            f"DC9={DC9_free:+.3f}",
        ),
        "VFD": FitOutcome(
            data["dataset"], p["decay"], "VFD_GREEN_600CELL", n, 1,
            chi2_vfd, a_vfd, b_vfd, a_vfd - a_free, A, DC9_eff,
            f"A={A:+.4f}, DC9_eff={DC9_eff:+.3f}",
        ),
    }


def bootstrap_amplitude(data: dict[str, Any], *, n_bootstrap: int = 500,
                        rng_seed: int = 12345) -> dict[str, float]:
    """Bin bootstrap (each bin = unit of all observables in that bin)."""
    p = _prepare(data)
    df = p["df"]
    bin_keys = list(zip(df["q2_lo"], df["q2_hi"]))
    unique_bins = sorted(set(bin_keys), key=lambda x: x[0])
    bin_to_rows = {b: [] for b in unique_bins}
    for i, k in enumerate(bin_keys):
        bin_to_rows[k].append(i)
    n_bins = len(unique_bins)

    rng = np.random.default_rng(rng_seed)
    amps = []
    for _ in range(n_bootstrap):
        sampled = rng.choice(n_bins, size=n_bins, replace=True)
        rows = []
        for s in sampled:
            rows.extend(bin_to_rows[unique_bins[s]])
        rows = np.array(rows, dtype=int)
        v_b = p["values"][rows]; e_b = p["errors"][rows]
        sm_b = p["sm_vec"][rows]; sl_b = p["slope_vec"][rows]
        ka_b = p["kappa"][rows]

        def loss(theta):
            a = float(theta[0])
            pred = sm_b + sl_b * (-a * ka_b)
            return chi2_fn(v_b, pred, errors=e_b)

        r = minimize(loss, x0=[1.5], method="Powell", bounds=[(-20.0, 20.0)],
                     options={"xtol": 1e-6, "ftol": 1e-8})
        amps.append(float(r.x[0]))
    arr = np.array(amps)
    return {
        "n_bootstrap": n_bootstrap,
        "amplitude_mean": float(np.mean(arr)),
        "amplitude_median": float(np.median(arr)),
        "amplitude_std": float(np.std(arr)),
        "amplitude_q05": float(np.percentile(arr, 5)),
        "amplitude_q95": float(np.percentile(arr, 95)),
        "fraction_negative": float(np.mean(arr < 0)),
    }


def region_splits(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Refit on low / central / high q^2 sub-windows. Bins are split by
    q^2 ranges that mirror WO-013 logic but adapted to whatever bins the
    dataset publishes."""
    df = data["observables"]
    bins_present = sorted(set(zip(df["q2_lo"], df["q2_hi"])), key=lambda x: x[0])

    def _filter(lo_bound, hi_bound):
        keep = [b for b in bins_present if lo_bound <= b[0] and b[1] <= hi_bound]
        return keep

    regions = {
        "low_q2":     _filter(0.0, 4.0),
        "central_q2": _filter(4.0, 8.5),
        "high_q2":    _filter(10.0, 20.0),
    }
    out = []
    for name, bins in regions.items():
        if not bins:
            out.append({"region": name, "n_bins": 0, "skipped": True})
            continue
        keep_mask = df.apply(lambda r: (float(r["q2_lo"]), float(r["q2_hi"])) in bins,
                             axis=1).to_numpy()
        sub_df = df[keep_mask].reset_index(drop=True)
        sub_cov = None
        if data.get("covariance") is not None:
            cov_idx = np.where(keep_mask)[0]
            sub_cov = np.asarray(data["covariance"])[np.ix_(cov_idx, cov_idx)]
        sub_data = {**data, "observables": sub_df,
                    "covariance": sub_cov,
                    "n_obs": len(sub_df)}
        fits = fit_dataset(sub_data)
        out.append({
            "region": name,
            "q2_range": (bins[0][0], bins[-1][1]),
            "n_bins": len(bins),
            "n_data": len(sub_df),
            "FREE_C9_chi2": fits["FREE_C9"].chi2,
            "FREE_C9_DC9": fits["FREE_C9"].A,
            "VFD_chi2": fits["VFD"].chi2,
            "VFD_A": fits["VFD"].A,
            "VFD_DC9_eff": fits["VFD"].DC9_eff_mean,
            "delta_aic": fits["VFD"].delta_aic_vs_FREE_C9,
        })
    return out


# ---------------------------------------------------------------------------
# Top-level run
# ---------------------------------------------------------------------------

LOADERS = {
    "Bs2phi-LHCb-2015": load_bs_phimumu_2015,
}


def run(*, n_bootstrap: int = 500,
        output_dir: Path | str = "reports") -> dict[str, Any]:
    main_rows = []
    bootstrap_summaries: dict[str, dict[str, float]] = {}
    region_results: dict[str, list[dict[str, Any]]] = {}
    for name, loader in LOADERS.items():
        print(f">> loading {name} ...")
        data = loader()
        n = data["n_obs"]
        bins = sorted(set(zip(data["observables"]["q2_lo"],
                              data["observables"]["q2_hi"])))
        print(f"   n_data={n}, bins={len(bins)}, decay={data['decay']}")
        print(f">> fitting {name} ...")
        fits = fit_dataset(data)
        for fr in fits.values():
            main_rows.append({
                "dataset": fr.dataset, "decay": fr.decay, "model": fr.model,
                "n_data": fr.n_data, "k_params": fr.k_params,
                "chi2": fr.chi2, "aic": fr.aic, "bic": fr.bic,
                "delta_aic_vs_FREE_C9": fr.delta_aic_vs_FREE_C9,
                "A_or_DC9": fr.A, "DC9_eff_mean": fr.DC9_eff_mean,
                "notes": fr.notes,
            })
        print(f">> bootstrap {name} ...")
        bootstrap_summaries[name] = bootstrap_amplitude(
            data, n_bootstrap=n_bootstrap)
        print(f">> region splits {name} ...")
        region_results[name] = region_splits(data)

    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    main_df = pd.DataFrame(main_rows)
    main_df.to_csv(out / "wo015_cross_channel.csv", index=False)
    boot_df = pd.DataFrame.from_dict(bootstrap_summaries, orient="index")
    boot_df.to_csv(out / "wo015_bootstrap.csv")
    region_rows = []
    for ds, regs in region_results.items():
        for r in regs:
            region_rows.append({"dataset": ds, **r})
    pd.DataFrame(region_rows).to_csv(out / "wo015_regions.csv", index=False)

    return {
        "main": main_df,
        "bootstrap": bootstrap_summaries,
        "regions": region_results,
        "provenance": PROVENANCE_VFD,
    }


def main() -> None:
    pd.set_option("display.width", 160)
    pd.set_option("display.max_colwidth", 70)
    res = run()

    print()
    print("=" * 100)
    print("WO-015 — Cross-channel frozen-kernel universality")
    print("=" * 100)
    print(res["main"].to_string(index=False))

    print()
    print("Bootstrap (per dataset):")
    for ds, s in res["bootstrap"].items():
        print(f"  {ds}: n={s['n_bootstrap']}, mean A = {s['amplitude_mean']:+.4f}, "
              f"std {s['amplitude_std']:.4f}, "
              f"90% CI [{s['amplitude_q05']:+.4f}, {s['amplitude_q95']:+.4f}], "
              f"frac A<0 = {s['fraction_negative']:.3f}")

    print()
    print("Region splits:")
    for ds, regs in res["regions"].items():
        print(f"  {ds}:")
        for r in regs:
            if r.get("skipped"):
                print(f"    {r['region']}: SKIPPED (no bins in window)")
                continue
            print(f"    {r['region']:11s}: q2={r['q2_range']}, "
                  f"n_data={r['n_data']}, FREE_C9 chi^2={r['FREE_C9_chi2']:.2f} "
                  f"DC9={r['FREE_C9_DC9']:+.3f}, "
                  f"VFD chi^2={r['VFD_chi2']:.2f} A={r['VFD_A']:+.3f}, "
                  f"dAIC={r['delta_aic']:+.3f}")

    # Universality verdict
    vfd = res["main"][res["main"]["model"] == "VFD_GREEN_600CELL"]
    print()
    print("Acceptance verdict:")
    print(f"  All A > 0?         {(vfd['A_or_DC9'] > 0).all()}")
    print(f"  All DC9_eff < 0?   {(vfd['DC9_eff_mean'] < 0).all()}")
    print(f"  All ΔAIC <= 0?     {(vfd['delta_aic_vs_FREE_C9'] <= 0).all()}")


if __name__ == "__main__":
    main()
