"""Ingestion of the LHCb HEPData supplementary archive (ins3094698, arXiv:2512.18053).

The archive layout (as extracted) provides one results JSON + one correlation
JSON per fit configuration:

    config_{i}_results.json
    config_{i}_correlation_{stat,syst,total}.json    (i in 1..6)

Schema (results JSON):
    {
      "independent_variables": {
        "name": "q2",
        "units": "GeV^2/c^4",
        "values": [[lo, hi], ... ]              # one entry per HEPData bin
      },
      "dependent_variables": {
        <obs_name>: {
          "units": <scalar multiplier, e.g. 0.1 for 10^-1>,
          "latex": <pretty label>,
          "values": [
            { "val": ..., "fit_val": ..., "fit_err": ...,
              "hi_total_err": ..., "lo_total_err": ...,
              "hi_syst_err": ..., "lo_syst_err": ...,
              "hi_stat_err": ..., "lo_stat_err": ...,
              "detector": ... }, ...
          ]
        }, ...
      }
    }

Schema (correlation JSON): a flat dict keyed by "{obs}_bin{N}" whose values
are dicts mapping the same {obs}_bin{N} identifiers to correlation coefficients.
Some entries are missing if HEPData treated a particular (obs, bin) pair as
exactly uncorrelated; we treat missing entries as 0.

This module follows the project's data_ingest contract: the returned dict
has 'observables' (DataFrame), 'metadata' (dict), and optionally 'covariance'
(np.ndarray aligned to the DataFrame row order).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .constants import (
    PROVENANCE_REPRODUCED,
    REQUIRED_OBSERVABLE_COLUMNS,
)
from .paper_tables import (
    DEFAULT_Q2_BIN_EDGES_GEV2,
    HEPDATA_BIN_INDEX_FOR_DEFAULT,
)

# HEPData observable name -> project schema observable name. Only a subset of
# the 27 HEPData observables map directly to our SM-baseline scaffold today.
HEPDATA_OBSERVABLE_ALIASES: dict[str, str] = {
    "P5p": "P5p",
    "P4p": "P4p",
    "P1": "P1",
    "P2": "P2",
    "dBF/dq2": "BR",
}

# Optional secondary scaling applied AFTER HEPData's intrinsic units multiplier,
# to bring the value into the project's preferred unit for that schema name.
# Example: HEPData stores dBF/dq2 with a 10^-8/GeV^2 multiplier (so values come
# out in absolute /GeV^2), but the project's SM baseline table for BR is in
# 10^-7/GeV^2 coefficient form. Multiplying by 1e7 here aligns them.
SCHEMA_VALUE_RESCALE: dict[str, float] = {
    "BR": 1.0e7,
}

DEFAULT_CONFIG_INDEX = 2  # Pi-basis, partially massive, standard binning


def _parse_units(raw) -> float:
    """HEPData encodes the units multiplier as either a float (e.g. 0.1) or a
    string like '10^-8'. Normalise to a positive float multiplier."""
    if isinstance(raw, (int, float)):
        return float(raw)
    if isinstance(raw, str):
        s = raw.strip()
        if "^" in s:
            base_str, exp_str = s.split("^", 1)
            try:
                return float(base_str) ** float(exp_str)
            except ValueError as exc:
                raise ValueError(f"cannot parse HEPData units string {raw!r}") from exc
        try:
            return float(s)
        except ValueError as exc:
            raise ValueError(f"cannot parse HEPData units string {raw!r}") from exc
    raise TypeError(f"unexpected HEPData units type {type(raw)!r}: {raw!r}")


def hepdata_archive_dir(extracted_root: Path | str) -> Path:
    """Resolve the directory containing the unzipped HEPData archive."""
    p = Path(extracted_root)
    if (p / "submission.yaml").exists():
        return p
    sub = p / "extracted"
    if (sub / "submission.yaml").exists():
        return sub
    raise FileNotFoundError(
        f"No submission.yaml found at {p} or {sub}; "
        "extract HEPData-ins3094698-v1.zip first."
    )


def _config_paths(extracted_dir: Path, config_index: int) -> dict[str, Path]:
    if not (1 <= config_index <= 6):
        raise ValueError(f"config_index must be in 1..6, got {config_index}")
    base = extracted_dir / f"config_{config_index}"
    return {
        "results": Path(f"{base}_results.json"),
        "correlation_total": Path(f"{base}_correlation_total.json"),
        "correlation_stat": Path(f"{base}_correlation_stat.json"),
        "correlation_syst": Path(f"{base}_correlation_syst.json"),
    }


def _load_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"missing HEPData file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def parse_config_results(
    results_path: Path | str,
    *,
    observables: Iterable[str] = ("P5p", "dBF/dq2"),
) -> tuple[pd.DataFrame, dict]:
    """Parse one config_{i}_results.json into the project schema DataFrame.

    Returns (observables_df, metadata).
    Only the LHCb exclusive bins (positions defined in HEPDATA_BIN_INDEX_FOR_DEFAULT)
    are included; the integrated [1.1, 6.0] / [15, 19] re-bins are excluded.
    """
    data = _load_json(Path(results_path))

    iv = data.get("independent_variables")
    dv = data.get("dependent_variables")
    if not isinstance(iv, dict) or iv.get("name") != "q2":
        raise ValueError(f"unexpected independent_variables structure in {results_path}")
    if not isinstance(dv, dict):
        raise ValueError(f"unexpected dependent_variables structure in {results_path}")

    hepdata_q2_bins: list[tuple[float, float]] = [
        (float(b[0]), float(b[1])) for b in iv["values"]
    ]

    requested = list(observables)
    missing = [o for o in requested if o not in dv]
    if missing:
        raise KeyError(
            f"observables {missing} not present in {results_path}; "
            f"available: {sorted(dv.keys())}"
        )

    rows: list[dict] = []
    # Ordering: project's DEFAULT_Q2_BIN_EDGES_GEV2 ordering, mapped to HEPData index.
    for proj_bin_idx, (q2_lo, q2_hi) in enumerate(DEFAULT_Q2_BIN_EDGES_GEV2):
        hep_idx = HEPDATA_BIN_INDEX_FOR_DEFAULT[proj_bin_idx]
        # Sanity: the HEPData bin edges at hep_idx must match the project edges.
        hb_lo, hb_hi = hepdata_q2_bins[hep_idx]
        if not (np.isclose(hb_lo, q2_lo) and np.isclose(hb_hi, q2_hi)):
            raise ValueError(
                f"bin-edge mismatch at proj_bin {proj_bin_idx}: "
                f"project=[{q2_lo}, {q2_hi}] vs hepdata=[{hb_lo}, {hb_hi}] (hep_idx={hep_idx})"
            )
        for hep_obs in requested:
            obs_block = dv[hep_obs]
            unit = _parse_units(obs_block.get("units", 1.0))
            entry = obs_block["values"][hep_idx]
            val = float(entry["val"]) * unit
            # Use the larger of |hi| and |lo| as a symmetric error proxy.
            hi_stat = abs(float(entry.get("hi_stat_err", 0.0))) * unit
            lo_stat = abs(float(entry.get("lo_stat_err", 0.0))) * unit
            hi_syst = abs(float(entry.get("hi_syst_err", 0.0))) * unit
            lo_syst = abs(float(entry.get("lo_syst_err", 0.0))) * unit
            stat_err = max(hi_stat, lo_stat)
            syst_err = max(hi_syst, lo_syst)
            if stat_err <= 0:
                # HEPData rows must have a non-zero stat error for chi^2.
                # If the supplementary table reports zero (exact theory point), skip.
                continue
            schema_obs = HEPDATA_OBSERVABLE_ALIASES.get(hep_obs, hep_obs)
            rescale = SCHEMA_VALUE_RESCALE.get(schema_obs, 1.0)
            rows.append(
                {
                    "q2_lo": float(q2_lo),
                    "q2_hi": float(q2_hi),
                    "observable": schema_obs,
                    "value": val * rescale,
                    "stat_err": stat_err * rescale,
                    "syst_err": syst_err * rescale,
                    "provenance": PROVENANCE_REPRODUCED,
                    # Internal-only column; preserved for covariance lookup.
                    "_hep_obs": hep_obs,
                    "_hep_bin_idx": int(hep_idx),
                }
            )

    if not rows:
        raise ValueError(f"no usable rows extracted from {results_path}")

    df = pd.DataFrame(
        rows,
        columns=list(REQUIRED_OBSERVABLE_COLUMNS) + ["_hep_obs", "_hep_bin_idx"],
    )

    metadata = {
        "source": f"hepdata:{Path(results_path).name}",
        "covariance_status": "absent_declared",  # populated by load_config when covariance is attached
        "date_accessed": "2026-04-28",
        "hepdata_record": "ins3094698",
        "hepdata_doi": "10.17182/hepdata.167733.v1",
        "hepdata_observables_requested": requested,
        "hepdata_q2_bins": hepdata_q2_bins,
        "project_q2_bins": DEFAULT_Q2_BIN_EDGES_GEV2,
        "notes": (
            "Real LHCb data from HEPData supplementary archive. Integrated re-bins "
            "[1.1, 6.0] and [15.0, 19.0] are intentionally excluded because they "
            "overlap their constituent bins and are not independent measurements."
        ),
    }
    return df, metadata


def parse_config_correlation(
    correlation_path: Path | str,
    df: pd.DataFrame,
) -> np.ndarray:
    """Build a covariance matrix aligned to df row ordering, from a HEPData correlation file.

    Requires `df` to carry the internal `_hep_obs` and `_hep_bin_idx` columns
    populated by parse_config_results.
    """
    cor = _load_json(Path(correlation_path))
    if not isinstance(cor, dict):
        raise ValueError(f"unexpected correlation file structure in {correlation_path}")

    required_cols = {"_hep_obs", "_hep_bin_idx", "stat_err", "syst_err"}
    if not required_cols.issubset(df.columns):
        raise ValueError(
            f"DataFrame must carry columns {required_cols}; got {set(df.columns)}"
        )

    n = len(df)
    keys = [f"{row['_hep_obs']}_bin{int(row['_hep_bin_idx'])}" for _, row in df.iterrows()]
    missing_rows = [k for k in keys if k not in cor]
    if missing_rows:
        raise KeyError(f"correlation matrix missing row keys: {missing_rows}")

    rho = np.eye(n, dtype=float)
    for i, ki in enumerate(keys):
        row_dict = cor[ki]
        for j, kj in enumerate(keys):
            if ki == kj:
                rho[i, j] = 1.0
                continue
            # Try both orderings: correlations are symmetric but the JSON
            # may only store one half.
            v = row_dict.get(kj)
            if v is None:
                v = cor.get(kj, {}).get(ki, 0.0)
            rho[i, j] = float(v)

    # Symmetrise to suppress floating drift / one-sided storage.
    rho = 0.5 * (rho + rho.T)

    # Diagonal sigma combines stat and syst in quadrature.
    sigma = np.sqrt(df["stat_err"].to_numpy(dtype=float) ** 2 + df["syst_err"].to_numpy(dtype=float) ** 2)
    cov = rho * np.outer(sigma, sigma)

    # Numerical safety: clip tiny negative eigenvalues by adding a small ridge
    # if needed (rare; happens when the covariance is borderline PD).
    try:
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        eps = 1e-10 * np.trace(cov) / n
        cov = cov + eps * np.eye(n)

    return cov


def load_config(
    archive_dir: Path | str,
    *,
    config_index: int = DEFAULT_CONFIG_INDEX,
    observables: Iterable[str] = ("P5p",),
    correlation_kind: str = "total",
) -> dict:
    """Load one fit configuration and return the project's standard data dict.

    correlation_kind: 'total' | 'stat' | 'syst' (or 'none' to skip).
    """
    extracted = hepdata_archive_dir(archive_dir)
    paths = _config_paths(extracted, config_index)

    df, metadata = parse_config_results(paths["results"], observables=observables)
    metadata["config_index"] = config_index
    metadata["correlation_kind"] = correlation_kind

    out: dict = {"observables": df, "metadata": metadata}

    if correlation_kind == "none":
        return out

    cor_key = f"correlation_{correlation_kind}"
    if cor_key not in paths:
        raise ValueError(f"unknown correlation_kind {correlation_kind!r}")
    covariance = parse_config_correlation(paths[cor_key], df)
    out["covariance"] = covariance
    metadata["covariance_status"] = "present"
    metadata["covariance_source"] = f"hepdata:{paths[cor_key].name}"
    return out


def public_dataframe(data: dict) -> pd.DataFrame:
    """Return df with ONLY the public schema columns (drop internal _hep_* helpers)."""
    df = data["observables"].copy()
    drop = [c for c in df.columns if c.startswith("_hep_")]
    return df.drop(columns=drop)
