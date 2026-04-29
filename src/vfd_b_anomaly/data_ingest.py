"""Schema-first data ingestion for LHCb B0 -> K*0 mu+ mu- observables.

The contract: data must be EXPLICIT about what is missing. There is no
silent fallback to zeros, no quiet swap from full covariance to diagonal,
no implicit synthesis. If a field is unavailable, the loader marks it as
such and a downstream consumer is responsible for deciding what to do.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .constants import (
    PROVENANCE_PLACEHOLDER,
    REQUIRED_OBSERVABLE_COLUMNS,
    SUPPORTED_OBSERVABLES,
)
from .paper_tables import DEFAULT_Q2_BIN_EDGES_GEV2

PROCESSED_DIR_DEFAULT = Path("data/processed")


def validate_observable_schema(data: dict[str, Any]) -> None:
    """Raise ValueError if `data` does not satisfy the observable contract.

    Required keys:
      - 'observables': pandas.DataFrame with REQUIRED_OBSERVABLE_COLUMNS
      - 'metadata': dict with at least 'source', 'covariance_status', 'date_accessed'
    """
    if not isinstance(data, dict):
        raise ValueError("data must be a dict")

    missing_keys = {"observables", "metadata"} - data.keys()
    if missing_keys:
        raise ValueError(f"missing required top-level keys: {sorted(missing_keys)}")

    df = data["observables"]
    if not isinstance(df, pd.DataFrame):
        raise ValueError("data['observables'] must be a pandas DataFrame")

    missing_cols = set(REQUIRED_OBSERVABLE_COLUMNS) - set(df.columns)
    if missing_cols:
        raise ValueError(f"observables missing columns: {sorted(missing_cols)}")

    if df.empty:
        raise ValueError("observables DataFrame is empty")

    # Bin edges must be ordered and form a proper bin (lo < hi).
    if not (df["q2_lo"] < df["q2_hi"]).all():
        raise ValueError("each row must satisfy q2_lo < q2_hi")

    # Errors must be strictly positive (zero error is undefined for chi^2).
    if (df["stat_err"] <= 0).any() or (df["syst_err"] < 0).any():
        raise ValueError("stat_err must be > 0 and syst_err must be >= 0 for every row")

    # Observable names must be supported.
    bad_obs = set(df["observable"]) - set(SUPPORTED_OBSERVABLES)
    if bad_obs:
        raise ValueError(f"unsupported observables present: {sorted(bad_obs)}")

    # Each (q2 bin, observable) row must have a provenance tag.
    if df["provenance"].isna().any() or (df["provenance"].astype(str).str.len() == 0).any():
        raise ValueError("every row must carry a non-empty provenance tag")

    # q^2 ordering within each observable must be monotonic.
    for obs, group in df.groupby("observable"):
        q2_lo = group["q2_lo"].to_numpy()
        if not np.all(np.diff(q2_lo) >= 0):
            raise ValueError(f"q2_lo not monotonic within observable {obs!r}")

    md = data["metadata"]
    if not isinstance(md, dict):
        raise ValueError("metadata must be a dict")
    for required_key in ("source", "covariance_status", "date_accessed"):
        if required_key not in md:
            raise ValueError(f"metadata missing key {required_key!r}")
    if md["covariance_status"] not in {"present", "absent_declared"}:
        raise ValueError(
            "metadata['covariance_status'] must be 'present' or 'absent_declared'; "
            "missing covariance must be declared explicitly, never implied"
        )


def load_lhcb_tables(source: str | Path) -> dict[str, Any]:
    """Load processed observables CSV + metadata JSON from a directory or file root.

    `source` can be:
      - a directory containing the standard processed filenames, or
      - a path stem (without extension) where {stem}.csv and {stem}.metadata.json exist.

    Returns the validated `data` dict (raises if invalid).
    """
    src = Path(source)
    if src.is_dir():
        csv_path = src / "lhcb_2512_18053_observables.csv"
        meta_path = src / "lhcb_2512_18053_metadata.json"
    else:
        csv_path = src.with_suffix(".csv")
        meta_path = src.with_suffix(".metadata.json")

    if not csv_path.exists():
        raise FileNotFoundError(f"observable CSV not found: {csv_path}")
    if not meta_path.exists():
        raise FileNotFoundError(f"metadata JSON not found: {meta_path}")

    df = pd.read_csv(csv_path)
    with meta_path.open("r", encoding="utf-8") as fh:
        metadata = json.load(fh)

    data = {"observables": df, "metadata": metadata}
    validate_observable_schema(data)
    return data


def load_covariance_matrix(source: str | Path | None) -> np.ndarray | None:
    """Load covariance matrix if available, else return None.

    None is the EXPLICIT signal of absence. Callers must handle it; they
    must not silently substitute a diagonal matrix.
    """
    if source is None:
        return None
    path = Path(source)
    if not path.exists():
        return None
    cov = np.load(path)
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1]:
        raise ValueError(f"covariance must be square 2-D, got shape {cov.shape}")
    if not np.allclose(cov, cov.T, atol=1e-10):
        raise ValueError("covariance matrix is not symmetric")
    return cov


def write_processed_dataset(
    data: dict[str, Any],
    out_path: str | Path,
    covariance: np.ndarray | None = None,
) -> dict[str, Path]:
    """Write observables CSV + metadata JSON (+ optional covariance .npy) under out_path.

    Returns a dict of the written paths. Validates schema before writing.
    """
    validate_observable_schema(data)
    out_dir = Path(out_path)
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "lhcb_2512_18053_observables.csv"
    meta_path = out_dir / "lhcb_2512_18053_metadata.json"

    data["observables"].to_csv(csv_path, index=False)
    with meta_path.open("w", encoding="utf-8") as fh:
        json.dump(data["metadata"], fh, indent=2, sort_keys=True)

    written = {"csv": csv_path, "metadata": meta_path}

    if covariance is not None:
        if data["metadata"].get("covariance_status") != "present":
            raise ValueError(
                "covariance supplied but metadata['covariance_status'] != 'present'; "
                "refusing to write inconsistent dataset"
            )
        n = len(data["observables"])
        if covariance.shape != (n, n):
            raise ValueError(
                f"covariance shape {covariance.shape} mismatched to observable vector length {n}"
            )
        cov_path = out_dir / "lhcb_2512_18053_covariance.npy"
        np.save(cov_path, covariance)
        written["covariance"] = cov_path

    return written


def make_placeholder_dataset(
    observables: tuple[str, ...] = SUPPORTED_OBSERVABLES,
    *,
    nominal_stat_err: float = 0.05,
    nominal_syst_err: float = 0.02,
) -> dict[str, Any]:
    """Build a placeholder dataset that conforms to the schema.

    Every row is tagged provenance=PROVENANCE_PLACEHOLDER. Values are set
    to NaN so any downstream consumer that treats them as real numbers
    will fail loudly. This is intentional: the test surface exists, but
    no real numerics leak in.
    """
    rows = []
    for obs in observables:
        if obs not in SUPPORTED_OBSERVABLES:
            raise ValueError(f"unsupported observable {obs!r}")
        for q2_lo, q2_hi in DEFAULT_Q2_BIN_EDGES_GEV2:
            rows.append(
                {
                    "q2_lo": float(q2_lo),
                    "q2_hi": float(q2_hi),
                    "observable": obs,
                    "value": float("nan"),
                    "stat_err": float(nominal_stat_err),
                    "syst_err": float(nominal_syst_err),
                    "provenance": PROVENANCE_PLACEHOLDER,
                }
            )
    df = pd.DataFrame(rows, columns=list(REQUIRED_OBSERVABLE_COLUMNS))
    metadata = {
        "source": "placeholder",
        "covariance_status": "absent_declared",
        "date_accessed": None,
        "notes": (
            "Placeholder schema only. Values are NaN by design; replace with digitised "
            "or supplementary-table data before any reproduction claim."
        ),
    }
    return {"observables": df, "metadata": metadata}
