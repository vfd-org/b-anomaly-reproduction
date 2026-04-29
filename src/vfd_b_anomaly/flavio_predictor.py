"""Generic flavio-on-demand SM and dO/dC9 predictor.

WO-014 needs SM predictions and Wilson-coefficient slopes for arbitrary q^2
bins across multiple decay processes (B^0->K*0, B^+->K*+, B_s->phi). The
project's hand-tabulated `sm_baseline` is fixed to the LHCb 2025 8-bin
grid — too rigid for cross-dataset validation. This module provides:

    sm_value(observable, q2_lo, q2_hi, decay)        -> float
    dO_dC9(observable, q2_lo, q2_hi, decay, delta=0.5) -> float
    predict(observable, q2_lo, q2_hi, decay, c9)     -> float

with a JSON-on-disk cache so a flavio call is made at most once per unique
(observable, q^2 bin, decay) combination across the project's lifetime.

Decay strings follow flavio's notation: "B0->K*mumu", "B+->K*+mumu",
"Bs->phimumu". Observable-prefix strings: "<P5p>", "<P4p>", "<P1>", "<P2>",
"<FL>", "<AFB>", etc.

The flavio QCDF interpolation issues a "do not trust above 6 GeV^2" warning
for B0->K*. We use the values anyway and note the caveat in the WO-014 report.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Iterable

warnings.simplefilter("ignore")

import numpy as np


_OBS_TO_FLAVIO = {
    "P5p": "<P5p>",
    "P4p": "<P4p>",
    "P1":  "<P1>",
    "P2":  "<P2>",
    "P3":  "<P3>",
    "P6p": "<P6p>",
    "P8p": "<P8p>",
    "FL":  "<FL>",
    "AFB": "<AFB>",
    "FH":  "<FH>",
    "S3":  "<S3>",
    "S4":  "<S4>",
    "S5":  "<S5>",
    "S7":  "<S7>",
    "S8":  "<S8>",
    "S9":  "<S9>",
    "BR":  "<dBR/dq2>",
}


def _decay_to_flavio(decay: str) -> str:
    norm = decay.lower().replace(" ", "").replace("μ", "mu")
    # flavio uses "B+->K*mumu" (no plus on K*) for the charged decay; isospin
    # partner of the neutral "B0->K*mumu". Keep the user-facing notation
    # explicit ("B+->K*+mumu") and translate on the way in.
    table = {
        "b0->k*mumu":   "B0->K*mumu",
        "b0->k*0mumu":  "B0->K*mumu",
        "b+->k*+mumu":  "B+->K*mumu",
        "b+->k*mumu":   "B+->K*mumu",
        "bs->phimumu":  "Bs->phimumu",
        "b0->kmumu":    "B0->Kmumu",
        "b+->k+mumu":   "B+->K+mumu",
        "b+->kmumu":    "B+->K+mumu",
    }
    if norm not in table:
        raise KeyError(f"Unknown decay {decay!r}; supported: {sorted(set(table.values()))}")
    return table[norm]


_DEFAULT_CACHE_PATH = Path(__file__).resolve().parents[2] / "data" / "processed" / "flavio_cache.json"


class FlavioPredictor:
    """Cached flavio SM and dO/dC9 predictor.

    Cache key is (decay, observable, q2_lo, q2_hi, delta_C9). Persisted as JSON
    so subsequent script runs do not re-call flavio.
    """

    def __init__(self, cache_path: Path | str = _DEFAULT_CACHE_PATH):
        self.cache_path = Path(cache_path)
        self._cache: dict[str, float] = {}
        if self.cache_path.exists():
            try:
                self._cache = json.loads(self.cache_path.read_text())
            except json.JSONDecodeError:
                self._cache = {}
        self._dirty = False
        self._flavio = None
        self._wilson = None

    def _ensure_flavio(self) -> None:
        if self._flavio is not None:
            return
        import flavio
        from wilson import Wilson
        self._flavio = flavio
        self._wilson = Wilson

    def _key(self, decay: str, observable: str, q2_lo: float, q2_hi: float,
             dC9: float = 0.0) -> str:
        return f"{decay}|{observable}|{q2_lo:.4f}|{q2_hi:.4f}|dC9={dC9:+.4f}"

    def _flavio_pred(self, decay: str, observable: str, q2_lo: float, q2_hi: float,
                     dC9: float) -> float:
        self._ensure_flavio()
        flavio = self._flavio
        Wilson = self._wilson
        decay_str = _decay_to_flavio(decay)
        obs_str = _OBS_TO_FLAVIO.get(observable, f"<{observable}>")
        full = f"{obs_str}({decay_str})"
        if abs(dC9) < 1e-12:
            return float(flavio.sm_prediction(full, q2min=q2_lo, q2max=q2_hi))
        wc = Wilson({"C9_bsmumu": dC9}, scale=4.8, eft="WET", basis="flavio")
        return float(flavio.np_prediction(full, wc, q2min=q2_lo, q2max=q2_hi))

    def _get(self, decay: str, observable: str, q2_lo: float, q2_hi: float,
             dC9: float = 0.0) -> float:
        key = self._key(decay, observable, q2_lo, q2_hi, dC9)
        if key not in self._cache:
            val = self._flavio_pred(decay, observable, q2_lo, q2_hi, dC9)
            self._cache[key] = val
            self._dirty = True
        return self._cache[key]

    def sm_value(self, observable: str, q2_lo: float, q2_hi: float,
                 decay: str = "B0->K*mumu") -> float:
        return self._get(decay, observable, q2_lo, q2_hi, 0.0)

    def dO_dC9(self, observable: str, q2_lo: float, q2_hi: float,
               decay: str = "B0->K*mumu", *, delta: float = 0.5) -> float:
        plus = self._get(decay, observable, q2_lo, q2_hi, +delta)
        minus = self._get(decay, observable, q2_lo, q2_hi, -delta)
        return (plus - minus) / (2 * delta)

    def predict(self, observable: str, q2_lo: float, q2_hi: float, c9_value: float,
                *, c9_sm: float = 4.27, decay: str = "B0->K*mumu") -> float:
        sm = self.sm_value(observable, q2_lo, q2_hi, decay)
        slope = self.dO_dC9(observable, q2_lo, q2_hi, decay)
        return sm + slope * (c9_value - c9_sm)

    def predict_vector(
        self,
        observables: Iterable[str],
        q2_los: Iterable[float],
        q2_his: Iterable[float],
        c9_values: np.ndarray | float,
        *,
        c9_sm: float = 4.27,
        decay: str = "B0->K*mumu",
    ) -> np.ndarray:
        observables = list(observables)
        q2_los = list(q2_los)
        q2_his = list(q2_his)
        if np.isscalar(c9_values):
            c9_arr = np.full(len(observables), float(c9_values))
        else:
            c9_arr = np.asarray(c9_values, dtype=float)
        out = np.empty(len(observables), dtype=float)
        for i, (obs, lo, hi) in enumerate(zip(observables, q2_los, q2_his)):
            out[i] = self.predict(obs, lo, hi, float(c9_arr[i]),
                                  c9_sm=c9_sm, decay=decay)
        return out

    def precompute(self, observables: Iterable[str],
                   q2_bins: Iterable[tuple[float, float]],
                   *, decay: str = "B0->K*mumu", delta: float = 0.5,
                   verbose: bool = False) -> None:
        """Warm the cache for a list of (observable, bin) pairs."""
        for obs in observables:
            for lo, hi in q2_bins:
                if verbose:
                    print(f"  precompute {decay} {obs} [{lo}, {hi}]")
                self.sm_value(obs, lo, hi, decay)
                self.dO_dC9(obs, lo, hi, decay, delta=delta)
        self.flush()

    def flush(self) -> None:
        if not self._dirty:
            return
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.cache_path.write_text(json.dumps(self._cache, indent=1, sort_keys=True))
        self._dirty = False


_DEFAULT = FlavioPredictor()

def default() -> FlavioPredictor:
    return _DEFAULT
