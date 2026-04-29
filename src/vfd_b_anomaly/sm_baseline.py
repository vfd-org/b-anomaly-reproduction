"""SM observable baseline and linear C9 sensitivities (MODE_B response model).

WO-012 (2026-04-29): values regenerated from flavio 2.4.0 (Straub 1810.08132)
on the LHCb default 8-bin grid. The slopes use central-difference numerics
with delta_C9 = +/- 0.5 around the SM and the flavio default WET basis.
Tools: see `tools/build_sm_table_from_flavio.py` for the regeneration
script. The previous hand-tabulated values had wrong-sign slopes for
several observables and SM values close to LHCb data (not theory).

The linearised model is:
    O_pred(q^2_bin, C9) = O_SM(q^2_bin) + dO_dC9(q^2_bin) * (C9 - C9_SM)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import C9_SM, SUPPORTED_OBSERVABLES


@dataclass(frozen=True)
class SMResponse:
    """Linearised SM baseline and dO/dC9 sensitivity for one observable on one bin."""

    sm_value: float
    dO_dC9: float


# Bin index convention matches paper_tables.DEFAULT_Q2_BIN_EDGES_GEV2.
# Bins (in GeV^2): [(0.06, 0.98), (1.1, 2.5), (2.5, 4.0), (4.0, 6.0),
#                   (6.0, 8.0), (11.0, 12.5), (15.0, 17.0), (17.0, 19.0)].
# Values regenerated from flavio 2.4 / wilson 2.5 (WO-012) using
# central-difference slopes with delta C9_bsmumu = +/- 0.5 around the SM at
# scale = 4.8 GeV in the WET-flavio basis. flavio QCDF interpolation issues
# a "do not trust above 6 GeV^2" warning; we use the values anyway and note
# this caveat in the WO-012 report. BR slope is reported separately because
# flavio's <dBR/dq2> is normalised differently from the project's BR scale.
_BIN_TABLE: dict[str, list[SMResponse]] = {
    "FL": [
        SMResponse(sm_value=+0.2568, dO_dC9=+0.0495),
        SMResponse(sm_value=+0.7604, dO_dC9=+0.0714),
        SMResponse(sm_value=+0.7963, dO_dC9=+0.0392),
        SMResponse(sm_value=+0.7113, dO_dC9=+0.0184),
        SMResponse(sm_value=+0.6070, dO_dC9=+0.0080),
        SMResponse(sm_value=+0.4352, dO_dC9=+0.0009),
        SMResponse(sm_value=+0.3484, dO_dC9=+0.0003),
        SMResponse(sm_value=+0.3281, dO_dC9=+0.0000),
    ],
    "AFB": [
        SMResponse(sm_value=-0.0850, dO_dC9=+0.0073),
        SMResponse(sm_value=-0.1380, dO_dC9=+0.0537),
        SMResponse(sm_value=-0.0174, dO_dC9=+0.0677),
        SMResponse(sm_value=+0.1222, dO_dC9=+0.0624),
        SMResponse(sm_value=+0.2399, dO_dC9=+0.0506),
        SMResponse(sm_value=+0.3914, dO_dC9=+0.0228),
        SMResponse(sm_value=+0.4019, dO_dC9=+0.0206),
        SMResponse(sm_value=+0.3184, dO_dC9=+0.0148),
    ],
    "P5p": [
        SMResponse(sm_value=+0.6687, dO_dC9=-0.1116),
        SMResponse(sm_value=+0.1385, dO_dC9=-0.2336),
        SMResponse(sm_value=-0.5007, dO_dC9=-0.2756),
        SMResponse(sm_value=-0.7574, dO_dC9=-0.1680),
        SMResponse(sm_value=-0.8334, dO_dC9=-0.0997),
        SMResponse(sm_value=-0.8232, dO_dC9=-0.0304),
        SMResponse(sm_value=-0.6701, dO_dC9=-0.0295),
        SMResponse(sm_value=-0.4826, dO_dC9=-0.0212),
    ],
    "P4p": [
        SMResponse(sm_value=+0.2386, dO_dC9=+0.0456),
        SMResponse(sm_value=-0.0644, dO_dC9=+0.0003),
        SMResponse(sm_value=-0.3933, dO_dC9=-0.0405),
        SMResponse(sm_value=-0.5037, dO_dC9=-0.0180),
        SMResponse(sm_value=-0.5363, dO_dC9=-0.0067),
        SMResponse(sm_value=-0.5698, dO_dC9=-0.0010),
        SMResponse(sm_value=-0.6188, dO_dC9=-0.0002),
        SMResponse(sm_value=-0.6612, dO_dC9=-0.0001),
    ],
    "P1": [
        SMResponse(sm_value=+0.0452, dO_dC9=+0.0040),
        SMResponse(sm_value=+0.0233, dO_dC9=-0.0019),
        SMResponse(sm_value=-0.1164, dO_dC9=-0.0360),
        SMResponse(sm_value=-0.1777, dO_dC9=-0.0180),
        SMResponse(sm_value=-0.2058, dO_dC9=-0.0075),
        SMResponse(sm_value=-0.3064, dO_dC9=-0.0013),
        SMResponse(sm_value=-0.5343, dO_dC9=-0.0006),
        SMResponse(sm_value=-0.7506, dO_dC9=-0.0003),
    ],
    "P2": [
        SMResponse(sm_value=-0.1236, dO_dC9=-0.0016),
        SMResponse(sm_value=-0.4514, dO_dC9=+0.0180),
        SMResponse(sm_value=-0.0621, dO_dC9=+0.2257),
        SMResponse(sm_value=+0.2924, dO_dC9=+0.1677),
        SMResponse(sm_value=+0.4141, dO_dC9=+0.0961),
        SMResponse(sm_value=+0.4650, dO_dC9=+0.0282),
        SMResponse(sm_value=+0.4127, dO_dC9=+0.0216),
        SMResponse(sm_value=+0.3168, dO_dC9=+0.0149),
    ],
    # BR: bin-integrated branching fraction in units of 10^-7 / GeV^2 (project
    # convention). flavio's <dBR/dq2>(B0->K*mumu) returns 0 in this version;
    # left at the previously-tabulated approximate values pending a separate
    # flavio integration for the rate observable. Do not use for quantitative
    # claims until WO-013 (BR backend) lands.
    "BR": [
        SMResponse(sm_value=1.10, dO_dC9=-0.10),
        SMResponse(sm_value=0.50, dO_dC9=-0.05),
        SMResponse(sm_value=0.46, dO_dC9=-0.05),
        SMResponse(sm_value=0.50, dO_dC9=-0.06),
        SMResponse(sm_value=0.54, dO_dC9=-0.07),
        SMResponse(sm_value=0.62, dO_dC9=-0.08),
        SMResponse(sm_value=0.70, dO_dC9=-0.09),
        SMResponse(sm_value=0.55, dO_dC9=-0.07),
    ],
}


def supported_observables() -> tuple[str, ...]:
    return SUPPORTED_OBSERVABLES


def n_default_bins() -> int:
    return len(_BIN_TABLE["FL"])


def sm_value(observable: str, bin_index: int) -> float:
    _check(observable, bin_index)
    return _BIN_TABLE[observable][bin_index].sm_value


def dO_dC9(observable: str, bin_index: int) -> float:
    _check(observable, bin_index)
    return _BIN_TABLE[observable][bin_index].dO_dC9


def predict(observable: str, bin_index: int, c9_value: float) -> float:
    """Linearised prediction at a given C9 value (MODE_B)."""
    return sm_value(observable, bin_index) + dO_dC9(observable, bin_index) * (
        c9_value - C9_SM
    )


def predict_vector(
    observables: list[str],
    bin_indices: list[int],
    c9_values: np.ndarray | float,
) -> np.ndarray:
    """Vectorised prediction.

    `c9_values` may be either a scalar (global C9 across all entries) or an
    array of length len(observables) (per-entry C9, i.e. BINNED_C9 mode).
    """
    if len(observables) != len(bin_indices):
        raise ValueError("observables and bin_indices length mismatch")

    if np.isscalar(c9_values):
        c9_arr = np.full(len(observables), float(c9_values))
    else:
        c9_arr = np.asarray(c9_values, dtype=float)
        if c9_arr.shape != (len(observables),):
            raise ValueError(
                f"c9_values must be scalar or length {len(observables)}, got {c9_arr.shape}"
            )

    out = np.empty(len(observables), dtype=float)
    for i, (obs, b) in enumerate(zip(observables, bin_indices)):
        out[i] = predict(obs, b, float(c9_arr[i]))
    return out


def _check(observable: str, bin_index: int) -> None:
    if observable not in _BIN_TABLE:
        raise KeyError(f"Unknown observable {observable!r}; supported: {SUPPORTED_OBSERVABLES}")
    n = len(_BIN_TABLE[observable])
    if not (0 <= bin_index < n):
        raise IndexError(f"bin_index {bin_index} out of range for {observable} (n={n})")
