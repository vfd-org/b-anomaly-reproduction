"""Placeholder for digitised paper tables.

Until machine-readable supplementary data from arXiv:2512.18053 is obtained,
this module exposes ONLY a placeholder bin grid and provenance marker.
It must NOT contain any paper-reported best-fit numbers presented as if
they were independently derived. See feedback_reproduce_first memory.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import PROVENANCE_PLACEHOLDER

# LHCb config_2 standard exclusive q^2 binning (arXiv:2512.18053).
# These edges are taken from the LHCb HEPData supplementary tables
# (config_2 = partially massive model, optimised Pi basis, standard
# binning). The two integrated re-bins [1.1, 6.0] and [15.0, 19.0]
# present in the HEPData are deliberately excluded here because they
# overlap their constituent bins and thus are NOT independent measurements.
# Edge provenance: reproduced_from_paper (HEPData ins3094698, config_2).
DEFAULT_Q2_BIN_EDGES_GEV2: list[tuple[float, float]] = [
    (0.06, 0.98),
    (1.10, 2.50),
    (2.50, 4.00),
    (4.00, 6.00),
    (6.00, 8.00),
    (11.00, 12.50),
    (15.00, 17.00),
    (17.00, 19.00),
]

# HEPData bin-index mapping: position in DEFAULT_Q2_BIN_EDGES_GEV2 above
# is NOT the same as the original HEPData bin index (HEPData orders
# 11.0-12.5 last among exclusive bins, after 15-17 and 17-19). This
# mapping is used by hepdata_ingest to look up the right HEPData entries.
HEPDATA_BIN_INDEX_FOR_DEFAULT: dict[int, int] = {
    0: 0,  # [0.06, 0.98]
    1: 1,  # [1.10, 2.50]
    2: 2,  # [2.50, 4.00]
    3: 3,  # [4.00, 6.00]
    4: 4,  # [6.00, 8.00]
    5: 7,  # [11.00, 12.50]  (HEPData index 7 in the original ordering)
    6: 5,  # [15.00, 17.00]
    7: 6,  # [17.00, 19.00]
}


@dataclass(frozen=True)
class PaperTableStatus:
    """Records whether digitised paper data is loaded.

    All paper-reported numerics live elsewhere or in CSV under data/processed.
    This dataclass exists so other modules can ask 'do we actually have paper
    data?' without reading hidden state.
    """

    have_observable_table: bool
    have_covariance: bool
    have_reported_delta_c9: bool
    notes: str


def status() -> PaperTableStatus:
    """Default status: nothing loaded. Update only when real digitised tables are added."""
    return PaperTableStatus(
        have_observable_table=False,
        have_covariance=False,
        have_reported_delta_c9=False,
        notes=(
            "No digitised paper-table contents committed. Project runs in "
            "placeholder/synthetic mode; see source_manifest.json."
        ),
    )


def default_bin_centers() -> np.ndarray:
    """Bin midpoints for the conventional anomaly grid (provenance: placeholder)."""
    edges = np.asarray(DEFAULT_Q2_BIN_EDGES_GEV2, dtype=float)
    return 0.5 * (edges[:, 0] + edges[:, 1])


def default_bin_provenance() -> str:
    return PROVENANCE_PLACEHOLDER
