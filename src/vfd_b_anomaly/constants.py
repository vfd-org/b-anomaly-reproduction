"""Physical and project-level constants.

Wilson-coefficient values are the standard low-scale (mu = m_b) SM benchmarks
quoted across the B -> K* mu mu literature. They are conventional reference
values, not measurements.
"""

from __future__ import annotations

import math

# Standard Model Wilson coefficients at mu = m_b (NNLL), conventional reference values.
C7_SM = -0.295
C9_SM = 4.27
C10_SM = -4.17

# q^2 reference scale used in log-phase residual (GeV^2).
Q2_REF_GEV2 = 4.0

# Charm-loop threshold scales in q^2 (GeV^2), narrow-resonance vetoes typical of LHCb.
J_PSI_Q2 = 9.59  # m_J/psi^2
PSI2S_Q2 = 13.59  # m_psi(2S)^2

# Standard charmonium veto windows excluded from anomaly fits.
CHARMONIUM_VETOES_GEV2 = [(8.0, 11.0), (12.5, 15.0)]

# Golden ratio - used as a default phase scale in log-phase closure residual.
PHI = (1.0 + math.sqrt(5.0)) / 2.0

# Project-level tags for content provenance in reports.
PROVENANCE_REPRODUCED = "reproduced_from_paper"
PROVENANCE_APPROXIMATED = "approximated_from_public_plot_or_table"
PROVENANCE_INFERRED = "inferred_by_model"
PROVENANCE_VFD = "vfd_hypothesis"
PROVENANCE_PLACEHOLDER = "placeholder_synthetic"

# Required column schema for processed observable tables.
REQUIRED_OBSERVABLE_COLUMNS = (
    "q2_lo",
    "q2_hi",
    "observable",
    "value",
    "stat_err",
    "syst_err",
    "provenance",
)

# Observables this project models. Adding more requires extending sm_baseline sensitivities.
SUPPORTED_OBSERVABLES = ("FL", "AFB", "P5p", "P4p", "P1", "P2", "BR")
