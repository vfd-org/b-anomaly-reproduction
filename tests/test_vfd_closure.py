from __future__ import annotations

import inspect

import numpy as np
import pytest

from vfd_b_anomaly import vfd_closure
from vfd_b_anomaly.vfd_closure import (
    ALL_MODES,
    MODE_CONSTANT,
    MODE_LOG_PHASE,
    MODE_THRESHOLD,
    PARAM_COUNT,
    closure_residual_q2,
    fit_vfd_closure,
)


@pytest.mark.parametrize("mode,expected_k", [(MODE_CONSTANT, 1), (MODE_LOG_PHASE, 2), (MODE_THRESHOLD, 2)])
def test_param_count_is_fixed_per_mode(mode, expected_k):
    assert PARAM_COUNT[mode] == expected_k


@pytest.mark.parametrize("mode", ALL_MODES)
def test_param_count_does_not_scale_with_grid_size(mode, synthetic_data_factory):
    """No hidden per-bin freedom: param count must be invariant under grid expansion."""
    data_small = synthetic_data_factory(delta_c9_true=-1.0)
    res_small = fit_vfd_closure(data_small, mode=mode)

    # Replicate the grid by duplicating observables (more rows -> still same n_params)
    data_big = synthetic_data_factory(delta_c9_true=-1.0,
                                      observables=("FL", "AFB", "P5p", "BR"))
    duplicated = data_big["observables"].copy()
    duplicated = duplicated.iloc[np.repeat(np.arange(len(duplicated)), 2)].reset_index(drop=True)
    data_big["observables"] = duplicated
    res_big = fit_vfd_closure(data_big, mode=mode)

    assert res_small.n_params == res_big.n_params == PARAM_COUNT[mode]


def test_residual_constant_returns_uniform_negative_delta_c9():
    q2 = np.array([1.0, 4.0, 8.0, 17.0])
    out = closure_residual_q2(q2, amplitude=0.5, mode=MODE_CONSTANT)
    assert out.shape == q2.shape
    np.testing.assert_allclose(out, -0.5)


def test_residual_log_phase_is_bounded():
    q2 = np.linspace(0.1, 19.0, 100)
    out = closure_residual_q2(q2, amplitude=1.0, mode=MODE_LOG_PHASE, theta=0.0)
    # |R(q^2)| <= 1 in log-phase mode, so |Δ| <= amplitude.
    assert np.max(np.abs(out)) <= 1.0 + 1e-12


def test_residual_threshold_localised_near_charm_resonances():
    q2 = np.linspace(0.1, 19.0, 200)
    out = closure_residual_q2(q2, amplitude=1.0, mode=MODE_THRESHOLD, sigma=1.0)
    # The residual should be largest in magnitude in the 9-14 GeV^2 region.
    abs_out = np.abs(out)
    in_window = (q2 > 8.0) & (q2 < 14.5)
    assert abs_out[in_window].max() > abs_out[~in_window].max()


def test_amplitude_negative_rejected():
    q2 = np.array([4.0])
    with pytest.raises(ValueError, match="amplitude must be >= 0"):
        closure_residual_q2(q2, amplitude=-0.1, mode=MODE_CONSTANT)


def test_unknown_mode_rejected():
    with pytest.raises(ValueError, match="unknown VFD mode"):
        closure_residual_q2(np.array([1.0]), amplitude=0.1, mode="quantum_woo")


def test_q2_must_be_positive():
    with pytest.raises(ValueError, match="strictly positive"):
        closure_residual_q2(np.array([0.0, 1.0]), amplitude=0.1, mode=MODE_CONSTANT)


@pytest.mark.parametrize("mode", ALL_MODES)
def test_fit_vfd_closure_recovers_negative_delta_c9_on_anomalous_synthetic(mode, synthetic_data_factory):
    data = synthetic_data_factory(delta_c9_true=-1.0)
    res = fit_vfd_closure(data, mode=mode)
    assert res.success
    # The closure model should produce a NEGATIVE effective ΔC9 (sign prior + anomaly direction).
    assert res.effective_delta_c9_mean < 0.0


def test_fit_vfd_closure_does_not_invent_anomaly_on_sm_data(synthetic_data_factory):
    data = synthetic_data_factory(delta_c9_true=0.0)
    res = fit_vfd_closure(data, mode=MODE_CONSTANT)
    # On true-SM data the amplitude should be near zero.
    assert res.params["amplitude"] < 0.25


def test_module_has_no_per_bin_free_fitting_helpers():
    """Source-level guard: nothing in vfd_closure should declare per-bin parameter arrays of length n_data."""
    src = inspect.getsource(vfd_closure)
    # We allow per-row evaluation arrays (delta_arr) but forbid optimiser bounds
    # whose length depends on grid size. The PARAM_COUNT dict is the only authority.
    for mode, k in PARAM_COUNT.items():
        # Only sanity: the seed function returns vectors of length k, not n_data.
        seed = vfd_closure._seed(mode)
        assert len(seed) == k
    # Forbid an obvious anti-pattern token that suggests per-bin freedom.
    assert "per_bin_amplitude" not in src
