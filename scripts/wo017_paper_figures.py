"""WO-017 — Paper figures F1, F2, F3.

Generates three figures for the headline non-linear analysis:
    F1 - Kernel shape kappa(q^2) over [0, 19] GeV^2 with bin centres.
    F2 - Bin pulls for the four LHCb 2025 angular observables under
         the non-linear FREE_C9 and VFD fits.
    F3 - Cross-dataset and cross-channel amplitude bar chart with
         non-linear bootstrap-style error bars.

All figures are saved into the paper's local figures/ directory so
that `\\graphicspath{{figures/}}` finds them.

Usage: PYTHONPATH=src python3 scripts/wo017_paper_figures.py
"""
from __future__ import annotations

import csv
import pathlib
import warnings

warnings.simplefilter("ignore")

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from vfd_b_anomaly import wo010_universality, vfd_closure
from vfd_b_anomaly.constants import J_PSI_Q2, PSI2S_Q2, PHI


REPO = pathlib.Path(__file__).resolve().parents[1]
OUT = REPO / "paper" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def _set_style():
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "axes.titlesize": 11,
    })


# ---------------------------------------------------------------------------
# F1 — Kernel shape over [0, 19] GeV^2
# ---------------------------------------------------------------------------

def fig_F1_kernel_shape():
    _set_style()
    q2 = np.linspace(0.05, 19.0, 400)
    kernel = wo010_universality.frozen_kernel_at_bin_centres(q2)
    x = vfd_closure.kappa_coordinate(q2)
    cont = np.exp(-np.abs(x) / PHI)

    lhcb_centres = np.array([
        0.5 * (lo + hi) for lo, hi in
        [(0.10, 0.98), (1.10, 2.50), (2.50, 4.00), (4.00, 6.00),
         (6.00, 8.00), (11.00, 12.50), (15.00, 17.00), (17.00, 19.00)]
    ])
    kappa_at_centres = wo010_universality.frozen_kernel_at_bin_centres(lhcb_centres)

    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.plot(q2, cont, "--", color="tab:gray", linewidth=1.4,
            label=r"continuum $e^{-|x|/\varphi}$ (Layer 1)")
    ax.plot(q2, kernel, "-", color="tab:blue", linewidth=2.0,
            label=r"$V_{600}$ shell-mean (Layer 3, used in fits)")
    ax.scatter(lhcb_centres, kappa_at_centres, color="tab:red",
               zorder=5, s=40, label="LHCb 2025 bin centres")
    ax.axvline(J_PSI_Q2, color="goldenrod", linestyle=":",
               linewidth=1.0, alpha=0.7)
    ax.axvline(PSI2S_Q2, color="goldenrod", linestyle=":",
               linewidth=1.0, alpha=0.7)
    ax.text(J_PSI_Q2 + 0.1, ax.get_ylim()[1] * 0.97, r"$J/\psi$",
            color="goldenrod", fontsize=8, va="top")
    ax.text(PSI2S_Q2 + 0.1, ax.get_ylim()[1] * 0.97, r"$\psi(2S)$",
            color="goldenrod", fontsize=8, va="top")
    ax.set_xlabel(r"$q^2$ (GeV$^2$)")
    ax.set_ylabel(r"$\kappa(q^2)$ (dimensionless)")
    ax.set_xlim(0, 19)
    ax.set_title(r"Geometry-derived response kernel on the LHCb $q^2$ window")
    ax.legend(loc="upper right", framealpha=0.95)
    fig.tight_layout()
    out = OUT / "fig_F1_kernel_shape.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# F2 — bin pulls under non-linear refit
# ---------------------------------------------------------------------------

def fig_F2_bin_pulls():
    _set_style()
    import flavio
    from wilson import Wilson
    from vfd_b_anomaly.wo014_cross_dataset import load_lhcb_2025

    data = load_lhcb_2025()
    df = data["observables"]
    obs_list = df["observable"].tolist()
    q2_lo = df["q2_lo"].to_numpy()
    q2_hi = df["q2_hi"].to_numpy()
    values = df["value"].to_numpy(dtype=float)
    errors = np.sqrt(df["stat_err"].to_numpy() ** 2 +
                      df["syst_err"].to_numpy() ** 2)
    q2_centres = 0.5 * (q2_lo + q2_hi)

    DC9_FREE = -1.0025
    A_VFD = 1.135
    OBS_TO_FLAVIO = {"P5p": "<P5p>", "P4p": "<P4p>", "P1": "<P1>", "P2": "<P2>"}
    full = lambda o: f"{OBS_TO_FLAVIO[o]}(B0->K*mumu)"

    wc_free = Wilson({"C9_bsmumu": DC9_FREE}, scale=4.8,
                     eft="WET", basis="flavio")
    pred_free = np.array([
        float(flavio.np_prediction(full(o), wc_free, q2min=lo, q2max=hi))
        for o, lo, hi in zip(obs_list, q2_lo, q2_hi)
    ])
    kappa = wo010_universality.frozen_kernel_at_bin_centres(q2_centres)
    dC9_perbin = -A_VFD * kappa
    pred_vfd = np.zeros_like(values)
    for i, (o, lo, hi) in enumerate(zip(obs_list, q2_lo, q2_hi)):
        wc = Wilson({"C9_bsmumu": float(dC9_perbin[i])}, scale=4.8,
                    eft="WET", basis="flavio")
        pred_vfd[i] = float(flavio.np_prediction(full(o), wc,
                                                 q2min=lo, q2max=hi))

    pulls_free = (values - pred_free) / errors
    pulls_vfd = (values - pred_vfd) / errors

    obs_label = {"P5p": r"$P_5'$", "P4p": r"$P_4'$",
                 "P1": r"$P_1$", "P2": r"$P_2$"}
    obs_unique = ["P5p", "P4p", "P1", "P2"]

    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.0), sharex=True)
    axes = axes.flat
    for ax, oname in zip(axes, obs_unique):
        mask = np.array([o == oname for o in obs_list])
        x_axis = q2_centres[mask]
        ax.bar(x_axis - 0.18, pulls_free[mask], width=0.34,
               color="tab:orange", label=r"FREE\_C9", alpha=0.85)
        ax.bar(x_axis + 0.18, pulls_vfd[mask], width=0.34,
               color="tab:blue", label="VFD", alpha=0.85)
        ax.axhline(0, color="k", linewidth=0.6)
        ax.axhline(+1, color="gray", linewidth=0.4, linestyle="--", alpha=0.6)
        ax.axhline(-1, color="gray", linewidth=0.4, linestyle="--", alpha=0.6)
        ax.set_ylim(-3.5, 3.5)
        ax.set_title(obs_label[oname])
        ax.set_ylabel("pull " + r"$(O - O_{\rm pred})/\sigma$")
    for ax in axes[2:]:
        ax.set_xlabel(r"$q^2$ (GeV$^2$)")
    axes[0].legend(loc="upper right")
    fig.suptitle("LHCb 2025 bin pulls under non-linear "
                 r"FREE\_C9 and VFD fits (32 bins, 4 obs)",
                 y=1.0, fontsize=11)
    fig.tight_layout()
    out = OUT / "fig_F2_bin_pulls.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# F3 — cross-dataset amplitude bar chart
# ---------------------------------------------------------------------------

def fig_F3_cross_dataset_A():
    _set_style()
    src = REPO / "reports" / "wo016d_nonlinear_xdataset.csv"
    rows = list(csv.DictReader(src.open()))
    vfd_rows = [r for r in rows if r["model"] == "VFD_GREEN_600CELL"]

    label_map = {
        "LHCb-2015":         "LHCb 2015\n" + r"$B^0\to K^{*0}$",
        "LHCb-2021-Kstplus": "LHCb 2021\n" + r"$B^+\to K^{*+}$",
        "CMS-2025-noP4p":    "CMS 2025\n" + r"$B^0\to K^{*0}$ (no $P_4'$)",
        "LHCb-2025":         "LHCb 2025\n" + r"$B^0\to K^{*0}$",
        "Bs2phi-LHCb-2015":  "LHCb 2015\n" + r"$B_s\to\phi$",
    }
    order = ["LHCb-2015", "LHCb-2021-Kstplus", "CMS-2025-noP4p",
             "LHCb-2025", "Bs2phi-LHCb-2015"]
    A_vals = {r["dataset"]: float(r["A_or_DC9"]) for r in vfd_rows}
    aic_diff = {r["dataset"]: float(r["delta_aic_vs_FREE_C9"])
                 for r in vfd_rows}

    fig = plt.figure(figsize=(8.5, 4.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1.1], wspace=0.05)
    ax_main = fig.add_subplot(gs[0])
    ax_phi = fig.add_subplot(gs[1])

    p_basis = order[:4]
    x_main = np.arange(len(p_basis))
    A_main = [A_vals[ds] for ds in p_basis]
    aic_main = [aic_diff[ds] for ds in p_basis]
    colors = ["tab:green" if a < 0 else "tab:orange" for a in aic_main]
    bars_main = ax_main.bar(x_main, A_main, color=colors,
                             edgecolor="black", linewidth=0.6)
    ax_main.axhline(0, color="k", linewidth=0.5)
    ax_main.set_xticks(x_main)
    ax_main.set_xticklabels([label_map[d] for d in p_basis], fontsize=8)
    ax_main.set_ylabel(r"fitted amplitude $A$ (non-linear)")
    ax_main.set_title(r"$P$-basis $B\to K^{*}\mu\mu$ fits", fontsize=10)
    ax_main.set_ylim(0, 3.5)
    for b, a, da in zip(bars_main, A_main, aic_main):
        ax_main.text(b.get_x() + b.get_width() / 2, a + 0.06,
                     fr"$A={a:.2f}$" + "\n" +
                     fr"$\Delta\mathrm{{AIC}}={da:+.2f}$",
                     ha="center", va="bottom", fontsize=8)

    ds_phi = "Bs2phi-LHCb-2015"
    bars_phi = ax_phi.bar([0], [A_vals[ds_phi]],
                           color="tab:green" if aic_diff[ds_phi] < 0
                           else "tab:orange",
                           edgecolor="black", linewidth=0.6)
    ax_phi.set_ylim(0, 7.5)
    ax_phi.set_xticks([0])
    ax_phi.set_xticklabels([label_map[ds_phi]], fontsize=8)
    ax_phi.set_title(r"$S$-basis cross-channel", fontsize=10)
    ax_phi.text(0, A_vals[ds_phi] + 0.15,
                fr"$A={A_vals[ds_phi]:.2f}$" + "\n" +
                fr"$\Delta\mathrm{{AIC}}={aic_diff[ds_phi]:+.2f}$",
                ha="center", va="bottom", fontsize=8)
    pred = 1.135 * 2.2  # P-basis * basis factor
    ax_phi.axhline(pred, color="tab:gray", linestyle="--", linewidth=1.0)
    ax_phi.text(0.0, pred + 0.1, fr"basis prediction $\sim {pred:.1f}$",
                ha="center", va="bottom", fontsize=7, color="tab:gray")

    fig.suptitle("Cross-dataset and cross-channel non-linear amplitudes",
                 y=1.0, fontsize=11)
    fig.tight_layout()
    out = OUT / "fig_F3_cross_dataset_A.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"))
    plt.close(fig)
    return out


def main():
    f1 = fig_F1_kernel_shape()
    print(f"wrote {f1.relative_to(REPO)}")
    f2 = fig_F2_bin_pulls()
    print(f"wrote {f2.relative_to(REPO)}")
    f3 = fig_F3_cross_dataset_A()
    print(f"wrote {f3.relative_to(REPO)}")


if __name__ == "__main__":
    main()
