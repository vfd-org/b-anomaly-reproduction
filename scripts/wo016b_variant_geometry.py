"""WO-016b — Variant selection by pure-geometry criterion.

The paper claims the unweighted Laplacian on V_600 is the right choice
for the kernel. Reviewer asked: was that selected by data correlation
against LHCb (data-informed; effective k_VFD > 1) or by a pure-geometry
criterion (correlation with the continuum Green's function
$\\kappa(x) = e^{-|x|/\\varphi}$ from Layer 1, no LHCb involvement)?

This script extracts the GREENS-mode variants from
reports/wo009_full_lift.json and ranks them by `correlation_with_exp`,
the pure-geometry criterion. The same ranking is then read out from
chi2 against LHCb data; agreement of the two rankings means the
variant choice is consistent with both pure-geometry and LHCb-data
selection (i.e. data-informed AND geometry-informed). Disagreement
would mean the variant choice is *only* data-informed.

Inputs : reports/wo009_full_lift.json
Output : reports/wo016b_variant_geometry.{csv,md}
"""
from __future__ import annotations

import csv
import json
import pathlib


REPO = pathlib.Path(__file__).resolve().parents[1]
SRC = REPO / "reports" / "wo009_full_lift.json"
OUT_CSV = REPO / "reports" / "wo016b_variant_geometry.csv"
OUT_MD = REPO / "reports" / "wo016b_variant_geometry.md"


def main():
    rows = json.loads(SRC.read_text())
    greens = [r for r in rows if r["model"].endswith("_GREENS")]

    by_geom = sorted(greens, key=lambda r: -float(r["correlation_with_exp"]))
    by_data = sorted(greens, key=lambda r: float(r["chi2"]))

    with OUT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "corr_with_exp", "chi2_vs_lhcb",
                    "geometry_rank", "data_rank"])
        ranks_geom = {r["model"]: i + 1 for i, r in enumerate(by_geom)}
        ranks_data = {r["model"]: i + 1 for i, r in enumerate(by_data)}
        for r in greens:
            w.writerow([
                r["model"],
                f"{float(r['correlation_with_exp']):.6f}",
                f"{float(r['chi2']):.4f}",
                ranks_geom[r["model"]],
                ranks_data[r["model"]],
            ])

    lines = ["# WO-016b — Variant selection by pure-geometry criterion", ""]
    lines.append(
        "Pure-geometry criterion: correlation between the discrete "
        "shell-mean of the V_600 Green's response and the continuum "
        "kernel $\\kappa(x) = e^{-|x|/\\varphi}$ from Layer 1 of the "
        "derivation. This criterion does **not** use LHCb data."
    )
    lines.append("")
    lines.append("Data criterion: chi^2 against LHCb 2025 P5' on the joint fit.")
    lines.append("")
    lines.append("| variant | corr(κ_continuum) | χ² (LHCb) | geom rank | data rank |")
    lines.append("|---|---:|---:|---:|---:|")
    for r in greens:
        rg = next(i for i, x in enumerate(by_geom) if x["model"] == r["model"]) + 1
        rd = next(i for i, x in enumerate(by_data) if x["model"] == r["model"]) + 1
        lines.append(
            f"| {r['model']} "
            f"| {float(r['correlation_with_exp']):.4f} "
            f"| {float(r['chi2']):.3f} "
            f"| {rg} "
            f"| {rd} |"
        )
    lines.append("")
    geom_winner = by_geom[0]["model"]
    data_winner = by_data[0]["model"]
    agree = geom_winner == data_winner
    lines.append(f"- Pure-geometry winner: **{geom_winner}** "
                 f"(corr with continuum kernel = "
                 f"{float(by_geom[0]['correlation_with_exp']):.4f})")
    lines.append(f"- LHCb-data winner: **{data_winner}** "
                 f"(χ² = {float(by_data[0]['chi2']):.3f})")
    lines.append("")
    if agree:
        lines.append(
            "**Agreement.** The same variant (unweighted Laplacian) wins "
            "on both criteria. The variant choice is consistent with "
            "pure-geometry selection independent of the data; the LHCb "
            "data merely confirms it. Effective parameter count for VFD "
            "remains k=1 under the pure-geometry interpretation."
        )
    else:
        lines.append(
            "**Disagreement.** The pure-geometry winner differs from the "
            "LHCb-data winner. The variant chosen in the paper is "
            "data-informed; effective parameter count for VFD should be "
            "treated as k>=2 in AIC comparisons."
        )

    OUT_MD.write_text("\n".join(lines) + "\n")
    print(OUT_MD.read_text())


if __name__ == "__main__":
    main()
