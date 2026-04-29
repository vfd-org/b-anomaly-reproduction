"""WO-016a — Akaike-weight stacking across cross-dataset and cross-channel fits.

Per dataset, w_M = exp(-Delta_M/2) / sum_M exp(-Delta_M/2)
where Delta_M = AIC_M - AIC_min over the model set {FREE_C9, VFD_GREEN_600CELL}.

Stacked weight (assuming dataset independence under the null) is the
product across datasets, renormalised over models. The relative
log-evidence is the sum of the AIC differences halved.

Inputs : reports/wo016d_nonlinear_xdataset.csv
Output : reports/wo016a_akaike_stack.{csv,md}
"""
from __future__ import annotations

import csv
import math
import pathlib
from collections import defaultdict


REPO = pathlib.Path(__file__).resolve().parents[1]
INPUTS = [REPO / "reports" / "wo016d_nonlinear_xdataset.csv"]
OUT_CSV = REPO / "reports" / "wo016a_akaike_stack.csv"
OUT_MD = REPO / "reports" / "wo016a_akaike_stack.md"

# Datasets in the headline cross-dataset claim. CMS-2025 is excluded in
# favour of CMS-2025-noP4p because the P4'-included CMS fit has known
# convention-mismatch chi^2 inflation; the noP4p fit is the clean comparison.
HEADLINE_DATASETS = [
    "LHCb-2015",
    "LHCb-2021-Kstplus",
    "CMS-2025-noP4p",
    "LHCb-2025",
    "Bs2phi-LHCb-2015",
]


def load_rows():
    rows = []
    for p in INPUTS:
        with p.open() as f:
            for r in csv.DictReader(f):
                rows.append(r)
    return rows


def main():
    rows = load_rows()
    by_dataset = defaultdict(dict)
    for r in rows:
        by_dataset[r["dataset"]][r["model"]] = r

    per_dataset = []
    for ds in HEADLINE_DATASETS:
        models = by_dataset[ds]
        aics = {m: float(v["aic"]) for m, v in models.items()}
        aic_min = min(aics.values())
        deltas = {m: a - aic_min for m, a in aics.items()}
        unnorm = {m: math.exp(-d / 2) for m, d in deltas.items()}
        Z = sum(unnorm.values())
        weights = {m: u / Z for m, u in unnorm.items()}
        per_dataset.append((ds, deltas, weights))

    log_w_total = defaultdict(float)
    for _, deltas, _ in per_dataset:
        for m, d in deltas.items():
            log_w_total[m] += -d / 2
    lmax = max(log_w_total.values())
    unnorm_stack = {m: math.exp(lw - lmax) for m, lw in log_w_total.items()}
    Z = sum(unnorm_stack.values())
    stacked = {m: u / Z for m, u in unnorm_stack.items()}

    with OUT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dataset", "model", "AIC_delta", "akaike_weight"])
        for ds, deltas, weights in per_dataset:
            for m in deltas:
                w.writerow([ds, m, f"{deltas[m]:.4f}", f"{weights[m]:.6f}"])
        w.writerow([])
        w.writerow(["STACKED", "model", "log_evidence_relative", "stacked_weight"])
        for m, lw in log_w_total.items():
            w.writerow(["STACKED", m, f"{lw:.4f}", f"{stacked[m]:.6f}"])

    lines = ["# WO-016a — Akaike-weight stack across five fits", ""]
    lines.append("Per-dataset AIC deltas and Akaike weights, plus stacked weight.")
    lines.append("Stacking assumes independence under the null hypothesis "
                 "(the five datasets share no observation-level information).")
    lines.append("")
    lines.append("| dataset | FREE_C9 ΔAIC | VFD ΔAIC | w(FREE_C9) | w(VFD) |")
    lines.append("|---|---:|---:|---:|---:|")
    for ds, deltas, weights in per_dataset:
        lines.append(
            f"| {ds} "
            f"| {deltas.get('FREE_C9', 0):.3f} "
            f"| {deltas.get('VFD_GREEN_600CELL', 0):.3f} "
            f"| {weights.get('FREE_C9', 0):.4f} "
            f"| {weights.get('VFD_GREEN_600CELL', 0):.4f} |"
        )
    lines.append("")
    lines.append("## Stacked")
    lines.append("")
    lines.append(f"- log-evidence(FREE_C9) − log-evidence(VFD) = "
                 f"{log_w_total['FREE_C9'] - log_w_total['VFD_GREEN_600CELL']:.3f}")
    lines.append(f"- Total ΔAIC sum (FREE_C9 vs VFD): "
                 f"{-2*(log_w_total['FREE_C9'] - log_w_total['VFD_GREEN_600CELL']):.3f}")
    lines.append("")
    lines.append("| model | stacked Akaike weight |")
    lines.append("|---|---:|")
    for m, w in sorted(stacked.items(), key=lambda kv: -kv[1]):
        lines.append(f"| {m} | {w:.4f} |")
    lines.append("")

    null_p_signs = 1 / 2 ** len(per_dataset)
    lines.append(f"Auxiliary check: under the null hypothesis "
                 f"P(VFD lower AIC on all {len(per_dataset)} fits) = "
                 f"$1/2^{{{len(per_dataset)}}}$ = {null_p_signs:.4f}.")

    OUT_MD.write_text("\n".join(lines) + "\n")
    print(OUT_MD.read_text())


if __name__ == "__main__":
    main()
