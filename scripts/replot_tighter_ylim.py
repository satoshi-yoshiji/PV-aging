"""
Re-render Section 2 KM and adjusted-cumulative-incidence figures with a
PER-PLOT data-driven y-axis (max cumulative incidence at xmax + 15% headroom,
snapped up to a tidy step).

The previous implementation derived a SHARED ymax from the raw KM curves and
re-used it for the adjusted-cuminc plots. Adjusted curves are usually tighter
than the raw stratified curves, so this left dead vertical space at the top of
the adjusted plots. This rewrite gives each output its own ymax.

Reads the already-saved source-data CSVs (results/source_data/*) and writes
both the PDF (results/figures/*) and the PNG copy in manuscript/figures/*.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/humgen/diabetes2/users/satoshi/misc/02.aging")
SRC = ROOT / "results" / "source_data"
PDF_DIR = ROOT / "results" / "figures"
PNG_DIR = ROOT / "manuscript" / "figures"

TITLES = {
    "all_cause_mortality": "All-cause mortality",
    "incident_HFRS5":      "HFRS-positive admission (>=5)",
    "zenin_composite":     "Zenin healthspan composite",
    "hip_fracture":        "Incident hip fracture",
    "falls":               "Incident falls",
    "delirium":            "Incident delirium",
    "pressure_ulcer":      "Incident pressure ulcer",
    "any_CVD_event":       "Any incident CVD event",
    "any_cancer_event":    "Any incident cancer",
    "incident_T2D":        "Incident T2D",
    "incident_HF":         "Incident HF",
    "incident_CKD":        "Incident CKD",
    "incident_dementia":   "Incident dementia",
}

PAL4 = ["#3b8bba", "#7fb069", "#f1a340", "#d7191c"]
GROUPS = ["Q1", "Q2", "Q3", "Q4"]


def tidy_ymax(y_observed: float) -> float:
    """Snap (y_observed * 1.15) up to a tidy step.

    Steps scale with magnitude so the axis carries 4-5 visible ticks
    regardless of endpoint rate.
    """
    target = max(y_observed * 1.15, 0.005)
    if target <= 0.015:   step = 0.0025
    elif target <= 0.025: step = 0.005
    elif target <= 0.05:  step = 0.01
    elif target <= 0.15:  step = 0.02
    elif target <= 0.40:  step = 0.05
    else:                 step = 0.10
    return math.ceil(target / step) * step


def derive_axes(df: pd.DataFrame) -> tuple[float, float]:
    """Per-plot xmax (99th-pct of times) and ymax from cuminc at that xmax."""
    xmax = float(np.nanpercentile(df["time"], 99))
    xmax = float(min(15.0, max(5.0, round(xmax))))
    per_group_max = []
    for g in GROUPS:
        sub = df[df["group"] == g]
        sub = sub[sub["time"] <= xmax]
        if len(sub):
            per_group_max.append(float(sub["cuminc"].max()))
    y_observed = max(per_group_max) if per_group_max else 0.005
    return xmax, tidy_ymax(y_observed)


def render_curves(df: pd.DataFrame, title: str, ymax: float, xmax: float,
                  is_adjusted: bool, out_stem: Path,
                  caption_suffix: str = "") -> None:
    """Render one figure (KM or adjusted-cuminc) with tight per-plot ylim."""
    fig, ax = plt.subplots(figsize=(7, 3.5))
    for i, g in enumerate(GROUPS):
        sub = df[df["group"] == g].sort_values("time")
        if len(sub) == 0:
            continue
        ax.plot(sub["time"], sub["cuminc"], color=PAL4[i % len(PAL4)],
                linewidth=1.5, label=g)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    brk = max(1, round(xmax / 5))
    ax.set_xticks(range(0, int(xmax) + 1, int(brk)))
    y_ticks = np.linspace(0, ymax, 5)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f"{v * 100:.1f}%" for v in y_ticks])
    ax.set_xlabel("Year", fontsize=10)
    ylabel = "Adjusted cumulative incidence (%)" if is_adjusted else "Cumulative incidence (%)"
    ax.set_ylabel(ylabel, fontsize=10)
    full_title = f"{title}{caption_suffix}"
    ax.set_title(full_title, fontsize=11)
    ax.legend(title="", fontsize=9)
    plt.tight_layout()
    fig.savefig(f"{out_stem}.pdf", dpi=300)
    fig.savefig(f"{out_stem}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    PNG_DIR.mkdir(parents=True, exist_ok=True)
    for label, title in TITLES.items():
        # Raw KM (model-independent)
        km_csv = SRC / f"section2_{label}_km.csv"
        if km_csv.exists():
            df_km = pd.read_csv(km_csv)
            xmax, ymax = derive_axes(df_km)
            render_curves(df_km, title, ymax, xmax, is_adjusted=False,
                          out_stem=PDF_DIR / f"section2_{label}_km")
            os.replace(str(PDF_DIR / f"section2_{label}_km.png"),
                       str(PNG_DIR / f"section2_{label}_km.png"))
            print(f"{label:24s}  km:                xmax={xmax:.0f}  ymax={ymax:.3f}")

        # Clinical-only adjcuminc
        clin_csv = SRC / f"section2_{label}_adjcuminc_clinical.csv"
        if clin_csv.exists():
            df_clin = pd.read_csv(clin_csv)
            xmax, ymax = derive_axes(df_clin)
            render_curves(df_clin, title, ymax, xmax, is_adjusted=True,
                          out_stem=PDF_DIR / f"section2_{label}_adjcuminc_clinical",
                          caption_suffix=" (Clinical-only)")
            os.replace(str(PDF_DIR / f"section2_{label}_adjcuminc_clinical.png"),
                       str(PNG_DIR / f"section2_{label}_adjcuminc_clinical.png"))
            print(f"{label:24s}  adjcuminc_clinical:  xmax={xmax:.0f}  ymax={ymax:.3f}")

        # Biomarkers (primary) adjcuminc
        bio_csv = SRC / f"section2_{label}_adjcuminc_biomarkers.csv"
        if bio_csv.exists():
            df_bio = pd.read_csv(bio_csv)
            xmax, ymax = derive_axes(df_bio)
            render_curves(df_bio, title, ymax, xmax, is_adjusted=True,
                          out_stem=PDF_DIR / f"section2_{label}_adjcuminc_biomarkers",
                          caption_suffix=" (Clinical + Biomarkers)")
            os.replace(str(PDF_DIR / f"section2_{label}_adjcuminc_biomarkers.png"),
                       str(PNG_DIR / f"section2_{label}_adjcuminc_biomarkers.png"))
            print(f"{label:24s}  adjcuminc_biomarkers: xmax={xmax:.0f}  ymax={ymax:.3f}")


if __name__ == "__main__":
    main()
