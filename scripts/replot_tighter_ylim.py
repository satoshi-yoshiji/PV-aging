"""
Re-render Section 2 KM and adjusted-cumulative-incidence figures with a
data-driven y-axis (max cumulative incidence at xmax + 15% headroom,
snapped up to a tidy step) instead of the previous fixed-floor heuristic.

Reads the already-saved source-data CSVs (results/source_data/*) and writes
both the PDF (results/figures/*) and the PNG copy in
manuscript/figures/*.
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

# Endpoint label -> nice plot title (matches run_section2_aging_prospective.ENDPOINTS)
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

# Standard 4-quartile palette consistent with PAL4 in ckd_shared
PAL4 = ["#3b8bba", "#7fb069", "#f1a340", "#d7191c"]

GROUPS = ["Q1", "Q2", "Q3", "Q4"]


def tidy_ymax(y_observed: float) -> float:
    """Snap y_observed (with 15% headroom) up to a tidy step.

    Uses 0.005 / 0.01 / 0.02 / 0.05 / 0.1 depending on magnitude so the
    axis carries 4-5 visible ticks regardless of endpoint rate.
    """
    target = max(y_observed * 1.15, 0.01)
    if target <= 0.025:   step = 0.005
    elif target <= 0.05:  step = 0.01
    elif target <= 0.15:  step = 0.02
    elif target <= 0.40:  step = 0.05
    else:                 step = 0.10
    return math.ceil(target / step) * step


def render_curves(df: pd.DataFrame, title: str, ymax: float, xmax: float,
                  is_adjusted: bool, out_stem: Path) -> None:
    """Render KM-style cumulative-incidence curves from a (time, group, cuminc) DF."""
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
    suffix = " (adjusted)" if is_adjusted else ""
    ax.set_title(f"{title}{suffix}", fontsize=11)
    ax.legend(title="", fontsize=9)
    plt.tight_layout()
    fig.savefig(f"{out_stem}.pdf", dpi=300)
    fig.savefig(f"{out_stem}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    PNG_DIR.mkdir(parents=True, exist_ok=True)
    # Pair each endpoint with both KM and adjcuminc; compute a SINGLE y-max per
    # endpoint shared between the two so they read consistently side-by-side.
    for label, title in TITLES.items():
        km_csv = SRC / f"section2_{label}_km.csv"
        adj_csv = SRC / f"section2_{label}_adjcuminc.csv"
        if not km_csv.exists():
            print(f"SKIP (no km source data): {label}")
            continue
        df_km = pd.read_csv(km_csv)
        df_adj = pd.read_csv(adj_csv) if adj_csv.exists() else None

        # x-axis cap: 99th percentile time across all rows
        xmax = float(np.nanpercentile(df_km["time"], 99))
        xmax = min(15.0, max(5.0, round(xmax)))

        # y-axis cap: max cumulative incidence at xmax across groups, +15%, snapped
        per_group_max = []
        for g in GROUPS:
            sub_km = df_km[df_km["group"] == g]
            sub_km = sub_km[sub_km["time"] <= xmax]
            if len(sub_km):
                per_group_max.append(float(sub_km["cuminc"].max()))
            if df_adj is not None:
                sub_adj = df_adj[df_adj["group"] == g]
                sub_adj = sub_adj[sub_adj["time"] <= xmax]
                if len(sub_adj):
                    per_group_max.append(float(sub_adj["cuminc"].max()))
        y_observed = max(per_group_max) if per_group_max else 0.05
        ymax = tidy_ymax(y_observed)
        print(f"{label:24s}  xmax={xmax:.0f}  observed_max={y_observed:.4f}  ymax={ymax:.3f}")

        render_curves(df_km, title, ymax, xmax, is_adjusted=False,
                      out_stem=PDF_DIR / f"section2_{label}_km")
        # Also write to manuscript PNG dir
        # (render_curves already wrote .png alongside .pdf next to the stem)
        # Move .png copy to manuscript
        os.replace(str(PDF_DIR / f"section2_{label}_km.png"),
                   str(PNG_DIR / f"section2_{label}_km.png"))

        if df_adj is not None and len(df_adj):
            render_curves(df_adj, title, ymax, xmax, is_adjusted=True,
                          out_stem=PDF_DIR / f"section2_{label}_adjcuminc")
            os.replace(str(PDF_DIR / f"section2_{label}_adjcuminc.png"),
                       str(PNG_DIR / f"section2_{label}_adjcuminc.png"))


if __name__ == "__main__":
    main()
