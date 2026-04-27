"""
Re-render Fig 1a (baseline HFRS distribution) as a stratified-by-ETP-quartile
view that reveals the dose-response. The previous histogram was dominated by
the zero-spike (most participants have no relevant ICD codes in the 2-year
baseline window) and made the rest of the distribution invisible.

This rewrite produces a 1x2 composite panel:
  Left  - 100% stacked bar of HFRS bin (0 / 0-1 / 1-3 / 3-5 / 5+) by ETP quartile.
  Right - bar of % of participants with HFRS >= 5 (Gilbert intermediate-risk
          threshold) by ETP quartile, with exact percentages annotated.

Both views directly illustrate the shift toward higher HFRS at higher
endotrophin and avoid the zero-spike legibility problem.

Output: section1_hfrs_baseline_hist.{pdf,png} (replaces the old single-panel
histogram so that the Fig 1 composite picks up the new layout automatically).
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

ROOT = Path("/humgen/diabetes2/users/satoshi/misc/02.aging")
sys.path.insert(0, str(ROOT))

from aging_shared import (
    load_full_df, prepare_aging_cohort, parse_icd_with_dates,
    cumulative_hfrs_baseline, make_quartile, save_figure,
)

INPUT = "/humgen/diabetes2/users/satoshi/misc/01.pv/curated_stats.tsv.gz"
PDF_DIR = ROOT / "results" / "figures"
PNG_DIR = ROOT / "manuscript" / "figures"

PAL4 = ["#3b8bba", "#7fb069", "#f1a340", "#d7191c"]
GROUPS = ["Q1", "Q2", "Q3", "Q4"]


def main():
    print("Loading data...")
    full = load_full_df(input_path=INPUT)
    full, df_aging, c6m, c6s = prepare_aging_cohort(full)

    long_df = parse_icd_with_dates(full)
    hfrs = cumulative_hfrs_baseline(long_df, window_days_back=730)

    df = df_aging.copy()
    df["eid"] = df["eid"].astype(str)
    df = df.merge(hfrs, on="eid", how="left")
    df["hfrs_baseline"] = df["hfrs_baseline"].fillna(0.0)
    df["quartile"] = make_quartile(df["col6a3_scaled"])

    # HFRS bin assignment
    bins = [-0.001, 0.001, 1.0, 3.0, 5.0, np.inf]
    bin_labels = ["0", "(0, 1]", "(1, 3]", "(3, 5]", ">5"]
    df["hfrs_bin"] = pd.cut(df["hfrs_baseline"], bins=bins, labels=bin_labels,
                            include_lowest=True)

    # 100%-stacked bar of HFRS bin by quartile
    counts = df.groupby(["quartile", "hfrs_bin"], observed=True).size().unstack(
        fill_value=0
    ).reindex(index=GROUPS, columns=bin_labels, fill_value=0)
    pcts = counts.div(counts.sum(axis=1), axis=0) * 100

    # % with HFRS >= 5 by quartile
    pct_high = df.groupby("quartile", observed=True).apply(
        lambda g: (g["hfrs_baseline"] >= 5).mean() * 100
    ).reindex(GROUPS)

    # Build the 1x2 composite
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.8),
                                    gridspec_kw={"width_ratios": [1.4, 1]})

    # Left: stacked bar (each bin gets a shade of grey -> red)
    bin_colors = ["#dcdcdc", "#f4cb91", "#f1a340", "#d7191c", "#7a0103"]
    bottom = np.zeros(len(GROUPS))
    x = np.arange(len(GROUPS))
    for j, (lab, col) in enumerate(zip(bin_labels, bin_colors)):
        vals = pcts[lab].values
        ax1.bar(x, vals, bottom=bottom, color=col, label=lab,
                edgecolor="white", linewidth=0.5)
        # annotate slices >= 4%
        for xi, v, b in zip(x, vals, bottom):
            if v >= 4:
                ax1.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                         color="black" if j <= 1 else "white", fontsize=8)
        bottom = bottom + vals
    ax1.set_xticks(x)
    ax1.set_xticklabels(GROUPS)
    ax1.set_xlabel("ETP quartile", fontsize=10)
    ax1.set_ylabel("Participants (%)", fontsize=10)
    ax1.set_ylim(0, 100)
    ax1.set_title("Baseline HFRS distribution by ETP quartile", fontsize=10)
    ax1.legend(title="Baseline HFRS", title_fontsize=8, fontsize=8,
               loc="lower right", framealpha=0.9, ncol=1)
    for spine in ("top", "right"):
        ax1.spines[spine].set_visible(False)

    # Right: % >= 5 by quartile
    bars = ax2.bar(x, pct_high.values, color=PAL4, edgecolor="white",
                   linewidth=0.5)
    for xi, v in zip(x, pct_high.values):
        ax2.text(xi, v + 0.05, f"{v:.1f}%", ha="center", va="bottom",
                 fontsize=9, color="#333")
    ax2.set_xticks(x)
    ax2.set_xticklabels(GROUPS)
    ax2.set_xlabel("ETP quartile", fontsize=10)
    ax2.set_ylabel("% with baseline HFRS >= 5", fontsize=10)
    top = max(pct_high.values) * 1.25
    ax2.set_ylim(0, top)
    ax2.set_title("High-frailty participants by ETP quartile", fontsize=10)
    for spine in ("top", "right"):
        ax2.spines[spine].set_visible(False)

    plt.tight_layout()
    save_figure(fig, str(PDF_DIR / "section1_hfrs_baseline_hist"))
    fig.savefig(PNG_DIR / "section1_hfrs_baseline_hist.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("Saved section1_hfrs_baseline_hist.{pdf,png}")
    print()
    print("Stacked-bar percentages (rows = ETP quartile):")
    print(pcts.round(1).to_string())
    print()
    print("% HFRS >= 5 by quartile:")
    print(pct_high.round(2).to_string())


if __name__ == "__main__":
    main()
