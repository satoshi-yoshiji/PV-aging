"""
Re-render Fig 1a as an alluvial / transition plot showing how the
cumulative Hospital Frailty Risk Score (HFRS) bin distribution evolves
across years 0, 4, 8 and 12 of follow-up.

For each participant and each time point t in {0, 4, 8, 12}, we compute
the cumulative HFRS as the sum of Gilbert 2018 weights over distinct
3-character ICD-10 codes recorded within (-2 years, t years] of baseline.
Each participant is then assigned to one of five HFRS bins:
  0  /  (0, 1]  /  (1, 3]  /  (3, 5]  /  >5

The alluvial draws a vertical bin-stack at each time point and fills the
flows between adjacent time points with cubic-Bezier curves, coloured by
the destination bin. The visual mirrors the R / ggalluvial reference at
https://longitudinalanalysis.com/visualizing-transitions-in-time-using-r-and-alluvial-graphs/
but uses our project's 5-tone HFRS palette and a Tufte-flavoured matplotlib
theme (no spines, minimal annotation).

Output: figures/section1_hfrs_alluvial.{pdf,png}
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
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

ROOT = Path("/humgen/diabetes2/users/satoshi/misc/02.aging")
sys.path.insert(0, str(ROOT))

from aging_shared import (
    load_full_df, prepare_aging_cohort, parse_icd_with_dates, HFRS_CODES,
    save_figure,
)

INPUT = "/humgen/diabetes2/users/satoshi/misc/01.pv/curated_stats.tsv.gz"
PDF_DIR = ROOT / "results" / "figures"
PNG_DIR = ROOT / "manuscript" / "figures"

# 5-tone palette consistent with the previous Fig 1a stacked bar
BIN_LABELS = ["0", "(0, 1]", "(1, 3]", "(3, 5]", ">5"]
BIN_COLORS = ["#dcdcdc", "#f4cb91", "#f1a340", "#d7191c", "#7a0103"]
BIN_EDGES = [-1e-6, 1e-6, 1.0, 3.0, 5.0, np.inf]

YEARS = [0, 4, 8, 12]


def hfrs_at_year(long_df: pd.DataFrame, year_t: float,
                 window_back_days: int = 730) -> pd.Series:
    """Cumulative HFRS per eid using distinct 3-char codes recorded in
    (-window_back_days, year_t * 365.25] days from baseline."""
    cutoff = year_t * 365.25
    df = long_df[(long_df["days_from_baseline"] >= -window_back_days)
                 & (long_df["days_from_baseline"] <= cutoff)].copy()
    df["w"] = df["code3"].map(HFRS_CODES).fillna(0.0)
    df = df[df["w"] > 0]
    if df.empty:
        return pd.Series(dtype=float)
    first = (df.sort_values("dxdate")
               .drop_duplicates(["eid", "code3"], keep="first"))
    return first.groupby("eid")["w"].sum()


def assign_bin(s: pd.Series) -> pd.Series:
    return pd.cut(s, bins=BIN_EDGES, labels=BIN_LABELS, include_lowest=True)


def alluvial(ax, columns_pct, transitions, x_positions,
             colors=BIN_COLORS, bin_labels=BIN_LABELS,
             year_labels=None, bar_width=0.10):
    """Draw the alluvial.

    columns_pct[i] : array of len(bin_labels) summing to 1.0 -- bin shares at x_positions[i]
    transitions[i]  : 2-D array, shape (n_bins, n_bins),
                      transitions[i][a, b] = fraction of total population that is
                      in bin a at x_positions[i] AND in bin b at x_positions[i+1].
    x_positions      : list of x coordinates (one per time point).
    """
    n_t = len(x_positions)
    n_bins = len(bin_labels)

    # Cumulative bin tops at each time point (for bar coordinates)
    # We stack from bin 0 at the bottom to bin n_bins-1 at the top.
    # bin_y_lo[i, b] = lower y of bar segment for bin b at time i
    bin_y_lo = np.zeros((n_t, n_bins))
    bin_y_hi = np.zeros((n_t, n_bins))
    for i in range(n_t):
        cum = 0.0
        for b in range(n_bins):
            bin_y_lo[i, b] = cum
            cum += columns_pct[i][b]
            bin_y_hi[i, b] = cum

    # Draw flow ribbons. For each adjacent pair (i, i+1):
    #   for each (a, b) cell of transitions[i],
    #     allocate a vertical sub-segment of bin a (at time i) AND of bin b (at i+1)
    #     proportional to transitions[i][a, b].
    for i in range(n_t - 1):
        x1 = x_positions[i] + bar_width / 2
        x2 = x_positions[i + 1] - bar_width / 2

        # Cursors per bin at t (proportions placed so far)
        used_at_t = np.zeros(n_bins)
        used_at_tp1 = np.zeros(n_bins)

        # Order: iterate destination-major (so destination bins get contiguous blocks)
        # We want the lower bins at the bottom of each side; iterate a in order, b in order.
        for a in range(n_bins):
            for b in range(n_bins):
                w = transitions[i][a, b]
                if w <= 0:
                    continue
                # Sub-range on the left (within bar a at time i)
                y1_lo = bin_y_lo[i, a] + used_at_t[a]
                y1_hi = y1_lo + w
                used_at_t[a] += w
                # Sub-range on the right (within bar b at time i+1)
                y2_lo = bin_y_lo[i + 1, b] + used_at_tp1[b]
                y2_hi = y2_lo + w
                used_at_tp1[b] += w
                # Ribbon coloured by SOURCE bin (origin), so a participant
                # moving from HFRS=0 (grey) to HFRS>5 (dark red) is shown
                # as a grey flow -- the colour communicates where the
                # participant came from.
                draw_flow(ax, x1, y1_lo, y1_hi, x2, y2_lo, y2_hi,
                          color=colors[a], alpha=0.55)

    # Draw the bin stacks (vertical bars at each time point)
    for i, x in enumerate(x_positions):
        for b in range(n_bins):
            h = bin_y_hi[i, b] - bin_y_lo[i, b]
            if h <= 0:
                continue
            rect = Rectangle((x - bar_width / 2, bin_y_lo[i, b]),
                             bar_width, h, facecolor=colors[b],
                             edgecolor="white", linewidth=0.5, zorder=3)
            ax.add_patch(rect)
            # Annotate the largest segment per bar with its percentage
            pct = h * 100
            if pct >= 4:
                ax.text(x, bin_y_lo[i, b] + h / 2, f"{pct:.0f}%",
                        ha="center", va="center", color="black" if b <= 1 else "white",
                        fontsize=8, zorder=4)

    # Axes
    if year_labels is None:
        year_labels = [f"Year {int(t)}" for t in x_positions]
    ax.set_xticks(x_positions)
    ax.set_xticklabels(year_labels, fontsize=11)
    ax.set_yticks([])
    ax.set_ylim(0, 1.001)
    ax.set_xlim(min(x_positions) - 1.0, max(x_positions) + 1.0)
    # Tufte-style: no spines, no grid
    for spine in ("top", "right", "left", "bottom"):
        ax.spines[spine].set_visible(False)


def draw_flow(ax, x1, y1_lo, y1_hi, x2, y2_lo, y2_hi, color, alpha,
              n_pts: int = 200, k: float = 8.0):
    """Draw a smooth alluvial ribbon between (x1, [y1_lo, y1_hi]) and
    (x2, [y2_lo, y2_hi]) using a high-resolution tanh-sigmoid S-curve
    on the top and bottom edges (closer to ggalluvial's geom_alluvium
    than a 4-point cubic-Bezier path).

    n_pts : number of sample points along each edge (200 -> very smooth).
    k     : sigmoid steepness in normalised x; 6-10 reads well at this scale.
    """
    t = np.linspace(0.0, 1.0, n_pts)
    s = 0.5 * (1 + np.tanh((t - 0.5) * k))
    xs = x1 + (x2 - x1) * t
    top_ys = y1_hi + (y2_hi - y1_hi) * s
    bot_ys = y1_lo + (y2_lo - y1_lo) * s
    # Build a closed polygon: top edge L->R, then bottom edge R->L
    polygon_x = np.concatenate([xs, xs[::-1]])
    polygon_y = np.concatenate([top_ys, bot_ys[::-1]])
    ax.fill(polygon_x, polygon_y, facecolor=color, edgecolor="none",
            alpha=alpha, antialiased=True, linewidth=0)


def main():
    print("Loading data...")
    full = load_full_df(input_path=INPUT)
    full, df_aging, _, _ = prepare_aging_cohort(full)
    long_df = parse_icd_with_dates(full)

    # Restrict to participants in the analytic cohort
    eids_keep = set(df_aging["eid"].astype(str))
    long_df = long_df[long_df["eid"].astype(str).isin(eids_keep)].copy()

    # Compute HFRS at each year
    print("Computing HFRS at year 0, 4, 8, 12 ...")
    hfrs_by_year = {}
    for t in YEARS:
        s = hfrs_at_year(long_df, t)
        # Reindex to ALL eids in the analytic cohort, fill missing with 0
        full_index = pd.Series(index=eids_keep, data=0.0, dtype=float, name="hfrs")
        full_index.update(s)
        hfrs_by_year[t] = full_index
        print(f"  year {t}: median HFRS = {full_index.median():.2f}; "
              f"% HFRS=0 = {(full_index == 0).mean()*100:.1f}%; "
              f"% HFRS>=5 = {(full_index >= 5).mean()*100:.2f}%")

    # Bin assignments per year
    bin_per_year = {t: assign_bin(hfrs_by_year[t]) for t in YEARS}

    # Total participants
    N = len(eids_keep)

    # Per-time bin shares (each row sums to 1)
    columns_pct = []
    for t in YEARS:
        counts = bin_per_year[t].value_counts(dropna=False).reindex(BIN_LABELS, fill_value=0)
        columns_pct.append(counts.values / N)

    # Adjacent-time transition matrices: trans[i][a, b] = P(bin@t_i = a AND bin@t_{i+1} = b)
    transitions = []
    for i in range(len(YEARS) - 1):
        a = bin_per_year[YEARS[i]]
        b = bin_per_year[YEARS[i + 1]]
        # Joint distribution: pd.crosstab returns a (n_bins, n_bins) df
        ct = pd.crosstab(a, b, dropna=False).reindex(index=BIN_LABELS,
                                                      columns=BIN_LABELS,
                                                      fill_value=0)
        transitions.append(ct.values / N)
        # Sanity check: transitions[i] sum to columns_pct[i]
        row_sums = transitions[-1].sum(axis=1)
        col_sums = transitions[-1].sum(axis=0)
        np.testing.assert_allclose(row_sums, columns_pct[i], atol=1e-6)
        np.testing.assert_allclose(col_sums, columns_pct[i + 1], atol=1e-6)
    print(f"All transition matrices consistency-checked (N = {N:,})")

    # ---- Render ----
    plt.rcParams.update({
        "font.family": "serif",  # Tufte-flavoured
        "font.size": 11,
    })
    fig, ax = plt.subplots(figsize=(8, 4))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    x_positions = [0, 4, 8, 12]
    alluvial(ax, columns_pct, transitions, x_positions,
             year_labels=["Year 0", "Year 4", "Year 8", "Year 12"])

    # Legend (HFRS bins)
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=BIN_COLORS[b], edgecolor="white",
                     label=BIN_LABELS[b]) for b in range(len(BIN_LABELS))]
    ax.legend(handles=handles, title="Cumulative HFRS",
              loc="upper left", bbox_to_anchor=(1.0, 1.0),
              frameon=False, fontsize=9, title_fontsize=9)

    ax.set_title(f"HFRS transitions across follow-up (n = {N:,})",
                 fontsize=12, pad=8)
    plt.tight_layout()

    pdf_path = PDF_DIR / "section1_hfrs_alluvial.pdf"
    png_path = PNG_DIR / "section1_hfrs_alluvial.png"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
