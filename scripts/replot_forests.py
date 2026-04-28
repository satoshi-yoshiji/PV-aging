"""Re-render the master forest (Fig 2) and the two supplementary forests
(S4, S5) in the thin-CI / no-caps / small-dot style of Nat. Cardiovasc. Res.
2024 Fig 3 (Yoshiji et al., 'Integrative proteogenomic analysis ...').

Source data: results/tables/section2_per_sd.csv,
             results/tables/section1_baseline_hfrs_OLS.csv,
             results/tables/section1_prevalent_syndromes_OR.csv.

Outputs: results/figures/section2_master_forest.{pdf,png}
         results/figures/section1_baseline_hfrs_forest.{pdf,png}
         results/figures/section1_prevalent_syndromes_forest.{pdf,png}
         + manuscript/figures/ PNG mirrors of each.
"""
from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import numpy as np
import pandas as pd


def _set_clean_ticks(ax, x_log: bool, xmin: float, xmax: float, ref: float):
    """Apply plain-number tick labels (no '10^0' notation) and a sensible
    set of major-tick positions for either log or linear axes."""
    if x_log:
        # Pick a tidy set of HR ticks within [xmin, xmax]
        candidates = [0.5, 0.7, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
        ticks = [t for t in candidates if xmin * 0.999 <= t <= xmax * 1.001]
        if ref is not None and ref not in ticks:
            ticks = sorted(ticks + [ref])
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_minor_locator(NullLocator())
    else:
        # Linear axis: 5 ticks spanning the range
        ticks = np.linspace(xmin, xmax, 5)
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _pos: (f"{v:.1f}".rstrip("0").rstrip(".") if abs(v) < 10 else f"{v:.0f}")
    ))

ROOT = Path("/humgen/diabetes2/users/satoshi/misc/02.aging")
PDF_DIR = ROOT / "results" / "figures"
PNG_DIR = ROOT / "manuscript" / "figures"

# Visual style: thin CI line, no caps, small filled dot, dashed reference line
DOT_KW = dict(s=22, color="#222222", zorder=3, edgecolors="none")
CI_KW = dict(color="#444444", linewidth=0.9, solid_capstyle="butt", zorder=2)
REF_KW = dict(color="#bbbbbb", linestyle="--", linewidth=0.8, zorder=1)


def forest_panel(ax, df, ratio_col="HR", lci_col="LCI", uci_col="UCI",
                 label_col="label", x_log=True, ref=1.0,
                 xlabel="HR per SD", title=None, xlim=None):
    """Draw a single thin-CI / no-caps / small-dot forest panel.

    df rows are drawn from bottom to top in the order given.
    """
    ys = np.arange(len(df))
    # CI line
    for y, lo, hi in zip(ys, df[lci_col].values, df[uci_col].values):
        ax.plot([lo, hi], [y, y], **CI_KW)
    # Dot
    ax.scatter(df[ratio_col].values, ys, **DOT_KW)
    # Reference line
    if ref is not None:
        ax.axvline(ref, **REF_KW)
    if x_log:
        ax.set_xscale("log")
    if xlim is not None:
        ax.set_xlim(*xlim)
    # Plain-number x-tick labels (no '10^0')
    xlo, xhi = ax.get_xlim()
    _set_clean_ticks(ax, x_log=x_log, xmin=xlo, xmax=xhi, ref=ref)
    ax.set_yticks(ys)
    ax.set_yticklabels(df[label_col].values, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=10)
    if title:
        ax.set_title(title, fontsize=11, pad=6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", labelsize=9)


# ---------- 1. Master forest (Fig 2): 3 panels x 13 endpoints ----------
def render_master_forest():
    df = pd.read_csv(ROOT / "results/tables/section2_per_sd.csv")
    # Use the +Biomarkers HRs to fix the y-axis order (sorted by HR)
    bio = df[df["model"] == "+Biomarkers"].copy().sort_values("HR_per_SD")
    order = list(bio["title"])

    # Panel labels (cosmetic)
    PANEL_LABELS = [("Base", "Base"),
                    ("+Clinical", "Clinical"),
                    ("+Biomarkers", "Clinical + Biomarkers")]

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.6), sharey=True)
    # Determine common x-limits across all panels for consistency
    all_lci = df["LCI_per_SD"].values
    all_uci = df["UCI_per_SD"].values
    xmin = max(0.5, np.nanmin(all_lci) * 0.95)
    xmax = max(np.nanmax(all_uci) * 1.05, 2.5)
    xlim = (xmin, xmax)

    for ax, (model_id, model_label) in zip(axes, PANEL_LABELS):
        sub = df[df["model"] == model_id].copy()
        # Sort to match the +Biomarkers order
        sub["title"] = pd.Categorical(sub["title"], categories=order, ordered=True)
        sub = sub.sort_values("title").reset_index(drop=True)
        sub = sub.rename(columns={
            "HR_per_SD": "HR", "LCI_per_SD": "LCI", "UCI_per_SD": "UCI",
            "title": "label",
        })
        forest_panel(ax, sub, ratio_col="HR", lci_col="LCI", uci_col="UCI",
                     label_col="label", x_log=True, ref=1.0,
                     xlabel="HR per SD endotrophin",
                     title=model_label, xlim=xlim)

    # Hide y-tick labels on the 2nd and 3rd panels (sharey)
    for ax in axes[1:]:
        ax.tick_params(axis="y", labelleft=False)

    fig.suptitle("Per-SD endotrophin associations across 13 prospective endpoints",
                 fontsize=12, y=0.995)
    plt.tight_layout()
    pdf = PDF_DIR / "section2_master_forest.pdf"
    png = PNG_DIR / "section2_master_forest.png"
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    # Also keep the Figure_2_master_forest.png mirror used by the manuscript
    import shutil
    shutil.copy(png, PNG_DIR / "Figure_2_master_forest.png")
    plt.close(fig)
    print(f"Saved master forest: {pdf}")


# ---------- 2. Suppl Fig S4: baseline HFRS beta-per-SD forest ----------
def render_baseline_hfrs_forest():
    df = pd.read_csv(ROOT / "results/tables/section1_baseline_hfrs_OLS.csv")
    sub = df[df["outcome"] == "HFRS (linear)"].copy()
    # Order: Base, +Clinical, +Biomarkers (top to bottom is conventional)
    order = ["Base", "+Clinical", "+Biomarkers"]
    label_map = {"Base": "Base", "+Clinical": "Clinical",
                 "+Biomarkers": "Clinical + Biomarkers"}
    sub["model_idx"] = sub["model"].map({m: i for i, m in enumerate(order)})
    sub = sub.dropna(subset=["model_idx"]).sort_values("model_idx", ascending=False)
    sub["label"] = sub["model"].map(label_map)
    sub = sub.rename(columns={"beta": "EST", "LCI": "LCI", "UCI": "UCI"})

    fig, ax = plt.subplots(figsize=(6.5, 2.8))
    forest_panel(ax, sub, ratio_col="EST", lci_col="LCI", uci_col="UCI",
                 label_col="label", x_log=False, ref=0.0,
                 xlabel="beta per SD endotrophin (HFRS units)",
                 title="Baseline HFRS (2-y window) ~ ETP")
    plt.tight_layout()
    pdf = PDF_DIR / "section1_baseline_hfrs_forest.pdf"
    png = PNG_DIR / "section1_baseline_hfrs_forest.png"
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    import shutil
    shutil.copy(png, PNG_DIR / "Figure_S4_baseline_HFRS_beta_forest.png")
    plt.close(fig)
    print(f"Saved S4: {pdf}")


# ---------- 3. Suppl Fig S5: prevalent-syndrome OR forest ----------
def render_prevalent_syndromes_forest():
    df = pd.read_csv(ROOT / "results/tables/section1_prevalent_syndromes_OR.csv")
    sub = df[(df["model"] == "+Biomarkers") & df["OR"].notna()].copy()
    # Order: hip fracture, falls, any frailty (skip too-sparse delirium / pressure ulcer)
    desired_order = ["Prevalent hip fracture", "Prevalent falls",
                     "Any prevalent frailty syndrome"]
    sub["__o"] = sub["outcome"].map({n: i for i, n in enumerate(desired_order)})
    sub = sub.dropna(subset=["__o"]).sort_values("__o", ascending=False)
    sub = sub.rename(columns={"outcome": "label"})

    fig, ax = plt.subplots(figsize=(6.5, 2.6))
    forest_panel(ax, sub, ratio_col="OR", lci_col="LCI", uci_col="UCI",
                 label_col="label", x_log=True, ref=1.0,
                 xlabel="OR per SD endotrophin (Clinical + Biomarkers)",
                 title="Prevalent frailty syndromes ~ ETP")
    plt.tight_layout()
    pdf = PDF_DIR / "section1_prevalent_syndromes_forest.pdf"
    png = PNG_DIR / "section1_prevalent_syndromes_forest.png"
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    import shutil
    shutil.copy(png, PNG_DIR / "Figure_S5_prevalent_syndromes_OR_forest.png")
    plt.close(fig)
    print(f"Saved S5: {pdf}")


if __name__ == "__main__":
    render_master_forest()
    render_baseline_hfrs_forest()
    render_prevalent_syndromes_forest()
