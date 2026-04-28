"""
Section 2 — Prospective Cox panel for ETP-aging.

13 endpoints x 3 adjustment models. Per-endpoint:
  - Per-SD Cox (HR for 1-SD col6a3)
  - Quartile Cox (Q2..Q4 vs Q1) plus p-trend (quartile-as-int)
  - KM cumulative-incidence curves by quartile (PDF + source data)
  - Adjusted cumulative-incidence curves by quartile (PDF + source data)
  - C-index on +Clinical only

Outputs:
  results/tables/section2_per_sd.csv
  results/tables/section2_quartiles.csv
  results/tables/section2_cohort_flow.csv
  results/tables/section2_hfrs_top_contributors.csv
  results/tables/section2_zenin_components.csv
  results/figures/section2_<endpoint>_km.pdf
  results/figures/section2_<endpoint>_adjcuminc.pdf
  results/figures/section2_master_forest.pdf
  results/source_data/section2_<endpoint>_*.csv
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Suppress lifelines convergence warnings (we look at p-values, not betas of dummies)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.insert(0, str(Path(__file__).parent))
from aging_shared import (
    LOG, PAL4, PAL5, XMAX_YEARS, YMAX_KM, YMAX_ADJ,
    load_full_df, prepare_aging_cohort,
    available_covars, build_cox_df, safe_coxph,
    extract_col6a3_per_sd, cox_hr_table,
    make_quartile,
    km_cumulative_incidence_plot, compute_adjusted_cuminc, plot_adjusted_cuminc,
    forest_plot, save_dataframe, save_figure,
    fmt_p, fmt_hr,
    ADJ_BASE, ADJ_CLIN, ADJ_BIO, ADJ_MODELS,
    HFRS_CODES, parse_icd_with_dates,
)

# (label, time_years_col, event_col, prevalent_col, plot_title, [caveat])
ENDPOINTS = [
    ("all_cause_mortality", "time_to_anydeath_years",        "incident_anydeath",         "prevalent_anydeath",        "All-cause mortality", None),
    ("incident_HFRS5",      "time_to_HFRS5_years",           "incident_HFRS5",            "prevalent_HFRS5",           "HFRS-positive admission (>=5)", None),
    ("zenin_composite",     "time_to_zenin_years",           "incident_zenin",            "prevalent_zenin",           "Zenin healthspan composite", None),
    ("hip_fracture",        "time_to_hip_fracture_years",    "incident_hip_fracture",     "prevalent_hip_fracture",    "Incident hip fracture", None),
    ("falls",               "time_to_falls_years",           "incident_falls",            "prevalent_falls",           "Incident falls", None),
    ("delirium",            "time_to_delirium_years",        "incident_delirium",         "prevalent_delirium",        "Incident delirium", None),
    ("pressure_ulcer",      "time_to_pressure_ulcer_years",  "incident_pressure_ulcer",   "prevalent_pressure_ulcer",  "Incident pressure ulcer", None),
    ("any_CVD_event",       "time_to_any_CVD_event_years",   "incident_any_CVD_event",    "prevalent_any_CVD_event",   "Any incident CVD event", "ICD-event substitute (no f.40001)"),
    ("any_cancer_event",    "time_to_any_cancer_event_years","incident_any_cancer_event", "prevalent_any_cancer_event","Any incident cancer", "ICD-event substitute (no f.40001)"),
    ("incident_T2D",        "time_to_T2D_years",             "incident_T2D_case_control", "prevalent_T2D",             "Incident T2D", None),
    ("incident_HF",         "time_to_HF_years",              "incident_HF",               "prevalent_HF",              "Incident HF", None),
    ("incident_CKD",        "time_to_CKD_years",             "incident_CKD",              "prevalent_CKD",             "Incident CKD", None),
    ("incident_dementia",   "time_to_dementia_years",        "incident_dementia",         "prevalent_dementia",        "Incident dementia", None),
]


def _tidy_ymax_from_curves(df, groups, xmax: float) -> float:
    """Compute a tight, tidy ymax for a (time, group, cuminc) DF.

    Takes the maximum cumulative-incidence value across groups at or before
    xmax, adds 15% headroom, and snaps up to a tidy step (0.0025 / 0.005 /
    0.01 / 0.02 / 0.05 / 0.1) so the axis carries 4-5 visible ticks.
    """
    import math
    import pandas as pd
    if not isinstance(df, pd.DataFrame) or len(df) == 0:
        return 0.05
    per_group_max = []
    for g in groups:
        sub = df[df["group"] == g]
        sub = sub[sub["time"] <= xmax]
        if len(sub):
            per_group_max.append(float(sub["cuminc"].max()))
    y_obs = max(per_group_max) if per_group_max else 0.005
    target = max(y_obs * 1.15, 0.005)
    if target <= 0.015:   step = 0.0025
    elif target <= 0.025: step = 0.005
    elif target <= 0.05:  step = 0.01
    elif target <= 0.15:  step = 0.02
    elif target <= 0.40:  step = 0.05
    else:                 step = 0.10
    return math.ceil(target / step) * step


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """BH-FDR adjusted q-values (no statsmodels dependency)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def run_one_endpoint(df_aging: pd.DataFrame, endpoint, model_name: str, covars: list[str],
                     fig_dir: Path, src_dir: Path) -> dict:
    label, t_col, ev_col, prev_col, title, caveat = endpoint
    LOG.info("== Endpoint: %s | model: %s ==", label, model_name)

    # 1. Drop prevalent cases of THIS endpoint
    df = df_aging.copy()
    if prev_col in df.columns:
        df = df[df[prev_col] != 1].copy()
    if t_col not in df.columns or ev_col not in df.columns:
        LOG.warning("Missing %s or %s — skipping", t_col, ev_col)
        return {"endpoint": label, "model": model_name, "skipped": True}
    df["time"] = pd.to_numeric(df[t_col], errors="coerce")
    df["event"] = pd.to_numeric(df[ev_col], errors="coerce").fillna(0).astype(int)
    df = df[(df["time"] > 0) & np.isfinite(df["time"])].copy()
    n_pre_cov = len(df)
    if len(df) == 0 or df["event"].sum() < 10:
        LOG.warning("Endpoint %s has too few rows/events; skipping", label)
        return {"endpoint": label, "model": model_name, "skipped": True}

    # 2. Per-SD Cox
    use_covars = available_covars(df, covars)
    cox_df = build_cox_df(df, "time", "event", "col6a3_scaled", use_covars).dropna()
    fit = safe_coxph(cox_df, "time", "event", penalizer=0.01, label=f"{label}_{model_name}_perSD")
    per_sd = extract_col6a3_per_sd(fit)
    if per_sd is None:
        return {"endpoint": label, "model": model_name, "skipped": True}
    n_used = len(cox_df)
    n_events_used = int(cox_df["event"].sum())

    # 3. Quartile Cox
    df["quartile"] = make_quartile(df["col6a3_scaled"])
    df_q = df.dropna(subset=["quartile"]).copy()
    df_q["q"] = df_q["quartile"]
    df_q_dummies = pd.get_dummies(df_q[["q"]], drop_first=True)
    df_q2 = pd.concat([df_q[["time", "event"] + use_covars], df_q_dummies], axis=1)
    cat_cols = [c for c in use_covars if not pd.api.types.is_numeric_dtype(df_q2[c])]
    num_cols = [c for c in use_covars if c not in cat_cols]
    parts = [df_q2[["time", "event"] + num_cols], df_q_dummies]
    if cat_cols:
        parts.append(pd.get_dummies(df_q2[cat_cols], drop_first=True))
    cox_q = pd.concat(parts, axis=1).dropna()
    fit_q = safe_coxph(cox_q, "time", "event", penalizer=0.01, label=f"{label}_{model_name}_quartile")
    quartile_rows = []
    if fit_q is not None:
        for q_lab in ["q_Q2", "q_Q3", "q_Q4"]:
            if q_lab in fit_q.summary.index:
                row = fit_q.summary.loc[q_lab]
                quartile_rows.append({
                    "term": q_lab.replace("q_", ""),
                    "HR": float(np.exp(row["coef"])),
                    "LCI": float(np.exp(row["coef lower 95%"])),
                    "UCI": float(np.exp(row["coef upper 95%"])),
                    "P": float(row["p"]),
                })

    # 4. p-trend (quartile-as-int)
    df_t = df_q.copy()
    qmap = {"Q1": 1, "Q2": 2, "Q3": 3, "Q4": 4}
    df_t["q_int"] = df_t["q"].map(qmap)
    cox_t = build_cox_df(df_t, "time", "event", "q_int", use_covars).dropna()
    fit_t = safe_coxph(cox_t, "time", "event", penalizer=0.01, label=f"{label}_{model_name}_trend")
    p_trend = float(fit_t.summary.loc["q_int", "p"]) if fit_t is not None and "q_int" in fit_t.summary.index else np.nan
    hr_trend = float(np.exp(fit_t.summary.loc["q_int", "coef"])) if fit_t is not None and "q_int" in fit_t.summary.index else np.nan

    # 5. C-index (only on +Clinical, the canonical stat)
    c_index = np.nan
    if model_name == "+Clinical" and fit is not None:
        try:
            c_index = float(fit.concordance_index_)
        except Exception:
            pass

    # 6a. Raw KM + +Clinical-adjusted cumulative incidence (run during +Clinical iteration)
    if model_name == "+Clinical":
        groups = ["Q1", "Q2", "Q3", "Q4"]
        try:
            xmax = int(min(XMAX_YEARS, np.nanpercentile(df_q["time"].values, 99)))
            xmax = max(5, xmax)
            # Data-driven ymax: pre-fit per-quartile KM to find the maximum
            # cumulative incidence at xmax across groups, then add 15% headroom
            # and snap up to a tidy step (avoids dead vertical space for rare
            # endpoints like hip fracture / delirium / pressure ulcer).
            from lifelines import KaplanMeierFitter
            import math
            cuminc_xmax = []
            for g in groups:
                sub = df_q[df_q["q"] == g]
                if len(sub) == 0:
                    continue
                k = KaplanMeierFitter().fit(sub["time"], sub["event"])
                sf = float(k.survival_function_at_times([xmax]).iloc[0])
                cuminc_xmax.append(1.0 - sf)
            y_obs = max(cuminc_xmax) if cuminc_xmax else 0.05
            target = max(y_obs * 1.15, 0.01)
            step = (0.005 if target <= 0.025 else
                    0.01 if target <= 0.05 else
                    0.02 if target <= 0.15 else
                    0.05 if target <= 0.40 else 0.10)
            ymax_km = math.ceil(target / step) * step
            fig, ax, kmf = km_cumulative_incidence_plot(
                df_q, "time", "event", "q", groups, PAL4,
                xmax=xmax, ymax=ymax_km, title=title,
            )
            save_figure(fig, str(fig_dir / f"section2_{label}_km"))
            plt.close(fig)
            # Source data
            t_grid = np.linspace(0, xmax, 200)
            rows = []
            for g, k in kmf.items():
                sf = k.survival_function_at_times(t_grid).values
                for tt, s_ in zip(t_grid, sf):
                    rows.append({"time": float(tt), "group": g, "cuminc": 1.0 - float(s_)})
            save_dataframe(pd.DataFrame(rows), str(src_dir / f"section2_{label}_km.csv"))
        except Exception as e:
            LOG.warning("KM failed for %s: %s", label, e)

        # Adjusted cuminc
        try:
            cox_q_dropped = cox_q.dropna()
            adj = compute_adjusted_cuminc(
                fit_q,
                df_design=df_q.loc[cox_q_dropped.index],
                cox_df_aligned=cox_q_dropped,
                groups=groups,
                group_var="q",
                covars_use=use_covars,
                time_var="time",
                event_var="event",
                xmax=xmax,
            )
            if adj is not None and len(adj) > 0:
                # Per-plot tight ylim derived from the ADJUSTED data, not raw KM
                ymax_adj = _tidy_ymax_from_curves(adj, groups, xmax)
                fig2 = plot_adjusted_cuminc(adj, groups, PAL4, xmax=xmax,
                                            ymax=ymax_adj,
                                            title=f"{title} (Clinical-only adjusted)")
                fig2.axes[0].set_ylabel(f"Adjusted cumulative incidence (%)", fontsize=10)
                save_figure(fig2, str(fig_dir / f"section2_{label}_adjcuminc_clinical"))
                plt.close(fig2)
                save_dataframe(adj, str(src_dir / f"section2_{label}_adjcuminc_clinical.csv"))
        except Exception as e:
            LOG.warning("Adjusted cuminc failed for %s: %s", label, e)

    # 6b. Clinical+Biomarkers-adjusted cumulative incidence (primary analysis)
    if model_name == "+Biomarkers":
        groups = ["Q1", "Q2", "Q3", "Q4"]
        try:
            xmax = int(min(XMAX_YEARS, np.nanpercentile(df_q["time"].values, 99)))
            xmax = max(5, xmax)
            cox_q_dropped = cox_q.dropna()
            adj = compute_adjusted_cuminc(
                fit_q,
                df_design=df_q.loc[cox_q_dropped.index],
                cox_df_aligned=cox_q_dropped,
                groups=groups,
                group_var="q",
                covars_use=use_covars,
                time_var="time",
                event_var="event",
                xmax=xmax,
            )
            if adj is not None and len(adj) > 0:
                # Per-plot tight ylim derived from the ADJUSTED data
                ymax_adj = _tidy_ymax_from_curves(adj, groups, xmax)
                fig2 = plot_adjusted_cuminc(adj, groups, PAL4, xmax=xmax,
                                            ymax=ymax_adj,
                                            title=f"{title} (Clinical + Biomarkers)")
                fig2.axes[0].set_ylabel(f"Adjusted cumulative incidence (%)", fontsize=10)
                save_figure(fig2, str(fig_dir / f"section2_{label}_adjcuminc_biomarkers"))
                plt.close(fig2)
                save_dataframe(adj, str(src_dir / f"section2_{label}_adjcuminc_biomarkers.csv"))
        except Exception as e:
            LOG.warning("+Biomarkers adjcuminc failed for %s: %s", label, e)

    return {
        "endpoint": label,
        "model": model_name,
        "title": title,
        "caveat": caveat or "",
        "n_total": n_pre_cov,
        "n_used": n_used,
        "n_events": n_events_used,
        "HR_per_SD": float(per_sd["HR"].iloc[0]),
        "LCI_per_SD": float(per_sd["LCI"].iloc[0]),
        "UCI_per_SD": float(per_sd["UCI"].iloc[0]),
        "P_per_SD": float(per_sd["P"].iloc[0]),
        "HR_trend_perQ": hr_trend,
        "P_trend": p_trend,
        "c_index": c_index,
        "_quartiles": quartile_rows,
    }


def section_2_cohort_flow(df_aging: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append({"step": "df_aging (analytic cohort)", "endpoint": "ALL", "N": len(df_aging), "events": np.nan})
    for ep in ENDPOINTS:
        label, t_col, ev_col, prev_col, _, _ = ep
        d = df_aging.copy()
        n_after_prev = int((d.get(prev_col, pd.Series(0, index=d.index)) != 1).sum())
        d = d[d.get(prev_col, 0) != 1]
        d["t"] = pd.to_numeric(d[t_col], errors="coerce")
        d["e"] = pd.to_numeric(d[ev_col], errors="coerce").fillna(0).astype(int)
        d = d[(d["t"] > 0) & np.isfinite(d["t"])]
        rows.append({"step": "drop prevalent + valid time",
                     "endpoint": label, "N": len(d), "events": int(d["e"].sum())})
    return pd.DataFrame(rows)


def hfrs_top_contributors(long_df: pd.DataFrame, top_n: int = 30) -> pd.DataFrame:
    df = long_df.copy()
    df["w"] = df["code3"].map(HFRS_CODES).fillna(0.0)
    df = df[df["w"] > 0]
    df["dx_first"] = df.groupby(["eid", "code3"])["dxdate"].transform("min")
    first = df[df["dxdate"] == df["dx_first"]].drop_duplicates(["eid", "code3"])
    by_code = (first.groupby("code3", as_index=False)
                    .agg(n_eid=("eid", "nunique"), total_weight=("w", "sum"),
                         per_code_weight=("w", "first")))
    by_code = by_code.sort_values("total_weight", ascending=False).head(top_n).reset_index(drop=True)
    return by_code


def zenin_component_breakdown(full_df: pd.DataFrame) -> pd.DataFrame:
    """Share of zenin events contributed by each first-component."""
    if "zenin_first_component" not in full_df.columns:
        return pd.DataFrame()
    inc = full_df[full_df.get("incident_zenin", 0) == 1]
    counts = inc["zenin_first_component"].value_counts(dropna=False).rename_axis("component").reset_index(name="n")
    counts["pct"] = 100 * counts["n"] / counts["n"].sum()
    return counts


def main(args):
    out = Path(args.outdir)
    fig_dir = out / "figures"; fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir = out / "tables";  tab_dir.mkdir(parents=True, exist_ok=True)
    src_dir = out / "source_data"; src_dir.mkdir(parents=True, exist_ok=True)

    LOG.info("Loading %s", args.input)
    full = load_full_df(input_path=args.input)
    full, df_aging, c6m, c6s = prepare_aging_cohort(full)
    LOG.info("Analytic cohort N=%d (col6a3 mean=%.3f sd=%.3f)", len(df_aging), c6m, c6s)

    # Diagnostics
    save_dataframe(section_2_cohort_flow(df_aging), str(tab_dir / "section2_cohort_flow.csv"))
    long_df = parse_icd_with_dates(full)
    save_dataframe(hfrs_top_contributors(long_df), str(tab_dir / "section2_hfrs_top_contributors.csv"))
    save_dataframe(zenin_component_breakdown(full), str(tab_dir / "section2_zenin_components.csv"))

    # Endpoint loop
    rows_per_sd = []
    rows_quartile = []
    for ep in ENDPOINTS:
        for model_name, covars in ADJ_MODELS:
            res = run_one_endpoint(df_aging, ep, model_name, covars, fig_dir, src_dir)
            if res.get("skipped"):
                continue
            qrows = res.pop("_quartiles", [])
            rows_per_sd.append(res)
            for qr in qrows:
                rows_quartile.append({"endpoint": res["endpoint"], "model": res["model"],
                                      "title": res["title"], "n_used": res["n_used"],
                                      **qr})

    df_sd = pd.DataFrame(rows_per_sd)
    df_q = pd.DataFrame(rows_quartile)

    # BH-FDR within each adjustment model on per-SD p
    if not df_sd.empty:
        df_sd["P_per_SD_q_BH"] = np.nan
        for model in df_sd["model"].unique():
            mask = df_sd["model"] == model
            df_sd.loc[mask, "P_per_SD_q_BH"] = _bh_fdr(df_sd.loc[mask, "P_per_SD"].values)
            mask2 = df_sd["model"] == model
            df_sd.loc[mask2, "P_trend_q_BH"] = _bh_fdr(df_sd.loc[mask2, "P_trend"].fillna(1.0).values)

    save_dataframe(df_sd, str(tab_dir / "section2_per_sd.csv"))
    save_dataframe(df_q, str(tab_dir / "section2_quartiles.csv"))

    # Master forest: per-SD HR across endpoints, faceted by adjustment model
    try:
        if not df_sd.empty:
            # Sort by Clinical+Biomarkers HR for the y-axis order
            bio = df_sd[df_sd["model"] == "+Biomarkers"].copy().sort_values("HR_per_SD")
            order = list(bio["title"])
            xmin = max(0.5, float(np.nanmin(df_sd["LCI_per_SD"])) * 0.95)
            xmax = max(float(np.nanmax(df_sd["UCI_per_SD"])) * 1.05, 2.5)
            panel_labels = [("Base", "Base"),
                            ("+Clinical", "Clinical"),
                            ("+Biomarkers", "Clinical + Biomarkers")]
            fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.6), sharey=True)
            for ax, (model_id, model_label) in zip(axes, panel_labels):
                sub = df_sd[df_sd["model"] == model_id].copy()
                if sub.empty:
                    continue
                sub["title"] = pd.Categorical(sub["title"], categories=order, ordered=True)
                sub = sub.sort_values("title").reset_index(drop=True)
                ys = np.arange(len(sub))
                # Thin CI line, no caps
                for y, lo, hi in zip(ys, sub["LCI_per_SD"].values, sub["UCI_per_SD"].values):
                    ax.plot([lo, hi], [y, y], color="#444444", linewidth=0.9,
                            solid_capstyle="butt", zorder=2)
                # Small filled dot
                ax.scatter(sub["HR_per_SD"].values, ys, s=22, color="#222222",
                           zorder=3, edgecolors="none")
                ax.axvline(1.0, color="#bbbbbb", linestyle="--", linewidth=0.8, zorder=1)
                ax.set_yticks(ys)
                ax.set_yticklabels(sub["title"], fontsize=9)
                ax.set_xscale("log")
                ax.set_xlim(xmin, xmax)
                ax.set_xlabel("HR per SD endotrophin", fontsize=10)
                ax.set_title(model_label, fontsize=11, pad=6)
                for spine in ("top", "right"):
                    ax.spines[spine].set_visible(False)
                ax.tick_params(axis="x", labelsize=9)
            for ax in axes[1:]:
                ax.tick_params(axis="y", labelleft=False)
            fig.suptitle("Per-SD endotrophin associations across 13 prospective endpoints",
                         fontsize=12, y=0.995)
            plt.tight_layout()
            save_figure(fig, str(fig_dir / "section2_master_forest"))
            plt.close(fig)
    except Exception as e:
        LOG.warning("Master forest failed: %s", e)

    LOG.info("Section 2 done. Tables in %s, figures in %s", tab_dir, fig_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="/humgen/diabetes2/users/satoshi/misc/01.pv/curated_stats.tsv.gz")
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    np.random.seed(args.seed)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    main(args)
