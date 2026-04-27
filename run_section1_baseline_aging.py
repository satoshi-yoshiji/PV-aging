"""
Section 1 — Cross-sectional baseline aging proxies.

The originally requested phenotypes (Williams 49-item FI, Fried, grip
strength, gait speed, hand-grip asymmetry) are not in curated_stats.tsv.gz.
This runner ships four feasible cross-sectional outcomes derived from
what IS in the file, plus a curation hook for the originally requested
analyses to fire automatically when the user later supplies a UKB extract.

Outputs (results/):
  tables/section1_baseline_hfrs_OLS.csv
  tables/section1_prevalent_syndromes_OR.csv
  tables/section1_multimorbidity_count.csv
  tables/section1_bioage_proxy.csv
  figures/section1_hfrs_baseline_hist.pdf
  figures/section1_baseline_hfrs_forest.pdf
  figures/section1_prevalent_syndromes_forest.pdf
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

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.insert(0, str(Path(__file__).parent))
from aging_shared import (
    LOG, PAL4,
    load_full_df, prepare_aging_cohort, parse_icd_with_dates,
    cumulative_hfrs_baseline, extend_with_baseline_phenotypes,
    selfrep_count, prevalent_disease_count,
    available_covars, save_dataframe, save_figure, fmt_p,
    forest_plot, forest_plot_OR,
    ADJ_BASE, ADJ_CLIN, ADJ_BIO, ADJ_MODELS,
    HFRS_CODES,
)


def _ols(df: pd.DataFrame, y: str, x: str, covars: list[str]) -> dict:
    """OLS of y ~ x + covars, returning beta(x), 95% CI, p, N."""
    import statsmodels.api as sm
    use = [c for c in covars if c in df.columns]
    Xcols = [x] + use
    sub = df[[y] + Xcols].dropna()
    cat = [c for c in use if not pd.api.types.is_numeric_dtype(sub[c])]
    num = [c for c in use if c not in cat]
    Xn = sub[[x] + num].astype(float)
    if cat:
        Xc = pd.get_dummies(sub[cat], drop_first=True).astype(float)
        X = pd.concat([Xn, Xc], axis=1)
    else:
        X = Xn
    X = sm.add_constant(X)
    try:
        m = sm.OLS(sub[y].astype(float), X).fit()
    except Exception as e:
        LOG.warning("OLS failed (%s ~ %s): %s", y, x, e)
        return {"beta": np.nan, "se": np.nan, "LCI": np.nan, "UCI": np.nan, "P": np.nan, "N": len(sub)}
    coef = m.params.get(x, np.nan)
    se = m.bse.get(x, np.nan)
    p = m.pvalues.get(x, np.nan)
    return {
        "beta": float(coef), "se": float(se),
        "LCI": float(coef - 1.96 * se), "UCI": float(coef + 1.96 * se),
        "P": float(p), "N": int(len(sub)),
    }


def _logit(df: pd.DataFrame, y: str, x: str, covars: list[str]) -> dict:
    """Logistic regression of y (0/1) ~ x + covars; returns OR for x."""
    import statsmodels.api as sm
    use = [c for c in covars if c in df.columns]
    Xcols = [x] + use
    sub = df[[y] + Xcols].dropna()
    sub = sub[sub[y].isin([0, 1])]
    if sub.empty or sub[y].sum() < 10 or (sub[y] == 0).sum() < 10:
        return {"OR": np.nan, "LCI": np.nan, "UCI": np.nan, "P": np.nan, "N": len(sub), "n_event": int(sub[y].sum())}
    cat = [c for c in use if not pd.api.types.is_numeric_dtype(sub[c])]
    num = [c for c in use if c not in cat]
    Xn = sub[[x] + num].astype(float)
    if cat:
        Xc = pd.get_dummies(sub[cat], drop_first=True).astype(float)
        X = pd.concat([Xn, Xc], axis=1)
    else:
        X = Xn
    X = sm.add_constant(X)
    try:
        m = sm.Logit(sub[y].astype(int), X).fit(disp=False, maxiter=100)
    except Exception as e:
        LOG.warning("Logit failed (%s ~ %s): %s", y, x, e)
        return {"OR": np.nan, "LCI": np.nan, "UCI": np.nan, "P": np.nan, "N": len(sub), "n_event": int(sub[y].sum())}
    coef = m.params.get(x, np.nan)
    se = m.bse.get(x, np.nan)
    p = m.pvalues.get(x, np.nan)
    return {
        "OR": float(np.exp(coef)),
        "LCI": float(np.exp(coef - 1.96 * se)),
        "UCI": float(np.exp(coef + 1.96 * se)),
        "P": float(p),
        "N": int(len(sub)),
        "n_event": int(sub[y].sum()),
    }


def _poisson(df: pd.DataFrame, y: str, x: str, covars: list[str]) -> dict:
    """Poisson regression of count y ~ x + covars; returns rate-ratio for x."""
    import statsmodels.api as sm
    use = [c for c in covars if c in df.columns]
    Xcols = [x] + use
    sub = df[[y] + Xcols].dropna()
    cat = [c for c in use if not pd.api.types.is_numeric_dtype(sub[c])]
    num = [c for c in use if c not in cat]
    Xn = sub[[x] + num].astype(float)
    if cat:
        Xc = pd.get_dummies(sub[cat], drop_first=True).astype(float)
        X = pd.concat([Xn, Xc], axis=1)
    else:
        X = Xn
    X = sm.add_constant(X)
    try:
        m = sm.GLM(sub[y].astype(float), X, family=sm.families.Poisson()).fit()
    except Exception as e:
        LOG.warning("Poisson failed (%s ~ %s): %s", y, x, e)
        return {"RR": np.nan, "LCI": np.nan, "UCI": np.nan, "P": np.nan, "N": len(sub)}
    coef = m.params.get(x, np.nan)
    se = m.bse.get(x, np.nan)
    p = m.pvalues.get(x, np.nan)
    return {
        "RR": float(np.exp(coef)),
        "LCI": float(np.exp(coef - 1.96 * se)),
        "UCI": float(np.exp(coef + 1.96 * se)),
        "P": float(p), "N": int(len(sub)),
    }


def section_1a_baseline_hfrs(full_df, df_aging, fig_dir: Path, tab_dir: Path):
    """OLS of baseline HFRS (2-y window) on col6a3_scaled."""
    LOG.info("Section 1a: baseline HFRS (2-y window) ~ col6a3_scaled")
    long_df = parse_icd_with_dates(full_df)
    hfrs = cumulative_hfrs_baseline(long_df, window_days_back=730)
    full_df["eid"] = full_df["eid"].astype(str)
    df = df_aging.copy()
    df["eid"] = df["eid"].astype(str)
    df = df.merge(hfrs, on="eid", how="left")
    df["hfrs_baseline"] = df["hfrs_baseline"].fillna(0.0)

    # Histogram
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.hist(df["hfrs_baseline"].clip(upper=20), bins=40, color="#444", alpha=0.85)
    ax.set_xlabel("Baseline HFRS (2-y window)")
    ax.set_ylabel("Number of participants")
    ax.set_title(f"Baseline HFRS distribution (median={df['hfrs_baseline'].median():.2f})")
    plt.tight_layout()
    save_figure(fig, str(fig_dir / "section1_hfrs_baseline_hist"))
    plt.close(fig)

    rows = []
    for model_name, covars in ADJ_MODELS:
        use = [c for c in covars if c in df.columns]
        for outcome, label in [
            ("hfrs_baseline", "HFRS (linear)"),
            ("log_hfrs_baseline", "log(HFRS+0.1)"),
        ]:
            if outcome == "log_hfrs_baseline":
                df["log_hfrs_baseline"] = np.log(df["hfrs_baseline"] + 0.1)
            r = _ols(df, outcome, "col6a3_scaled", use)
            r.update({"outcome": label, "model": model_name})
            rows.append(r)

    out = pd.DataFrame(rows)
    save_dataframe(out, str(tab_dir / "section1_baseline_hfrs_OLS.csv"))

    # Forest of beta-per-SD across models for the linear outcome
    sub = out[out["outcome"] == "HFRS (linear)"].copy()
    if not sub.empty:
        sub["label"] = sub["model"]
        sub = sub.rename(columns={"beta": "HR"})
        sub["HR_e"] = sub["HR"]
        try:
            fig, ax = plt.subplots(figsize=(6, 2.6))
            ys = np.arange(len(sub))
            ax.errorbar(sub["HR"], ys,
                        xerr=[sub["HR"] - sub["LCI"], sub["UCI"] - sub["HR"]],
                        fmt="o", color="#222", ecolor="#666", capsize=2)
            ax.axvline(0.0, color="#bbb", linestyle="--")
            ax.set_yticks(ys); ax.set_yticklabels(sub["label"], fontsize=10)
            ax.set_xlabel("Beta per SD ETP (HFRS units)")
            ax.set_title("Baseline HFRS (2-y) ~ ETP")
            plt.tight_layout()
            save_figure(fig, str(fig_dir / "section1_baseline_hfrs_forest"))
            plt.close(fig)
        except Exception as e:
            LOG.warning("HFRS forest failed: %s", e)
    return out


def section_1b_prevalent_syndromes(df_aging, fig_dir: Path, tab_dir: Path):
    """Logistic regression of prevalent frailty-syndrome flags on col6a3_scaled."""
    LOG.info("Section 1b: prevalent frailty syndromes ~ col6a3_scaled")
    df = df_aging.copy()
    df["prevalent_any_frailty"] = (
        df.get("prevalent_hip_fracture", 0).astype(int)
        | df.get("prevalent_falls", 0).astype(int)
        | df.get("prevalent_delirium", 0).astype(int)
        | df.get("prevalent_pressure_ulcer", 0).astype(int)
    ).astype(int)

    rows = []
    for model_name, covars in ADJ_MODELS:
        for outcome, title in [
            ("prevalent_hip_fracture", "Prevalent hip fracture"),
            ("prevalent_falls", "Prevalent falls"),
            ("prevalent_delirium", "Prevalent delirium"),
            ("prevalent_pressure_ulcer", "Prevalent pressure ulcer"),
            ("prevalent_any_frailty", "Any prevalent frailty syndrome"),
        ]:
            r = _logit(df, outcome, "col6a3_scaled", covars)
            r.update({"outcome": title, "model": model_name})
            rows.append(r)

    out = pd.DataFrame(rows)
    save_dataframe(out, str(tab_dir / "section1_prevalent_syndromes_OR.csv"))

    # Forest of OR for +Clinical model
    sub = out[(out["model"] == "+Clinical") & out["OR"].notna()].copy()
    if not sub.empty:
        sub["label"] = sub["outcome"]
        try:
            fig, ax = plt.subplots(figsize=(7, 3.2))
            ys = np.arange(len(sub))
            ax.errorbar(sub["OR"], ys,
                        xerr=[sub["OR"] - sub["LCI"], sub["UCI"] - sub["OR"]],
                        fmt="o", color="#222", ecolor="#666", capsize=2)
            ax.axvline(1.0, color="#bbb", linestyle="--")
            ax.set_yticks(ys); ax.set_yticklabels(sub["label"], fontsize=10)
            ax.set_xscale("log")
            ax.set_xlabel("OR per SD ETP (+Clinical)")
            ax.set_title("Prevalent frailty syndromes ~ ETP")
            plt.tight_layout()
            save_figure(fig, str(fig_dir / "section1_prevalent_syndromes_forest"))
            plt.close(fig)
        except Exception as e:
            LOG.warning("Prevalent syndromes forest failed: %s", e)
    return out


def section_1c_multimorbidity(full_df, df_aging, tab_dir: Path):
    """Poisson regression of multimorbidity counts on col6a3_scaled."""
    LOG.info("Section 1c: multimorbidity ~ col6a3_scaled")
    df = df_aging.copy()
    df["selfrep_count"] = selfrep_count(full_df).reindex(df.index).fillna(0).astype(int) \
        if "self_reported_conditions" in full_df.columns else 0
    # Recompute on df_aging order:
    df["selfrep_count"] = selfrep_count(df).astype(int)
    df["prevalent_disease_count"] = prevalent_disease_count(df).astype(int)

    rows = []
    for model_name, covars in ADJ_MODELS:
        for outcome, title in [
            ("selfrep_count", "Self-reported condition count"),
            ("prevalent_disease_count", "Prevalent disease count (8 conditions)"),
        ]:
            r = _poisson(df, outcome, "col6a3_scaled", covars)
            r.update({"outcome": title, "model": model_name})
            rows.append(r)

    out = pd.DataFrame(rows)
    save_dataframe(out, str(tab_dir / "section1_multimorbidity_count.csv"))
    return out


def section_1d_bioage(df_aging, tab_dir: Path):
    """Lab-based bio-age proxy: principal component of standardised
    age, eGFR, CRP, HbA1c, HDL-C, total-C, NT-proBNP. Then OLS of
    composite and (composite - age) on col6a3_scaled."""
    LOG.info("Section 1d: bio-age proxy ~ col6a3_scaled")
    biomarkers = ["age", "egfr", "CRP", "hba1c_prevent",
                  "hdl_c_prevent", "total_c_prevent", "ntprobnp"]
    use = [c for c in biomarkers if c in df_aging.columns]
    if len(use) < 4:
        LOG.warning("Bio-age proxy: insufficient biomarker columns; skipping")
        return pd.DataFrame()
    M = df_aging[use].apply(pd.to_numeric, errors="coerce")
    # Sign convention: higher value -> older. eGFR and HDL go down with age, so flip them.
    flip = {"egfr": -1.0, "hdl_c_prevent": -1.0}
    for c in use:
        if c in flip:
            M[c] = M[c] * flip[c]
    M_z = (M - M.mean(skipna=True)) / M.std(skipna=True)
    bio = M_z.mean(axis=1, skipna=True) * np.sqrt(len(use))  # composite z-sum
    df = df_aging.copy()
    df["bioage_z"] = bio
    df["delta_bioage"] = bio - ((df["age"] - df["age"].mean()) / df["age"].std())
    rows = []
    for model_name, covars in ADJ_MODELS:
        for outcome, title in [
            ("bioage_z", "Bio-age composite (z)"),
            ("delta_bioage", "Bio-age − age (delta)"),
        ]:
            # Drop 'age' and 'egfr' from covars when these are inside the composite
            cov_use = [c for c in covars if c not in {"age", "egfr",
                                                       "hba1c_prevent",
                                                       "hdl_c_prevent",
                                                       "total_c_prevent",
                                                       "CRP", "ntprobnp"}]
            r = _ols(df, outcome, "col6a3_scaled", cov_use)
            r.update({"outcome": title, "model": model_name})
            rows.append(r)
    out = pd.DataFrame(rows)
    save_dataframe(out, str(tab_dir / "section1_bioage_proxy.csv"))
    return out


def main(args):
    out = Path(args.outdir)
    fig_dir = out / "figures"; fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir = out / "tables";  tab_dir.mkdir(parents=True, exist_ok=True)

    LOG.info("Loading %s", args.input)
    full = load_full_df(input_path=args.input)
    full = extend_with_baseline_phenotypes(full, args.ukb_extract)
    full, df_aging, c6m, c6s = prepare_aging_cohort(full)

    section_1a_baseline_hfrs(full, df_aging, fig_dir, tab_dir)
    section_1b_prevalent_syndromes(df_aging, fig_dir, tab_dir)
    section_1c_multimorbidity(full, df_aging, tab_dir)
    section_1d_bioage(df_aging, tab_dir)

    # Optional grip/walk/FI block (no-op unless ukb_extract supplied real columns)
    optional_cols = ["grip_max", "grip_asymmetry", "walk_pace",
                     "williams_FI", "fried_score", "frail_fried"]
    avail = [c for c in optional_cols if c in df_aging.columns]
    if avail:
        LOG.info("Section 1e: ETP vs UKB raw aging phenotypes: %s", avail)
        rows = []
        for model_name, covars in ADJ_MODELS:
            for outcome in avail:
                if outcome == "frail_fried":
                    r = _logit(df_aging, outcome, "col6a3_scaled", covars)
                    r.update({"outcome": outcome, "model": model_name, "kind": "OR"})
                else:
                    r = _ols(df_aging, outcome, "col6a3_scaled", covars)
                    r.update({"outcome": outcome, "model": model_name, "kind": "beta"})
                rows.append(r)
        save_dataframe(pd.DataFrame(rows), str(tab_dir / "section1e_grip_walk_FI.csv"))
    else:
        LOG.info("Section 1e: no UKB extract columns found; skipping grip/walk/FI block")

    LOG.info("Section 1 done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="/humgen/diabetes2/users/satoshi/misc/01.pv/curated_stats.tsv.gz")
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--ukb-extract", default=None,
                        help="Optional path to a UKB extract with grip/walk/FI fields by eid")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    np.random.seed(args.seed)
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    main(args)
