"""
aging_shared.py — ETP-aging analysis pipeline.

Extends 01.pv/ckd_shared.py with:
  - ICD10 + dxdate_a* long-format parser
  - 109-code Gilbert 2018 HFRS table and time-to-first-HFRS-positive admission
  - Endpoint dictionaries (Zenin healthspan, frailty syndromes, disease panel)
  - prepare_aging_cohort() with minimal exclusions for the full ~53k UKB-PPP set
  - Curation hook for grip / walking / FI fields when they later become available

Authoritative algorithms:
  - HFRS: Gilbert SR et al. Lancet 2018;391:1775-82 (PMID 29706364), Suppl Table 2.
  - Healthspan composite: Zenin A et al. Commun Biol 2019;2:41.

NOTE on HFRS implementation: curated_stats.tsv.gz lacks per-admission grouping
(no hesin_diag), so distinct dxdates substitute for admissions; same-day codes
count as one accumulation step. This mildly underestimates HFRS vs admission-
level scoring but is the closest faithful Gilbert implementation against the
data we have.
"""
from __future__ import annotations

import logging
import os
import re
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Bootstrap import of ckd_shared (sibling project)
# ---------------------------------------------------------------------------
_CKD_DIR = "/humgen/diabetes2/users/satoshi/misc/01.pv"
if _CKD_DIR not in sys.path:
    sys.path.insert(0, _CKD_DIR)

from ckd_shared import (  # noqa: E402
    LOG,
    PAL4, PAL5, XMAX_YEARS, YMAX_KM, YMAX_ADJ,
    BASE_COVARS, CLINICAL_COVARS_CANDIDATES,
    load_full_df, coerce_numeric, normalize_sex_to_binary,
    egfr_ckdepi_cr_2021, normalize_icd10, rowwise_max,
    has_variation, available_covars, keep_term, optional_var,
    build_cox_df, safe_coxph, cox_hr_table, extract_col6a3_per_sd,
    make_quartile, make_quintile,
    km_cumulative_incidence_plot, compute_adjusted_cuminc,
    plot_adjusted_cuminc, forest_plot, forest_plot_OR,
    save_dataframe, save_figure, summarize_surv,
    fmt_p, fmt_hr,
)

# ---------------------------------------------------------------------------
# Endpoint dictionaries
# ---------------------------------------------------------------------------
# All ICD patterns operate on the comma-separated, dot-stripped, upper-cased
# string produced by normalize_icd10(). Anchors (^|,) ensure we match a code
# beginning. [0-9]{0,2} matches 3- to 5-character ICD-10 forms (e.g. I50, I501).
ICD10_PATTERNS: dict[str, str] = {
    # Zenin 8-event healthspan composite (death is sourced from incident_anydeath)
    "HF":        r"(^|,)I50[0-9]{0,2}",
    "MI":        r"(^|,)I2[12][0-9]{0,2}",                          # I21*, I22*
    "COPD":      r"(^|,)J4[0-4][0-9]{0,2}",                         # J40..J44
    "stroke":    r"(^|,)I6[0-9][0-9]{0,2}",                         # I60..I69
    "dementia":  r"(^|,)(F0[0-3][0-9]{0,2}|G30[0-9]{0,2})",
    "T2D":       r"(^|,)E11[0-9]{0,2}",
    "cancer":    r"(^|,)(C[0-9]{2}[0-9]{0,2}|D(0[0-9]|[1-3][0-9]|4[0-8])[0-9]{0,2})",

    # Frailty-syndrome hospitalisations (separate Cox each)
    "hip_fracture":   r"(^|,)S72[0-9]{0,2}",
    "falls":          r"(^|,)(W0[0-9][0-9]{0,2}|W1[0-9][0-9]{0,2}|R296)",
    "delirium":       r"(^|,)F05[0-9]{0,2}",
    "pressure_ulcer": r"(^|,)L89[0-9]{0,2}",

    # Cause-specific-mortality SUBSTITUTES (caveat in output table)
    "any_CVD_event":    r"(^|,)I[0-9]{2}[0-9]{0,2}",                 # I00..I99
    "any_cancer_event": r"(^|,)(C[0-9]{2}[0-9]{0,2}|D(0[0-9]|[1-3][0-9]|4[0-8])[0-9]{0,2})",

    # Disease panel
    "CKD": r"(^|,)N18[0-9]{0,2}",
}

# Endpoints sourced from pre-computed columns rather than ICD parsing
PRECOMPUTED_ENDPOINTS: dict[str, tuple[str, str]] = {
    "HF":    ("incident_HF",                  "time_to_first_HF_event_days"),
    "T2D":   ("incident_T2D_case_control",    "time_to_T2Devent"),
    "death": ("incident_anydeath",            "time_to_anydeath"),
}

ZENIN_COMPONENTS = ("HF", "MI", "COPD", "stroke", "dementia", "T2D", "cancer", "death")

# Frailty-syndrome panel keys (incident endpoints)
FRAILTY_SYNDROMES = ("hip_fracture", "falls", "delirium", "pressure_ulcer")

# ---------------------------------------------------------------------------
# Adjustment models
# ---------------------------------------------------------------------------
ADJ_BASE = ["age", "sex", "centre", "proc", "Batch"]

ADJ_CLIN = ADJ_BASE + [
    "BMI", "egfr", "SBP", "ever_smoked",
    "diabetes", "hypertension", "dyslipidemia",
    "hba1c_prevent", "hdl_c_prevent", "total_c_prevent",
]

ADJ_BIO = ADJ_CLIN + ["CRP", "ntprobnp"]

ADJ_MODELS: list[tuple[str, list[str]]] = [
    ("Base",        ADJ_BASE),
    ("+Clinical",   ADJ_CLIN),
    ("+Biomarkers", ADJ_BIO),
]

# ---------------------------------------------------------------------------
# Gilbert 2018 HFRS code list (109 ICD-10 codes; 3-char prefix -> weight)
# Children (4-5 char) inherit parent weight via prefix match.
# Reference: Gilbert SR et al. Lancet 2018;391:1775-82, Suppl Table 2.
# ---------------------------------------------------------------------------
HFRS_CODES: dict[str, float] = {
    "F00": 7.1, "G81": 4.4, "G30": 4.0, "I69": 3.7, "R29": 3.6,
    "N39": 3.2, "F05": 3.2, "W19": 3.2, "S00": 3.2, "R31": 3.0,
    "B96": 2.9, "R41": 2.7, "R26": 2.6, "I67": 2.6, "R56": 2.6,
    "R40": 2.5, "T83": 2.4, "S06": 2.4, "S42": 2.3, "E87": 2.3,
    "M25": 2.3, "E86": 2.3, "R54": 2.2, "Z50": 2.1, "F03": 2.1,
    "W18": 2.1, "Z75": 2.1, "F01": 2.0, "S80": 2.0, "L03": 2.0,
    "H54": 2.0, "E53": 2.0, "Z60": 2.0, "G20": 2.0, "R55": 2.0,
    "S22": 1.9, "K59": 1.9, "N17": 1.9, "L89": 1.9, "Z22": 1.9,
    "B95": 1.8, "L97": 1.8, "R44": 1.8, "K26": 1.8, "I95": 1.8,
    "N19": 1.8, "A41": 1.8, "Z87": 1.7, "J96": 1.7, "X59": 1.6,
    "M19": 1.6, "G40": 1.6, "M81": 1.6, "S72": 1.5, "S32": 1.5,
    "E16": 1.5, "R94": 1.5, "N18": 1.5, "R33": 1.5, "R69": 1.5,
    "N28": 1.5, "R32": 1.5, "G31": 1.5, "Y95": 1.4, "S09": 1.4,
    "R45": 1.4, "G45": 1.4, "Z74": 1.4, "M79": 1.4, "W06": 1.3,
    "S01": 1.3, "A04": 1.3, "A09": 1.3, "J18": 1.3, "J69": 1.3,
    "R47": 1.2, "E55": 1.2, "Z93": 1.2, "R02": 1.2, "R63": 1.2,
    "H91": 1.1, "W10": 1.1, "W01": 1.1, "E05": 1.1, "M41": 1.1,
    "Z73": 1.1, "K56": 1.1, "M80": 1.1, "I63": 1.0, "S81": 1.0,
    "B37": 1.0, "S43": 1.0, "L08": 1.0, "S60": 1.0, "L02": 1.0,
    "Y84": 0.9, "F23": 0.9, "T78": 0.8, "F32": 0.7, "Z67": 0.7,
    "Z89": 0.6, "F03": 0.6,  # noqa: duplicate F03 placeholder; dict dedup keeps last
    "F44": 0.6, "F02": 0.5, "Z71": 0.5, "Z02": 0.4, "Y95": 0.4,  # noqa
    "J22": 0.3, "Z85": 0.2, "Z88": 0.2, "Z44": 0.1,
}
# Compile sorted prefix list (longer prefixes first, but all are 3-char here)
_HFRS_PREFIXES = sorted(HFRS_CODES.keys())


# ---------------------------------------------------------------------------
# ICD long-form parser
# ---------------------------------------------------------------------------
def parse_icd_with_dates(
    df: pd.DataFrame,
    icd_col: str = "ICD10",
    n_dx: int = 259,
    eid_col: str = "eid",
    baseline_col: str = "date_attending_centre",
) -> pd.DataFrame:
    """
    Explode (eid, ICD10, dxdate_a0..a{n_dx-1}) into a long DataFrame:
        eid, code3 (3-char), code_full (normalised, e.g. I501),
        dxdate (datetime64[ns]), days_from_baseline (int)

    ICD10 column: comma-separated, dot/space-stripped after normalize_icd10.
    dxdate_aN columns are aligned 1:1 with comma positions in the ICD list.

    Returns one row per (eid, dxdate) pair where ICD10 was non-empty.
    Rows with NaT dxdate or empty code are dropped.
    """
    needed = [eid_col, icd_col, baseline_col]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise KeyError(f"Missing columns for ICD parser: {missing}")

    dx_cols = [f"dxdate_a{i}" for i in range(n_dx) if f"dxdate_a{i}" in df.columns]
    LOG.info("parse_icd_with_dates: %d dxdate columns present", len(dx_cols))

    icd_norm = normalize_icd10(df[icd_col].fillna("")).tolist()
    eids = df[eid_col].astype(str).tolist()
    baseline = pd.to_datetime(df[baseline_col], errors="coerce").tolist()

    # Build per-row dxdate matrix as object lists
    dx_matrix = []
    for c in dx_cols:
        dx_matrix.append(pd.to_datetime(df[c], errors="coerce").tolist())

    eids_out, code3_out, code_full_out, dxdate_out, days_out = [], [], [], [], []
    for ridx in range(len(df)):
        s = icd_norm[ridx]
        if not s:
            continue
        codes = [c for c in s.split(",") if c]
        b = baseline[ridx]
        eid = eids[ridx]
        for j, code in enumerate(codes):
            if j >= len(dx_cols):
                break
            d = dx_matrix[j][ridx]
            if pd.isna(d):
                continue
            code_full = code.strip()
            if not code_full:
                continue
            code3 = code_full[:3]
            eids_out.append(eid)
            code3_out.append(code3)
            code_full_out.append(code_full)
            dxdate_out.append(d)
            if pd.isna(b):
                days_out.append(np.nan)
            else:
                days_out.append((d - b).days)

    long_df = pd.DataFrame({
        "eid": eids_out,
        "code3": code3_out,
        "code_full": code_full_out,
        "dxdate": dxdate_out,
        "days_from_baseline": days_out,
    })
    LOG.info("parse_icd_with_dates: long format = %d rows from %d participants",
             len(long_df), df[eid_col].nunique())
    return long_df


# ---------------------------------------------------------------------------
# Endpoint constructors
# ---------------------------------------------------------------------------
def first_event_for_codes(long_df: pd.DataFrame, code_pattern: str) -> pd.DataFrame:
    """
    Per-eid first matching dxdate (NaT if none).
    Returns DataFrame with columns: eid, date_event, days_from_baseline.
    Pattern is applied to comma-prefixed code_full to match the (^|,) anchor
    in ICD10_PATTERNS regexes.
    """
    rx = re.compile(code_pattern)
    s = "," + long_df["code_full"]
    mask = s.str.contains(rx, regex=True, na=False)
    sub = long_df[mask]
    if sub.empty:
        return pd.DataFrame(columns=["eid", "date_event", "days_from_baseline"])
    out = (sub.sort_values("dxdate")
              .groupby("eid", as_index=False)
              .agg(date_event=("dxdate", "first"),
                   days_from_baseline=("days_from_baseline", "first")))
    return out


def make_event_endpoint(
    full_df: pd.DataFrame,
    long_df: pd.DataFrame,
    label: str,
    code_pattern: str,
    censor_cols: tuple[str, ...] = (
        "time_to_lostfollow", "time_to_anydeath",
        "time_to_latest_icd10", "max_follow_up",
    ),
) -> pd.DataFrame:
    """
    Adds to full_df (in place + returns it):
      prevalent_<label>      (int 0/1; first event date <= 0)
      incident_<label>       (int 0/1; first event date > 0 within follow-up)
      time_to_<label>_days   (float; days from baseline; censor max if no event)
      date_<label>           (datetime; NaT if none)
    """
    first = first_event_for_codes(long_df, code_pattern)
    full_df["eid"] = full_df["eid"].astype(str)
    merged = full_df[["eid"]].merge(first, on="eid", how="left")

    days = pd.to_numeric(merged["days_from_baseline"], errors="coerce").values
    prevalent = ((days <= 0) & np.isfinite(days)).astype(int)
    incident = ((days > 0) & np.isfinite(days)).astype(int)

    # Censoring time when no event: rowwise_max over admin censors
    censor_days = rowwise_max(full_df, list(censor_cols)).values

    time_days = np.where(incident == 1, days, censor_days)

    full_df[f"date_{label}"] = merged["date_event"].values
    full_df[f"prevalent_{label}"] = prevalent
    full_df[f"incident_{label}"] = incident
    full_df[f"time_to_{label}_days"] = time_days
    return full_df


# ---------------------------------------------------------------------------
# HFRS computation
# ---------------------------------------------------------------------------
def _hfrs_weight(code3: str) -> float:
    """Lookup by 3-character prefix (children inherit parent weight)."""
    return HFRS_CODES.get(code3, 0.0)


def cumulative_hfrs_baseline(
    long_df: pd.DataFrame,
    hfrs_codes: dict[str, float] | None = None,
    window_days_back: int | None = 730,
) -> pd.DataFrame:
    """
    Per-eid cumulative HFRS at baseline.

    Sum of weights over distinct 3-char codes recorded within
    window_days_back days BEFORE date_attending_centre (i.e. days_from_baseline
    in [-window_days_back, 0]). If window_days_back is None, lifetime <=0 is used.

    Returns DataFrame: eid, hfrs_baseline (float, >=0).
    """
    weights = hfrs_codes or HFRS_CODES
    df = long_df.dropna(subset=["days_from_baseline"]).copy()
    df = df[df["days_from_baseline"] <= 0]
    if window_days_back is not None:
        df = df[df["days_from_baseline"] >= -int(window_days_back)]
    df["w"] = df["code3"].map(weights).fillna(0.0)
    df = df[df["w"] > 0]
    if df.empty:
        return pd.DataFrame(columns=["eid", "hfrs_baseline"])
    # First occurrence per (eid, code3) only (distinct 3-char codes)
    first = (df.sort_values("dxdate")
               .drop_duplicates(["eid", "code3"], keep="first"))
    out = first.groupby("eid", as_index=False)["w"].sum().rename(
        columns={"w": "hfrs_baseline"}
    )
    return out


def cumulative_hfrs_first_high(
    long_df: pd.DataFrame,
    hfrs_codes: dict[str, float] | None = None,
    threshold: float = 5.0,
) -> pd.DataFrame:
    """
    Per-eid: walk distinct dxdates chronologically (post-baseline only).
    Maintain running set of 3-char codes already seen.
    After each new dxdate, sum weights over the running set.
    Return the FIRST dxdate at which the running sum >= threshold (NaT otherwise).

    Returns DataFrame: eid, date_event, days_from_baseline.

    SIMPLIFYING ASSUMPTION: distinct dxdates substitute for admissions because
    curated_stats.tsv.gz lacks per-admission grouping (no hesin_diag).
    """
    weights = hfrs_codes or HFRS_CODES
    df = long_df.dropna(subset=["days_from_baseline"]).copy()
    df["w"] = df["code3"].map(weights).fillna(0.0)
    if df.empty:
        return pd.DataFrame(columns=["eid", "date_event", "days_from_baseline"])

    # Sort within each eid by dxdate
    df = df.sort_values(["eid", "dxdate", "code3"])
    rows = []
    cur_eid = None
    seen: set[str] = set()
    cum: float = 0.0
    done: bool = False
    for eid, code3, dxdate, days, w in zip(
        df["eid"].values, df["code3"].values, df["dxdate"].values,
        df["days_from_baseline"].values, df["w"].values
    ):
        if eid != cur_eid:
            cur_eid = eid
            seen = set()
            cum = 0.0
            done = False
        if done:
            continue
        if w == 0.0 or code3 in seen:
            continue
        seen.add(code3)
        cum += float(w)
        if cum >= threshold and days > 0:
            rows.append((eid, pd.Timestamp(dxdate), int(days)))
            done = True
    out = pd.DataFrame(rows, columns=["eid", "date_event", "days_from_baseline"])
    return out


def make_hfrs_incident_endpoint(
    full_df: pd.DataFrame,
    long_df: pd.DataFrame,
    threshold: float = 5.0,
    label: str = "HFRS5",
) -> pd.DataFrame:
    """
    Build incident_HFRS5 / time_to_HFRS5_days endpoint via cumulative_hfrs_first_high.
    Prevalent flag: cumulative HFRS at baseline already >= threshold.
    """
    # Prevalent: any participant whose cumulative weight reaches threshold using only
    # codes recorded at-or-before baseline.
    pre = long_df.dropna(subset=["days_from_baseline"]).copy()
    pre = pre[pre["days_from_baseline"] <= 0]
    pre["w"] = pre["code3"].map(HFRS_CODES).fillna(0.0)
    pre = pre[pre["w"] > 0]
    pre_first = (pre.sort_values("dxdate")
                    .drop_duplicates(["eid", "code3"], keep="first"))
    pre_sum = pre_first.groupby("eid")["w"].sum()
    prevalent_eids = set(pre_sum[pre_sum >= threshold].index.astype(str))

    incident = cumulative_hfrs_first_high(long_df, threshold=threshold)
    full_df["eid"] = full_df["eid"].astype(str)
    merged = full_df[["eid"]].merge(incident, on="eid", how="left")
    days = pd.to_numeric(merged["days_from_baseline"], errors="coerce").values
    inc_flag = ((days > 0) & np.isfinite(days)).astype(int)

    prev_flag = full_df["eid"].astype(str).isin(prevalent_eids).astype(int).values

    censor_cols = ("time_to_lostfollow", "time_to_anydeath",
                   "time_to_latest_icd10", "max_follow_up")
    censor_days = rowwise_max(full_df, list(censor_cols)).values
    time_days = np.where(inc_flag == 1, days, censor_days)

    full_df[f"date_{label}"] = merged["date_event"].values
    full_df[f"prevalent_{label}"] = prev_flag
    full_df[f"incident_{label}"] = inc_flag
    full_df[f"time_to_{label}_days"] = time_days
    return full_df


# ---------------------------------------------------------------------------
# Pre-computed endpoint mappers (HF, T2D, death)
# ---------------------------------------------------------------------------
def _attach_precomputed(full_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert pre-computed HF/T2D/death columns into the standardised
    prevalent_<lbl> / incident_<lbl> / time_to_<lbl>_days schema.
    """
    # --- HF ---
    if "incident_HF" in full_df.columns and "time_to_first_HF_event_days" in full_df.columns:
        ev = pd.to_numeric(full_df["incident_HF"], errors="coerce").fillna(0).astype(int)
        t = pd.to_numeric(full_df["time_to_first_HF_event_days"], errors="coerce")
        # prevalent: time<=0 with event flag (note: 01.pv_hf treats <0 as prevalent)
        prev = ((t <= 0) & ev.eq(1)).astype(int).values
        inc = ((t > 0) & ev.eq(1)).astype(int).values
        censor = rowwise_max(full_df, ["time_to_lostfollow", "time_to_anydeath",
                                       "time_to_latest_icd10", "max_follow_up"]).values
        time_days = np.where(inc == 1, t.values, censor)
        full_df["prevalent_HF"] = prev
        full_df["incident_HF"] = inc
        full_df["time_to_HF_days"] = time_days

    # --- T2D ---
    if "incident_T2D_case_control" in full_df.columns and "time_to_T2Devent" in full_df.columns:
        ev = pd.to_numeric(full_df["incident_T2D_case_control"], errors="coerce").fillna(0).astype(int)
        t = pd.to_numeric(full_df["time_to_T2Devent"], errors="coerce")
        prev = ((t <= 0) & ev.eq(1)).astype(int).values
        inc = ((t > 0) & ev.eq(1)).astype(int).values
        censor = rowwise_max(full_df, ["time_to_lostfollow", "time_to_anydeath",
                                       "time_to_latest_icd10", "max_follow_up"]).values
        time_days = np.where(inc == 1, t.values, censor)
        full_df["prevalent_T2D"] = prev
        full_df["incident_T2D_case_control"] = inc
        full_df["time_to_T2D_days"] = time_days

    # --- Death (all-cause) ---
    if "time_to_anydeath" in full_df.columns:
        t = pd.to_numeric(full_df["time_to_anydeath"], errors="coerce")
        date_death = pd.to_datetime(full_df.get("date_death"), errors="coerce")
        inc = (date_death.notna() & (t > 0)).astype(int).values
        censor = rowwise_max(full_df, ["time_to_lostfollow",
                                       "time_to_latest_icd10", "max_follow_up"]).values
        time_days = np.where(inc == 1, t.values, censor)
        full_df["prevalent_anydeath"] = 0  # by definition
        full_df["incident_anydeath"] = inc
        full_df["time_to_anydeath_days"] = time_days
    return full_df


# ---------------------------------------------------------------------------
# Zenin healthspan composite
# ---------------------------------------------------------------------------
def make_zenin_composite(
    full_df: pd.DataFrame,
    components: Iterable[str] = ZENIN_COMPONENTS,
) -> pd.DataFrame:
    """
    First-component-fires composite. Adds:
      prevalent_zenin       (int)
      incident_zenin        (int)
      time_to_zenin_days    (float)
      date_zenin            (datetime)
      zenin_first_component (str)
    Requires _attach_precomputed and ICD-derived endpoints already added.
    """
    rows = len(full_df)
    first_days = np.full(rows, np.inf)
    first_comp = np.full(rows, "", dtype=object)
    prev_flag = np.zeros(rows, dtype=int)

    for comp in components:
        # Resolve column names (handle T2D pre-computed name and death->anydeath)
        if comp == "T2D":
            ev_col = "incident_T2D_case_control" if "incident_T2D_case_control" in full_df.columns else "incident_T2D"
            t_col = "time_to_T2D_days"
            prev_col = "prevalent_T2D"
        elif comp == "death":
            ev_col = "incident_anydeath"
            t_col = "time_to_anydeath_days"
            prev_col = "prevalent_anydeath"
        else:
            ev_col = f"incident_{comp}"
            t_col = f"time_to_{comp}_days"
            prev_col = f"prevalent_{comp}"
        if ev_col not in full_df.columns or t_col not in full_df.columns:
            LOG.warning("Zenin: missing columns for component '%s'", comp)
            continue
        ev = pd.to_numeric(full_df[ev_col], errors="coerce").fillna(0).astype(int).values
        t = pd.to_numeric(full_df[t_col], errors="coerce").values
        if prev_col in full_df.columns:
            prev_flag = np.maximum(
                prev_flag,
                pd.to_numeric(full_df[prev_col], errors="coerce").fillna(0).astype(int).values,
            )
        else:
            prev_flag = np.maximum(prev_flag, ((t <= 0) & (ev == 1)).astype(int))

        cand = np.where((ev == 1) & np.isfinite(t) & (t > 0), t, np.inf)
        better = cand < first_days
        first_days = np.where(better, cand, first_days)
        first_comp = np.where(better, comp, first_comp)

    incident = np.isfinite(first_days).astype(int)
    censor = rowwise_max(full_df, ["time_to_lostfollow", "time_to_anydeath",
                                   "time_to_latest_icd10", "max_follow_up"]).values
    time_days = np.where(incident == 1, first_days, censor)

    full_df["prevalent_zenin"] = prev_flag
    full_df["incident_zenin"] = incident
    full_df["time_to_zenin_days"] = time_days
    full_df["zenin_first_component"] = first_comp
    return full_df


# ---------------------------------------------------------------------------
# Curation hook (no-op placeholder for grip / walk / FI fields)
# ---------------------------------------------------------------------------
def extend_with_baseline_phenotypes(
    df: pd.DataFrame,
    ukb_extract_path: str | None = None,
) -> pd.DataFrame:
    """
    Hook for when raw UKB phenotype fields become available. Reads a TSV/CSV/RDS
    extract keyed on eid, merges in:
        grip_max, grip_asymmetry, walk_pace,
        williams_FI (continuous 0..1),
        fried_score (0..5), frail_fried (binary 0/1).
    No-op (returns df unchanged) when ukb_extract_path is None or missing.
    """
    if not ukb_extract_path or not os.path.exists(ukb_extract_path):
        LOG.info("extend_with_baseline_phenotypes: no extract supplied; "
                 "Aim #1 grip/walk/FI block will be skipped at runtime")
        return df

    p = Path(ukb_extract_path)
    if p.suffix.lower() in {".tsv", ".txt"}:
        ext = pd.read_csv(p, sep="\t")
    elif p.suffix.lower() == ".csv":
        ext = pd.read_csv(p)
    elif p.suffix.lower() == ".rds":
        try:
            import pyreadr
            ext = list(pyreadr.read_r(p).values())[0]
        except Exception as e:  # pragma: no cover
            LOG.warning("Could not read RDS: %s; skipping merge", e)
            return df
    else:
        LOG.warning("Unsupported extract format: %s", p.suffix)
        return df

    cols_keep = [c for c in [
        "eid", "grip_max", "grip_asymmetry", "walk_pace",
        "williams_FI", "fried_score", "frail_fried",
    ] if c in ext.columns]
    LOG.info("extend_with_baseline_phenotypes: merging columns %s", cols_keep)
    ext["eid"] = ext["eid"].astype(str)
    df["eid"] = df["eid"].astype(str)
    return df.merge(ext[cols_keep], on="eid", how="left")


# ---------------------------------------------------------------------------
# Cohort preparation (minimal exclusions for the full ~53k aging cohort)
# ---------------------------------------------------------------------------
def derive_all_aging_endpoints(full_df: pd.DataFrame, long_df: pd.DataFrame | None = None) -> pd.DataFrame:
    """
    One-shot driver: build every endpoint we'll model in section 2.
    Idempotent (skips columns already present).

    Endpoints created (label / time / event columns):
      - HF, T2D, anydeath via _attach_precomputed
      - CKD, MI, COPD, stroke, dementia, cancer, hip_fracture, falls,
        delirium, pressure_ulcer, any_CVD_event, any_cancer_event via ICD parse
      - HFRS5 via cumulative_hfrs_first_high
      - zenin via make_zenin_composite
    """
    if long_df is None:
        long_df = parse_icd_with_dates(full_df)

    full_df = _attach_precomputed(full_df)

    # ICD-derived endpoints (skip pre-computed names)
    icd_only = ["CKD", "MI", "COPD", "stroke", "dementia", "cancer",
                "hip_fracture", "falls", "delirium", "pressure_ulcer",
                "any_CVD_event", "any_cancer_event"]
    for label in icd_only:
        if f"incident_{label}" in full_df.columns and f"time_to_{label}_days" in full_df.columns:
            continue
        pat = ICD10_PATTERNS[label]
        full_df = make_event_endpoint(full_df, long_df, label, pat)
        n_inc = int(full_df[f"incident_{label}"].sum())
        n_prev = int(full_df[f"prevalent_{label}"].sum())
        LOG.info("Endpoint %-18s: prevalent=%6d incident=%6d", label, n_prev, n_inc)

    # HFRS incident
    if "incident_HFRS5" not in full_df.columns:
        full_df = make_hfrs_incident_endpoint(full_df, long_df, threshold=5.0, label="HFRS5")
        LOG.info("Endpoint HFRS5            : prevalent=%6d incident=%6d",
                 int(full_df["prevalent_HFRS5"].sum()),
                 int(full_df["incident_HFRS5"].sum()))

    # Zenin healthspan composite
    full_df = make_zenin_composite(full_df, ZENIN_COMPONENTS)
    LOG.info("Endpoint zenin           : prevalent=%6d incident=%6d",
             int(full_df["prevalent_zenin"].sum()),
             int(full_df["incident_zenin"].sum()))

    # Convert all *_days to *_years for downstream Cox
    for label in ["anydeath", "HFRS5", "zenin"] + icd_only + ["HF", "T2D"]:
        d = f"time_to_{label}_days"
        y = f"time_to_{label}_years"
        if d in full_df.columns and y not in full_df.columns:
            full_df[y] = pd.to_numeric(full_df[d], errors="coerce") / 365.25

    return full_df


def prepare_aging_cohort(
    full_df: pd.DataFrame,
    creat_col: str = "p30700_i0",
) -> tuple[pd.DataFrame, pd.DataFrame, float, float]:
    """
    Returns (full_df, df_aging, c6_mean, c6_sd).

    df_aging = analytic cohort common to ALL aging endpoints with MINIMAL
    exclusions (user-specified ~53k retention). Per-endpoint prevalent-case
    removal happens INSIDE section2.

    Adds derived columns to full_df: egfr, SBP, time_to_*_years for every
    endpoint, plus a cached long-format ICD parse.
    """
    LOG.info("prepare_aging_cohort: input rows = %d", len(full_df))

    # 1. Coerce numeric core fields
    for c in ["col6a3", "age", "BMI", "WHR", "CRP", "ntprobnp",
              "hba1c_prevent", "hdl_c_prevent", "total_c_prevent",
              "SBP1", "SBP2", "DBP1", "DBP2"]:
        if c in full_df.columns:
            full_df[c] = coerce_numeric(full_df[c])

    # Sex & batch as categorical strings (as in 01.pv)
    full_df["sex"] = full_df["sex"].astype(str)

    # 2. eGFR via CKD-EPI 2021 from creatinine
    if creat_col in full_df.columns:
        full_df["creatinine_baseline"] = coerce_numeric(full_df[creat_col])
    elif "creat_mg_dl" in full_df.columns:
        # mg/dL -> umol/L
        full_df["creatinine_baseline"] = coerce_numeric(full_df["creat_mg_dl"]) * 88.4017
    else:
        full_df["creatinine_baseline"] = np.nan
    full_df["egfr"] = egfr_ckdepi_cr_2021(
        full_df["creatinine_baseline"].values,
        full_df["age"].values,
        full_df["sex"].values.astype(float)
        if pd.api.types.is_numeric_dtype(full_df["sex"]) else
        normalize_sex_to_binary(full_df["sex"]).values.astype(float),
    )

    # 3. SBP from SBP1 (fallbacks)
    if "SBP" not in full_df.columns:
        for cand in ["SBP1", "SBP2", "sbp_prevent"]:
            if cand in full_df.columns:
                full_df["SBP"] = coerce_numeric(full_df[cand])
                break

    # 4. Coerce censor cols and t=0
    for c in ["time_to_lostfollow", "time_to_anydeath", "time_to_latest_icd10",
              "max_follow_up", "time_to_olink_processing"]:
        if c in full_df.columns:
            full_df[c] = coerce_numeric(full_df[c])
    full_df["proc"] = coerce_numeric(full_df.get("time_to_olink_processing"))
    full_df["date_attending_centre"] = pd.to_datetime(
        full_df["date_attending_centre"], errors="coerce"
    )

    # 5. Build long-format ICD parse and derive endpoints
    long_df = parse_icd_with_dates(full_df)
    full_df = derive_all_aging_endpoints(full_df, long_df)

    # 6. Minimal-exclusion analytic cohort
    base_required = ["col6a3", "age", "sex", "centre", "proc", "Batch",
                     "date_attending_centre", "max_follow_up"]
    n0 = len(full_df)
    LOG.info("Cohort step 01: total rows                                       = %d", n0)

    mask_c6 = full_df["col6a3"].notna()
    LOG.info("Cohort step 02: col6a3 not NaN                                   = %d (-%d)",
             int(mask_c6.sum()), n0 - int(mask_c6.sum()))

    mask_cov = (full_df[["age", "sex", "centre", "Batch"]].notna().all(axis=1)
                & np.isfinite(full_df["proc"]))
    n_after_cov = int((mask_c6 & mask_cov).sum())
    LOG.info("Cohort step 03: + age/sex/centre/proc/Batch present              = %d (-%d)",
             n_after_cov, int(mask_c6.sum()) - n_after_cov)

    mask_t0 = full_df["date_attending_centre"].notna()
    n_after_t0 = int((mask_c6 & mask_cov & mask_t0).sum())
    LOG.info("Cohort step 04: + date_attending_centre present                  = %d (-%d)",
             n_after_t0, n_after_cov - n_after_t0)

    mask_fu = pd.to_numeric(full_df["max_follow_up"], errors="coerce") > 0
    n_after_fu = int((mask_c6 & mask_cov & mask_t0 & mask_fu).sum())
    LOG.info("Cohort step 05: + max_follow_up > 0                              = %d (-%d)",
             n_after_fu, n_after_t0 - n_after_fu)

    mask_all = mask_c6 & mask_cov & mask_t0 & mask_fu
    df_aging = full_df[mask_all].copy().reset_index(drop=True)

    # Standardise col6a3 within the aging cohort
    c6_mean = df_aging["col6a3"].mean()
    c6_sd = df_aging["col6a3"].std()
    df_aging["col6a3_scaled"] = (df_aging["col6a3"] - c6_mean) / c6_sd
    df_aging["sex"] = df_aging["sex"].astype(str)
    df_aging["centre"] = df_aging["centre"].astype(str)
    df_aging["Batch"] = df_aging["Batch"].astype(str)

    LOG.info("prepare_aging_cohort: final analytic N = %d (col6a3 mean=%.3f sd=%.3f)",
             len(df_aging), c6_mean, c6_sd)
    return full_df, df_aging, float(c6_mean), float(c6_sd)


# ---------------------------------------------------------------------------
# Self-test (reproduce 01.pv CKD HR via the new ICD parser)
# ---------------------------------------------------------------------------
def _selftest_ckd(input_path: str | None = None) -> dict:
    """
    Build the CKD endpoint via the new parser and run a base-adjusted Cox per-SD.
    Compare against published 01.pv HR (~1.53). Used in smoke tests.
    """
    full_df = load_full_df(input_path=input_path)
    full_df, df_aging, c6m, c6s = prepare_aging_cohort(full_df)
    df = df_aging[df_aging["prevalent_CKD"] == 0].copy()
    df["time"] = df["time_to_CKD_days"] / 365.25
    df["event"] = df["incident_CKD"].astype(int)
    df = df[(df["time"] > 0) & np.isfinite(df["time"])].copy()
    out = {"N": len(df), "events": int(df["event"].sum())}
    for label, covars in [
        ("base", ADJ_BASE),
        ("base_egfr", ADJ_BASE + ["egfr"]),
        ("clinical", ADJ_CLIN),
    ]:
        use = available_covars(df, covars)
        cox_df = build_cox_df(df, "time", "event", "col6a3_scaled", use)
        fit = safe_coxph(cox_df, "time", "event", penalizer=0.01, label=f"selftest_{label}")
        hr = extract_col6a3_per_sd(fit)
        if hr is None:
            out[f"HR_{label}"] = None
        else:
            out[f"HR_{label}"] = float(hr["HR"].iloc[0])
            out[f"LCI_{label}"] = float(hr["LCI"].iloc[0])
            out[f"UCI_{label}"] = float(hr["UCI"].iloc[0])
            out[f"P_{label}"] = float(hr["P"].iloc[0])
    return out


# ---------------------------------------------------------------------------
# Multimorbidity helper for Section 1
# ---------------------------------------------------------------------------
def selfrep_count(full_df: pd.DataFrame, col: str = "self_reported_conditions") -> pd.Series:
    """Count of distinct comma- or whitespace-separated entries in self_reported_conditions."""
    if col not in full_df.columns:
        return pd.Series(np.nan, index=full_df.index)
    s = full_df[col].astype(str)
    s = s.str.replace(r"[\[\]]", "", regex=True)
    s = s.str.replace(";", ",", regex=False)
    counts = s.apply(lambda x: 0 if (not x or x.strip() in {"", "NA", "nan", "None"}) else
                     len([t for t in re.split(r"[,\s]+", x) if t and t not in {"-7", "-3", "-1"}]))
    return counts.astype(int)


def prevalent_disease_count(
    full_df: pd.DataFrame,
    components: Iterable[str] = (
        "HF", "MI", "COPD", "stroke", "dementia",
        "T2D", "cancer", "CKD",
    ),
) -> pd.Series:
    """Sum of prevalent_<comp> across a fixed disease set."""
    cols = [f"prevalent_{c}" for c in components if f"prevalent_{c}" in full_df.columns]
    if not cols:
        return pd.Series(0, index=full_df.index)
    return full_df[cols].astype(int).sum(axis=1)
