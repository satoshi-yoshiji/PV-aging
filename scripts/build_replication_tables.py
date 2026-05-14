"""Generate replication-specific supplementary tables: ICD-10 endpoint
definitions (Table S9), the full 109-code Gilbert 2018 HFRS table
(Table S10), and the adjustment-model covariate sets (Table S11).

Outputs CSVs into manuscript/supplementary_tables/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/humgen/diabetes2/users/satoshi/misc/02.aging")
sys.path.insert(0, str(ROOT))

from aging_shared import (
    ICD10_PATTERNS, HFRS_CODES, ADJ_BASE, ADJ_CLIN, ADJ_BIO,
)

OUT = ROOT / "manuscript" / "supplementary_tables"
OUT.mkdir(parents=True, exist_ok=True)


# ---- Table S9: ICD-10 endpoint definitions ---------------------------------
# Display-friendly description of each pattern + the actual regex used by the
# parser (which operates on a comma-prefixed normalised ICD string).
ENDPOINT_DESCRIPTIONS = {
    "HF":               ("Heart failure (Zenin component)",
                         "I50*",
                         "All four-character children of I50 (e.g. I500, I501, I509)."),
    "MI":               ("Myocardial infarction",
                         "I21*, I22*",
                         "I21 (acute MI) and I22 (subsequent MI)."),
    "COPD":             ("Chronic obstructive pulmonary disease",
                         "J40 - J44",
                         "Bronchitis and emphysema codes used by the Zenin healthspan composite."),
    "stroke":           ("Stroke",
                         "I60-I64",
                         "Acute cerebrovascular events: subarachnoid haemorrhage, intracerebral haemorrhage, "
                         "other nontraumatic intracranial haemorrhage, cerebral infarction, and stroke not "
                         "specified; excludes I69 sequelae."),
    "dementia":         ("Dementia (any subtype)",
                         "F00 - F03, G30*",
                         "F00-F03 cover Alzheimer / vascular / other dementias; G30 = Alzheimer's disease."),
    "T2D":              ("Type-2 diabetes (Zenin component)",
                         "E11*",
                         "Note: incident T2D is also taken from the pre-computed `incident_T2D_case_control` column."),
    "cancer":           ("Cancer (Zenin component)",
                         "C00-C97, D00-D09, D37-D48",
                         "Malignant neoplasms (C00-C97), in-situ neoplasms (D00-D09), and neoplasms of "
                         "uncertain/unknown behaviour (D37-D48); excludes benign neoplasms (D10-D36)."),
    "hip_fracture":     ("Frailty syndrome: hip fracture",
                         "S72.0-S72.2",
                         "Proximal femur / hip fracture codes: femoral neck, pertrochanteric/intertrochanteric, "
                         "and subtrochanteric fracture; excludes shaft/lower-end femur fractures."),
    "falls":            ("Frailty syndrome: falls",
                         "W00 - W19, R29.6",
                         "External-cause fall codes plus R29.6 (history of falling)."),
    "delirium":         ("Frailty syndrome: delirium",
                         "F05*",
                         "Delirium not induced by alcohol / other substances."),
    "pressure_ulcer":   ("Frailty syndrome: pressure ulcer",
                         "L89*",
                         "Decubitus / pressure-ulcer codes."),
    "any_CVD_event":    ("Any incident CVD event (cause-specific-mortality substitute)",
                         "I00 - I99",
                         "Any first occurrence in the ICD-10 circulatory chapter; reported with caveat in the manuscript."),
    "any_cancer_event": ("Any incident cancer (cause-specific-mortality substitute)",
                         "C00-C97, D00-D09, D37-D48",
                         "Same code set as the cancer Zenin component, used as a cause-specific-mortality substitute."),
    "CKD":              ("Chronic kidney disease",
                         "N18*",
                         "ICD-10 CKD codes; matches the published 01.pv (CKD) pipeline definition."),
}

rows = []
for ep_id, regex in ICD10_PATTERNS.items():
    name, code_summary, notes = ENDPOINT_DESCRIPTIONS.get(
        ep_id, (ep_id, "", ""))
    rows.append({
        "endpoint_id": ep_id,
        "endpoint_label": name,
        "icd10_codes": code_summary,
        "regex_used": regex,
        "notes": notes,
    })
df_s9 = pd.DataFrame(rows)
df_s9.to_csv(OUT / "Table_S9_ICD10_endpoint_definitions.csv", index=False)
print(f"Saved Table S9 with {len(df_s9)} endpoints")


# ---- Table S10: Gilbert 2018 HFRS code table ------------------------------
# 109 ICD-10 three-character prefixes + their Gilbert 2018 weights, with the
# inheritance rule (longer 4-5-character child codes inherit the parent weight)
# documented in the file's docstring column.
hfrs_rows = []
for code3, weight in sorted(HFRS_CODES.items(), key=lambda kv: -kv[1]):
    hfrs_rows.append({"icd10_3char_prefix": code3, "weight": weight})
df_s10 = pd.DataFrame(hfrs_rows)
df_s10["weight_threshold_intermediate_risk"] = ">= 5.0"
df_s10["weight_threshold_high_risk"] = ">= 15.0"
df_s10["inheritance_rule"] = (
    "4- and 5-character child codes inherit the parent's weight via 3-char prefix match."
)
df_s10["reference"] = (
    "Gilbert SR et al. Lancet 2018;391:1775-82, Supplementary Table 2."
)
df_s10.to_csv(OUT / "Table_S10_Gilbert_HFRS_codes.csv", index=False)
print(f"Saved Table S10 with {len(df_s10)} HFRS codes "
      f"(total weight = {df_s10['weight'].sum():.1f})")


# ---- Table S11: adjustment-model covariate sets ---------------------------
def model_row(name, label, covars, role):
    return {
        "model_id": name,
        "manuscript_label": label,
        "role": role,
        "n_covariates": len(covars),
        "covariates": "; ".join(covars),
    }

s11_rows = [
    model_row("Base", "Base", ADJ_BASE,
              "minimal adjustment for technical and demographic covariates"),
    model_row("+Clinical", "Clinical",
              [c for c in ADJ_CLIN if c not in ADJ_BASE],
              "added clinical-risk-factor covariates (cumulative on top of Base)"),
    model_row("+Biomarkers", "Clinical + Biomarkers (PRIMARY)",
              [c for c in ADJ_BIO if c not in ADJ_CLIN],
              "added CRP and NT-proBNP (cumulative on top of Clinical); pre-specified primary analysis"),
]
df_s11 = pd.DataFrame(s11_rows)
df_s11.to_csv(OUT / "Table_S11_adjustment_models.csv", index=False)
print(f"Saved Table S11 with {len(df_s11)} models")
print(f"   Base       has {len(ADJ_BASE)} covariates")
print(f"   +Clinical  has {len(ADJ_CLIN)} cumulative covariates")
print(f"   +Biomarkers has {len(ADJ_BIO)} cumulative covariates")
