# PV-aging: ETP and aging-related outcomes (UKB-PPP)

Analysis pipeline for circulating endotrophin (ETP, encoded by `col6a3` in
this dataset) and aging-related outcomes in the UK Biobank Pharma Proteomics
Project (UKB-PPP, n ≈ 53k baseline; n ≈ 44.6k after minimal exclusions).

This repository extends the previously published CKD (`01.pv`) and HF
(`01.pv_hf`) ETP analyses to a panel of aging endpoints:

- All-cause mortality
- Hospital Frailty Risk Score (Gilbert 2018) — incident time-to-first
  cumulative-score-≥5 admission
- Healthspan composite (Zenin 2019; first of HF/MI/COPD/stroke/dementia/T2D/
  cancer/death)
- Frailty-syndrome hospitalisations: hip fracture, falls, delirium, pressure
  ulcer (separate Cox each)
- Cardiovascular and cancer event composites (cause-specific-mortality
  substitutes; flagged with a caveat)
- Disease panel: incident T2D, HF, CKD, dementia
- Cross-sectional: baseline HFRS, prevalent frailty syndromes,
  multimorbidity counts, lab-based bio-age proxy

## Pipeline

```
aging_shared.py                    # ICD/dxdate parser, HFRS table, endpoint helpers
run_section1_baseline_aging.py     # Cross-sectional aging proxies
run_section2_aging_prospective.py  # Cox panel (13 endpoints x 3 adjustment models)
run_all_aging.py                   # Top-level wrapper
```

Reuses helpers from `/humgen/diabetes2/users/satoshi/misc/01.pv/ckd_shared.py`
(Cox / KM / forest / cumulative-incidence machinery).

## Data

Input: `/humgen/diabetes2/users/satoshi/misc/01.pv/curated_stats.tsv.gz`
(45,479 rows × 3,525 columns; UKB-PPP Olink Explore 3072 + curated phenotypes).

Endpoints are derived on the fly:

- `ICD10` column is a comma-separated lifetime list of HES-linked ICD-10
  codes per participant; `dxdate_a0..a258` are the paired dates (1:1 with
  comma positions).
- The Gilbert 2018 HFRS table (109 codes; 3-character prefixes inheriting
  weight to 4–5-char children) is encoded in `aging_shared.HFRS_CODES`.

## Adjustment models

| Model        | Covariates                                                                 |
| ------------ | --------------------------------------------------------------------------- |
| Base         | age, sex, centre, time_to_olink_processing, Batch                           |
| +Clinical    | + BMI, eGFR, SBP, ever_smoked, diabetes, hypertension, dyslipidemia, HbA1c, HDL-C, total-C |
| +Biomarkers  | + CRP, NT-proBNP                                                            |

Time scale: years from `date_attending_centre`. Censoring follows the
`rowwise_max(time_to_lostfollow, time_to_anydeath, time_to_latest_icd10,
max_follow_up)` recipe used in `01.pv`. Cox penalizer = 0.01 (matches
`01.pv_hf`).

## Outputs

```
results/
  tables/
    section1_baseline_hfrs_OLS.csv
    section1_prevalent_syndromes_OR.csv
    section1_multimorbidity_count.csv
    section1_bioage_proxy.csv
    section2_per_sd.csv
    section2_quartiles.csv
    section2_cohort_flow.csv
    section2_hfrs_top_contributors.csv
    section2_zenin_components.csv
  figures/
    section1_*.pdf
    section2_<endpoint>_km.pdf, section2_<endpoint>_adjcuminc.pdf
    section2_master_forest.pdf
  source_data/
    section2_<endpoint>_*.csv
Results_and_Methods.md             # manuscript-style write-up
```

## Running

```bash
conda activate ckd_analysis
cd /humgen/diabetes2/users/satoshi/misc/02.aging
python run_all_aging.py \
  --input /humgen/diabetes2/users/satoshi/misc/01.pv/curated_stats.tsv.gz \
  --outdir results --seed 42
```

For UGER cluster submission, see `uger/run_section2.qsub` (Section 2 is the
heavy step; Section 1 finishes in a few minutes locally).

## Notes

- HFRS is implemented with the Gilbert 2018 weights using distinct ICD
  dxdates as the unit of accumulation (`curated_stats.tsv.gz` lacks
  `hesin_diag` admission grouping). This mildly underestimates HFRS vs
  admission-level scoring; documented in the module docstring.
- Cause-specific (CVD, cancer) mortality is **not** computed because the
  underlying cause-of-death code (UKB f.40001) is absent from
  `curated_stats.tsv.gz`. Reported as `any_CVD_event` / `any_cancer_event`
  composites with a `caveat="ICD-event substitute (no f.40001)"` tag.
- Williams 49-item FI, Fried phenotype, grip strength, walking pace,
  hand-grip asymmetry: a `extend_with_baseline_phenotypes(df, path)` hook
  is in place to merge a future UKB extract by eid; the Section 1 runner
  automatically activates an extra block when those columns become
  available.

## Citations

- Gilbert SR et al. _Lancet_ 2018;391:1775–82.
- Zenin A et al. _Commun Biol_ 2019;2:41.
- Inker LA et al. _N Engl J Med_ 2021;385:1737–49 (CKD-EPI 2021).
