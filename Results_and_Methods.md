# Circulating endotrophin and the prospective burden of frailty, multimorbidity and mortality in the UK Biobank Pharma Proteomics Project

Yoshiji *et al.*

---

## Abstract

Endotrophin, the C-terminal cleavage fragment of the α3 chain of type VI collagen (COL6A3), marks active extracellular-matrix remodelling and has been linked to incident heart failure and chronic kidney disease. Whether circulating endotrophin tracks the broader programmes of biological aging — frailty, healthspan compression, and the spectrum of age-related multimorbidity — has been unclear. Using Olink Explore 3072 measurements in 44,642 UK Biobank participants with up to 13 years of hospital and mortality follow-up, we show that one standard-deviation higher baseline endotrophin (col6a3, NPX) is associated, after adjustment for age, sex, body-mass index, smoking, blood pressure, eGFR, glycaemia and lipids, with **42 % higher all-cause mortality** (HR 1.42, 95 % CI 1.38–1.47; *P* = 4.7 × 10⁻¹⁰²), **26 % higher rate of incident frailty** (Hospital Frailty Risk Score ≥ 5; HR 1.26, 1.23–1.29; *P* = 5 × 10⁻⁶⁸), and **graded compression of healthspan** (Zenin 8-event composite, HR 1.15, 1.12–1.17). The signal extends to the constituent frailty syndromes (hip fracture, falls, delirium, pressure ulcer; HR 1.15–1.31), to the disease panel (HF 1.27, CKD 1.48, T2D 1.09, dementia 1.08), and to cumulative cardiovascular and cancer events. A cross-sectional analysis at baseline shows the same pattern for prior frailty syndromes, multimorbidity counts, and a lab-based biological-age proxy. Together these data position endotrophin as a single circulating index of accelerated aging that integrates frailty, mortality risk and the leading age-related diseases.

---

## Results

### Cohort

The analytic cohort comprised 44,642 UK Biobank participants with valid baseline endotrophin (col6a3, Olink Explore 3072 NPX) and core covariates (Methods; **Extended Data Table 1**). Median follow-up was 13.7 years through linked Hospital Episode Statistics (HES) records and the National Death Register. Endotrophin was approximately normally distributed within the cohort (mean 0.027 NPX, SD 0.42); per-SD effects below correspond to one SD increment within this cohort. Endpoints were derived from a long-format parsing of the per-participant ICD-10 lifetime list and aligned diagnosis-date fields (1:1 with the comma-separated codes; **Methods, ICD-10 endpoint construction**), with pre-computed columns used where available for heart failure, type 2 diabetes and all-cause mortality.

### Endotrophin tracks baseline frailty and multimorbidity

We first examined whether elevated endotrophin already accompanies frailty and multimorbidity at the time of the baseline assessment (**Fig. 1a–c**, **Extended Data Tables 2–4**).

A continuous baseline Hospital Frailty Risk Score (HFRS), computed as the sum of Gilbert *et al.* 2018 weights over distinct three-character ICD-10 codes recorded in the two years before the assessment-centre visit, increased by **β = +0.14 HFRS units per 1-SD endotrophin** in the +Clinical model (95 % CI 0.13–0.15; *P* = 7.1 × 10⁻¹²⁰; *n* = 34,248). The signal was insensitive to the choice of linear vs. log-transformed outcome (β = +0.13 on the log scale) and survived further adjustment for C-reactive protein and NT-proBNP (β = +0.13; *P* = 2.4 × 10⁻⁹⁶).

Prevalent frailty-syndrome ICD codes recorded at or before baseline showed the same pattern. Per-SD endotrophin was associated with **OR 2.28 (1.83–2.85)** for prior hip fracture, **1.34 (1.23–1.45)** for prior falls and **1.38 (1.27–1.50)** for any prior frailty syndrome (all +Clinical). Prevalent delirium (n = 5) and pressure ulcer (n = 15) were too sparse for stable estimation.

Two parallel multimorbidity indices reinforced the result: per-SD endotrophin was associated with a **2.8 % increase** in self-reported chronic-condition count (RR 1.028, 1.008–1.040; *P* = 2.7 × 10⁻⁴) and a **22 % increase** in prevalent ICD-defined disease count across eight common conditions (HF, MI, COPD, stroke, dementia, T2D, cancer, CKD; RR 1.22, 1.18–1.26; *P* = 2.6 × 10⁻³⁴).

A composite biological-age proxy — a z-sum of age, eGFR (CKD-EPI 2021), CRP, HbA1c, HDL-C, total cholesterol and NT-proBNP, with eGFR and HDL-C sign-flipped so that higher values mean older — also rose monotonically with endotrophin (**β = +0.43 SD per 1-SD endotrophin**, *P* < 1 × 10⁻³⁰⁰), as did the chronological-age-adjusted residual ("delta-bioage", β = +0.22). The cross-sectional results indicate that endotrophin is already enriched in older, more frail and more multimorbid participants at the moment of measurement.

### Endotrophin predicts all-cause mortality and incident frailty

Across 44,642 participants and up to 13.7 years of prospective follow-up, baseline endotrophin was associated with all 13 pre-specified time-to-event outcomes after multivariable adjustment (**Fig. 2**, **Table 1**). For the +Clinical model:

- **All-cause mortality** (4,293 deaths in 34,248 with complete covariates): per-SD HR 1.42 (1.38–1.47), *P* = 4.7 × 10⁻¹⁰²; concordance index 0.748. The dose-response was monotonic across endotrophin quartiles, with **Q4 vs Q1 HR 1.74 (1.59–1.90)**, *P* = 3.1 × 10⁻³⁴. Adjustment for CRP and NT-proBNP attenuated the per-SD HR to 1.32 (1.28–1.36), confirming a partly inflammation- and cardiac-stress-mediated channel but a substantial residual signal.

- **Incident frailty** (HFRS-positive admission, cumulative score ≥ 5; 8,287 events): per-SD HR 1.26 (1.23–1.29), *P* = 5.0 × 10⁻⁶⁸; C-index 0.706. The signal persisted after adjustment for CRP and NT-proBNP (HR 1.22, *P* = 1.6 × 10⁻⁴⁸), arguing for an endotrophin-specific contribution to frailty trajectories beyond systemic inflammation. The Hospital Frailty Risk Score was implemented faithful to Gilbert *et al.* 2018 using distinct ICD-10 dxdates as the unit of accumulation (Methods); the participants who eventually crossed the ≥ 5 threshold contributed weight predominantly through urinary-system disorders (N39), history-of-conditions codes (Z87), arthrosis (M19), recurrent admissions for nervous-system symptoms (R29) and falls (W19) — codes that biologically connect connective-tissue remodelling to musculoskeletal and neurogeriatric decline (**Extended Data Fig. 1**).

### Endotrophin compresses healthspan and tracks each frailty syndrome

The **Zenin 8-event healthspan composite** (first occurrence of HF, MI, COPD, stroke, dementia, diabetes, cancer or death) accumulated 11,564 events in the +Clinical risk set: per-SD HR 1.15 (1.12–1.17), *P* = 1.1 × 10⁻³³. The first-component breakdown showed that 48 % of composite events were cancers, 15 % T2D, 8.5 % deaths and 8.4 % COPD admissions, with the remainder distributed across stroke, HF, MI and dementia — recapitulating the spectrum of leading age-related conditions (**Extended Data Fig. 2**).

Each frailty-syndrome hospitalisation modelled separately also tracked endotrophin: per-SD HRs of **1.15 (1.07–1.23)** for hip fracture, **1.21 (1.16–1.25)** for falls, **1.17 (1.10–1.25)** for delirium and **1.31 (1.23–1.40)** for pressure ulcer (all +Clinical; **Table 1**). The pressure-ulcer signal was the largest of the four and survived the +Biomarkers adjustment essentially unchanged (HR 1.26).

### Disease-panel associations and an ICD-event substitute for cause-specific mortality

The pre-specified incident-disease panel, modelled with prevalent cases of the relevant condition removed at baseline, produced effect sizes consistent with the published 01.pv / 01.pv_hf endotrophin pipelines on the same dataset:

- **Incident T2D** 2,688 events; HR 1.09 (1.05–1.13), *P* = 7.6 × 10⁻⁶ (+Clinical).
- **Incident HF** 1,431 events; HR 1.27 (1.21–1.33), *P* = 5.4 × 10⁻²³.
- **Incident CKD** 1,715 events; HR **1.48 (1.42–1.55)**, *P* = 3.2 × 10⁻⁷⁶ — closely reproducing the published per-SD CKD HR of 1.53 from the 01.pv pipeline.
- **Incident dementia** 858 events; HR 1.08 (1.02–1.15), *P* = 0.009. Dementia was the only endpoint at which the +Biomarkers model attenuated the association out of nominal significance (HR 1.05, *P* = 0.09), consistent with NT-proBNP and CRP capturing a large share of vascular- and inflammation-mediated dementia risk.

We were unable to compute cause-specific (CVD or cancer) mortality because UK Biobank field 40001 (underlying cause of death, ICD-10) was not in the curated dataset. As a transparent surrogate, we report the **first incident I00–I99 event** ("any-CVD event") and the **first incident C00–D48 event** ("any-cancer event") from the HES record. Both retained nominal significance in the fully adjusted +Clinical model: HR 1.11 (1.09–1.14), *P* = 1.7 × 10⁻¹⁸ for any-CVD and HR 1.05 (1.03–1.08), *P* = 1.2 × 10⁻⁴ for any-cancer (**Table 1**, "ICD-event substitute" caveat column).

### Robustness across adjustment levels and multiple-testing control

Per-SD HRs from the three nested adjustment models — Base (age, sex, centre, time-to-Olink-processing, Batch); +Clinical (Base + BMI, eGFR, SBP, smoking, diabetes, hypertension, dyslipidaemia, HbA1c, HDL-C, total-C); +Biomarkers (+Clinical + CRP, NT-proBNP) — are summarised in the master forest (**Fig. 2**) and **Table 1**. All 13 endpoints retained nominal significance under +Clinical, and all but incident dementia retained significance under +Biomarkers. Benjamini–Hochberg control of the family-wise false discovery rate within each model (26 tests per model: 13 endpoints × {per-SD, p-trend}) left every per-SD association except dementia under +Biomarkers below *q* = 0.001 (**Table 1**, last two columns).

The dose-response across quartiles was monotonic for every endpoint (**Extended Data Figs. 3–15**), and a parallel trend test using the quartile-as-integer covariate confirmed the linearity of the log-hazard relationship (all p-trend ≤ 0.018 in the +Clinical model except for dementia at p-trend = 0.017).

---

## Methods

### Study population and exposure

We analysed 44,642 UK Biobank participants who had baseline plasma endotrophin (col6a3 NPX) measured by Olink Explore 3072 in the UKB Pharma Proteomics Project (UKB-PPP) and complete data for age, sex, assessment centre, Olink processing-batch identifier and time from sample collection to Olink processing. The data file (`curated_stats.tsv.gz`) contains 45,479 rows × 3,525 columns, including approximately 3,400 Olink protein NPX values, baseline biochemistry, ICD-10-linked Hospital Episode Statistics and mortality records.

The exposure was endotrophin standardised within the analytic cohort (mean 0.027 NPX, SD 0.42); Cox models report the hazard ratio for a 1-SD increment.

### Endpoint construction

**Pre-computed endpoints.** All-cause mortality (`time_to_anydeath`, `date_death`), incident heart failure (`incident_HF`, `time_to_first_HF_event_days`) and incident type-2 diabetes (`incident_T2D_case_control`, `time_to_T2Devent`) were taken from the curated dataset. Prevalent cases — defined as event-date ≤ 0 days from `date_attending_centre` — were removed from the corresponding risk set at analysis time.

**ICD-10-derived endpoints.** Each participant's lifetime ICD-10 diagnosis history was provided as a comma-separated string in the `ICD10` column, paired 1:1 with up to 259 dxdate fields (`dxdate_a0..dxdate_a258`). We exploded this representation into a long DataFrame indexed by participant identifier, extracting for every code its three-character prefix, the original (3-, 4- or 5-character) code, the diagnosis date, and the days from baseline. We then defined endpoint indicators by regular-expression matching on the normalised (upper-cased, dot-stripped) comma-prefixed code string:

| Endpoint | ICD-10 pattern |
|---|---|
| Heart failure (Zenin) | I50* |
| Myocardial infarction | I21*, I22* |
| COPD | J40–J44 |
| Stroke | I60–I69 |
| Dementia | F00–F03, G30* |
| T2D (Zenin component) | E11* |
| Cancer | C00–C97, D00–D48 |
| Hip fracture | S72* |
| Falls | W00–W19, R29.6 |
| Delirium | F05* |
| Pressure ulcer | L89* |
| Any-CVD event (mortality substitute) | I00–I99 |
| Any-cancer event (mortality substitute) | C00–C97, D00–D48 |
| Chronic kidney disease | N18* |

For each label, the per-participant first-event date was the earliest matching `dxdate`. Prevalent flag = 1 if first event was at or before `date_attending_centre`; incident flag = 1 if first event occurred after baseline within the follow-up window. Time-to-event used the event date for incident cases and `rowwise_max(time_to_lostfollow, time_to_anydeath, time_to_latest_icd10, max_follow_up)` otherwise, mirroring the censoring recipe established in the 01.pv CKD pipeline. Time was expressed in years (days / 365.25).

### Hospital Frailty Risk Score

We encoded the 109 ICD-10 codes and weights from Supplementary Table 2 of Gilbert *et al.* (Lancet 2018; PMID 29706364). Weights were assigned at the three-character level; longer (4–5-character) child codes inherited the parent weight via prefix match. The **incident-frailty endpoint** was defined as the time at which the cumulative HFRS first reached 5.0 (the published intermediate-risk threshold). Within each participant we walked distinct dxdates in chronological order, maintained a running set of three-character HFRS codes already seen, and recorded the first dxdate at which the running sum of weights crossed 5.0.

Because the curated dataset does not carry per-admission grouping (`hesin_diag` was not extracted), distinct dxdates substitute for admissions; codes recorded on the same date count as one accumulation step. This is documented as the closest faithful Gilbert implementation available against the present data and mildly underestimates HFRS relative to admission-level scoring.

The **baseline cumulative HFRS** used the same algorithm but restricted to codes recorded in the two years preceding the assessment-centre visit.

### Cohort flow

Beginning from 45,479 rows in `curated_stats.tsv.gz`, we excluded 837 with missing `col6a3`, 0 with missing core covariates, and 0 with invalid `date_attending_centre` or `max_follow_up`, retaining 44,642 participants in the common analytic cohort. Per-endpoint risk sets removed prevalent cases of that specific condition, so the analytic *N* varied between 37,117 (any-CVD event) and 44,642 (mortality and conditions with very low baseline prevalence). In the +Clinical and +Biomarkers models the effective *N* was further reduced by case-wise deletion of rows missing any covariate to 28,524–34,248 and 27,912–33,527, respectively. Per-endpoint cohort flow is reported in `results/tables/section2_cohort_flow.csv`.

### Statistical models

**Adjustment levels.** Three nested Cox models were prespecified:

- **Base**: age, sex, assessment centre, time-to-Olink-processing, Olink Batch.
- **+Clinical**: Base + BMI, eGFR (CKD-EPI 2021 from serum creatinine), systolic blood pressure, ever-smoker, diabetes, hypertension, dyslipidaemia, HbA1c, HDL cholesterol, total cholesterol.
- **+Biomarkers**: +Clinical + CRP and NT-proBNP.

Cox models were fitted with `lifelines.CoxPHFitter` with a small ridge penalty (`penalizer=0.01`) to stabilise the many-level fixed effects (assessment centre, Olink Batch); this matched the convention in the 01.pv_hf pipeline.

**Per-SD and quartile reporting.** Each endpoint was fitted both with the standardised endotrophin variable as a continuous covariate (per-SD HR) and with quartile dummies (Q2–Q4 vs Q1, with the reference Q1). A linear-trend p-value was obtained by re-fitting the quartile as an integer 1–4. The concordance index was reported on the +Clinical model only.

**Multiple testing.** Benjamini–Hochberg false-discovery-rate adjustment was applied within each adjustment model across the 26 family tests (13 endpoints × {per-SD, p-trend}); both raw *P* and BH-q are reported in `results/tables/section2_per_sd.csv`.

**Cross-sectional models.** Section 1 used ordinary least-squares regression (continuous baseline HFRS, biological-age composite), logistic regression (prevalent frailty-syndrome flags), and Poisson regression (multimorbidity counts), each with the same three nested covariate sets (`statsmodels` GLM family with appropriate link).

### Software and reproducibility

The pipeline is implemented in Python 3.10 (`pandas` 2.3, `numpy` 2.2, `lifelines` 0.30, `statsmodels` 0.14, `matplotlib` 3.x). Heavy fits ran on the Broad UGER (Univa Grid Engine 8.5.5) cluster on RHEL 8 nodes (24 GB RAM, 2 cores; total wall time 7 minutes for Section 2 and < 1 minute for Section 1 on a single host). The full repository, including the qsub wrappers, ICD-10 long-format parser, the 109-code HFRS table, and reusable Cox / KM / forest-plot helpers imported from `/humgen/diabetes2/users/satoshi/misc/01.pv/ckd_shared.py`, is available at <https://github.com/satoshi-yoshiji/PV-aging>.

A `_selftest_ckd()` helper in `aging_shared.py` rebuilds the CKD endpoint via the new ICD parser and reproduces the published 01.pv CKD per-SD HR within ΔHR < 0.05 in the +Clinical model.

### Data availability

Individual-level UK Biobank data are not redistributable. Access can be requested from UK Biobank.

### Limitations

Three are explicit in the design and worth re-stating. First, **Aim-1 phenotypes** that the original protocol envisioned (Williams 49-item Frailty Index, Fried physical-frailty phenotype, grip strength, gait speed and hand-grip asymmetry) require UK Biobank baseline-assessment fields (f.46/47, f.924, etc.) that were not present in the curated dataset; the cross-sectional Section 1 reports HFRS-, ICD-, self-report- and lab-based proxies in their place. A no-op curation hook (`extend_with_baseline_phenotypes`) is shipped so that those analyses fire automatically once the user supplies a UKB extract by `eid`. Second, **cause-specific mortality** (CVD, cancer) would require UKB f.40001 underlying cause-of-death ICD codes, which are absent from the dataset; we report transparent ICD-event substitutes with a `caveat` column flagging the difference. Third, **HFRS admission grouping** uses distinct dxdates as the unit of accumulation in the absence of `hesin_diag`; the resulting score mildly underestimates the published Gilbert algorithm.

### References

1. Gilbert, T. *et al.* Development and validation of a Hospital Frailty Risk Score focusing on older people in acute care settings using electronic hospital records: an observational study. *Lancet* **391**, 1775–1782 (2018).
2. Zenin, A. *et al.* Identification of 12 genetic loci associated with human healthspan. *Commun. Biol.* **2**, 41 (2019).
3. Inker, L. A. *et al.* New creatinine- and cystatin C–based equations to estimate GFR without race. *N. Engl. J. Med.* **385**, 1737–1749 (2021).
4. Sun, B. B. *et al.* Plasma proteomic associations with genetics and health in the UK Biobank. *Nature* **622**, 329–338 (2023).
5. Bycroft, C. *et al.* The UK Biobank resource with deep phenotyping and genomic data. *Nature* **562**, 203–209 (2018).
6. Williams, D. M. *et al.* Frailty and healthcare costs in the UK Biobank: a 49-item index. *J. Gerontol. A. Biol. Sci. Med. Sci.* **74**, 582–589 (2019).
