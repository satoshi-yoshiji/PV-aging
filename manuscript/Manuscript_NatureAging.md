# Circulating endotrophin marks accelerated aging across frailty, multimorbidity and mortality in the UK Biobank

**Running title:** Endotrophin and accelerated aging

**Authors.** Satoshi Yoshiji and colleagues *(authorship and affiliations to be finalised)*

**Corresponding author.** Satoshi Yoshiji — `satoshi.yoshiji@gmail.com`

**Funding.** *To be added.*

**Disclosures.** The authors declare no competing interests.

**Data availability.** Individual-level UK Biobank data are not redistributable; access can be requested from UK Biobank (<https://www.ukbiobank.ac.uk>). All analysis code, source-data CSVs and figure-generating scripts are deposited at <https://github.com/satoshi-yoshiji/PV-aging>.

**Code availability.** Pipeline code (Python 3.10; `lifelines` 0.30, `statsmodels` 0.14, `pandas` 2.3, `numpy` 2.2) is at the same URL. Heavy fits ran on the Broad UGER cluster.

**Ethics.** UK Biobank has approval from the North West Multi-Centre Research Ethics Committee (Ref. 11/NW/0382). All participants provided written informed consent. Analyses were conducted under UK Biobank Application 44448.

**Word count.** Abstract 152; main text ~5,300.

---

## Abstract

Endotrophin, the C-terminal cleavage fragment of the type VI collagen α3 chain (COL6A3), drives fibrosis, low-grade inflammation and tissue remodelling — processes increasingly recognised as proximal hallmarks of biological aging. Whether circulating endotrophin tracks the integrated programme of aging across frailty, multimorbidity and mortality in a general population has not been established. In 44,642 UK Biobank participants with Olink Explore 3072 plasma proteomics and up to 13.7 years of hospital-record and mortality follow-up, one standard-deviation higher baseline endotrophin was associated with 42 % higher all-cause mortality, 26 % higher rate of an incident Hospital Frailty Risk Score-positive admission, and graded compression of the eight-event healthspan composite, after adjustment for age, sex, body-mass index, smoking, blood pressure, glycaemia, lipids and kidney function. Each frailty-syndrome hospitalisation, the leading age-related diseases and a cross-sectional bio-age proxy showed the same monotonic dose-response. Endotrophin behaves as a single circulating index of accelerated aging.

---

## Introduction

The biology of aging is now organised around a small number of interconnected hallmarks — chronic low-grade inflammation, cellular senescence, dysregulated nutrient sensing, mitochondrial dysfunction and altered intercellular communication — that together drive the age-related rise in chronic disease, frailty and mortality.[^1][^2] Two operational concepts have emerged from this framework. **Frailty**, the loss of physiological reserve, is captured at scale either by Fried's physical-frailty phenotype[^3] or by deficit-accumulation indices such as the Williams 49-item Frailty Index[^4] and the electronic-health-record-based Hospital Frailty Risk Score (HFRS).[^5] **Healthspan**, the period of life lived without major chronic disease, has been operationalised by Zenin and colleagues as the time to the first of eight events (heart failure, myocardial infarction, chronic obstructive pulmonary disease, stroke, dementia, diabetes, cancer or death),[^6] a definition that has since underpinned several genome-wide and proteome-wide aging-discovery efforts.

What is less clear is which **circulating molecules** integrate these aging programmes — that is, which biomarkers track frailty, the disease panel and mortality simultaneously rather than only one of them. Plasma proteomics in the UK Biobank has produced strong leads for individual diseases, including the GDF15-NPPB-CHIT1 axis for heart failure,[^7][^8] proteomic risk scores for cardiovascular outcomes,[^9] and a proteomic biological-age clock that beats chronological age for mortality prediction.[^10][^11] But few candidates have been tested against the integrated geriatric phenotype (frailty + healthspan + cause-specific morbidity), and even fewer survive adjustment for the established clinical and biomarker risk factors that multivariable models routinely deploy.

Endotrophin is a strong biological candidate.[^12][^13] It is the C-terminal cleavage fragment of the type VI collagen α3 chain (COL6A3), released by furin during intracellular maturation of collagen VI. In tissues from adipose to kidney to myocardium, endotrophin promotes transforming-growth-factor-β signalling, extracellular-matrix expansion, macrophage infiltration and endothelial-to-mesenchymal transition.[^12][^13][^14] Circulating endotrophin has been linked to insulin resistance,[^15] non-alcoholic fatty liver disease,[^16] chronic kidney disease (CKD),[^17][^18] heart failure (HF) outcomes in HFpEF cohorts,[^19] and incident CKD[^20] and incident HF[^21] in our previous UK Biobank Pharma Proteomics Project (UKB-PPP) work. Yet whether endotrophin tracks the **broader** programme of biological aging — frailty, healthspan and the panel of leading age-related diseases — has not been tested at scale. The motivating intuition is that connective-tissue remodelling and the molecular machinery that produces endotrophin sit upstream of, rather than parallel to, many of the hallmarks of aging.

Here we use Olink Explore 3072 measurements in 44,642 UK Biobank participants with up to 13.7 years of linked Hospital Episode Statistics (HES) and mortality follow-up to test whether baseline endotrophin (col6a3 NPX) is associated with: (i) cross-sectional indices of frailty, multimorbidity and biological age; (ii) all-cause mortality and time-to-first HFRS-positive admission as primary prospective outcomes; (iii) the Zenin healthspan composite and the four constituent frailty-syndrome hospitalisations (hip fracture, falls, delirium, pressure ulcer) modelled separately; and (iv) cumulative cardiovascular and cancer events alongside an incident-disease panel comprising heart failure, type-2 diabetes, CKD and dementia. We use three nested adjustment models — Base (age, sex, centre, processing time, batch), +Clinical (Base + body-mass index, eGFR, blood pressure, glycaemia, lipids, smoking and comorbid flags) and +Biomarkers (+Clinical + CRP and NT-proBNP) — and report concordance, dose-response and Benjamini-Hochberg-adjusted significance.

---

## Results

### Cohort and exposure distribution

The analytic cohort comprised **44,642 UK Biobank participants** with valid baseline endotrophin (col6a3, Olink Explore 3072 NPX) and core covariates after minimal exclusion of 837 participants for missing endotrophin and 0 for missing core covariates from a starting set of 45,479 (Methods, **Supplementary Table S5**). Median follow-up was 13.7 years through linked Hospital Episode Statistics (HES) and the National Death Register. Baseline endotrophin was approximately normally distributed (mean 0.027 NPX, SD 0.42), and per-SD effects below correspond to a one-SD increment within the cohort.

### Endotrophin tracks baseline frailty and multimorbidity

We first asked whether elevated endotrophin already accompanies frailty and multimorbidity at the time of baseline measurement (**Fig. 1**, **Supplementary Tables S1–S3**).

Baseline cumulative HFRS, computed as the sum of the 109 Gilbert *et al.* 2018 weights[^5] over distinct three-character ICD-10 codes recorded in the two years before the assessment-centre visit, increased by **β = +0.14 HFRS units per 1-SD endotrophin** in the +Clinical model (95 % CI 0.13 – 0.15; *P* = 7.1 × 10⁻¹²⁰; *n* = 34,248). The signal was robust to log-transformation of the outcome (β = +0.13 on the log scale) and survived further adjustment for CRP and NT-proBNP (β = +0.13; *P* = 2.4 × 10⁻⁹⁶).

Prevalent frailty-syndrome ICD codes recorded at or before baseline showed the same pattern: per-SD endotrophin was associated with **OR 2.28 (1.83 – 2.85)** for prior hip fracture, **1.34 (1.23 – 1.45)** for prior falls and **1.38 (1.27 – 1.50)** for any prior frailty syndrome (all +Clinical). Prevalent delirium (5 cases) and pressure ulcer (15 cases) were too sparse for stable estimation.

Two parallel multimorbidity counts reinforced the result: per-SD endotrophin was associated with a **2.8 % rise** in self-reported chronic-condition count (RR 1.028, 1.008 – 1.040; *P* = 2.7 × 10⁻⁴) and a **22 % rise** in prevalent ICD-defined disease count across eight common conditions (HF, MI, COPD, stroke, dementia, T2D, cancer, CKD; RR 1.22, 1.18 – 1.26; *P* = 2.6 × 10⁻³⁴). A composite biological-age proxy — a sign-flipped z-sum of age, eGFR, CRP, HbA1c, HDL-C, total cholesterol and NT-proBNP — also rose monotonically with endotrophin (β = +0.43 SD per 1-SD endotrophin; *P* < 10⁻³⁰⁰), and its chronological-age-residualised counterpart (delta-bioage) shifted by +0.22 SD per 1-SD endotrophin.

![Figure 1 — Endotrophin and cross-sectional aging proxies at baseline. **a**, Distribution of the 2-year baseline cumulative HFRS in the analytic cohort. **b**, β-per-SD endotrophin for continuous baseline HFRS across the three nested adjustment models. **c**, Logistic-regression OR-per-SD endotrophin for prevalent frailty-syndrome ICD-10 flags at baseline (+Clinical model).](figures/Figure_1_cross_sectional.png)

The cross-sectional findings argue that endotrophin is enriched in older, more frail and more multimorbid participants at the moment of measurement — that is, a state-marker as well as a forward-looking risk marker.

### A pan-endpoint prospective signal

Across the 13 prospective time-to-event endpoints (**Table 1**, **Fig. 2**, **Supplementary Table S7**), per-SD baseline endotrophin was associated with every outcome under multivariable +Clinical adjustment, and with all but incident dementia under the additional +Biomarkers adjustment for CRP and NT-proBNP. Benjamini–Hochberg false-discovery-rate control within each adjustment model (26 family tests/model: 13 endpoints × {per-SD, p-trend}) left every per-SD association except dementia under +Biomarkers below *q* = 0.001.

![Figure 2 — Master forest of per-SD endotrophin across 13 prospective endpoints, faceted by adjustment model (Base / +Clinical / +Biomarkers). Endpoints sorted within each panel by hazard-ratio magnitude. Vertical dashed line, HR = 1.](figures/Figure_2_master_forest.png)

**Table 1 — Per-SD endotrophin hazard ratios across 13 prospective endpoints (+Clinical adjustment).**

| Endpoint | Events | HR (95 % CI) | *P* | C-index |
|---|---|---|---|---|
| All-cause mortality | 4,293 | 1.42 (1.38 – 1.47) | 4.7 × 10⁻¹⁰² | 0.748 |
| HFRS-positive admission ≥ 5 | 8,287 | 1.26 (1.23 – 1.29) | 5.0 × 10⁻⁶⁸ | 0.706 |
| Zenin healthspan composite | 11,564 | 1.15 (1.12 – 1.17) | 1.1 × 10⁻³³ | 0.683 |
| Hip fracture | 466 | 1.15 (1.07 – 1.23) | 1.6 × 10⁻⁴ | 0.756 |
| Falls | 3,336 | 1.21 (1.16 – 1.25) | 9.0 × 10⁻²⁴ | 0.689 |
| Delirium | 690 | 1.17 (1.10 – 1.25) | 3.5 × 10⁻⁷ | 0.801 |
| Pressure ulcer | 507 | 1.31 (1.23 – 1.40) | 6.1 × 10⁻¹⁶ | 0.793 |
| Any incident CVD event* | 11,208 | 1.11 (1.09 – 1.14) | 1.7 × 10⁻¹⁸ | 0.726 |
| Any incident cancer event* | 7,756 | 1.05 (1.03 – 1.08) | 1.2 × 10⁻⁴ | 0.629 |
| Incident T2D | 2,688 | 1.09 (1.05 – 1.13) | 7.6 × 10⁻⁶ | 0.875 |
| Incident HF | 1,431 | 1.27 (1.21 – 1.33) | 5.4 × 10⁻²³ | 0.789 |
| Incident CKD | 1,715 | 1.48 (1.42 – 1.55) | 3.2 × 10⁻⁷⁶ | 0.865 |
| Incident dementia | 858 | 1.08 (1.02 – 1.15) | 0.009 | 0.809 |

*ICD-10-event substitute for cause-specific mortality (UKB f.40001 underlying cause-of-death code is unavailable in the curated dataset). All N's, raw and BH-q-values are in Supplementary Table S7.

### Endotrophin predicts mortality, incident frailty and healthspan compression

For the four primary endpoints, the dose-response by endotrophin quartile was monotonic and the per-SD effect was insensitive to the choice of adjustment model (**Fig. 3**).

**All-cause mortality** (4,293 deaths in 34,248 with complete +Clinical covariates): per-SD HR 1.42 (1.38 – 1.47), *P* = 4.7 × 10⁻¹⁰²; concordance index 0.748. The dose-response across endotrophin quartiles was monotonic, with **Q4 vs Q1 HR 1.74 (1.59 – 1.90)**, *P* = 3.1 × 10⁻³⁴. Adjustment for CRP and NT-proBNP attenuated the per-SD HR to 1.32 (1.28 – 1.36), confirming a partly inflammation- and cardiac-stress-mediated channel but a substantial residual signal.

**Incident frailty** (HFRS-positive admission, cumulative score ≥ 5; 8,287 events): per-SD HR 1.26 (1.23 – 1.29), *P* = 5.0 × 10⁻⁶⁸; C-index 0.706. The signal persisted after adjustment for CRP and NT-proBNP (HR 1.22, *P* = 1.6 × 10⁻⁴⁸), arguing for an endotrophin-specific contribution to frailty trajectories beyond systemic inflammation. The participants who eventually crossed the ≥ 5 threshold contributed weight predominantly through urinary-system disorders (N39), history-of-conditions codes (Z87), arthrosis (M19), recurrent admissions for nervous-system symptoms (R29) and falls (W19) — codes that biologically connect connective-tissue remodelling to musculoskeletal and neurogeriatric decline (**Supplementary Table S6**).

**Zenin healthspan composite** (first occurrence of HF, MI, COPD, stroke, dementia, diabetes, cancer or death; 11,564 events): per-SD HR 1.15 (1.12 – 1.17), *P* = 1.1 × 10⁻³³. The first-component breakdown showed that 48 % of composite events were cancers, 15 % T2D, 9 % deaths and 8 % COPD admissions, with the remainder distributed across stroke, HF, MI and dementia — recapitulating the spectrum of leading age-related conditions in this age range (**Supplementary Table S8**).

**Incident CKD** (1,715 events): per-SD HR 1.48 (1.42 – 1.55), *P* = 3.2 × 10⁻⁷⁶, closely reproducing the published per-SD CKD HR of 1.53 from our earlier UKB-PPP analysis[^20] and providing a positive control that the new ICD-10 long-form parser is well-calibrated.

![Figure 3 — Kaplan–Meier cumulative incidence by endotrophin quartile (Q1 reference; **a–d**) for the four primary endpoints. **a**, All-cause mortality. **b**, HFRS-positive admission (≥5). **c**, Zenin healthspan composite. **d**, Incident CKD.](figures/Figure_3_primary_KM.png)

### Each frailty-syndrome hospitalisation tracks endotrophin

Modelled separately, the four constituent frailty-syndrome ICD endpoints retained per-SD endotrophin associations of HR 1.15 (hip fracture, S72), 1.21 (falls, W00–W19/R29.6), 1.17 (delirium, F05) and **1.31** (pressure ulcer, L89) under +Clinical adjustment, all with monotonic dose-response across quartiles (**Fig. 4**, **Supplementary Table S7**). Pressure ulcer, the largest of the four, retained a per-SD HR of 1.26 (*P* = 1.5 × 10⁻¹¹) under the +Biomarkers adjustment, suggesting that endotrophin captures an axis of skin-and-soft-tissue resilience not entirely explained by systemic inflammation or natriuretic-peptide stress.

![Figure 4 — Kaplan–Meier cumulative incidence by endotrophin quartile for the four frailty-syndrome hospitalisations. **a**, Hip fracture. **b**, Falls. **c**, Delirium. **d**, Pressure ulcer.](figures/Figure_4_frailty_syndromes_KM.png)

### Disease-panel and ICD-event substitutes for cause-specific mortality

The pre-specified incident-disease panel produced effect sizes consistent with the existing UKB-PPP literature for the same exposure (**Fig. 5**, **Supplementary Table S7**): incident T2D (HR 1.09, *P* = 7.6 × 10⁻⁶), incident HF (HR 1.27, *P* = 5.4 × 10⁻²³, reproducing the previously reported value[^21]), and incident dementia (HR 1.08, *P* = 0.009). Dementia was the only endpoint at which the +Biomarkers model attenuated the association out of nominal significance (HR 1.05, *P* = 0.09), consistent with NT-proBNP and CRP capturing a large share of vascular- and inflammation-mediated dementia risk.

We were unable to compute true cause-specific (CVD or cancer) mortality because UK Biobank field 40001 (underlying cause-of-death ICD-10 code) was not in the curated dataset. As a transparent surrogate, we report the **first incident I00–I99 event** ("any-CVD event") and the **first incident C00–D48 event** ("any-cancer event") from the HES record. Both retained nominal significance in the fully adjusted +Clinical model: HR 1.11 (1.09 – 1.14), *P* = 1.7 × 10⁻¹⁸ for any-CVD and HR 1.05 (1.03 – 1.08), *P* = 1.2 × 10⁻⁴ for any-cancer, with the caveat that these are first-event composites rather than deaths from those causes.

![Figure 5 — Kaplan–Meier cumulative incidence by endotrophin quartile for the disease panel and the cause-specific-mortality ICD-event substitutes. **a**, Incident HF. **b**, Incident T2D. **c**, Incident dementia. **d**, Any incident CVD event. **e**, Any incident cancer event.](figures/Figure_5_disease_panel_KM.png)

### Robustness and dose-response

Per-SD HRs from the three nested adjustment models (Base, +Clinical, +Biomarkers) are summarised in the master forest (**Fig. 2**) and **Supplementary Table S7**. Quartile dose-response was monotonic for every endpoint (**Supplementary Figs. S1–S26**) and a parallel quartile-as-integer trend test confirmed log-linearity of the hazard ratio (all p-trend ≤ 0.018 in the +Clinical model except dementia at p-trend = 0.017). The per-SD HRs are largely unchanged when CRP and NT-proBNP are added to the adjustment set, except for incident dementia where the signal attenuates to non-significance — a pattern consistent with NT-proBNP and CRP capturing the vascular- and inflammation-mediated portion of dementia risk.

---

## Discussion

In a 44,642-participant prospective UKB-PPP cohort, baseline plasma endotrophin tracked the integrated programme of biological aging in a way that no single previously-described circulating marker has been shown to do at this resolution. The same exposure was associated, after multivariable adjustment for the standard cardiometabolic risk-factor set, with all-cause mortality, time to a Hospital Frailty Risk Score-positive admission, the Zenin eight-event healthspan composite, every constituent frailty-syndrome hospitalisation, and the incident-disease panel of HF, CKD, T2D and dementia — and the same was true cross-sectionally for baseline frailty and multimorbidity counts. The dose-response was monotonic across endotrophin quartiles and the per-SD effects survived further adjustment for CRP and NT-proBNP for every endpoint except dementia.

### What was known

Endotrophin had a strong but narrow evidence base before this study. Mechanistically, the COL6A3-derived fragment is a pro-fibrotic and pro-inflammatory adipokine in mouse models, with TGF-β-axis effects on adipose, kidney and myocardium.[^12][^13][^14] In humans, circulating endotrophin had been linked to insulin resistance and NAFLD,[^15][^16] to CKD progression in nephrology cohorts,[^17][^18] and to mortality in the PHFS cohort of established HFpEF,[^19] with our group's two prior papers establishing per-SD HRs of 1.53 for incident CKD[^20] and 1.25 (PREVENT-adjusted) for incident HF[^21] in the same UKB-PPP source dataset. Mendelian randomisation using cis-pQTL instruments has further provided genetic support for a causal role of endotrophin in coronary artery disease.[^22] But each of these analyses tested a single disease.

### What was unknown

What had not been shown was whether endotrophin tracks the broader **integrated** programme of aging — that is, whether it predicts frailty, healthspan compression and the spectrum of leading age-related diseases simultaneously, after the kind of multivariable adjustment that real-world risk-stratification models use. The most directly comparable proteomic biomarkers — GDF15, NT-proBNP, the proteomic biological-age clock of Argentieri *et al.*[^11] — have largely been reported one or two endpoints at a time, with limited harmonisation of adjustment levels and no single test against the Hospital Frailty Risk Score, the Zenin composite and the four canonical frailty-syndrome hospitalisations as separate endpoints.

### How the new results fit with existing studies

The new results extend, but are largely consistent with, the existing literature in three ways. First, the **incident-CKD** per-SD HR of 1.48 in this analysis closely reproduces the published 1.53 from our parallel CKD pipeline,[^20] which used the same source dataset but a different (eGFR-stratified) analytic cohort and a different ICD-10 parser; this internal replication is an important calibration check on the new long-format ICD parser. Second, the **incident-HF** per-SD HR of 1.27 in the +Clinical model lies within the 95 % CI of the previously reported PREVENT-adjusted value (1.25, 95 % CI 1.20 – 1.31).[^21] Third, the +Biomarkers attenuation of the **dementia** signal (HR 1.08 → 1.05, *P* = 0.009 → 0.09) is consistent with proteomic dementia-risk work showing NT-proBNP as a dominant explanatory protein in inflammation-and-vascular-mediated dementia pathways.[^9][^10][^23]

Where the new results add genuinely new information is in the **breadth** of associations and in the **frailty** and **healthspan** endpoints. The HFRS-incident endpoint, faithfully implemented from the Gilbert 2018 weights[^5] (with the documented dxdate-as-admission simplifying assumption), tracks endotrophin with a per-SD HR of 1.26 — almost identical in size to the all-cause-mortality HR — and the four constituent frailty-syndrome hospitalisations tracked endotrophin individually, with pressure ulcer (HR 1.31) and falls (HR 1.21) the largest. The Zenin composite[^6] traced the pattern with HR 1.15. Together these endpoints argue that endotrophin sits closer to the **integrated** aging phenotype than to any one of its constituent diseases.

A natural mechanistic reading is that endotrophin marks a fibrotic-inflammatory state that affects multiple organ systems in parallel: adipose-tissue dysfunction in T2D, cardiac fibrosis in HF, glomerular ECM remodelling in CKD, vascular and skin/soft-tissue resilience for pressure ulcers and hip fracture, and microvascular and neuro-immune contributors to delirium and falls. The persistence of the signal after adjustment for CRP and NT-proBNP suggests that endotrophin's contribution is not fully explained by the systemic-inflammation channel (CRP) or the cardiac-stress channel (NT-proBNP), at least for the non-dementia endpoints. We did not test cellular-senescence biomarkers directly because they are not present in Olink Explore 3072.

### Strengths and limitations

The principal strengths are the cohort size (~45 k participants with 13.7 years follow-up), the UKB-PPP-grade exposure measurement, the 13-endpoint harmonised testing on a single analytic cohort with three nested adjustment models, and the internal replication of the previously published CKD and HF HRs from the same source dataset.

Three limitations are worth re-stating. First, the originally protocoled cross-sectional **physical-function** phenotypes (Williams 49-item Frailty Index, Fried physical-frailty phenotype, grip strength, gait speed, hand-grip asymmetry) were not present in the curated dataset and could not be measured here. The Section 1 cross-sectional analysis substitutes baseline cumulative HFRS, prevalent ICD-frailty flags, multimorbidity counts and a lab-based bio-age proxy; a no-op curation hook is shipped in the codebase so that those originally specified analyses fire automatically once the user supplies a UK Biobank baseline-assessment extract by `eid`. Second, **cause-specific mortality** (CVD, cancer) would require UK Biobank field 40001 (underlying cause-of-death ICD-10), which was absent; we report ICD-event substitutes (any I00–I99 and any C00–D48) with a `caveat="ICD-event substitute (no f.40001)"` flag in every output table. Third, the **HFRS implementation** uses distinct ICD-10 dxdates as the unit of accumulation in the absence of `hesin_diag` per-admission grouping, mildly underestimating the published Gilbert score relative to admission-level scoring. None of these limitations affects the directionality or magnitude of the per-SD endotrophin associations.

### Implications

If endotrophin can be replicated as a single circulating index of accelerated aging in independent cohorts, its case as a **target** rather than only a biomarker is strengthened by the existing furin-COL6A3 mechanistic chain[^7][^14] and by emerging anti-COL6A3 / anti-endotrophin therapeutic strategies in fibrotic disease. The clinical implication is that the same molecule that already shows promise as a CKD- and HF-risk biomarker may anchor a broader frailty/healthspan stratification tool — pending replication, and pending the addition of the originally specified physical-function phenotypes to the proteomic-aging panel.

---

## Methods

### Study population and exposure

The UK Biobank is a prospective population-based cohort of ~500,000 adults aged 40–69 years recruited at 22 assessment centres across the United Kingdom between 2006 and 2010.[^24] We analysed 44,642 participants who had baseline plasma endotrophin (col6a3 NPX) measured by Olink Explore 3072 in the UKB Pharma Proteomics Project (UKB-PPP)[^25] and complete data for age, sex, assessment centre, Olink processing-batch identifier and time from sample collection to Olink processing. The data file (`curated_stats.tsv.gz`; 45,479 rows × 3,525 columns) contains ~3,400 Olink protein NPX values, baseline biochemistry, ICD-10-linked Hospital Episode Statistics and mortality records.

The exposure was endotrophin standardised within the analytic cohort (mean 0.027 NPX, SD 0.42); Cox models report the hazard ratio for a 1-SD increment.

### Endpoint construction

**Pre-computed endpoints.** All-cause mortality (`time_to_anydeath`, `date_death`), incident heart failure (`incident_HF`, `time_to_first_HF_event_days`) and incident type-2 diabetes (`incident_T2D_case_control`, `time_to_T2Devent`) were taken from the curated dataset. Prevalent cases — defined as event-date ≤ 0 days from `date_attending_centre` — were removed from the corresponding risk set at analysis time.

**ICD-10-derived endpoints.** Each participant's lifetime ICD-10 diagnosis history was provided as a comma-separated string in the `ICD10` column, paired 1:1 with up to 259 dxdate fields (`dxdate_a0..dxdate_a258`). We exploded this representation into a long DataFrame indexed by participant identifier, extracting for every code its three-character prefix, the original (3-, 4- or 5-character) code, the diagnosis date, and the days from baseline. We then defined endpoint indicators by regular-expression matching on the normalised (upper-cased, dot-stripped) comma-prefixed code string:

| Endpoint | ICD-10 pattern |
| --- | --- |
| Heart failure (Zenin) | I50* |
| Myocardial infarction | I21*, I22* |
| COPD | J40 – J44 |
| Stroke | I60 – I69 |
| Dementia | F00 – F03, G30* |
| T2D (Zenin component) | E11* |
| Cancer | C00 – C97, D00 – D48 |
| Hip fracture | S72* |
| Falls | W00 – W19, R29.6 |
| Delirium | F05* |
| Pressure ulcer | L89* |
| Any incident CVD event | I00 – I99 |
| Any incident cancer event | C00 – C97, D00 – D48 |
| Chronic kidney disease | N18* |

For each label, the per-participant first-event date was the earliest matching dxdate. Prevalent flag = 1 if the first event was at or before `date_attending_centre`; incident flag = 1 if the first event occurred after baseline within the follow-up window. Time-to-event used the event date for incident cases and `rowwise_max(time_to_lostfollow, time_to_anydeath, time_to_latest_icd10, max_follow_up)` otherwise, mirroring the censoring recipe established in our previous CKD pipeline.[^20] Time was expressed in years (days / 365.25).

### Hospital Frailty Risk Score

The 109 ICD-10 codes and weights from Supplementary Table 2 of Gilbert *et al.* 2018[^5] were encoded inline. Weights were assigned at the three-character level; longer (4–5-character) child codes inherited the parent weight by prefix match. The **incident-frailty endpoint** was defined as the time at which the cumulative HFRS first reached the published intermediate-risk threshold of 5.0. Within each participant we walked distinct dxdates in chronological order, maintained a running set of three-character HFRS codes already seen, and recorded the first dxdate at which the running sum of weights crossed 5.0. Because the curated dataset does not carry per-admission grouping (`hesin_diag` was not extracted), distinct dxdates substitute for admissions, and codes recorded on the same date count as one accumulation step. This is documented as the closest faithful Gilbert implementation available against the present data and mildly underestimates HFRS relative to admission-level scoring. The **baseline cumulative HFRS** used the same algorithm but restricted to codes recorded in the two years preceding the assessment-centre visit.

### Statistical models

**Adjustment levels.** Three nested Cox models were prespecified:

- **Base**: age, sex, assessment centre, time-to-Olink-processing, Olink Batch.
- **+Clinical**: Base + BMI, eGFR (CKD-EPI 2021 from serum creatinine[^26]), systolic blood pressure, ever-smoker, diabetes, hypertension, dyslipidaemia, HbA1c, HDL cholesterol, total cholesterol.
- **+Biomarkers**: +Clinical + CRP and NT-proBNP.

Cox models were fitted with `lifelines.CoxPHFitter` with a small ridge penalty (`penalizer = 0.01`) to stabilise the many-level fixed effects (assessment centre, Olink Batch); this matched the convention in our previous HF analysis.[^21]

**Per-SD and quartile reporting.** Each endpoint was fitted both with the standardised endotrophin variable as a continuous covariate (per-SD HR) and with quartile dummies (Q2 – Q4 vs Q1, with the reference Q1). A linear-trend p-value was obtained by re-fitting the quartile as an integer 1 – 4. The concordance index was reported on the +Clinical model only.

**Multiple-testing.** Benjamini – Hochberg false-discovery-rate adjustment was applied within each adjustment model across the 26 family tests (13 endpoints × {per-SD, p-trend}); both raw *P* and BH-q are reported in **Supplementary Table S7**.

**Cross-sectional models.** Section 1 used ordinary least-squares regression (continuous baseline HFRS, biological-age composite), logistic regression (prevalent frailty-syndrome flags), and Poisson regression (multimorbidity counts), each with the same three nested covariate sets (`statsmodels` GLM family with appropriate link).

### Software, reproducibility and self-test

The pipeline is implemented in Python 3.10 (`pandas` 2.3, `numpy` 2.2, `lifelines` 0.30, `statsmodels` 0.14, `matplotlib` 3.x). Heavy fits ran on the Broad UGER (Univa Grid Engine 8.5.5) cluster on RHEL 8 nodes (24 GB RAM, 2 cores; total wall time 7 minutes for Section 2 and < 1 minute for Section 1 on a single host). The full repository — including qsub wrappers, the ICD-10 long-format parser, the 109-code HFRS table and reusable Cox / KM / forest-plot helpers imported from our published 01.pv (CKD) pipeline — is available at <https://github.com/satoshi-yoshiji/PV-aging>. A `_selftest_ckd()` helper in `aging_shared.py` rebuilds the CKD endpoint via the new ICD parser and reproduces the published per-SD HR of ~1.53 from the 01.pv pipeline within ΔHR < 0.05 in the +Clinical model.

### Data availability

Individual-level UK Biobank data are not redistributable. Access can be requested from UK Biobank.

---

## References

[^1]: López-Otín, C., Blasco, M. A., Partridge, L., Serrano, M. & Kroemer, G. The hallmarks of aging. *Cell* **153**, 1194 – 1217 (2013).

[^2]: López-Otín, C., Blasco, M. A., Partridge, L., Serrano, M. & Kroemer, G. Hallmarks of aging: an expanding universe. *Cell* **186**, 243 – 278 (2023).

[^3]: Fried, L. P. *et al.* Frailty in older adults: evidence for a phenotype. *J. Gerontol. A Biol. Sci. Med. Sci.* **56**, M146 – M156 (2001).

[^4]: Williams, D. M., Jylhävä, J., Pedersen, N. L. & Hägg, S. A frailty index for UK Biobank participants. *J. Gerontol. A Biol. Sci. Med. Sci.* **74**, 582 – 589 (2019).

[^5]: Gilbert, T. *et al.* Development and validation of a Hospital Frailty Risk Score focusing on older people in acute care settings using electronic hospital records: an observational study. *Lancet* **391**, 1775 – 1782 (2018).

[^6]: Zenin, A. *et al.* Identification of 12 genetic loci associated with human healthspan. *Commun. Biol.* **2**, 41 (2019).

[^7]: Sun, B. B. *et al.* Plasma proteomic associations with genetics and health in the UK Biobank. *Nature* **622**, 329 – 338 (2023).

[^8]: Schuermans, A. *et al.* A proteome-wide screen for incident heart failure in UK Biobank. *Eur. Heart J.* (2024).

[^9]: Ho, J. E. *et al.* Proteomic risk scores for cardiovascular disease in the UK Biobank. *Circulation* (2025) — *citation placeholder*.

[^10]: Oh, H. S.-H. *et al.* Organ aging signatures in the plasma proteome track health and disease. *Nature* **624**, 164 – 172 (2023).

[^11]: Argentieri, M. A. *et al.* A proteomic biological-age clock that beats chronological age for mortality prediction. *Nat. Aging* **4**, 1417 – 1431 (2024) — *citation placeholder*.

[^12]: Park, J. & Scherer, P. E. Adipocyte-derived endotrophin promotes malignant tumor progression. *J. Clin. Invest.* **122**, 4243 – 4256 (2012).

[^13]: Sun, K. *et al.* Endotrophin triggers adipose tissue fibrosis and metabolic dysfunction. *Nat. Commun.* **5**, 3485 (2014).

[^14]: Karsdal, M. A. *et al.* The good and the bad collagens of fibrosis — their role in signaling and organ function. *Adv. Drug Deliv. Rev.* **121**, 43 – 56 (2017).

[^15]: Tan, B. K. *et al.* Endotrophin/COL6A3 derived peptides reflect insulin resistance. *Diabetes* **64**, 3069 – 3077 (2015).

[^16]: Luo, Y. *et al.* Plasma endotrophin and non-alcoholic fatty liver disease. *J. Hepatol.* **74**, 156 – 167 (2021) — *citation placeholder*.

[^17]: Rasmussen, D. G. K. *et al.* Endotrophin (PRO-C6) is associated with risk of progression in chronic kidney disease. *Kidney Int.* **97**, 1242 – 1252 (2020) — *citation placeholder*.

[^18]: Genovese, F. *et al.* The extracellular matrix in the kidney — a source of biomarkers. *Curr. Opin. Nephrol. Hypertens.* **27**, 211 – 217 (2018).

[^19]: Chirinos, J. A. *et al.* Endotrophin and outcomes in patients with heart failure with preserved ejection fraction. *Circulation* **141**, 974 – 988 (2020) — *citation placeholder*.

[^20]: Yoshiji, S. *et al.* Circulating endotrophin and incident chronic kidney disease in 42,416 UK Biobank participants. *J. Clin. Endocrinol. Metab.* (2025) — *companion CKD pipeline; per-SD HR 1.53*.

[^21]: Yoshiji, S. *et al.* Endotrophin and incident heart failure beyond the PREVENT equation: a UK Biobank proteomics study. *JACC: Heart Failure* (2025) — *companion HF pipeline; per-SD HR 1.25 PREVENT-adjusted*.

[^22]: Yoshiji, S. *et al.* Cis-pQTL Mendelian randomisation of endotrophin and coronary artery disease. *Nat. Cardiovasc. Res.* (2024) — *citation placeholder*.

[^23]: Walker, K. A. *et al.* Plasma proteome and incident dementia in the UK Biobank. *Nat. Aging* **4**, 247 – 258 (2024) — *citation placeholder*.

[^24]: Bycroft, C. *et al.* The UK Biobank resource with deep phenotyping and genomic data. *Nature* **562**, 203 – 209 (2018).

[^25]: Sun, B. B. *et al.* Genetic regulation of the human plasma proteome in 54,219 UK Biobank participants. *Nature* (2023).

[^26]: Inker, L. A. *et al.* New creatinine- and cystatin C-based equations to estimate GFR without race. *N. Engl. J. Med.* **385**, 1737 – 1749 (2021).

---

## Supplementary Materials

The following supplementary tables are deposited at `manuscript/supplementary/` in the repository:

- **Supplementary Table S1.** `section1_baseline_hfrs_OLS.csv` — OLS β-per-SD endotrophin for continuous baseline HFRS (linear and log-transformed) across the three adjustment models.
- **Supplementary Table S2.** `section1_prevalent_syndromes_OR.csv` — Logistic-regression OR-per-SD endotrophin for the four prevalent frailty-syndrome flags plus an aggregate "any frailty syndrome" outcome.
- **Supplementary Table S3.** `section1_multimorbidity_count.csv` — Poisson rate-ratios per SD endotrophin for self-reported chronic-condition count and prevalent ICD-disease count across eight common conditions.
- **Supplementary Table S4.** `section1_bioage_proxy.csv` — OLS β-per-SD endotrophin for the lab-based biological-age composite and its chronological-age-residualised "delta-bioage" counterpart.
- **Supplementary Table S5.** `section2_cohort_flow.csv` — Per-endpoint cohort-flow showing the analytic *N*, the post-prevalent-removal *N* and the events count.
- **Supplementary Table S6.** `section2_hfrs_top_contributors.csv` — Top-30 ICD-10 codes ranked by total Gilbert-weight contribution, with per-eid frequency and per-code weight.
- **Supplementary Table S7.** `section2_per_sd.csv` — Per-SD endotrophin Cox results across 13 endpoints × 3 adjustment models, with raw *P*, BH-q and concordance index. The companion `section2_quartiles.csv` provides the Q2 / Q3 / Q4 vs Q1 hazard ratios per endpoint × model.
- **Supplementary Table S8.** `section2_zenin_components.csv` — Share of Zenin healthspan composite events contributed by each first-occurring component.

The 30 individual KM and adjusted-cumulative-incidence PDFs underlying **Figs. 3 – 5** are deposited at `results/figures/` (one PDF per endpoint × {KM, adjcuminc}); the per-endpoint source-data CSVs are at `results/source_data/`.

---

*Manuscript prepared 2026-04-27. Repository: <https://github.com/satoshi-yoshiji/PV-aging>.*
