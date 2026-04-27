# Manuscript — PV-aging

Manuscript draft for the ETP-aging analysis.

```
manuscript/
├── Manuscript_ETP_aging.v1.md  # full manuscript with embedded figures
├── figures/
│   ├── Figure_1_cross_sectional.png        # composite, 3 panels
│   ├── Figure_2_master_forest.png          # 13-endpoint forest, faceted by adjustment model
│   ├── Figure_3_primary_KM.png             # composite, 4 primary endpoints (mortality, HFRS, Zenin, CKD)
│   ├── Figure_4_frailty_syndromes_KM.png   # composite, 4 frailty syndromes
│   ├── Figure_5_disease_panel_KM.png       # composite, disease panel + ICD-event substitutes
│   ├── section1_*.png                      # raw per-figure PNGs (Section 1)
│   └── section2_*.png                      # raw per-endpoint KM and adjusted-cuminc PNGs
└── supplementary/
    ├── section1_baseline_hfrs_OLS.csv          # Suppl. Table S1
    ├── section1_prevalent_syndromes_OR.csv     # Suppl. Table S2
    ├── section1_multimorbidity_count.csv       # Suppl. Table S3
    ├── section1_bioage_proxy.csv               # Suppl. Table S4
    ├── section2_cohort_flow.csv                # Suppl. Table S5
    ├── section2_hfrs_top_contributors.csv      # Suppl. Table S6
    ├── section2_per_sd.csv                     # Suppl. Table S7 (per-SD HRs)
    ├── section2_quartiles.csv                  # Suppl. Table S7 (quartile HRs)
    └── section2_zenin_components.csv           # Suppl. Table S8
```

The original PDF figures (300-DPI vector) live in `../results/figures/`; the
KM and adjusted-cuminc source-data CSVs live in `../results/source_data/`.
The `figures/Figure_*.png` files in this directory are 200-DPI PNG composites
built for inline rendering of the markdown manuscript.

To regenerate the composite figures after rerunning the pipeline:

```bash
conda activate ckd_analysis
cd /humgen/diabetes2/users/satoshi/misc/02.aging
python <<'EOF'
import fitz, os, shutil
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
from pathlib import Path
src_pdf = Path("results/figures"); dst = Path("manuscript/figures"); dst.mkdir(exist_ok=True)
for f in sorted(os.listdir(src_pdf)):
    if f.endswith(".pdf"):
        doc = fitz.open(src_pdf/f); doc[0].get_pixmap(dpi=200).save(dst/f.replace(".pdf",".png")); doc.close()
# (rebuild Figure_1..5 composites — see scripts/ in commit history)
EOF
```
