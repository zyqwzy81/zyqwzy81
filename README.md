# Intravenous Phosphate Supplementation in ICU Patients with Hypophosphatemia: A Target Trial Emulation Study

This repository contains the analysis code for a target trial emulation study evaluating the effect of intravenous phosphate supplementation on 28-day mortality in critically ill patients with hypophosphatemia, using the MIMIC-IV database (2008–2022).

## Study Overview

- **Design**: Sequential target trial emulation (Cloning, Censoring, Weighting framework)
- **Population**: ICU patients with mild-to-moderate hypophosphatemia (phosphate 0.32–0.80 mmol/L)
- **Exposure**: Intravenous phosphate supplementation within 24 hours of hypophosphatemia detection
- **Primary Outcome**: 28-day all-cause mortality
- **Analysis Method**: ITT analysis via SEQTaRget R package, with PP and CCW sensitivity analyses
- **Key Findings**: No significant overall mortality benefit (HR=0.987; 95% CI, 0.893–1.092; p=0.802), with a potential protective effect in the lower end of the mild-to-moderate range (baseline phosphate <0.50 mmol/L; HR=0.806; p=0.022)

## Repository Structure

```
.
├── README.md                          # This file
├── README_CN.md                       # Chinese version of README
├── data_extraction/
│   ├── python/                        # Python scripts for data extraction (PostgreSQL → Parquet)
│   └── sql/pure_control/              # SQL for the post-hoc pure-control analyses (M0_40–M0_48)
├── analysis/
│   └── r/
│       ├── v28.4/                     # Revision code (current): multi-enrollment SEQTaRget pipeline
│       │   ├── pure_control/          # Post-hoc IV vs no-supplement mechanism analyses (PC01–PC12b)
│       │   └── v28.4_*.R              # Main/PP/subgroup/piecewise/secondary/continuous/90-day/Figure 2
│       ├── config/                    # Configuration files (covariates, paths, parameters)
│       ├── main/                      # Primary ITT analysis (original submission)
│       ├── pp_analysis/               # Per-protocol analysis
│       ├── sensitivity/               # Sensitivity analyses (CCA, grace period, truncation, MNAR, LASSO vs DAG)
│       ├── secondary/                 # Secondary outcomes analysis
│       ├── subgroup/                  # Subgroup analyses
│       ├── interaction/               # Interaction tests
│       ├── external_validation/       # eICU external validation
│       ├── ccw_engine/                # Batch CCW engine (memory-efficient SEQTaRget alternative)
│       ├── utils/                     # Helper functions and diagnostics
│       └── tests/                     # Unit tests
└── docs/                              # Documentation
```

## Code Versions

- **v28.4 (`analysis/r/v28.4/`, current — journal revision)**: rebuilt multi-enrollment (trial-period level) sequential target trial emulation with a natural-calendar daily timeline, including ITT/PP main analyses, subgroup and piecewise hazards, 90-day mortality, secondary binary/continuous outcomes, the full-horizon weighted Kaplan–Meier analysis, and the post-hoc pure-control analyses (`pure_control/`). Common modules: `v28.4_common.R`, `v28.4_multi_common.R`, `v28.4_stay_common.R`.
- The stage-organized directories (`main/`, `pp_analysis/`, `sensitivity/`, …) contain the original v28.3 submission code retained for reference.

**Database credential**: scripts never embed the PostgreSQL password. Set it via the `MIMIC_PG_PASSWORD` environment variable before running any script that connects to the local MIMIC-IV database (e.g. `MIMIC_PG_PASSWORD=xxx Rscript v28.4_secondary_outcomes.R`).

## Dependencies

### R Packages
- **SEQTaRget** (≥0.1.0): Sequential target trial emulation framework
- **survival**: Cox proportional hazards models
- **mice** (≥3.16.0): Multiple imputation
- **glmnet**: LASSO variable selection
- **ggplot2**, **survminer**, **forestplot**: Publication-grade visualizations
- **cobalt**: Covariate balance assessment
- **EValue**: E-value sensitivity analysis
- **arrow**: Parquet data reading
- **sandwich**, **lmtest**: Cluster-robust standard errors

### Python Packages
- **pandas**, **pyarrow**: Parquet data handling
- **SQLAlchemy**: PostgreSQL database connection

## Data Source

This study uses the **MIMIC-IV v3.1** database (PhysioNet, https://physionet.org/content/mimiciv/3.1/). Access requires:
1. PhysioNet credentialed access
2. CITI Data or Specimens Only Research completion
3. Signed data use agreement

The analysis dataset is built from 9 MIMIC-IV tables via the Python scripts in `data_extraction/`.

## Reproducing the Analysis

### Step 1: Data Extraction
```bash
cd data_extraction/python
python build_sequential_trials.py
python build_analysis_dataset.py
```

### Step 2: Run Analysis
```bash
cd analysis/r
Rscript misc/run_all_analyses.R
```

## Registration

The study protocol was retrospectively registered in the Open Science Framework:
https://osf.io/yc6kb/

## Citation

If you use this code, please cite:

[Citation information to be added upon publication]

## License

This code is available for research purposes. The MIMIC-IV database is subject to separate data use agreements.