# Intravenous Phosphate Supplementation in ICU Patients with Hypophosphatemia: A Target Trial Emulation Study

This repository contains the analysis code for a target trial emulation study evaluating the effect of intravenous phosphate supplementation on 28-day mortality in critically ill patients with hypophosphatemia, using the MIMIC-IV database (2008–2022).

## Study Overview

- **Design**: Sequential target trial emulation (Cloning, Censoring, Weighting framework)
- **Population**: ICU patients with hypophosphatemia (phosphate <0.80 mmol/L)
- **Exposure**: Intravenous phosphate supplementation within 24 hours of hypophosphatemia detection
- **Primary Outcome**: 28-day all-cause mortality
- **Analysis Method**: ITT analysis via SEQTaRget R package, with PP and CCW sensitivity analyses
- **Key Findings**: No significant overall mortality benefit (HR=0.987; 95% CI, 0.893–1.092; p=0.802), with potential protective effect in severe hypophosphatemia (phosphate <0.50 mmol/L; HR=0.806; p=0.022)

## Repository Structure

```
.
├── README.md                          # This file
├── README_CN.md                       # Chinese version of README
├── data_extraction/
│   └── python/                        # Python scripts for data extraction (PostgreSQL → Parquet)
├── analysis/
│   └── r/
│       ├── config/                    # Configuration files (covariates, paths, parameters)
│       ├── main/                      # Primary ITT analysis
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