# Intravenous Phosphate Supplementation in ICU Patients with Hypophosphatemia: A Target Trial Emulation Study

This repository contains the analysis code for a target trial emulation study evaluating the effect of intravenous phosphate supplementation on 28-day mortality in critically ill patients with hypophosphatemia, using the MIMIC-IV database (2008–2022).

## Study Overview

- **Design**: Sequential target trial emulation (Cloning, Censoring, Weighting framework)
- **Population**: ICU patients with mild-to-moderate hypophosphatemia (phosphate 0.32–0.80 mmol/L)
- **Exposure**: Intravenous phosphate supplementation within 24 hours of hypophosphatemia detection
- **Primary Outcome**: 28-day all-cause mortality
- **Analysis Method**: ITT analysis via SEQTaRget R package, with PP and CCW sensitivity analyses
- **Key Findings**: No significant overall mortality benefit (HR=0.982; 95% CI, 0.897–1.076; p=0.695). A lower phosphate-stratum exploratory result (baseline phosphate <0.50 mmol/L; HR=0.765; p=0.0097 for interaction) was downgraded to hypothesis-generating because of unit-dependence and lack of multiplicity correction; stay-level sensitivity was non-significant (HR=0.746; p=0.117).

## Repository Structure

```
.
├── README.md                              # This file
├── README_CN.md                           # Chinese version of README
├── data_extraction/
│   ├── python/                            # Python scripts for data extraction (PostgreSQL → Parquet)
│   │   ├── extract_event_times.py         #   KDIGO stages + PRBC infusion event times for v28.4 secondary outcomes
│   │   ├── build_sequential_trials.py     #   Sequential trial origin construction (predecessor of all analyses)
│   │   ├── augment_full_dataset.py        #   Augmented trial-level dataset with baseline covariates
│   │   └── …                              #   Other pre-analysis data builders
│   └── sql/pure_control/                  # SQL for the post-hoc pure-control analyses (M0_40–M0_48)
├── analysis/
│   └── r/
│       ├── main_multi.R                   # Main ITT (trial-period level sequential target trial)
│       ├── piecewise_multi.R              # Piecewise hazards + per-protocol with IPCW
│       ├── subgroup_multi.R               # Subgroup heterogeneity tests (trial-period level)
│       ├── secondary_multi.R              # Secondary outcomes (binary + continuous, MI-bootstrap)
│       ├── death90_multi.R                # 90-day mortality (patients.dod source-verified)
│       ├── fig2_fullkm.R                  # Full-horizon weighted Kaplan–Meier (no discharge truncation)
│       ├── binary_continuous.R            # Binary + continuous secondary outcomes (stay-level sensitivity)
│       ├── daily_threshold.R              # Daily threshold sweep (stay-level)
│       ├── stay_common.R / multi_common.R / common.R   # Shared modules (stay-level / trial-period / generic)
│       ├── generate_all_figures.R         # Publication figures (consumes result tables)
│       └── pure_control/                  # Post-hoc IV vs no-supplement mechanism analyses (PC01–PC12b)
└── .gitignore
```

## Analysis Code

All committed R scripts implement the **revision** analysis (journal major-revision round, 2026-10-10), using a rebuilt **trial-period level sequential target trial emulation** with a natural-calendar daily timeline. This is the definitive analysis for the manuscript. Common modules: `common.R`, `multi_common.R`, `stay_common.R`. Source scripts still embed absolute paths (`/data/小米云盘/文档/论文/血磷/TTE/v28.4/code/r/...`) to their local working locations — the version prefix has been stripped from filenames inside the repo, but the local path prefix intentionally remains so that scripts continue to run on the author's workstation.

**Database credential**: scripts never embed the PostgreSQL password. Set it via the `MIMIC_PG_PASSWORD` environment variable before running any script that connects to the local MIMIC-IV database (e.g. `MIMIC_PG_PASSWORD=xxx Rscript secondary_outcomes.R`).

**Patient-level data (parquet, csv) is not committed** — the analysis tables they produce are written locally; authors with PhysioNet credentialing can rebuild them from MIMIC-IV using the `data_extraction/python/` scripts.

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