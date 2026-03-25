# UW–RTSL Cardiovascular Health Investment Case Assessment

> **Projected health impact of Resolve to Save Lives cardiovascular health programs (2017–2024) through 2050**

---

## Research Question

What is the expected health impact of Resolve to Save Lives (RTSL) cardiovascular health programs from 2017 to 2024 in the populations targeted by these programs, projected to the year 2050?

---

## Overview

This repository contains the simulation model and analysis code used to estimate deaths and DALYs averted through RTSL cardiovascular health (CVH) programs under maintained 2024 implementation levels, compared to a counterfactual in which 2017 conditions persist into the future.

The model supports four CVH intervention modules:

- **Hypertension control** — antihypertensive treatment scale-up in Global Hearts program catchment populations
- **Sodium reduction** — dietary sodium reduction or replacement policies
- **Trans-fat elimination (TFA)** — elimination of industrially produced trans-fatty acids through national policies
- **Statin treatment** — lipid-lowering therapy for cardiovascular disease prevention

Results are intended to demonstrate the return on donor investment in RTSL programs and to support an interactive online simulation tool co-produced by the University of Washington (UW) and RTSL.

---

## Repository Structure

```
uwrtsl-assessment/
├── code/           # R scripts for data preparation, calibration, modeling, and analysis
├── data/           # Input data (raw and adjusted); not tracked by Git — see Data section
├── docs/           # Documentation, methods notes, and supplementary materials
├── output/         # Model outputs (RDS files per country, summary tables); not tracked by Git
├── scenarios/      # Scenario definition files
├── .gitignore
├── .Rhistory
├── README.md
└── uwrtsl-assessment.Rproj
```

---

## Country / Population Scope

| Intervention | Population |
|---|---|
| Hypertension control | Catchment population covered by Global Hearts HTN program as of end of 2024 |
| Sodium reduction | To be defined (TBD) |
| Trans-fat elimination | National populations covered by passed or implemented TFA policies by end of 2024 |
| Statin treatment | Countries with available baseline coverage data |


---

## Methods

### Model Design

The model uses a **multiplicative, illness-death (multi-state) framework** that tracks individuals through healthy, diseased, and dead states for five CVD causes:

| Abbreviation | Cause |
|---|---|
| `ihd` | Ischemic heart disease |
| `istroke` | Ischemic stroke |
| `hstroke` | Intracerebral hemorrhage |
| `hhd` | Hypertensive heart disease |


### Data Sources

| Data | Source |
|---|---|
| Population projections | UN World Population Prospects (UNWPP) 2024 |
| Baseline mortality & incidence rates | GBD 2023 (calibrated) |
| Blood pressure distributions | Country-specific BP data (`bp_data6.csv`) |
| COVID-19 excess mortality | UNWPP excess mortality estimates |
| Relative risks (hypertension) | GBD 2019 relative risks; Ettehad *et al.* meta-analysis |
| HTN control scale-up trajectories | `covfxn2.csv` (RTSL program data) |
| Statin scenarios | Country-level statin coverage scenarios |
| TFA scenarios | Country-level TFA policy scenarios |

### Intervention Modules

#### 1. Antihypertensive Treatment
- Linear scale-up from baseline to a target control level (default: 50% controlled by 2040)
- BP-category-specific relative risks from GBD 2019 and Ettehad *et al.*
- Aggregate weighted coverage computed across hypertensive BP bins

#### 2. Sodium Reduction
- Percent or absolute reduction in population mean sodium intake
- Scale-up from `saltyear1` to `saltyear2`
- BP shift derived from sodium–BP dose-response

#### 3. Trans-Fat Elimination
- Target TFA % energy → 0% (full elimination)
- Policy effect lagged two years from `tfa_policy_start_year`
- Applied to countries with implemented or passed TFA policies by end of 2024

#### 4. Statin Treatment
- Linear scale-up to target coverage (default: 60% by 2050)
- Adherence adjustment for incidence (`adherence_ir`) and case fatality (`adherence_cf`) states
- Cause-specific relative risk reductions applied multiplicatively

### Counterfactual

2017 conditions (BP distributions, coverage levels, mortality rates) are assumed to persist unchanged into the future. The intervention scenario assumes 2024 implementation levels are maintained through 2050, with no further expansion beyond enrolled countries/populations.

### Trend Adjustments (Optional Flags)

| Flag | Description |
|---|---|
| `run_adjustment_model` | Apply GBD 2023 age-sex-cause–specific IR/CF adjustment factors |
| `run_bgmx_trend` | Apply forecasted downward trend to background mortality |
| `run_CF_trend` | Apply forecasted downward trend to case fatality rates |
| `run_CF_trend_ihme` | Use IHME-based (vs. GBD) CVD case fatality trends |
| `run_CF_trend_80` | Apply 80% of forecasted CFR trend |

### Parallelization

Country-level simulations are run in parallel using `doParallel` / `foreach`. Each country is processed independently, with results saved as individual `.rds` files in `output/out_model/` and per-country log files written for diagnostic review.

---

## Scenarios

The model supports running combinations of interventions simultaneously. Default scenarios:

| Scenario | Interventions Active |
|---|---|
| `baseline` | None (counterfactual) |
| `bp_only` | Antihypertensive treatment |
| `sodium_only` | Sodium reduction |
| `tfa_only` | Trans-fat elimination |
| `statins_only` | Statin treatment |
| `bp_sodium` | Antihypertensive + sodium |
| `bp_sodium_tfa` | Antihypertensive + sodium + TFA |
| `all_four` | All four interventions combined |

Scenario parameters are fully configurable to support sensitivity analyses and to generate updated estimates when country policies change or new countries are enrolled.

---

## Key Parameters (Default Values)

| Parameter | Value | Description |
|---|---|---|
| `target_control` | 0.50 | HTN control target (50%) |
| `control_start_year` | 2026 | HTN scale-up start |
| `control_target_year` | 2040 | HTN target achievement year |
| `saltmet` | `"percent"` | Sodium reduction method |
| `salteff` | 0.30 | 30% sodium reduction |
| `saltyear1` | 2026 | Sodium scale-up start |
| `saltyear2` | 2030 | Sodium scale-up end |
| `tfa_target_tfa` | 0 | TFA elimination target |
| `tfa_policy_start_year` | 2028 | TFA policy effect start |
| `statin_target_coverage` | 0.60 | Statin coverage target (60%) |
| `statin_start_year` | 2026 | Statin scale-up start |
| `statin_target_year` | 2050 | Statin target achievement year |
| `adherence_ir` | 0.575 | Statin adherence (incidence) |
| `adherence_cf` | 0.664 | Statin adherence (case fatality) |
| `ncores` | 6 | Parallel CPU cores |

---

## Setup & Usage

### Prerequisites

R (≥ 4.2 recommended) with the following packages:

```r
install.packages(c("data.table", "dplyr", "readxl", "ggplot2", "doParallel", "foreach"))
```

A helper functions file (`functions_review_6_dm.R`) must be sourced before running the model.

### Directory Configuration

Set the following directory paths before running:

```r
wd_raw   <- "/path/to/raw/"        # Raw GBD/RTSL input files
wd       <- "/path/to/data/"       # Processed data files
wd_data  <- "/path/to/adjusted/"   # Calibrated baseline rate RDS files
wd_outp  <- "/path/to/output/"     # Model output destination
```

### Required Input Files

| File | Location | Description |
|---|---|---|
| `wpp.adj.Rda` | `wd` | UNWPP-adjusted mortality including COVID excess |
| `PopulationsAge20_2050.csv` | `wd_data/UN2024/` | UNWPP 2024 population projections |
| `PopulationsSingleAge0050.rds` | `wd_data/UN2024/` | Single-age population file |
| `bp_data6.csv` | `wd` | Country blood pressure distributions |
| `covfxn2.csv` | `wd_data` | HTN program coverage trajectories |
| `adjusted*.rds` | `wd_data` | Calibrated baseline transition rate files |
| `ettehad_rr_bp_reduction_10mmHg.csv` | `wd` | Ettehad relative risks |
| `ettehad_rr_bp_reduction_effects.xlsx` | `wd` | Ettehad cumulative BP bin effects |
| `IHME_GBD_2019_RELATIVE_RISKS_*.xlsx` | `wd_raw/GBD/` | GBD 2019 HTN relative risks |

### Running the Model

```r
# 1. Source helper functions
source("functions_review_6_dm.R")

# 2. Set working directories (see above)

# 3. Configure intervention parameters

# 4. Open and run the main model script:
#    code/7__model_100_multiplicative_func_new3_gbd23_noaroc.R
```

Country-level outputs are saved to `output/out_model/model_output_part_<Country>.rds`.  
Per-country logs are saved to `output/out_model/log_<Country>.txt`.

---

## Expected Outputs

- **Deaths averted** (2017–2024 and 2017–2050) by country, intervention, sex, age, and CVD cause
- **DALYs averted** corresponding to the above
- Scenario comparisons across single and combined interventions
- Validation summaries for each country run

---

## Deliverables

1. **Interim results** for RTSL advocacy and fundraising communications
2. **UW "CVH Investment Case" interactive simulation tool** — an online, user-friendly tool for generating country-level outcomes under different RTSL CVH interventions (hypertension treatment, statin treatment, sodium reduction/replacement, and trans-fat elimination), co-produced by UW and RTSL teams

---

## Team & Collaboration

This project is a collaboration between the **University of Washington (UW)** and **Resolve to Save Lives (RTSL)**.

---

## License

*To be specified.*

---

## Citation

*To be added upon publication.*
