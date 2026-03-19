# Project Skills Overview for UW–RTSL Cardiovascular Health Investment Case Assessment

This document summarizes the technical, analytical, and domain-specific skills required to contribute effectively to the **UW–RTSL Cardiovascular Health Investment Case Assessment**, a multi-country simulation project estimating the long-term health impact of Resolve to Save Lives (RTSL) cardiovascular health programs.

---

## Core Competencies

### Epidemiologic Modeling
- Experience with **multi-state (illness–death) models** for chronic disease.
- Familiarity with **multiplicative risk frameworks** and cause-specific hazard adjustments.
- Ability to integrate **relative risks**, **exposure distributions**, and **intervention effect sizes** into population-level projections.

### Data Engineering & Management
- Proficiency handling large, multi-country datasets.
- Experience with **UNWPP population projections**, **GBD incidence/mortality**, and **risk factor distributions**.
- Ability to manage structured directories for raw, processed, and adjusted data.

### R Programming & Workflow Automation
- Strong command of **R**, including:
  - `data.table`, `dplyr`, `readxl`, `ggplot2`
  - Parallelization using `doParallel` and `foreach`
  - Writing modular, reproducible scripts
- Ability to read, modify, and extend complex modeling scripts.

### Intervention Modeling
- Understanding of cardiovascular health interventions:
  - **Hypertension control** (BP distributions, treatment effects, control trajectories)
  - **Sodium reduction** (dose–response BP shifts, scale-up timelines)
  - **Trans-fat elimination** (policy timing, lag structures)
  - **Statin treatment** (coverage, adherence, incidence vs. case-fatality effects)
- Ability to parameterize and run **single** and **combined** intervention scenarios.

### Calibration & Trend Adjustment
- Familiarity with calibration of baseline incidence and mortality rates.
- Understanding of optional trend flags:
  - Background mortality trends
  - Case fatality trends (GBD - IHME)
  - Adjustment models for age–sex–cause transitions

### Scenario Design & Sensitivity Analysis
- Ability to define, modify, and document scenario files.
- Experience running multi-country simulations with varying assumptions.
- Skill in interpreting scenario outputs for policy relevance.

### Output Interpretation & Communication
- Ability to synthesize results into:
  - **Deaths and DALYs averted**
  - Cause-, age-, and sex-specific summaries
  - Cross-scenario comparisons
- Experience preparing materials for **policy audiences**, **donor communications**, and **interactive tools**.

---

## Technical Environment

### Required Tools
- R (≥ 4.2)
- RStudio project structure
- Git/GitHub for version control

### Key Input Files
- UNWPP 2024 population projections
- GBD 2023 calibrated baseline rates
- Country-specific BP distributions
- RTSL program coverage trajectories
- Relative risk datasets (GBD 2019, Ettehad et al.)

### Directory Structure
- `code/` for modeling scripts
- `data/` and `adjusted/` for inputs
- `output/` for per-country `.rds` results
- `scenarios/` for intervention definitions

---

## Collaboration Expectations
- Ability to work within a **multi-institutional team** (UW + RTSL).
- Clear documentation of code changes and assumptions.
- Responsiveness to iterative model updates and new data inputs.
- Commitment to reproducibility and transparent analytical workflows.

---

## Summary
Contributors to this project must combine **epidemiologic modeling expertise**, **strong R programming skills**, and **deep familiarity with cardiovascular health interventions**. The work requires comfort with large datasets, multi-country simulation workflows, and translating technical outputs into policy-relevant insights.

---
