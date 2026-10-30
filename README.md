# Joint exposure to wildfire smoke and extreme heat and wearable-measured sleep: analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23700934.svg)](https://doi.org/10.5281/zenodo.23700934)

Code, software-environment records, and the dated analysis plan for the manuscript

> Zhang C, Trevathan E, Qian Z. *Joint Exposure to Wildfire Smoke and Extreme Heat and Wearable-Measured Sleep Characteristics Among US Adults.* (manuscript under submission, 2026)

The study is a time-stratified self-matched (case-crossover) analysis of 1,633,414 person-nights from 30,268 adults with
consumer-wearable sleep records in the All of Us Research Program (Controlled Tier curated data repository version 9,
C2025Q4R6), 2016 through 2023. Exposures are daily wildfire-specific PM2.5 (smoke nights: > 20 µg/m³) and daily mean
temperature (hot nights: at or above the location- and calendar-month-specific 90th percentile), assigned at the
3-digit ZIP code level. The primary outcome is total sleep time in the main sleep episode; secondary outcomes cover
sleep continuity, timing, and architecture.

**This repository contains no data.** No participant-level data, no All of Us data of any kind, and no aggregate
results are included. All of Us data are available to authorized users through the Researcher Workbench
(https://www.researchallofus.org/). The exposure inputs are public (see below).

## Repository layout

| Directory | Contents | Where it runs |
|---|---|---|
| `analysis_plan/` | The dated analysis plan (2026-08-26), the revision plan (2026-08-28), and two prespecification addenda (2026-10-12). Checksums and caveats are in `PROVENANCE.md`. | n/a |
| `local_exposure_prep/` | Builds the daily 3-digit ZIP exposure panels from public data and quantifies ZCTA-to-ZIP3 aggregation error. `requirements-lock.txt` is the `pip freeze` of the environment used. | Local machine |
| `workbench/` | Cohort extraction, the main case-crossover analysis, and all revision and sensitivity analyses. `ENVIRONMENT.md` describes the Workbench application; `requirements-vm-freeze.txt` is the `pip freeze` captured on it. | All of Us Researcher Workbench |
| `figures/` | Scripts that draw the main figures and supplementary eFigures from aggregate outputs. `requirements-figures.txt` lists the package versions used. | Local machine |

Script comments and docstrings are in Chinese (the working language of the analyst); every docstring states the
script's inputs, outputs, and run command. English summaries follow.

## Pipeline

### 1. Exposure panels (`local_exposure_prep/`, public data only)

| Script | What it does |
|---|---|
| `download_prism.sh`, `download_prism_tmin.sh` | Download PRISM 4-km daily mean and daily minimum temperature GeoTIFFs for 2016-2023 (https://prism.oregonstate.edu/). |
| `03_smokepm_zip3_prep.py` | Three stages (`smoke`, `temp`, `merge`). Aggregates daily wildfire-specific PM2.5 at ZCTA resolution (Stanford ECHO Lab smokePM v2.0; Childs et al., *Environ Sci Technol* 2022) and PRISM daily mean temperature sampled at ZCTA centroids to 3-digit ZIP code-days using 2020 Census ZCTA population weights; writes `out/smoke_temp_zip3_daily.csv`. The panel is restricted to the contiguous United States (PRISM coverage). |
| `05_tmin_zip3_prep.py` | Same construction for PRISM daily minimum temperature (used in the night-time temperature sensitivity models). |
| `07_zip3_aggregation_error.py` | Treats ZCTA-day exposure as the reference and ZIP3-day exposure as the assigned value, and reports population-weighted sensitivity, positive predictive value, specificity, and kappa for smoke, hot, and joint days; also computes the per-ZIP3 smoke-classification heterogeneity index used by `workbench/13_zip3_heterogeneity_restrict.py`. The restriction rules were written before any All of Us result was seen (2026-09-22). |

Inputs expected under `local_exposure_prep/raw/` (not included): the smokePM v2.0 daily ZCTA predictions
(`smokePM2pt5_predictions_daily_zcta_20060101-20231231.rds`), 2020 Census ZCTA population (DHC table P1), the
Census ZCTA gazetteer (centroids), and the PRISM daily files.

### 2. Cohort, exposure linkage, and models (`workbench/`, inside the All of Us Researcher Workbench)

| Script | What it does |
|---|---|
| `01_cdr_recon.py` | Read-only reconnaissance of the Controlled Tier tables and columns (wearables sleep tables, ZIP3, dates). Run cell by cell in a notebook. |
| `02_sleep_cohort_extract.py` | Extracts main-sleep-episode nights and participant covariates into a workspace BigQuery dataset. Run cell by cell in a notebook. Nothing participant-level leaves the Workbench. |
| `04_casecrossover_analysis.py` | Links nights to the uploaded ZIP3-day exposure panel and fits the prespecified time-stratified case-crossover models (strata = person × calendar month × day of week; within-stratum fixed-effects OLS), including distributed lags (0-7), continuous smoke, the no-temperature-adjustment model, and the next-day (lead) negative control. |
| `05_sensitivity.py` | Alternative smoke thresholds (10 and 35 µg/m³) and stratification by device data source. |
| `06_phaseE.py` | Machine-learning companions: counterfactual expected-sleep baseline (gradient boosting), heterogeneous-effect quintiles with GATES calibration, prespecified subgroups, and the first smoke × heat joint model. |
| `07_remaining.py` | Inverse-probability weighting for nightly wear, time negative control, exclusion of hot days, the June 2023 East Coast smoke episode, and West versus East comparisons. |
| `08_tables.py` | Attrition waterfall and the descriptive tables (participant and night level, by year). |
| `09_final.py` | Time negative control (redefined), wider sleep-window sensitivity (60-960 minutes), exclusion of 2023, and E-values. |
| `10_negctrl.py` | Time negative control restricted to low-smoke ZIP3-years in fire-season months. |
| `11_revision.py` | Joint (compound) exposure analyses: four-category smoke × heat models, continuous joint response curves and surface, the 3 × 3 threshold robustness matrix, formal subgroup interaction tests, secondary outcomes, Table 1, and Figure 1 count exports. |
| `12_reconcile.py` | Reconciles the joint-exposure night counts between the pre-filter and regression frames and writes the final four-category and interaction models. |
| `12_recluster_tempflex.py` | Variance and temperature-specification ladder: standard errors clustered by person, ZIP3, ZIP3-day, and two-way, plus ZIP3-level wild cluster bootstrap *P* values; temperature adjustment from none through linear, quadratic, natural splines of mean temperature, and night-time minimum temperature; delta-method confidence intervals for the joint-minus-additive contrast. |
| `13_ipw.py` | Inverse-probability-of-wear weighting for the joint-exposure interaction, with the wear model trained on a sampled person-day panel. |
| `13_zip3_heterogeneity_restrict.py` | Sensitivity analyses restricted to ZIP3 areas with low within-area smoke heterogeneity (rules prespecified in `07_zip3_aggregation_error.py`). |
| `14_region_by_heterogeneity.py` | Post hoc regional decomposition of the heterogeneity-restricted result. |
| `15_subgroup_characteristics.py` | Participant characteristics by subgroup and the remaining marginal subgroup estimates. |
| `19_sleep_domains_rebuild.py`, `20_sleep_domains_joint.py` | Sleep-domain outcomes (continuity, timing, architecture, and composite sleep-health indicators) rebuilt from the wearables summary tables and fitted with the same model specification; prespecified in the 2026-10-12 addenda. |
| `21_table1_quartiles.py` | Table 1 overall and by quartile of participant-mean wildfire-specific PM2.5. |
| `22_table1_smoking.py` | Smoking status (never, former, current, not reported) from the All of Us Lifestyle survey for Table 1, overall and by quartile; same sample and quartile cut points as `21_table1_quartiles.py`. Added 2026-10-30 (v1.1.0) at a co-author's request. |

Run order and the execution environment are described in `workbench/ENVIRONMENT.md`. Every export from the Workbench is an aggregate table in which any cell corresponding to 1-20 participants is suppressed, in line with the
All of Us Data and Statistics Dissemination Policy.

### 3. Figures (`figures/`)

| Script | Output |
|---|---|
| `make_fig_compound_zcta_map.py` | Figure 1: hot-and-smoky days per ZCTA, 2016-2023 (public exposure data only). |
| `make_fig_jno_v9.py` | Figure 2 (joint exposure categories) and eFigure 7 (continuous joint response curves). |
| `make_efig_jno_v9.py` | eFigures 1-6 (flowchart, subgroups, distributed lag, causal-forest calibration, specification robustness, exposure nights over time and by region). |
| `lu_palette.py` | Shared palette used by the map figure. |
| `make_figures.py`, `make_fig1_overview.py` to `make_fig4_interaction_robustness.py` | Figures of earlier versions of the manuscript, kept for completeness. |

The figure scripts read the aggregate CSVs exported from the Workbench (and, for some eFigures, the supplement
tables); those files are not distributed here. The numbers they contain are the ones reported in the manuscript and
its supplement.

## Software environments

* Workbench: see `workbench/ENVIRONMENT.md` and `workbench/requirements-vm-freeze.txt`.
* Local exposure preparation: `local_exposure_prep/requirements-lock.txt` (Python 3.14.7).
* Figures: `figures/requirements-figures.txt` (Python 3.9.6).

## Prespecification and provenance

`analysis_plan/` holds the plan written before data access (2026-08-26), the revision plan written after the first
full run (2026-08-28), and two addenda (2026-10-12) that prespecified the composite sleep-health and sleep-domain
outcomes after coauthor review and before those models were fitted. `PROVENANCE.md` gives each document's checksum,
the checksum recorded in the project log at finalisation where one exists, and one known discrepancy.

## Citing

Please cite the manuscript and this repository (`CITATION.cff`). All released versions are archived on Zenodo under the concept DOI https://doi.org/10.5281/zenodo.23700934, which
resolves to the latest version. Version DOIs: v1.0.0 (2026-10-27) https://doi.org/10.5281/zenodo.23700935;
v1.1.0 (2026-10-30, adds `22_table1_smoking.py`) is listed on the Zenodo page.

## License and contact

MIT License (see `LICENSE`). Cai Zhang, Department of Health and Kinesiology, University of Illinois
Urbana-Champaign (caiz2@illinois.edu).
