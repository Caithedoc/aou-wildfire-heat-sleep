# Workbench execution environment

The scripts in this directory were run inside the All of Us Researcher Workbench (Verily Workbench), in the
workspace `aou-wildfire-sleep`, on a Compute Engine Jupyter application created from the program's `aou-jupyter`
app configuration on 2026-08-28 (`AoU_Jupyter_ComputeEngine_20260828`). Data were read from the Controlled Tier
curated data repository version 9 (C2025Q4R6) through BigQuery; intermediate tables were written to a
workspace-scoped BigQuery dataset (`wildfire_sleep`) and a Parquet cache on the VM. Nothing participant-level
left the Workbench.

`requirements-vm-freeze.txt` is the output of `pip list --format=freeze` captured on the same application on
2026-10-27, after all analyses reported in the manuscript had been run; no packages were added or upgraded on the VM
between the analyses and the capture. The interpreter was Python 3.12.13 (header line of the file).

Each script states in its docstring the inputs it expects, the outputs it writes, and the command used to run it.
Run order (2026-08-28 unless stated): 01 (schema reconnaissance, read-only; notebook cells) -> 02 (cohort and
outcome extraction; notebook cells) -> upload of the local exposure panels (see `../local_exposure_prep/`) ->
04 (main case-crossover analysis) -> 05, 06, 07, 08, 09, 10 (sensitivity, machine-learning companions, descriptive
tables, negative controls) -> 11 (joint-exposure revision), 12_reconcile, 13_ipw -> 12_recluster_tempflex
(2026-08-31) -> 13_zip3_heterogeneity_restrict, 14 (2026-09-22) -> 15 (2026-09-24) -> 19 then 20 (sleep-domain
outcomes, 2026-10-12) -> 21 (Table 1 by exposure quartile, 2026-10-13).

The scripts in this directory are the copies that were on the Workbench application on 2026-10-27, except 01 and
02, which were run as notebook cells and kept only locally. `04_casecrossover_analysis.py` and `07_remaining.py`
are the as-run versions after two in-place edits made on the VM (a Parquet cache for the analysis-night table in
04; a syntax fix to a date filter in 07); `11_revision.py` differs from the local copy only by a trailing newline.
