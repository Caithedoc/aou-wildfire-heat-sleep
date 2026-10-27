# Workbench execution environment

The scripts in this directory were run inside the All of Us Researcher Workbench (Verily Workbench), in the
workspace `aou-wildfire-sleep`, on a Compute Engine Jupyter application created from the program's `aou-jupyter`
app configuration on 2026-08-28 (`AoU_Jupyter_ComputeEngine_20260828`). Data were read from the Controlled Tier
curated data repository version 9 (C2025Q4R6) through BigQuery; intermediate tables were written to a
workspace-scoped BigQuery dataset (`wildfire_sleep`) and a Parquet cache on the VM. Nothing participant-level
left the Workbench.

`requirements-vm-freeze.txt` will be the output of `pip freeze` captured on the same application after all analyses
reported in the manuscript had been run (no packages were added or upgraded on the VM after the analyses). It is
committed before the tagged release; if this file is absent in the commit you are reading, the capture had not yet
been committed.

Each script states in its docstring the inputs it expects, the outputs it writes, and the command used to run it.
Run order: 01 (schema reconnaissance, read-only) -> 02 (cohort and outcome extraction) -> upload of the local
exposure panels (see `../local_exposure_prep/`) -> 04 (main case-crossover analysis) -> 06 (machine-learning
companions) -> 11, 12, 13, 14, 15, 16, 17, 18 (revision and sensitivity analyses, in numerical order) ->
19 then 20 (sleep-domain outcomes) -> 21 (Table 1 by exposure quartile).
