# Actual-case prediction checks

Run the maintained check with a completed real case and a separate output folder:

```sh
Rscript tests/real-data/prediction.R completed_case_directory QA_directory
```

The check uses every original hold-out location, saved full-data models and
actual reporting-scenario climate values. It compares paired execution with
ordinary engine predictions for total, driver-only and driver-removed
combinations, and checks the response API under both clamp settings. No models
are refitted, observations generated, trees reduced or spatial folds changed.

## Recorded check: 16 November 2026

The case has 135 fitted predictors; 6 climate variables change in the reporting
scenario. All 3,542 hold-out locations are used. There are 84 engine/policy/
combination checks. Maximum differences: RF 0, MaxEnt 0, GAM 1.11e-15, including
actual future values outside training ranges.

Windows, R 4.4.3, castSDM 0.40.0, Intel Xeon w7-3455 (24 physical / 48 logical
cores); prediction is limited to 2 threads. For 13 matched combinations,
ordinary predictions took 25.86 / 27.15 s and paired execution 7.86 / 7.85 s
(unclamped / clamped). Each is one execution-order-dependent measurement on
this machine, not a repeated cold/warm benchmark or an improvement estimate
against an earlier castSDM release. Response API setup is checked separately
and is not included in those paired/ordinary timings.

CSV results and original input content hashes are under `results/2026-11-16`.
They verify mathematical and implementation consistency, not independent
spatial explanation stability, causal identification or whole-workflow speed.
A publication comparison still needs matched outputs and domains, repetitions,
setup and I/O costs, memory measurements and another environment.
