# Environmental drivers of habitat change: development and validation

Scientific question: which environmental changes improve or reduce predicted habitat suitability, where do they act, and how reproducible are those spatial conclusions?

## Implementation contract

- One paired-scenario API accepts current/future predictor tables or aligned rasters and an existing fitted SDM. Keep one model set, one set of weights and one suitability scale throughout all comparisons.
- Total change: future prediction minus current prediction.
- Driver-only effect: replace that driver's current values with its future values, holding other predictors at current values.
- Driver-retention contrast: full future prediction minus the future prediction with that driver held at current values. These are scenario comparisons in a fitted model, not identified ecological causal effects.
- Single-driver effects need not add to total change. Report the remainder explicitly; it includes joint response and any changes in omitted drivers, so call it an interaction remainder only when all changing predictors are included.
- Preserve source data, existing APIs and `inst` cartographic assets. Remove only demonstrably unused implementation or obsolete active workflow code. Archive historical evidence rather than treating it as the current benchmark.

## Advantage experiment

Generate independent known-response datasets, with heterogeneous spatial responses, nonlinear response, joint effects, null drivers, correlated predictors and out-of-range futures. Prespecify 20 independently generated datasets per scenario, 1,500 training observations and a 40 x 40 prediction grid. Compare local actual-future replacement with a global average response contrast and a uniform mean-shift map using the same fitted model and target variable. Also compare with the exact known response surface. Score map RMSE, signed-effect agreement, hotspot recovery and predictive loss. Show all scenarios and failures; superiority is a result to measure, not an assumption.

## Reliability experiment

For each independent dataset, refit under a second sampling realization and a second forest seed; evaluate map correlation, error and driver-ranking agreement separately. Compare RF and GAM on the same inputs without changing the target. Record extrapolation coverage and accuracy inside/outside training ranges. Check table/raster equivalence, block-size invariance, static-driver zero effects, additive cancellation and weighted-model consistency with independent numerical oracles. Pixels are not independent experimental replicates.

## Figures

R only. Figure 1: known-truth spatial plate (truth / local scenario / global response / error), with one prespecified first replicate rather than a selected example. Figure 2: paired method-advantage results across independent datasets, including joint and out-of-range scenarios. Figure 3: resampling, model-seed and model-family reliability, plus support/error diagnostics. Use Arial, 183 mm width, vector PDF/SVG and 600 dpi TIFF/PNG, source CSVs and final-size visual QA. Single-species figures must use package plot methods.

## Case and delivery

Run Schizopygopsis (genus-level SchiSOFT collection, explicitly labelled) through occurrence QC, thinning, background generation, selection, four engines, nested spatial CV, held-out evaluation, full-data refit, current projection, all 18 configured futures, factor effects, scenario attribution, source plotting and ODMAP reporting. Record limitations of presence/background suitability and national background. Deliver output manifests, all configured scenario/engine statuses, source hashes and session information. Complete package tests and a short R CMD check, then commit and push the reviewed source to the existing GitHub remote.
