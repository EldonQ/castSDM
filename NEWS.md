# castSDM 0.8.1

## Breaking changes

* `cast_select(method = "cpi")` gains a **spatial-block inference layer**,
  now the default (`inference = "block"`): per-observation log-loss
  differences are averaged within a fine spatial grid and the one-sided
  t-test runs on the block means (`df = n_blocks - 1`, cluster-robust), so
  spatially autocorrelated observations are no longer treated as
  independent *and* the fold-level test's power loss is recovered. Requires
  `lon`/`lat` columns; falls back to `inference = "fold"` (with a warning)
  when coordinates are absent. New argument `n_blocks` controls the block
  count (auto ~ `min(50, max(20, n/20))`). `"fold"` and `"observation"`
  remain for comparison. CPI significance values will differ from 0.8.0.

# castSDM 0.8.0

Design-review release (3 critical / 19 major findings addressed, each with
regression tests): spatially honest CPI inference, contiguous default
spatial blocks, segfault-free raster batching, and consistent NA / weight
semantics across all prediction paths.

## Breaking changes

* `cast_select(method = "cpi")` now runs an in-package CPI estimator
  (`.cast_cpi_core`) with **fold-aggregated inference by default**
  (`inference = "fold"`): one log-loss difference per cross-fitting fold is
  tested (df = folds - 1) instead of treating n spatially autocorrelated
  observations as i.i.d., so p-values and FDR control stay calibrated on
  spatial data. The old per-observation t-test remains available as
  `inference = "observation"`. Cross-fitting folds are now stratified by
  the response, and rare species get automatic fold reduction with
  actionable error messages instead of ranger crashes. The `cpi` package
  is no longer required; `knockoff` (already used for the knockoff matrix)
  replaces it in Suggests. CPI significance values will differ from 0.7.0.
* `block_method = "grid"` in `cast_cv()` / `make_spatial_folds()` /
  `cast_prepare()` now builds **spatially contiguous blocks** (connected
  grid cells) instead of count-balanced random packing, and gains a
  `buffer` argument that excludes training points within `buffer` distance
  of held-out points. The legacy behaviour is available as
  `block_method = "grid_random"`. CV metrics, ensemble weights, and
  consensus frequencies will shift relative to 0.7.0 (generally downward,
  removing spatial leakage).
* AUC is now computed with a fixed `direction = "<"`: a model that ranks
  presences *below* backgrounds reports AUC < 0.5 instead of being flipped
  to 1 - AUC. Ensemble composite scores (`2*AUC-1 + TSS + CBI`) can now be
  negative, so reversed models are correctly down-weighted.
* Default `dml_folds` / `select_dml_folds` is now **10** (was 5) in
  `cast_select()`, `cast()`, and `cast_batch()`, and the YAML template
  follows suit.
* `cast_ensemble()` / `cast_ensemble_raster()`: zero-weight models
  (composite score < 0.5) no longer count towards the cross-model
  `hss_sd`; with fewer than two contributing models `hss_sd` is `NA`
  instead of a misleading zero or single-model value.
* `cast_rep()` prediction surfaces are now aggregated per model as
  `hss_mean_<model>` / `hss_sd_<model>` column pairs; the previous
  first-model-only `hss_mean` / `hss_sd` columns are gone.
* `cast_batch()` per-species seeds are now derived from `seed` plus a
  stable hash of the species name, so a species gets the same seed
  regardless of which subset of the batch it runs in. Batch results are
  bit-identical between interrupted+resumed and uninterrupted runs, but
  differ from 0.7.0 runs at the same `seed`.
* `cast_background()` validates study-area/raster geometry consistency
  (`terra::compareGeom`) and checks all layers for NA when defining valid
  sampling cells; additional occurrence columns are preserved in the
  output. Misaligned masks now error instead of silently mis-indexing.

## Bug fixes

* `cast_batch()` no longer segfaults when `raster_stack` / `future_rasters`
  carry in-memory `SpatRaster` objects: cache signatures hash a metadata
  fingerprint (`.cast_cfg_fingerprint`) instead of serialising external
  pointers, and rasters cross to PSOCK workers via `terra::wrap()` /
  `unwrap()`. Step-cache signatures now include `env_data`, so changing
  the prediction grid no longer replays stale predictions (M13).
* `cast_batch()` result files are written atomically (temp file + rename),
  and resume validates the RDS instead of trusting `file.exists`, so a
  truncated write can no longer mark a species permanently done.
  `cast_batch(parallel = TRUE)` now actually sets a `future` plan (and
  restores the previous plan and `CASTSDM_ROOT` on exit).
* `cast_predict_tiled()` writes `NA` back to cells with NA covariates
  instead of fabricating suitability via training-median imputation
  (matching `cast_ensemble_raster()`), restores the user's `future` plan
  on exit, and sizes workers from the intra-model thread budget.
* `cast_ensemble_raster()` renormalises the weights of contributing
  models per block after excluding non-finite predictions, fixing a
  systematic downward compression of HSS towards zero.
* `cast_project()` validates current/future environmental row alignment
  (`.cast_align_future_env`) before differencing — length or order
  mismatches now abort instead of silently recycling; each scenario runs
  under its own `tryCatch`; `.save_prediction_tif()` honours the input
  raster's `res`/`crs` instead of hard-coding 0.05°/EPSG:4326.
* `print.cast_batch()` displays the metric columns again (field-name drift
  between assembly `auc/tss/cbi` and printing `auc_mean/...` fixed).
* Response columns are validated as 0/1 at all entry points
  (`.cast_check_response`) and non-numeric predictors are rejected
  uniformly across fit/evaluate/predict/CV paths
  (`.cast_check_numeric_predictors`) instead of silent `as.numeric()`
  coercion.
* MESS computation (`.cast_mess`) is vectorised (`findInterval`/ecdf), so
  raster-scale extrapolation flagging is now feasible.

## New features

* `cast_select(method = "cpi")` gains `knockoff_reps` (default 3): the
  knockoff matrix is redrawn per replicate and per-variable statistics are
  replicate medians, stabilising Monte-Carlo p-values.
* `cast_ensemble_raster()` gains `clamp`, `extrapolation`, and
  `max_memory_mb` arguments and writes a block-computed `mess.tif`
  extrapolation layer alongside `hss.tif` / `hss_sd.tif`, bringing the
  raster path up to the documented extrapolation-control behaviour.
* DML selector runs `n_rep = 3` repetitions (median reported), following
  the DoubleML recommendation instead of relying on a single fit.

# castSDM 0.7.0

## OpenCodeReview hardening pass

Full-repository scan (alibaba/open-code-review) + fix round:

* `cast_predict_tiled()`: fixed a critical double-skip bug (the NA prototype
  written at init caused the write-back loop to skip prediction on first
  runs, leaving all-NA outputs) and a Windows file-lock issue (final write
  no longer re-opens the output file).
* `cast_esm()`: the validation split now requires >= 4 observations per
  class; below that the rare-species fallback fits on all rows and weights
  by in-sample AUC (`val_auc = NA`).
* `cast_run_step()`: errors raised inside the step expression now propagate
  instead of being swallowed by the peakRAM wrapper.
* `cast_batch()`: per-species errors inside PSOCK workers no longer abort
  the whole batch; the metrics table now normalises columns when a batch
  mixes CV and hold-out species (rbind mismatch fixed).
* `cast_ensemble_raster()`: cross-model SD now counts contributors per
  block, the mask geometry is validated before use, and empty
  model intersections abort with a clear error.
* `cast_project_raster()`: per-scenario statistics are now computed also
  when the change raster already exists (`overwrite = FALSE`).
* `make_spatial_folds()`: degenerate-safe binning (constant coordinates,
  fewer distinct values than folds, kmeans guard).
* `cast_vif()`, `cast_fit()`: non-finite / non-numeric predictors are
  rejected with targeted errors instead of silent coercion.
* `cast_background()`: the background set is topped back up when cells are
  NA in layers other than the first.
* `cast_report_odmap()`: null-safe screen diagnostics.
* `cast_study_area()`: robust CRS handling (unconditional idempotent
  transform; non-EPSG representations supported).
* `plot.cast_screen_comparison()`: named retention counts (fixed NA
  subscripting).
* `.cast_digest()`: the no-`digest` fallback no longer overflows 32-bit
  integers.
* `bioclim-io`: exact model-token matching and deterministic (sorted)
  recursive file matching.
* Deduplicated the N-SDM score/weight computation into shared helpers
  (`.cast_ensemble_scores()`, `.cast_ensemble_weights()`).

## Consensus selection across spatial folds

* `cast_cv()` now stores `selection_freq`: each predictor's fold-level
  selection frequency, the package's spatial-stability diagnostic.
* New `cast_consensus()` aggregates the fold selections into a consensus
  variable set (predictors retained in at least `threshold` of folds,
  default 0.5) and returns a `cast_select` object that plugs directly into
  `cast_fit()`. Validated on real China fish data: consensus sets are
  ~2x more stable across CV reruns than raw fold selections (Jaccard 0.66
  vs 0.30) and transfer as well as the all-predictor model with a third
  of the variables.

## Replicate design, uncertainty layers, ODMAP reporting, rare-species routing

* `cast_rep()` repeats background sampling `n_reps` times and aggregates
  evaluation metrics (mean +/- SD), selection frequencies, and per-cell
  prediction mean/SD across replicates.
* `cast_ensemble()` now returns a cross-model `hss_sd` uncertainty column;
  `cast_ensemble_raster()` additionally writes an `hss_sd.tif` raster.
  `cast_counterfactual()` returns `delta_sd` (cross-model SD of the
  change) and a `mean_abs_delta_sd` summary.
* `cast_report_odmap()` renders an ODMAP-aligned (Zurell et al. 2020)
  Markdown report of the model, data, assessment, and prediction settings.
* `cast_batch()` gains occurrence-count gates: species with fewer than
  `esm_min` presences are skipped with a warning; species between
  `esm_min` and `min_occ` are routed to the ESM pipeline automatically
  (`esm_used = TRUE` in the result).
* `cast()` and `cast_batch()` expose `select_dml_folds` and `num_threads`;
  the YAML template covers the full current argument set.

## Defect fixes and engineering hardening

* `plot.cast_select()` labels the DML axis `|DML statistic|` instead of the
  wrong `|CPI statistic|`; fold maps extend their palette beyond 10 folds.
* `cast_project()` marks cells with non-binary predictions as `NA` instead of
  an empty-string change class; summary counts ignore them explicitly.
* `cast_vif()` rejects non-finite predictors with a targeted error and fails
  loudly instead of looping when VIF cannot be computed.
* `cast_ensemble()` excludes models with non-finite predictions (with a
  warning) and renormalises weights, instead of silently zeroing them.
* Removed the dead `threshold_method` argument from `cast_ensemble()`,
  `cast_project()`, `cast_ensemble_raster()`, `cast_project_raster()`, and
  `cast_save_future_projection()`; binary thresholds always follow the
  max-TSS rule.
* `load_basemap()` restores the user's global `sf` s2 setting after reading.
* `cast_cv()` checks for `pROC` with a friendly error and isolates per-fold
  metric failures.
* `cast_select()` reports the CPI/DML 200-tree cap instead of truncating
  silently; the batch pipeline falls back to `"cpi"` (not `"dml"`).
* `cast_predict()` errors with the missing column names when a prediction
  grid lacks fitted predictors.
* CPI confidence intervals are clamped at zero (CPI is non-negative).
* Batch resource logging writes one CSV per species and merges them after
  the run, removing the parallel-append race.
* Step-level checkpoints are keyed on a run signature (data, config, seed,
  models) and tile-level checkpoints on a per-fit fingerprint, so stale
  caches are never replayed after a configuration or data change;
  `overwrite = TRUE` now refreshes tile checkpoints too.
* `cast()` and `cast_batch()` refit final models on the full data set
  before spatial prediction (`refit_full = TRUE`), reuse every record for
  published maps, and report explicitly when the ensemble is skipped.
* Selection scores now flag `fallback` (kept through the `min_vars` floor)
  and `forced` (kept through `force_include`) predictors; the RF baseline
  in `cast_screen_comparison()` is budget-matched to CPI; `cast_esm()`
  accepts a `screen` so the rare-species route can use conditional
  selection instead of univariate ranking.
* Documentation alignment: README, package overview, and NEWS now describe
  the CPI default with the DML effect-reporting layer; N-SDM is cited as
  Adde et al. (2020).

# castSDM 0.6.0

## New core

* `cast_screen_comparison()` contrasts the castSDM conditional predictive-impact
  screen against the associational baselines researchers commonly use
  (correlation filter, stepwise VIF, univariate marginal screen, and
  random-forest permutation importance) on the same shared candidate pool, with
  a `plot()` tile matrix that makes the conditional-vs-associational contrast
  explicit for the Results narrative.
* `cast_select(method = "dml")` is the new default selector. It fits a
  double machine-learning partially linear model (`DoubleML` + `mlr3`
  random-forest nuisance) per candidate predictor and keeps those whose
  Neyman-orthogonal effect survives Benjamini-Hochberg FDR control. The FDR
  level (`alpha`) is the only ecological tuning choice.
* `method = "rf"` is retained as a conventional permutation-importance
  benchmark for comparison.
* Added `cast_effect()`: a tidy DML effect table with FDR-adjusted
  significance and confidence intervals, plus a `plot()` coefficient (forest)
  plot.
* Added `cast_counterfactual()`: a g-computation what-if on the current
  climate that shifts one predictor, holds the rest fixed, and maps the change
  in habitat suitability, plus a diverging-map `plot()` method.
* `cast_cv()` re-runs selection inside every outer training fold and stores
  fold-specific selections, preventing feature-selection leakage.

## Removed

* Removed the `"stable"` selector and its invariance machinery (data-derived
  environments, binomial-GLM invariance tests, greedy minimum stable-set
  search), the `stable_no_invariance` ablation, `cast_domains()`, and the
  `glmnet` dependency. Variable choice is now the single FDR-controlled DML
  criterion.

## Existing workflow

* Retains data preparation, background sampling, RF/BRT/MaxEnt/GAM and ESM,
  tiled raster prediction, ensembles, CMIP6 projection, batch checkpointing,
  YAML configuration, and plotting.

## Plotting and stability

* China maps always show the complete, standard national frame plus the South
  China Sea inset; a range-restricted species is never cropped to its coloured
  cells. Seamless colour is achieved instead by filling the whole basemap with
  each map's null/background colour, so cells outside the accessible area
  continue the surface rather than tearing it at an artificial edge.
* Suitability maps (`plot.cast_predict()`, `plot.cast_ensemble()`) use a pale
  grey to deep teal-blue sequential ramp whose low anchor equals the basemap
  fill; range-change maps (`plot.cast_project()`) use a refined gain/loss/
  stable-present palette over a neutral `stable_absent` background; what-if maps
  (`plot.cast_counterfactual()`) diverge around a neutral grey midpoint.
* `plot.cast_select()` gains a `top` argument and now colours predictors by
  ecological class with readable labels; `plot.cast_cv()` fold maps share the
  full-China frame and carry a spatial-block narrative.
* Ensemble/counterfactual raster prediction is chunked by rows to bound peak
  memory on national 1 km grids, preventing out-of-memory segfaults.

# castSDM 0.2.0--0.5.1

Earlier experimental releases explored DAG, role-prior, CATE, RF-only, and
stable/invariance screening designs. These interfaces are retired in 0.6.0.
