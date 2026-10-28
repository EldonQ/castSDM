# castSDM 0.23.0

- New `cast_factorial_attribution()`: a 2 x 2 factorial causal scenario
  attribution for climate and land use. It computes climate and land-use main
  effects, the joint change, the interaction and the two conditional effects
  from four aligned current/climate-only/LULC-only/joint suitability surfaces.
  Continuous and binary range summaries, optional support and cross-model SD
  layers, a contrast GeoTIFF and an explicit claim contract are returned in a
  `cast_factorial` object.
- The factorial design fixes the other scenario component, so no DAG or
  adjustment set is required. The contrasts are model-based causal scenario
  attributions under stated assumptions, not identified real-world causal
  effects.
- `cast_attribution()` remains the within-scenario driver-only/retained
  decomposition; the two APIs answer different questions.
- Package title and description now foreground factorial causal scenario
  attribution while retaining the predictive engine and attribution contracts.
- `cast_ensemble_raster()` reports block progress at quartiles when
  `verbose = TRUE`, so long national projections are no longer silent between
  the grid message and the final raster.

# castSDM 0.22.2

- `cast_project_raster()` projection statistics now name centroid coordinates
  `*_centroid_x` / `*_centroid_y`, matching their projected CRS units instead
  of the misleading `lon` / `lat` labels.

# castSDM 0.22.1

- `cast_attribution()` now reports a `remainder` row in `summary` (total change
  minus the sum of driver-only effects), matching the existing `remainder`
  raster layer and documenting non-additivity explicitly.
- Attribution progress remains available through `verbose`; long-running
  callers should checkpoint driver subsets rather than relying on a single
  monolithic surface.

# castSDM 0.22.0

- `cast_attribution()` gains `clamp`; set it to `TRUE` to apply the same
  training-range clamp as `cast_project_raster()`.
- Projection statistics gain `n_current_present_total`, `n_unpaired_current`
  and `pct_turnover`; `pct_change` now uses the paired current range, matching
  the raster implementation.
- The screen message reports both the post-filter candidate count and the
  final elastic-net/ensure selection count.

# castSDM 0.21.0

- Variable selection: `method = "screen"` is an independent elastic-net
  screen implemented in the package (variance/unique-value filter, optional
  construct de-duplication, correlation pruning, `glmnet` ranking, and a
  `ncov`/`maxncov` cap). It also supports `keep` and `ensure` coverage. No
  external selection library is vendored.
- Engine set: RF, MaxEnt and GAM. The previous boosting engine is no longer
  available. GAM is fitted unweighted and calibrated with the same stored
  logit-scale prevalence offset as MaxEnt.
- `cast_select()` diagnostics record `final_ranks`, `coefficients`, `lambda`
  and `cap_binding`; forced predictors are never reported as
  collinear-thinned.

# castSDM 0.20.0

Breaking release: one production pipeline with an explicit evaluation protocol.

- Hold-out diagnostics (2026-10-20): `cast_prepare()` returns `split_type`
  (`spatial`, `stratified`, `random`) alongside the descriptive `split`
  string, so callers can distinguish a spatial hold-out from a fallback.

- Raster maintenance (2026-10-19): ensemble, single-shift maps, tiled prediction
  and paired attribution share row-block I/O and recoverable file publication.
  Full-grid output accumulators are removed. Mask geometry/read failures abort;
  masks now require one aligned layer and finite non-zero inclusion values.
- Tiled checkpoints fingerprint input content and callback dependencies, verify
  stored prediction digests and require matching completion receipts before
  skipping outputs. Prediction errors stop without publishing partial rasters.
  Predictor names `x`/`y` are preserved rather than mistaken for coordinates.
- Thinning puts unknown years last; a sole valid geographic record is sampled
  correctly. Recency alone is not evidence of coordinate precision.
- Effect maps are backed by streamed double-precision GeoTIFFs, including when
  no filename is supplied. Existing filenames require explicit overwrite.

- New `cast_protocol()` and `cast_run()` replace `cast()`. The protocol is the
  single declaration point for engines, selection, CV design, thresholds,
  prevalence calibration, ensemble policy and extrapolation handling.
  `cast_run()` echoes the effective values, and with `output_dir` writes
  per-stage artifacts with content-digest receipts into `manifest.rds`
  (stages/, tables/, odmap_report.md).
- Evaluation protocol: selection, thresholds and ensemble weights are
  estimated on the training split; the hold-out split is evaluated exactly
  once with pooled out-of-fold thresholds (`cast_evaluate(thresholds = )`);
  published predictions use the full-data refit; the ensemble threshold is
  optimised on the out-of-fold ensemble surface that is actually thresholded.
  Per-fold TSS/SEDI/kappa now use the pooled out-of-fold threshold instead of
  each fold's in-sample training cut, which had collapsed the held-out TSS of
  flexible learners (RF).
- AOA threshold follows Meyer & Pebesma (2021): 75th percentile + 1.5 x IQR of
  the pooled training-side nearest-neighbour DI; held-out DI and its share
  above the threshold are validation evidence (`share_outside_oof`).
- Prevalence calibration: with `prevalence_target`, every engine -- including
  MaxEnt, whose fit accepts no row weights -- is calibrated to the same
  training-sample prevalence through a stored logit-scale offset, so ensemble
  members share one scale. `cast_protocol()` defaults to 0.5.
- Spatial blocking: `block_method = "blockcv"` (default) builds square blocks
  with blockCV, sized from the spatial-autocorrelation range
  (`block_size` overrides); dependency-free `"grid"` fallback retained,
  `"env"` kept. `grid_random` and `cluster` removed.
- RNG hygiene: every stage seeds locally (`withr`) and restores the caller's
  RNG state; global `set.seed()` side effects removed; `ranger` prediction no
  longer advances the global RNG.
- Removed dead code: `cast()`, `cast_load_future_envs()`,
  `cast_save_future_projection()`, `grid_random`, `cluster` blocking.

# castSDM 0.19.0

- Accept external fitted models through `cast_as_fit()` with explicit prediction
  callbacks and training references; accept external `fold_ids` in `cast_cv()`.
- Preserve explicit RF mtry during tuning and refitting; forward engine settings
  through spatial CV and validate integer parameter ranges before fitting.
- Reject missing/one-class fitting responses and entirely nonfinite predictors;
  report imputation counts without claiming observed environmental support.
- Required pipeline CV/refit/ensemble stages now stop on failure. GAM smooth
  failures warn when a linear-only fallback is used.
- Reduce the maintained case report to four main figures; single-shift audits
  and ensemble permutation are explicit optional diagnostics. Add nearest-neighbour
  display resampling for binary/categorical maps. Retain `inst` assets and examples.
- Add paired current/future spatial attribution with frozen model weights,
  driver-only and retained-driver contrasts, range flags and model disagreement.
- Add source-native attribution plots and share cartographic rendering with
  the single-species workflow. Preserve all cartographic assets.
- Correct unavailable-CV label overlap and binary range legends; the case
  replot entry records original computation and revised plotting sources separately.
- Allow single-shift effects to use the same model weights as projection.
- Reframe the package around environmental drivers of habitat change.
- Allow independent future raster scenarios to use the caller's future plan,
  with explicit raster transport and a per-worker memory budget. Reuse exactly
  equivalent driver comparisons when only one predictor changes.
- Uniform background sampling now proposes mask cells without replacement and
  accepts only complete finite environments, avoiding a full multi-layer scan.
  The design is unchanged but seeded sample sequences differ from earlier versions.
- Single-shift maps evaluate only rows retained by their existing hard range
  mask, avoiding predictions that would immediately be discarded.

# castSDM 0.18.0 (unreleased)

## Fixed

### Silent-wrong-result class (the audit's highest-severity findings)

* `cast_cv()` no longer swallows a per-fold prediction failure. A model that
  failed to predict in a fold produced an all-`NA` vector with **no warning**,
  and that surface flowed into the out-of-fold predictions and thence into the
  ensemble threshold — a model that never predicted looked like a model that
  predicted poorly. Prediction, training-set prediction and threshold
  calibration failures inside a fold now each warn with the fold index, the
  engine name and the underlying message.
* `cast_cv()` gained a **fold-coverage floor** on the threshold guard. The old
  guard only checked that at least 10 out-of-fold values existed in total, so a
  surface that was complete *within the folds that worked* but missing every
  fold where the model failed still produced a "valid" threshold. At least 80%
  of folds (never fewer than 2) must now produce predictions; otherwise a
  warning names the shortfall. The per-model count is returned as
  `cv$pred_folds`.
* `cast_ensemble()` applies the same coverage test before thresholding the
  ensemble mean, and falls back to the weighted mean of per-model thresholds
  (with a warning) instead of deriving an operating point from a partial
  surface.
* `cast_ensemble_importance()` now records which models fail to predict, and
  **refuses to score** when a model fails on the unpermuted baseline. It used
  to let `.ensemble_rowmean()` silently reallocate that model's weight, so a
  "4-model ensemble" was scored on 3 models and still reported as 4.
  Permutations are additionally restricted to the same model set as the
  baseline so the comparison keeps a fixed denominator.

### Correctness and robustness

* `cast_predict_tiled()` re-read the entire predictor stack from disk once per
  model, because the model loop wrapped the tile loop. The loops are now
  ordered tile-outer / model-inner, so each tile is read exactly once
  (~4x less I/O for a four-model ensemble) while still holding one tile at a
  time in memory.
* `cast_predict_tiled()` no longer writes an all-`NA` prototype GeoTIFF before
  the sweep. An interrupted run used to leave that prototype behind, and the
  next `overwrite = FALSE` run would report "exists; skipping" and hand back a
  raster of pure `NA`.
* Raster output is staged in a temporary sibling file before publication
  (`.cast_write_raster_atomic()`) across `cast_predict_tiled()`,
  `cast_ensemble_raster()`, `cast_effect_map()`, `cast_project()` and
  `cast_project_raster()`. Resume logic keys off `file.exists()`, which
  previously meant a half-written GeoTIFF was accepted as complete. Rename
  is the normal publication path; the copy fallback is not atomic, so this
  is not a guarantee for every filesystem or interruption condition.
* `cast_fit()` and `cast_evaluate()` now reject a non-`data.frame` input.
  `cast_evaluate(fit, NULL)` used to pass validation silently (`NULL[["x"]]`
  is `NULL`, `all(logical(0) %in% c(0, 1))` is `TRUE`) and return all-`NA`
  metrics as though the evaluation had merely been inconclusive.
* `drivers = character(0)` and `models = character(0)` no longer silently
  expand to "every fitted predictor"/"every model". These used to be rewritten
  by the package's own `%||%`, whose zero-length-is-absent rule made the
  subsequent "must be non-empty" check unreachable. Both now use `is.null()`,
  so an explicit empty vector is reported as the error it is.

### Argument validation

Documented ranges that were not enforced at runtime now abort instead of
producing a wrong or empty result:

* `train_fraction` in `cast()` / `cast_prepare()` must lie in the open
  interval (0, 1). `0` built an empty training set and the whole pipeline ran
  on it; `1` produced all-`NA` evaluation metrics.
* `buffer` in `cast_cv()` must be finite. `Inf` passed the non-negative test
  and then excluded every training row.
* `tile_size` in `cast_predict_tiled()`, `buffer_km` / `min_cells` in
  `cast_study_area()`, `n_blocks` in `cast_prepare()`, `threshold` in
  `cast_vif()`, `num_threads` in `cast_fit()`, and `prefix` / `fig_dpi` /
  `min_score` / `min_metric_folds` in `cast_save_future_projection()`.

### Filename safety

Scenario names and output prefixes are pasted into output filenames. A
scenario named `../../etc/cron.d/x` escaped the output directory. Both are now
sanitised by `.cast_safe_name()` in `cast_project()`, `cast_project_raster()`,
`cast_ensemble_raster()` and `cast_save_future_projection()`, with a warning
naming every scenario that was changed and an error on post-sanitisation
collisions.

## Performance

* The cross-model SD and the `method = "median"` combination were computed with
  `apply(mat, 1, <closure>)` — one R-level callback per cell. At product scale
  (~3e7 cells) that is tens of millions of interpreted calls, which dominated
  the cost of the ensemble itself. Both are now vectorised
  (`.row_sd_finite()`, `.ensemble_rowmedian()`), using `matrixStats` when
  available. Measured on a 2e6 x 4 block: **0.27 s versus 19.1 s (71x)**,
  with values identical to the previous implementation.
* `matrixStats` is a new (optional) Suggests dependency; the package falls back
  to the original `apply()` form when it is absent, so results are unchanged
  either way.

## Internal structure

* `R/engines.R` (new) holds a learner registry: engine key -> backing package,
  label, plot colour, and the mapping from `tune$best` fields to the
  `cast_fit()` formals. Adding a fifth learner previously required edits in
  **four files and eleven dispatch sites**; a forgotten site silently produced
  an `NA` model or an invisible plot colour rather than an error.
* `R/ensemble.R` (1399 lines) was split into `ensemble.R` (core: scoring,
  weights, thresholds, aggregation, grouped and permutation entry points) and
  `ensemble-raster.R` (block engine and GeoTIFF output). The raster half is now
  the only ensemble module performing low-level `terra` I/O.
* TSS had three independent implementations that had already drifted in their
  NA handling; they now share `.cast_tss()`. The two AUC implementations are
  deliberately kept separate — `compute_auc()` serves the selector's inner
  loop and `pROC` serves evaluation's documented `direction = "<"` fix — but
  are now verified equal and documented as such.
* `cast_ensemble_by()`'s example referenced a `"glm"` engine that does not
  exist; corrected to the four supported engines.

## Infrastructure

* GitHub Actions are pinned to commit SHAs instead of floating tags (a tag can
  be re-pointed upstream). `permissions: read-all` was already in place.
* Version floors added for `cli (>= 3.6.0)`, `ggplot2 (>= 3.4.0)`,
  `sf (>= 1.0.0)` and `terra (>= 1.7.0)`, whose APIs this package depends on.

## Fixed (previous release notes retained below)

* `cast_ensemble()` / `cast_ensemble_raster()` no longer mislabel a
  weighted-min_score fallback as a deliberate fit. The fallback branch used
  `fallback_used <<- TRUE` inside the `switch()` block, which skips the
  function's local frame in R and stored the flags in the enclosing
  environment instead: the weights were correct, but `weights_fallback`
  stayed `FALSE` and `fallback_reason` stayed `NULL`, hiding the fact that
  the ensemble had collapsed to a single model (found in a production case
  where all four engine composite scores were negative and only GAM was
  silently retained). Plain assignment now records the fallback and the
  `cli_warn` emitted at fallback time is unchanged.

## AOA applicability, grouped/decay ensembles, ensemble permutation importance, kappa/omission/MPA metrics, environmental blocking (A2 deep-review round)

* **Area of Applicability (AOA; Meyer & Pebesma 2021)**: `cast_cv()` gains
  `aoa = TRUE`, which calibrates the dissimilarity-index (DI) threshold from
  the pooled held-out-fold DI values (95th percentile) in the standardized
  predictor space (model-independent, uniform weights; requires `FNN`).
  The calibration is stored in the `aoa` component of the `cast_cv` object
  together with the per-row held-out DI (`di_oof`).
* `cast_predict()` gains `cv`, `aoa` and `vi`: with `aoa = TRUE` each
  prediction row is flagged against the calibrated threshold (`aoa` logical
  column) and receives its raw dissimilarity index (`aoa_di`). The DI is
  computed on the unclamped input (clamping never softens the flag, same
  convention as MESS). Optional `vi` weights the standardized dimensions;
  a warning is issued when the supplied `cv` threshold was calibrated
  uniform.
* `cast_ensemble_raster()` gains `aoa_cv`: writes `<prefix>_aoa.tif`
  (flag) and `<prefix>_aoa_di.tif` (raw DI) block by block, on the
  unclamped block input.
* **Score-power ensemble weighting**: `cast_ensemble()` and
  `cast_ensemble_raster()` gain `decay` (default 1). With
  `method = "weighted"`, weights become `score^decay` (negative scores
  clamped at 0) renormalised; `decay = 0` degenerates to equal weights and
  `decay > 1` sharpens towards the best-scoring models. Raw scores are
  still reported in `model_scores` and `min_score` filters on the raw
  scores. The used exponent is stored in the `decay` field of the
  `cast_ensemble` object.
* `cast_ensemble_by()`: builds one ensemble per caller-specified group of
  models (the explicit castSDM analogue of biomod2's `em.by` grouping),
  returning a `cast_ensemble_grouped` object whose `predictions` carry
  `hss_<group>`, `hss_sd_<group>` and `binary_<group>` columns; groups may
  overlap. Print method included.
* `cast_ensemble_importance()`: permutation importance on the ensemble
  layer - each predictor is permuted `n_perm` times (default 25), the
  ensemble is recomputed with frozen per-model weights, and importance is
  `1 - cor(original, permuted)` averaged over permutations (Pearson or
  Spearman). Returns a `cast_ensemble_importance` object sorted descending.
* **Threshold-free evaluation metrics**: `evaluate_model_full()` and the
  whole CV chain (`cast_cv()` fold metrics, `cast_evaluate()`) now also
  report `kappa` (Cohen's kappa at the fold threshold, paired with TSS),
  and the threshold-free `omission_5` / `omission_10` (proportion of
  presence records below the 1-E quantile of all predictions; E = 5% / 10%)
  and `mpa` (maximum predicted area at 90% training sensitivity; Pearson
  2007-style). Aggregated columns `*_mean`/`*_sd`/`*_n_folds` follow the
  existing convention.
* **Environmental blocking**: `block_method = "env"` in `cast_cv()` and
  `cast_prepare()` - k-means partitions the scaled numeric predictor
  matrix (medians imputed, zero-variance columns dropped), so folds are
  contiguous in environmental space rather than geographic space
  (Roberts et al. 2017-style). Split labels distinguish
  `"environmental block (env)"` from spatial blocking.

# castSDM 0.17.0 (unreleased)

## Change-map coding, range-change statistics, binning and buffer units (P2 fast-win round)

* `cast_project()` and `cast_project_raster()` gain `coding`: change maps
  can be written with the package's internal encoding (`"cast"`: gain=1,
  loss=-1, stable_present=2, stable_absent=0) or with biomod2's
  `bm_RangeSize()` convention (`"biomod2"`: Gain=1, Loss=-2,
  Stable_Pres=-1, Stable_Abs=0) so saved GeoTIFFs feed
  BIOMOD_RangeSize visualisation scripts without remapping.
* Projection statistics gain `n_current_range` / `n_future_range` pixel
  counts and the biomod2-convention `pct_loss` (loss over the current
  range) and `pct_gain` (gain over the future range), in both
  `cast_project()` and `cast_project_raster()`.
* `cast_background(strategy = "environmental")` gains `bin_method`
  (`"pca2d"` PC1 x PC2 grid as before, or `"kmeans"` partitioning of the
  scaled environment space, the sdm convention) and `bin_k` (bin/cluster
  count, default 20).
* `cast_background(cell_thin = TRUE)` now warns when occurrences fall
  outside the raster or on NA-valued cells, reporting the lost share
  (previously they silently merged into the generic duplicate count).
* `cast_cv()` gains `buffer_unit`: the spatial exclusion band can be
  `"deg"` (planar degrees, as before) or `"km"` (great-circle Haversine
  kilometres). `"km"` requires decimal-degree coordinates and aborts on
  projected inputs.
* `print.cast_effect_table()` reports the masking rate (number of masked
  drivers and the lowest observed support) so the `min_support` default's
  effect on correlated drivers is visible.
* Documentation: the ensemble Score is described accurately as the mean of
  the finite components instead of a fixed /3 formula; the fixed 101-bin CBI
  is flagged as an internal diagnostic
  only, not comparable with ecospat-derived Boyce; `cast_vif()` documents
  the fewer-than-five-variables stopping guard.

# castSDM 0.16.0 (unreleased)

## Prevalence weighting, SRE/disk background, thinning, response curves, replication (P1 review round)

* `cast_fit()` gains `prevalence_target`: presence/background rows are
  weighted (weight = target/n per class) so the fitted prevalence equals the
  target, the standard presence-background correction. RF accepts the
  weights; GAM and MaxEnt are fitted unweighted and calibrated with a stored
  logit-scale offset.
  Evaluation metrics remain unweighted. `cast_cv()` and `cast()` pass the
  target into every fold and the full-data refit.
* `cast_background()` gains two sampling strategies (Barbet-Massin et al.
  2012): `strategy = "sre"` draws cells whose environments fall inside the
  per-variable presence quantile envelope (`sre_quantile`, default 0.025),
  and `strategy = "disk"` draws cells in a distance band (`disk_min` /
  `disk_max`, map units) around presences. Both restrict the NA top-up pool
  to the same design; an empty pool aborts with a hint instead of falling
  back to random.
* New `cast_thin()`: greedy great-circle (Haversine) thinning of occurrence
  records to a minimum distance (`min_dist_km`), with a seeded random
  visiting order for reproducibility; rows are returned in original order.
* New `cast_response_curves()`: single-variable response curves with the
  other predictors fixed at training medians or marginalized as partial
  dependence (`type = "pdp"`, Friedman 2001), plus optional two-variable
  interaction grids (`pair`); tidy long-format output with a ggplot-based
  `plot()` method for both layouts.
* `cast_background()` gains `n_rep`: independent pseudo-absence replicate
  sets, each drawn under its own sub-seed (a named `cast_background_reps`
  list; `n_rep = 1` keeps the plain data.frame return).
* `cast_cv()` gains `n_repeat`: independent spatial-blocking repeats.
  Metrics are averaged within each repeat and reported with the
  **between-repeat** standard deviation (blocking randomness rather than
  fold noise); `fold_metrics` gains a `replicate` column; selection
  frequencies pool all repeats; out-of-fold predictions and thresholds
  refer to the first repeat. `cast()` passes `n_repeat` through.

# castSDM 0.15.0 (unreleased)

## Hyperparameter tuning, user background, reference fixes (P0 review round)

* `cast_fit()` gains `tune = TRUE`: a small per-engine grid search (RF `mtry`
  scored by out-of-bag TSS; MaxEnt feature classes x regularisation multiplier
  scored by stratified random-fold TSS), informed by the biomod2/flexsdm grids and
  trimmed for cost. Combinations that fail to fit score `NA`; failing grids
  keep the fixed defaults with a warning. `cast_cv(tune = TRUE)` runs the
  grid inside every outer training fold (nested, leak-free);
  `cast()` passes `tune` / `tune_folds` through and reuses the training-split
  hyperparameters for the full-data refit.
* `cast_fit()` exposes `maxent_classes`, `maxent_regmult` and `rf_mtry` for
  explicit hyperparameter control.
* `cast_background()` gains `user_table` (user-defined background from exact
  points) and `bias_raster` (target-group background: cells drawn with
  probability proportional to a bias surface; Phillips et al. 2009), closing
  the gap between the package's own reference list and its feature set.
* The moving-window Boyce index no longer drops ratio-0 windows: only
  undefined 0/0 windows are excluded, aligning with ecospat.boyce and
  removing an optimistic bias that propagated into ensemble weights.
* `predict_single_model()` passes `clamp` explicitly to maxnet: the MaxEnt
  engine no longer clamps predictions silently to the training range while
  the other engines extrapolate. Clamping is now solely the package-level
  `clamp` argument of `cast_predict()` / `cast_ensemble_raster()`.
* Aligned the ensemble-score documentation with the implemented
  >=2-components rule.
* Removed the dead internal helper `find_tss_threshold()`.

# castSDM 0.14.0 (unreleased)

## Model defaults follow the SDM literature; 0.12.0 stubs removed

* BREAKING: `cast_fit()` / `cast_cv()` default `rf_ntree` changed from 300 to
  500. Fitted models and cached CV outputs must be recomputed.
* `cast_select()` bounds its search cost with a response-stratified subsample
  (every presence row kept, only the background evenly thinned) instead of a
  plain systematic subsample, so a large background (e.g. ~25 presences in
  10,000 rows) no longer starves the inner folds of the minority class.
* BREAKING: the 0.12.0 defunct stubs are removed outright --
  `cast_sensitivity()`, `cast_dose_response()`, `cast_effect_support()`,
  `cast_necessity()` and their `plot()`/`print()` methods no longer exist.
* Documentation: `plot.cast_importance()` describes the forward-selection
  path (not permutation importance); ensemble fallbacks and empty-selection
  wording updated.

# castSDM 0.13.0 (unreleased)

## Stage 2 is now a spatial forward search

* BREAKING: `cast_select()` stage 2 no longer calibrates a shift-effect
  statistic against a permutation null. It now runs a forward search on the
  inner cross-validated loss (Brier by default, `"auc"` optional) of a
  probability random forest: spatial folds when coordinates exist, starting
  from the prespecified set (or the best pair), admitting a candidate only
  when its paired improvement exceeds a 2-SE guard, and stopping by itself.
  There is no predictor-count cap (`ncov`/`maxncov` removed), no p-values,
  and no fallback set: an empty selection honestly means nothing improved
  on the intercept-only model. `n_perm` and `shift_size` are removed;
  `cast()` gains `select_metric` / `select_tolerance` / `select_n_folds`.
* `scores` now carries the admission `step_added` (`0` for prespecified) and
  `loss_gain` instead of `interventional_effect` / `p_value` /
  `null_threshold` / `passed_null` / `fallback`; `diagnostics` carries the
  forward `path`, `inner_method`, `null_loss`, `final_loss` and `n_fits`
  instead of null-distribution summaries. `cast_importance()` reports the
  admission path. Existing selection caches and downstream model/CV outputs
  must be recomputed.

# castSDM 0.12.0 (unreleased)

## One shift, one map, hard masking

* BREAKING: `cast_effect_table()` and `cast_effect_map()` now take a single
  raw-unit shift per driver (`shift`, default `1` in `shift_type = "raw"`;
  `"sd"` converts through the stored training SDs) instead of a symmetric
  shift set. Outputs are `mean_abs_dHSS` / `mean_signed_dHSS` over
  supported rows or cells only, with `n`, `n_supported`, `support` and a
  `masked` flag. Units outside the training quantile box are masked (`NA`
  estimates), not ranked. `outside_range_fraction`, `pct_gain` / `pct_loss`,
  `median_signed_dHSS`, `max_gain` / `max_loss` and `n_shifts` are removed;
  existing effect outputs must be recomputed.
* BREAKING: `cast_dose_response()`, `cast_effect_support()`,
  `cast_necessity()` and `cast_sensitivity()` (including the
  `marginaleffects` backend) are removed; calls abort with a forwarding
  error. `cast_select(method = "tramicp")` and its `environment` /
  `icp_*` arguments are removed, as are the `select_environment` /
  `select_icp_*` arguments of `cast()` and `cast_cv()`.
* `cast_select()` stage 2 keeps the shift-effect statistic but calibrates it
  against a within-stratum permutation null: each survivor is shuffled only
  among rows with similar values on the other survivors (one k-means
  stratification per predictor), so the null respects the observed joint
  support. The second attribution column (`perm_importance`), `p_adjusted`
  and `importance_agreement` are removed; `cast_importance()` reports the
  single conditional-effect column. Existing selection caches and downstream
  model/CV outputs must be recomputed.

# castSDM 0.11.0 (unreleased)

## Optional invariant selection and scenario comparisons

* `cast_select(method = "tramicp", environment = ...)` delegates exhaustive
  binary-logistic invariant selection to the optional tramicp package. It uses
  supplied environment groups, not automatic spatial partitions, and preserves
  empty/no-accepted-set results without fallback selection. The output does not
  establish ecological causality or identify an adjustment set.
* `cast()` forwards the environment column and ICP settings to each training
  screen. `cast_cv()` protects training-data arguments and retains empty/failed
  fold statuses, including diagnostics on all-fold failure conditions.
  For invariant selection, `cast()` propagates all-fold failure rather than
  discarding its diagnostics and substituting hold-out evaluation.
* Retired selection arguments and unknown methods now raise errors instead
  of being ignored or replaced with `two_stage`. Removed an unused binning
  helper and unused internal arguments; AUC and TSS now share one ROC
  calculation, and GAM scenario comparisons avoid an unused data copy.
* `cast_sensitivity(backend = "marginaleffects", model = "gam")` provides
  signed rowwise scenario contrasts through marginaleffects without changing
  the native default. It reports suitability contrasts, not causal effects or
  calibrated occurrence probabilities; no standard errors are computed.

## Correct support evaluation and causal interpretation

* Single-predictor MaxEnt fits preserve data-frame dimensions and support the
  linear-only feature set without changing the default background augmentation.
  Nested spatial CV now reports fold-level fitting errors rather than silently
  discarding failed folds.
* `cast_select(keep = ...)` protects prespecified predictors during thinning,
  null screening and the complexity cap; they count toward `ncov` and cannot
  be silently dropped. `kept_by_design` and `selected_reason = "prespecified"`
  distinguish this decision from statistical evidence. `cast(select_keep = ...)`
  passes the same set to the training screen and every nested CV fold. This
  safeguard does not discover or verify a causal adjustment set.
* `cast_effect_map()` adds `support_<driver>` layers and `support_probs`.
  Each layer reports the fraction of requested shifts inside the training
  quantile box per cell. Attached-table support remains the minimum per-shift
  coverage across cells, not the mean of the new layer. Effects remain unmasked.
* Dose-response curves now average paired engine deltas, matching tables and
  maps when a model prediction is missing for only one member of a pair.
* Corrected the Hooker, Mentch and Zhou (2021) journal and DOI.
* Support fractions now evaluate `newdata`, not just the training population,
  and require both baseline and shifted rows to lie inside all predictors'
  training quantile bounds. Table and curve diagnostics use the same complete
  rows as their effect estimates; `max_rows` and `support_probs` are validated.
* `range_supported` replaces the misleading `estimable` curve column. Coverage
  of at least 0.5 is a descriptive range screen, not causal estimability.
  Missing coverage is displayed as unknown, not as supported.
* Correction to the 0.10.0 description below: a quantile box does not establish
  joint positivity, and additive shifts can leave joint support. Contrasts
  outside the box are still returned, not filtered out. Disagreement with
  permutation importance does not identify proxies; selection does not find
  a valid causal adjustment set. Existing support outputs must be recomputed.

## Stage-1 ranking sees curvature; stage-2 output is capped by sample size

* `cast_select()` stage 1 now ranks by a univariate quadratic-logistic
  signal (the smaller p-value of the two `poly(x, 2)` terms), so U-shaped
  drivers outrank noise that a linear correlation would tie with.
  Predictors whose GLM fails fall back to the marginal-correlation order.
* The retained set is capped at `ncov` predictors (new argument, default
  `ceiling(log2(n_presence))`, at most `maxncov = 12`), ordered by the
  interventional effect. When nothing clears the null, the fallback keeps
  the top `ncov` by effect instead of the whole stage-1 set.
* `scores` gains `stage1_p`, `stage1_rank`, `effect_rank` and
  `selected_reason` (`"null"`, `"null+top-ncov"`, `"fallback-top-ncov"`,
  `"full"`, `"excluded"`); `plot.cast_select()` now draws the
  interventional effect. `cast()` gains `select_ncov` / `select_maxncov`.
* `cast_effect_table()` gains a `support` column (worst 1-99% support
  fraction over the driver's shift set, the range diagnostic
  `cast_effect_support()` reports per shift), and two paired plots:
  `plot.cast_effect_table()` (magnitude, direction, hollow = low support)
  and `plot.cast_necessity()` (held-out AUC cost with fold range).
* Fixed: predicting from a fit object reloaded in a fresh session (e.g. via
  `readRDS`) silently returned all-`NA` contrasts, because the engine
  namespaces that register the `predict` S3 methods were never loaded.
  `predict_single_model()` now loads the engine namespace before predicting.
* Existing selection caches and downstream model/CV outputs must be
  recomputed.

# castSDM 0.10.0

## Stage-2 selection is now an interventional contrast

* `cast_select()` selects on the **interventional effect**: the change a
  predictor causes in the fitted probability when it is shifted while every
  other predictor stays at its observed value, averaged over `+shift_size` and
  `-shift_size` training SD. The permutation null is retained, so calibration
  still compares like with like. `scores` gains
  `interventional_effect`; `perm_importance` stays as a diagnostic.
* Why: permutation breaks a predictor's correlation with the others, so for
  collinear predictors the permuted rows leave the observed data support and
  the score is governed by the model's extrapolation behaviour (Hooker, Mentch
  & Zhou 2021). An interventional shift keeps the rest of the vector observed.
* `diagnostics$importance_agreement` reports the Spearman agreement between the
  two statistics; large disagreement marks predictors the forest leans on but
  barely responds to when changed, which is the signature of a collinear
  stand-in.
* New arguments `shift_size` and `max_rows`. Existing selection caches and
  downstream model/CV outputs must be recomputed.

## Interventional effect products gain a positivity diagnostic

* `cast_dose_response()` sweeps the size of an additive shift and reports the
  whole curve, so a saturating or asymmetric response is no longer collapsed to
  one number. `cast_effect_support()` reports, per driver and shift size, the
  fraction of observed predictor vectors still inside the training support.
* The support rule makes positivity operational: a shift answered mostly by
  extrapolation is flagged `estimable = FALSE`, never silently answered anyway.
* A shift that leaves the training range everywhere is refused rather than
  answered. `plot()` methods for both objects.
* `cast_importance()` now reports both attribution columns and the rank
  agreement, and its plot shows the interventional effect.

# castSDM 0.9.6

* Effect maps now average paired changes across models before taking their
  magnitude, matching effect tables even when model responses cancel. Both
  products use the same prediction implementation and per-cell summaries.
  `n` counts complete cells/rows; `n_shifts` reports the intervention count.
* Effect tables accept predictor-only frames, reject invalid shifts explicitly,
  and report `outside_range_fraction` (a marginal support warning, not a joint
  positivity test). Table and map omit non-finite predictor rows consistently.
* Stage-2 importance references are now feature-specific, rather than pooled
  across features. `passed_null` and `fallback` distinguish null exceedance
  from retaining the stage-1 set when nothing passes. Existing selection and
  downstream model/CV caches must be recomputed; old results are historical.
* Documentation separates model responses, RF predictive necessity and causal
  identification. Agreement of the two diagnostics is not an identification
  test. Citation metadata now follows the installed package version.

# castSDM 0.9.5

* Fixed: `cast_project_raster()` reported a meaningless `centroid_shift_km` on
  any projected grid. It fed the current and future centroids to a haversine
  formula, which assumes decimal degrees, so a raster in metres (e.g. an Albers
  or Web Mercator grid — the usual case for a national 1 km stack) produced
  shifts inflated by orders of magnitude: a real 24.7 km range shift was
  reported as 15915.1 km. The call site now dispatches on the raster's CRS,
  using haversine for lon/lat grids and planar distance scaled by
  `terra::linearUnits()` otherwise. `cast_project()`'s data-frame path is
  unchanged; its documented contract is EPSG:4326. No other projection
  statistic was affected, so existing suitability and range-size outputs stand.

# castSDM 0.9.4

* Fixed: `cast_effect_map()` aborted with `'data' must be a data.frame, not a
  matrix or an array` for every engine except `rf`. It handed each raster block
  to the engine predictors as a matrix, which `predict.gbm()` and `mgcv`'s
  predict method reject because they route through `model.frame()`. Effect maps
  from a `gam` fit or an ensemble containing one were therefore unreachable. `cast_effect_table()`, `cast_predict()` and `cast_predict_tiled()`
  were never affected.

# castSDM 0.9.3

* Fixed: `cast_necessity()` aborted with `length(object) == 1 is not TRUE`
  whenever a spatial fold had to be skipped. The skip logic was correct, but
  the warning reporting it used two `cli` pluralisation markers with no
  length-one quantity in scope, so the warning itself errored. This hit any
  species sparse enough that a fold contains no presences, which is exactly
  when the skip was needed.

# castSDM 0.9.2

* `cast_necessity()` gains a `num_threads=` argument passed through to
  `ranger`. It defaults to 1 (bit-for-bit reproducible) and lets callers
  parallelise the many refits — one full model plus one per driver, per fold —
  which matters when knocking out a full predictor stack across many species.

# castSDM 0.9.1

* `cast_necessity()` gains an optional `folds=` argument: an integer vector
  assigning each row of `data` to a fold. When supplied it overrides internal
  fold construction and `block_method`, so a caller holding pre-computed or
  frozen spatial folds (e.g. a pre-registered validation protocol) can run the
  knockout on exactly those folds. The returned object reports
  `block_method = "custom"` in that case.

# castSDM 0.9.0

The package's claim is narrowed to what it can defend: **attribution as an
audited pair**. `cast_effect_table()` / `cast_effect_map()` report the
interventional response to each driver, and the new `cast_necessity()`
reports what dropping that driver costs in held-out AUC. Where the two
disagree, the driver is substitutable by a collinear partner and attribution
to it is not identified. Variable selection is no longer presented as causal
at all: it is a two-stage screen scoped to parsimony and projection
robustness.

## Breaking changes

* `cast_select()` now takes `method = c("two_stage", "full")` only. The
  conditional (CPI) and DML selectors, the `"rf"` benchmark, and the
  `alpha`, `min_vars`, `max_candidates`, `dml_folds`, `force_include` and
  `n_folds` arguments are removed. Retired arguments land in `...` and warn
  rather than silently changing the screen.
  - Stage 1: rank by absolute marginal association, greedily keep predictors
    with pairwise |r| <= 0.7 (Dormann et al. 2013).
  - Stage 2: keep predictors whose random-forest permutation importance
    exceeds the 95th percentile of a null built by refitting on `n_perm`
    permuted responses (Altmann et al. 2010). A bare `importance > 0` rule
    retains roughly half of all uninformative predictors, because the null is
    centred on zero rather than bounded by it.
* `cast()` gains `select_n_perm` and drops `select_alpha`,
  `select_min_vars`, `select_max_candidates` and `select_n_folds`;
  `select_method` defaults to `"two_stage"`.
* `cast_effect_table()` / `cast_effect_map()` average over a frozen symmetric
  shift set (`shifts = c(-2, -1, 1, 2)` SD) instead of one arbitrary step,
  and report magnitude (`mean_abs_dHSS`) separately from direction
  (`mean_signed_dHSS`). Effect estimates no longer carry confidence
  intervals; a single signed step at one magnitude was not a stable ranking.
* `cast_screen_comparison()` and its `new_*` / `plot` methods are removed.
* `cast_cv()` now returns `oof`, the out-of-fold prediction surface, which
  `cast_ensemble()` needs to threshold the ensemble itself.

## New

* `cast_necessity()`: per-driver knockout diagnostic. Builds its own spatial
  folds and refits both the full and the knocked-out model inside every
  training fold, so `mean_dAUC` is an out-of-sample cost rather than an
  in-sample importance score. `pct_folds_positive` is reported so a stricter
  rule can be applied without refitting.

## Bug fixes

* `cast_ensemble()`: the ensemble score no longer silently substitutes 0 for
  a missing CBI or drops components with `na.rm`; a model with an incomplete
  metric row is excluded and warned about, and an all-NA score set aborts
  instead of producing uniform weights.
* `cast_ensemble()`: the ensemble threshold is now optimised on the ensemble
  out-of-fold surface rather than averaging per-model thresholds (which
  optimised nothing). Older `cast_cv` objects without `oof` fall back to the
  weighted mean of per-model thresholds, with a warning.
* `cast_ensemble()`: non-finite per-model predictions renormalise the weights
  per cell instead of poisoning the row or excluding the whole model, and a
  multi-layer mask no longer misaligns per-block cell indexing.
* `cast_ensemble()`: the early return when outputs already exist now carries
  the weights and threshold instead of `NULL`.
* `cast_cv()`: the spatial-buffer exclusion no longer builds a full n x n
  distance matrix, which exhausted memory on national data sets.
* `cast_predict_tiled()`: tiles are folded into the output as they are
  produced instead of being accumulated in a list, restoring the documented
  one-tile-at-a-time memory bound.
* `cast_project()`: change-class tallies use `terra::freq()` instead of
  reading the whole change raster into memory.
* Basemaps are read, reprojected and validated once per session, and a CRS
  carrying no EPSG code no longer errors the reprojection guard.
* `plot.cast_predict()` switches to `geom_raster()` on large grids, matching
  `plot.cast_ensemble()`.
* `summary.cast_result()` no longer emits an empty "best model" line when
  every AUC is NA.

# castSDM 0.8.0

Release refactor: a single, focused niche — conditional variable selection
that separates true ecological drivers from collinear bystanders. The
double-machine-learning (DML) selector, the ensemble-of-small-models (ESM)
rare-species route, and the self-built multi-species batch/replicate
infrastructure were removed; the remaining conditional (CPI) selector was
re-audited for hand-tuned truncation parameters.

## Breaking changes

* `cast_select(method = "cpi")` is now the only conditional selector
  (`method = "rf"` remains as the associational benchmark). The DML selector
  and its `dml_folds` / `select_dml_folds` arguments are removed; use
  `n_folds` / `select_n_folds`.
* Causal naming is dropped in favour of conditional naming: the importance
  reporter is now `cast_importance()` (was `cast_effect()`) and the what-if
  mapper is `cast_sensitivity()` (was `cast_counterfactual()`).
* Removed `cast_esm()`, `cast_rep()`, `cast_batch()`, `cast_batch_resume()`,
  `cast_run_from_config()`, `cast_config_template()`, and
  `cast_worker_budget()` (plus their `new_*`/`print`/`plot` methods). The
  package now targets single-species workflows; orchestration is left to the
  caller.

## Parameter audit (hand-tuned truncation removed)

* `max_candidates = NULL` (was a hidden 30-candidate RF pre-screen): every
  non-constant predictor is tested conditionally by default; a positive
  integer still enables an optional RF pre-screen for computational
  feasibility.
* `min_vars = 0` (was 3): an empty selection is allowed when nothing passes
  FDR. Ecological priors now enter only through `force_include`; topped-up
  predictors are flagged `fallback`.
* Knockoff estimator: predictors are Gaussian-quantile-transformed
  (rank -> `qnorm`) before `create.second_order()` so the multivariate-Gaussian
  exchangeability assumption holds for skewed bioclimatic predictors, and a
  **single** knockoff draw is used (the previous median-of-`knockoff_reps`
  p-value was not a valid pooled p-value). `knockoff_reps` is removed.
* The silent 200-tree cap on the RF nuisance/benchmark forests is removed;
  `num_trees` now passes through unchanged (default 300).
* `n_blocks` for `inference = "block"` is data-adaptive
  (~`min(50, max(20, n/20))`) and no longer hard-coded.

## Bug fixes

* Degenerate CPI significance layers (constant aggregated log-loss
  differences) now report `p = 1` instead of `NA`.
* Removed stale references to the deleted batch/ESM/replicate/DML code paths
  from printing, plotting, and the ODMAP report.

# castSDM 0.7.0

## OpenCodeReview hardening pass

Full-repository scan (alibaba/open-code-review) + fix round:

* `cast_predict_tiled()`: fixed a critical double-skip bug (the NA prototype
  written at init caused the write-back loop to skip prediction on first
  runs, leaving all-NA outputs) and a Windows file-lock issue (final write
  no longer re-opens the output file).
* `cast_ensemble_raster()`: cross-model SD now counts contributors per
  block, the mask geometry is validated before use, and empty
  model intersections abort with a clear error.
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
* Deduplicated the ensemble score/weight computation into shared helpers
  (`.cast_ensemble_scores()`, `.cast_ensemble_weights()`).

## Consensus selection across spatial folds

* `cast_cv()` now stores `selection_freq`: each predictor's fold-level
  selection frequency, the package's spatial-stability diagnostic.
* New `cast_consensus()` aggregates the fold selections into a consensus
  variable set (predictors retained in at least `threshold` of folds,
  default 0.5) and returns a `cast_select` object that plugs directly into
  `cast_fit()`.

## Defect fixes

* `cast_project()` marks cells with non-binary predictions as `NA` instead of
  an empty-string change class; summary counts ignore them explicitly.
* `cast_ensemble()` excludes models with non-finite predictions (with a
  warning) and renormalises weights, instead of silently zeroing them.
* `cast_predict()` errors with the missing column names when a prediction
  grid lacks fitted predictors.
* CPI confidence intervals are clamped at zero (CPI is non-negative).
* Selection scores flag `fallback` (kept through the `min_vars` floor)
  and `forced` (kept through `force_include`) predictors; the RF baseline
  in `cast_screen_comparison()` is budget-matched to CPI.
* Spatial-block CPI inference (cluster-robust t on block means) added as the
  default `inference = "block"`, with `"fold"` and `"observation"` kept for
  comparison.

# castSDM 0.6.0

## New core

* `cast_screen_comparison()` contrasts the castSDM conditional predictive-impact
  screen against the associational baselines researchers commonly use
  (correlation filter, stepwise VIF, univariate marginal screen, and
  random-forest permutation importance) on the same shared candidate pool.
* `method = "rf"` is retained as a conventional permutation-importance
  benchmark for comparison.
* `cast_cv()` re-runs selection inside every outer training fold and stores
  fold-specific selections, preventing feature-selection leakage.

## Plotting and stability

* China maps always show the complete, standard national frame plus the South
  China Sea inset; a range-restricted species is never cropped to its coloured
  cells.
* Suitability maps use a pale grey to deep teal-blue sequential ramp whose
  low anchor equals the basemap fill; range-change maps use a refined
  gain/loss/stable-present palette; what-if maps diverge around a neutral
  grey midpoint.
* `plot.cast_select()` gains a `top` argument and colours predictors by
  ecological class; `plot.cast_cv()` fold maps share the full-China frame.

# castSDM 0.2.0--0.5.1

Earlier experimental releases explored DAG, role-prior, CATE, RF-only, DML,
ESM, and stable/invariance screening designs. These interfaces are retired in
0.8.0.
