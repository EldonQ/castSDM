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
  whole curve, so a saturating or threshold response is no longer collapsed to
  one number. `cast_effect_heatmap()` bins the observed data by the intervened
  predictor and an effect modifier. `cast_effect_support()` reports, per driver
  and shift size, the fraction of observed predictor vectors still inside the
  training support.
* The support rule makes positivity operational: a shift answered mostly by
  extrapolation is reported and drawn as an empty cell, never silently coloured
  in. Bins and points below the support threshold are flagged `supported =
  FALSE` / `estimable = FALSE`.
* A shift that leaves the training range everywhere is refused rather than
  answered. `plot()` methods for all three objects.
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
  from a `brt` or `gam` fit — or any ensemble containing one — were therefore
  unreachable. `cast_effect_table()`, `cast_predict()` and `cast_predict_tiled()`
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

* `cast_ensemble()`: the N-SDM score no longer silently substitutes 0 for a
  missing CBI or drops components with `na.rm`; a model with an incomplete
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
* Deduplicated the N-SDM score/weight computation into shared helpers
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
