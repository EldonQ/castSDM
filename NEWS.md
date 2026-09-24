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
