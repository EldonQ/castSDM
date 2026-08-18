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
