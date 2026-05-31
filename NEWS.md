# castSDM 0.3.0

## Breaking changes

This is a major refactoring release that transforms castSDM from a
"causal-inference-powered SDM" to a **"Causal-Aware SDM"** focused on
DAG-guided variable selection and standard multi-algorithm ensemble
prediction.

### Removed

* **ATE module**: `cast_ate()`, `new_cast_ate()`, `plot.cast_ate()`,
  `print.cast_ate()` are removed. The DML-based Average Treatment Effect
  estimation is no longer part of the pipeline.
* **CI-MLP neural network**: The `cast`, `mlp_ate`, and `mlp` model types
  are removed from `cast_fit()`. The pure-R MLP backend
  (`cast_mlp_fit()` / `cast_mlp_predict()`) and all torch-based code are
  removed. `torch` is no longer a suggested dependency.
* **Feature engineering**: `cast_features()` is removed. DAG interaction
  features and ATE-weighted inputs are no longer constructed.
* **E-value sensitivity**: `cast_evalue()` is removed (depended on ATE).
* **NOTEARS**: The `notears_linear` DAG structure method is removed
  (depended on torch).
* **Old screening**: `cast_screen()` (the ATE+DAG+RF adaptive screening
  function) is replaced by `cast_select()`.

### New

* **`cast_select()`**: DAG-guided variable selection using Markov Blanket
  extraction + RF permutation importance. MB variables are auto-selected;
  non-MB variables must exceed the median MB importance. Returns a
  `cast_select` object with `$selected`, `$scores`, and `$roles`
  (parent/child/co_parent/predictive).
* **`cast_ensemble()`**: Performance-weighted, best-model, or equal-weight
  ensemble prediction from CV metrics.
* **`cast_project()`**: Future climate range projection with gain/loss/stable
  change maps and centroid shift statistics.
* **`cast_dag()` now includes presence**: The response variable is included
  as a node in DAG learning (converted to factor for bnlearn's Conditional
  Gaussian network). This enables Markov Blanket extraction. Default
  structure method changed from `bootstrap_hc` to `pc` with `mi-cg` test.
* **DAG-guided CATE confounders**: `cast_cate()` now accepts a `dag`
  argument and uses `MB(presence) \ {treatment}` as confounders.

### Changed

* Default `models` in `cast()` and `cast_batch()` changed from
  `c("cast", "rf", "brt")` to `c("rf", "brt", "maxent", "gam")`.
* Default `dag_structure_method` changed from `"bootstrap_hc"` to `"pc"`.
* `cast_roles()` now reads MB-based roles (parent/child/co_parent/predictive)
  from `cast_select()` output, with fallback to degree-based heuristic.
* `cast_shap_fit()` now only supports `which = "rf"` (the `"cast"` path
  is removed).
* `cast_tune()` now only supports RF and BRT grids (CI-MLP grid removed).
* `cast_fit()` signature simplified: removed `ate`, `backend`, `n_epochs`,
  `n_runs`, `patience`, `val_fraction`, focal loss parameters, and all
  NN-related arguments.
* `cast_cv()` signature simplified: removed `ate`, `n_epochs`, `n_runs`.
* `cast_batch()` signature simplified: removed all ATE, NN, and NOTEARS
  parameters.
* SHAP panel layout changed from 2x3 (XGB/RF/CAST) to 2x2 (XGB/RF).

# castSDM 0.2.0

## Windows platform improvements

### Performance and reliability (P0)

* **CI-MLP pure R backend** (`cast_mlp_fit()` / `cast_mlp_predict()`): the
  default CI-MLP backend is now a pure R implementation (BLAS-accelerated
  matrix operations with LayerNorm, SiLU, AdamW, focal loss, and early
  stopping). The ~1.5 GB `torch` runtime is no longer required for
  installation or parallel execution. `cast_fit(backend = "torch")` is
  still available as an opt-in alternative.
* **Tile-based raster prediction** (`cast_predict_tiled()`): reads and
  predicts covariate rasters tile by tile with per-tile checkpoints and
  optional `future.apply` parallelism. Peak RAM is bounded regardless of
  grid size.
* **Multi-layer worker budget** (`cast_worker_budget()`): distributes a
  fixed pool of workers across species, intra-species steps, and CV folds
  to fill available cores without oversubscription on Windows PSOCK
  clusters.
* **Checkpoint and crash-safe resume** (`cast_batch_resume()`): every
  pipeline step writes a `.steps/<name>.rds` cache; `cast_batch_resume()`
  scans output directories and re-runs only incomplete species. Elapsed
  wall-clock time and peak RAM are logged to `resource_log.csv`.

### Extensibility and method completeness (P1)

* **YAML config-driven workflow** (`cast_run_from_config()`,
  `cast_config_template()`): a single YAML file specifies all species,
  environmental data, hyperparameters, and parallel settings for
  reproducible, version-controlled batch runs.
* **GAM algorithm** added to `cast_fit()` (via `mgcv::gam`).
* **Ensemble of Small Models** (`cast_esm()`): bivariate GLM/GAM
  sub-models with AUC-weighted averaging — the recommended fallback for
  rare species where BRT would overfit.

### Evaluation and observability (P3)

* **Unified hyperparameter tuning** (`cast_tune()`, `cast_default_grid()`):
  per-algorithm grid search scored on `AUC + TSS - 1` across stratified
  splits. `cast_fit(tune_result = ...)` consumes the result directly.
* **Resource profiling vignette** (`vignettes/perf-profile.Rmd`):
  ggplot2-based visualisation of per-species per-step wall time and peak
  RAM from `resource_log.csv`.

## New features

* 32 bundled species datasets (ISEA3H Res9 hexagonal grid, China) plus
  `china_env_grid` for spatial prediction.

# castSDM 0.1.0

## New features

* Core pipeline: `cast()` for one-step causal SDM modeling.
* Modular functions: `cast_dag()`, `cast_ate()`, `cast_screen()`,
  `cast_roles()`, `cast_features()`, `cast_fit()`, `cast_evaluate()`,
  `cast_predict()`, `cast_cate()`.
* S3 classes with `print()`, `summary()`, and `plot()` methods.
* VIF-based collinearity screening via `cast_vif()`.
* CI-MLP neural network with residual blocks, focal loss, and cosine
  annealing LR schedule.
* Four DAG structure learning methods: bootstrap HC, PC, BiDAG BGe,
  and NOTEARS linear.
* E-value sensitivity analysis via `cast_ate()` (`compute_evalue = TRUE`).
* Backdoor criterion check via `cast_roles()`.
* SHAP value computation for XGBoost models via `cast_shap_xgb()`.
* Spatial CV support in `cast_cv()`.
* Batch multi-species workflow via `cast_batch()`.
* Publication-quality spatial heatmaps (HSS & CATE) via helper functions
  in `inst/examples/cast_spatial_heatmap_helpers.R`.

## Bug fixes

* CI-MLP training loop: fixed torch dataloader iterator returning a
  symbol instead of a list on exhaustion.
* ATE multiple testing correction changed from Bonferroni to FDR (BH).
* CATE heatmap: HSS mask now applied after spatial interpolation to
  prevent nearest-neighbor infill into non-habitat areas.
* Consistency matrix: NA values displayed as "NA" on grey background
  instead of "0.000".
* `cast_dag()`: suppressed benign bnlearn v-structure warnings.
