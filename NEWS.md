# castSDM 0.3.0

## Focus

This release frames castSDM as a causal-aware SDM toolkit focused on
response-focused Markov Blanket screening, standard multi-algorithm SDM
ensembles, spatial validation, and future projection.

## Breaking changes

* Removed the previous ATE, neural-network, feature-engineering, E-value, and
  NOTEARS components from the active workflow.
* Replaced the older adaptive screening path with `cast_select()`.
* Removed old role terminology from new outputs. New runs report screening
  roles rather than parent/child causal labels.

## New and changed

* `cast_dag()` now supports `response_as_sink = TRUE` by default, forbidding
  directed `presence -> environmental predictor` edges.
* `cast_select()` combines response-focused Markov Blanket membership with RF
  permutation importance.
* New screening roles:
  * `mb_direct`
  * `mb_associated`
  * `importance_added`
  * `importance_screened`
* Default DAG structure learning is `mb_first`, a two-stage Markov Blanket
  discovery plus local PC workflow.
* `cast_cate()` uses the same response-focused MB definition as
  `cast_select()` when a `cast_dag` object is supplied.
* `cast_ensemble()` supports weighted, best-model, and equal-weight ensemble
  prediction.
* `cast_project()` supports future climate range projection with
  gain/loss/stable change maps and centroid shift statistics.
* `cast_batch()` exposes prediction and ensemble controls aligned with
  `cast()`: `do_predict`, `do_ensemble`, and `ensemble_method`.
* Raster-native helpers were added for background sampling, bioclim loading,
  ensemble prediction, and future projection.

## Interpretation

castSDM outputs should be read as data-informed screening structures, not as
confirmed ecological causal mechanisms. The package is intended to make
predictor screening explicit, sparse, reproducible, and easier to audit across
single-species and multi-species SDM workflows.

# castSDM 0.2.0

## Windows platform improvements

* Tile-based raster prediction keeps peak RAM bounded on large rasters.
* Multi-layer worker budgeting distributes cores across species and
  intra-species work on Windows PSOCK clusters.
* Checkpoint and crash-safe resume support per-species recovery.
* YAML config-driven workflows support reproducible batch runs.
* GAM and Ensemble of Small Models support were added.
* Resource profiling vignette was added.

## Data

* Bundled species datasets and `china_env_grid` were added for examples and
  spatial prediction.

# castSDM 0.1.0

## Initial features

* Core modular SDM pipeline with DAG learning, screening, fitting,
  evaluation, prediction, CATE surfaces, spatial CV, and batch workflows.
* S3 classes with `print()`, `summary()`, and `plot()` methods.
* VIF-based collinearity screening.
* Four DAG structure-learning routes available at the time.
