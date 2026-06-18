# castSDM <img src="man/figures/logo.png" align="right" height="139" alt="castSDM logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml)
[![License: GPL (>= 3)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

castSDM (Causal-Aware Species Distribution Modeling) integrates spatially invariant causal screening with standard species distribution modeling algorithms. It is designed for reproducible SDM workflows where variable screening, spatial validation, ensemble prediction, future projection, and optional CATE surfaces are kept in one auditable pipeline.

The package should be interpreted as a causal-aware screening toolkit, not as proof of causal mechanisms from observational ecological data. Its default selector searches for predictors whose effects are predictive, spatially stable, direction-consistent, and non-redundant. DAG and Markov Blanket outputs remain available as audit structures and legacy comparators.

## What castSDM Does

- Learns a response-focused DAG or MB screening graph with `presence` included as a node.
- Defaults to `response_as_sink = TRUE`, forbidding `presence -> environmental predictor` edges.
- Selects predictors with invariant screening: RF importance, bootstrap stability, spatial-block effect consistency, and correlation-cluster redundancy control.
- Assigns screening roles: `invariant_driver`, `stable_predictive`, `predictive_rescue`, and `redundant_proxy`.
- Fits RF, BRT, MaxEnt, and GAM models, then evaluates them with AUC, TSS, and CBI.
- Supports spatial block cross-validation, performance-weighted ensembles, future projections, raster prediction, and optional CATE estimation.
- Scales from one-species workflows with `cast()` to multi-species workflows with `cast_batch()` and checkpointed resume.

## Installation

```r
# install.packages("pak")
pak::pak("EldonQ/castSDM")

# or
# install.packages("devtools")
devtools::install_github("EldonQ/castSDM")
```

For a full workflow, install the optional backends you plan to use:

```r
install.packages(c(
  "car", "bnlearn", "BiDAG", "ranger", "gbm", "maxnet", "mgcv",
  "pROC", "ggplot2", "patchwork", "sf", "terra", "future",
  "future.apply", "grf", "ggraph", "igraph", "viridisLite",
  "yaml", "peakRAM"
))
```

## Pipeline

```r
species_data + env_data
       |
  cast_prepare()          # train/test split
       |
  cast_dag()              # response-focused DAG / MB graph
       |                    response_as_sink = TRUE by default
  cast_select()           # invariant causal screening
       |
  screen$roles            # invariant_driver / stable_predictive / rescue
       |
  cast_fit()              # RF, BRT, MaxEnt, GAM
       |
  cast_cv()               # spatial block cross-validation
       |
  cast_evaluate()         # hold-out AUC, TSS, CBI
       |
  cast_predict()          # habitat suitability maps
       |
  cast_ensemble()         # weighted / best / equal ensemble
       |
  cast_project()          # future climate range change
       |
  cast_cate() [optional]  # CATE via causal forests
```

## Workflow at a Glance

| Stage | Main functions | Purpose |
|-------|----------------|---------|
| Utility | `get_env_vars()`, `cast_vif()` | Detect environmental predictors and optionally remove severe collinearity |
| Preparation | `cast_prepare()` | Validate data and create train/test splits |
| Screening graph | `cast_dag()` | Learn response-focused DAG or MB screening graph |
| Variable screening | `cast_select()` | Spatial invariance + RF/stability + redundancy control; assign screening roles |
| Fitting and evaluation | `cast_fit()`, `cast_evaluate()`, `cast_cv()` | Fit RF/BRT/MaxEnt/GAM and evaluate transferability |
| Rare species | `cast_esm()` | Ensemble of Small Models fallback for low-presence species |
| Spatial outputs | `cast_predict()`, `cast_predict_tiled()` | Map habitat suitability in memory or tile by tile |
| Ensemble and projection | `cast_ensemble()`, `cast_project()` | Ensemble suitability and future range change |
| Raster workflows | `cast_load_bioclim()`, `cast_background()`, `cast_ensemble_raster()`, `cast_project_raster()` | Raster-native sampling, prediction, and projection |
| CATE | `cast_cate()` | Spatially heterogeneous treatment effects with MB-guided covariates |
| Wrappers | `cast()`, `cast_batch()`, `cast_batch_resume()` | One-species or multi-species end-to-end workflows |
| Config-driven | `cast_run_from_config()`, `cast_config_template()` | Drive a batch from a YAML file |
| Parallelism | `cast_worker_budget()` | Allocate workers across species and intra-species stages |

## Model and Screening Options

### Fitted models

- `rf`: Random Forest via `ranger`.
- `brt`: Boosted Regression Trees via `gbm`.
- `maxent`: MaxEnt via `maxnet`.
- `gam`: Generalized Additive Model via `mgcv`.
- `esm`: Ensemble of Small Models for rare species.

### Variable screening

- `invariant_screen`: default. Fast spatial-block invariant screening with an adaptive score break and correlation-cluster sparsity; it does not target a fixed number of predictors unless `max_vars` is set.
- `mb_rf`: legacy Markov Blanket + RF importance fusion.
- `rf`: pure RF/stability screening with redundancy control.

### DAG and MB learners

- `mb_first`: two-stage MB discovery plus local PC; available for audit and legacy screening.
- `pc`: constraint-based PC algorithm.
- `bootstrap_hc`: bootstrap-aggregated score-based DAG learning.
- `bidag_bge`: Bayesian MAP search via BiDAG.

## Dependency Notes

Hard imports are intentionally small: `cli` and `grid`. Extended functionality is activated through suggested packages.

| Purpose | Packages |
|---------|----------|
| Collinearity | `car` |
| DAG / graph tools | `bnlearn`, `BiDAG`, `igraph`, `ggraph`, `pcalg` |
| Model fitting | `ranger`, `gbm`, `mgcv`, `maxnet` |
| Evaluation | `pROC` |
| CATE | `grf` |
| Spatial / plotting | `ggplot2`, `patchwork`, `sf`, `terra`, `viridisLite` |
| Parallel execution | `future`, `future.apply`, `pkgload` |
| Config | `yaml` |
| Resource monitoring | `peakRAM` |
| Vignettes | `rmarkdown`, `knitr` |

Functions that rely on optional packages emit informative errors if the required backend is not installed.

## Interpretation Guidance

- DAG outputs should be treated as data-informed structural screening hypotheses.
- Invariant screening is a testable proxy for transferable causal signal: variables should remain predictive and direction-consistent across spatial/environmental blocks.
- Markov Blanket screening is retained for audit, but dense MB outputs should be treated as non-selective diagnostics.
- Screening roles describe how variables entered the selected set; they are not confirmed ecological mechanisms.
- CATE estimates remain conditional on measured covariates and model assumptions.

## Citation

```bibtex
@software{castSDM,
  author  = {Liqiang Q},
  title   = {castSDM: Causal-Aware Species Distribution Modeling},
  year    = {2026},
  url     = {https://github.com/EldonQ/castSDM},
  version = {0.3.0}
}
```

## License

GPL (>= 3)
