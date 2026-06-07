# castSDM <img src="man/figures/logo.png" align="right" height="139" alt="castSDM logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml)
[![License: GPL (>= 3)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

castSDM (Causal-Aware Species Distribution Modeling) integrates response-focused Markov Blanket screening with standard species distribution modeling algorithms. It is designed for reproducible SDM workflows where variable screening, spatial validation, ensemble prediction, future projection, and optional CATE surfaces are kept in one auditable pipeline.

The package should be interpreted as a causal-aware screening toolkit, not as proof of causal mechanisms from observational ecological data. DAG and Markov Blanket outputs are data-informed screening structures that help make predictor selection more explicit, sparse, and reproducible.

## What castSDM Does

- Learns a response-focused DAG or MB screening graph with `presence` included as a node.
- Defaults to `response_as_sink = TRUE`, forbidding `presence -> environmental predictor` edges.
- Selects predictors by combining Markov Blanket membership with RF permutation importance.
- Assigns screening roles: `mb_direct`, `mb_associated`, `importance_added`, and `importance_screened`.
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
  cast_select()           # Markov Blanket + RF importance screening
       |
  screen$roles            # mb_direct / mb_associated / importance_*
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
| Variable screening | `cast_select()` | Markov Blanket extraction + RF importance; assign screening roles |
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

### DAG and MB learners

- `mb_first`: two-stage MB discovery plus local PC; the default.
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
- Markov Blanket screening is useful for sparse, auditable variable selection, but its causal interpretation depends on observational assumptions.
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
