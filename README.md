# castSDM <img src="man/figures/logo.png" align="right" height="139" alt="castSDM logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml)
[![License: GPL (>= 3)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

castSDM (Causal-Aware Species Distribution Modeling) is an R package that integrates causal DAG-guided variable selection with standard species distribution modeling algorithms. The core innovation is using Markov Blanket extraction from a learned DAG to select ecologically relevant predictors --- a theoretically grounded approach that is novel in the SDM literature.

The package supports both end-to-end and modular use. You can run a full one-species workflow with `cast()`, process many species with `cast_batch()`, drive reproducible runs from a YAML config via `cast_run_from_config()`, or call each stage separately depending on how much control you need.

castSDM is optimised for Windows workstations: all parallel paths use `future::multisession` (PSOCK), tile-based prediction keeps RAM bounded on large rasters, and per-species checkpoints enable crash-safe batch runs without HPC infrastructure.

## What castSDM does

- Learns a species-specific causal DAG among environmental predictors (with presence as a node) via PC, bootstrap HC, or BiDAG algorithms.
- Extracts the Markov Blanket of the response variable and fuses it with RF permutation importance for principled variable selection.
- Assigns causal roles (parent / child / co-parent / predictive) to selected variables based on DAG topology.
- Fits an ensemble of standard SDM algorithms (RF, BRT, MaxEnt, GAM) with performance-weighted, best-model, or equal-weight ensemble prediction.
- Projects species ranges under future climate scenarios and quantifies gain/loss/stable range cells with centroid shift.
- Produces habitat suitability maps, spatial cross-validation summaries, inter-model consistency diagnostics, SHAP explanations, and optional CATE surfaces with DAG-guided confounders.
- Scales from focal species analysis to multi-species batch workflows with per-species checkpoints, crash-safe resume, tile-based prediction, and unified hyperparameter tuning.

## Installation

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("EldonQ/castSDM")

# or
# install.packages("devtools")
devtools::install_github("EldonQ/castSDM")
```

castSDM keeps hard imports lightweight and installs many advanced capabilities on demand through suggested packages. For a typical full workflow, the following packages are commonly needed:

```r
install.packages(c(
  "car", "bnlearn", "BiDAG", "dagitty", "ranger", "gbm", "maxnet",
  "mgcv", "pROC", "ggplot2", "patchwork", "sf", "terra", "future",
  "future.apply", "fastshap", "xgboost", "grf", "viridisLite",
  "yaml", "peakRAM"
))
```

## Pipeline

```
species_data + env_data
       |
  cast_prepare()          # train/test split
       |
  cast_dag()              # DAG learning (PC / bootstrap HC / BiDAG)
       |                    with presence as a node
  cast_select()           # Markov Blanket + RF importance fusion
       |
  cast_roles()            # parent / child / co_parent / predictive
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
  cast_cate() [optional]  # DAG-guided CATE via causal forests
       |
  cast_shap_xgb() [opt.]  # SHAP explanations
```

## Workflow at a glance

| Stage | Main functions | Purpose |
|-------|----------------|---------|
| Utility | `get_env_vars()`, `cast_vif()` | Detect environmental predictors and optionally remove severe collinearity |
| Preparation | `cast_prepare()` | Validate data and create train/test splits |
| Causal structure | `cast_dag()`, `cast_backdoor()` | Learn predictor DAG with presence as a node; check adjustment sets |
| Variable selection | `cast_select()`, `cast_roles()` | Markov Blanket extraction + RF importance; assign causal roles |
| Fitting and evaluation | `cast_fit()`, `cast_evaluate()`, `cast_cv()`, `cast_tune()` | Fit RF/BRT/MaxEnt/GAM, evaluate, and tune hyperparameters |
| Rare species | `cast_esm()` | Ensemble of Small Models (bivariate GLM/GAM) fallback for low-presence species |
| Spatial outputs | `cast_predict()`, `cast_predict_tiled()`, `cast_consistency()` | Map habitat suitability (in-memory or tile-based for large rasters), compare models |
| Ensemble & projection | `cast_ensemble()`, `cast_project()` | Performance-weighted ensemble; future climate range change |
| CATE | `cast_cate()` | Spatially heterogeneous treatment effects with DAG-guided confounders |
| Interpretation | `cast_shap_xgb()`, `cast_shap_fit()`, `cast_report()` | Explain fitted models and export HTML summary |
| Wrappers | `cast()`, `cast_batch()`, `cast_batch_resume()` | Run one-species or multi-species end-to-end workflows, with crash-safe resume |
| Config-driven | `cast_run_from_config()`, `cast_config_template()` | Drive a full batch from a single YAML file |
| Parallelism | `cast_worker_budget()` | Allocate workers across species and intra-species stages |

## Model and interpretation options

### Fitted models

- `rf`: Random Forest (via `ranger`).
- `brt`: Boosted Regression Trees (via `gbm`).
- `maxent`: MaxEnt (via `maxnet`).
- `gam`: Generalized Additive Model (via `mgcv`).
- `esm`: Ensemble of Small Models for rare species.

### DAG structure learners

- `pc`: Constraint-based PC algorithm (default; supports mixed data via `mi-cg`).
- `bootstrap_hc`: Bootstrap-aggregated score-based DAG learning.
- `bidag_bge`: Bayesian MAP search via BiDAG.

### SHAP routes

- `cast_shap_xgb()`: XGBoost surrogate on raw environmental predictors for a standardized explanation layer.
- `cast_shap_fit()`: Direct SHAP computation for fitted RF models via fastshap.

## Dependency notes

Hard imports are currently limited to `cli`, `grid`, and `scales`.

Extended functionality is activated through suggested packages:

| Purpose | Packages |
|---------|----------|
| Collinearity | `car` |
| DAG / graph tools | `bnlearn`, `BiDAG`, `dagitty`, `igraph`, `ggraph`, `pcalg` |
| Model fitting | `ranger`, `gbm`, `mgcv`, `maxnet` |
| Evaluation | `pROC` |
| CATE | `grf` |
| SHAP | `xgboost`, `fastshap` |
| Spatial / plotting | `ggplot2`, `patchwork`, `sf`, `terra`, `viridisLite` |
| Future climate | `geodata` |
| Parallel execution | `future`, `future.apply`, `pkgload` |
| Config / tuning | `yaml` |
| Resource monitoring | `peakRAM` |
| Reporting / vignettes | `rmarkdown`, `knitr` |

Functions that rely on optional packages will emit informative errors if the required backend is not installed.

## Interpretation guidance

castSDM is best understood as a causal-aware SDM toolkit, not as a guarantee of causal identification from observational ecological data. In practice:

- DAG outputs should be treated as data-informed structural hypotheses, not confirmed mechanism.
- The Markov Blanket provides a theoretically minimal sufficient set for prediction, which improves variable selection over purely correlation-based methods.
- CATE estimates remain conditional on measured covariates and modelling assumptions.
- The main strength of the package is operational: it makes causal structure and principled variable selection available inside one reproducible SDM workflow.

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
