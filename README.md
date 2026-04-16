# castSDM <img src="man/figures/logo.png" align="right" height="139" alt="castSDM logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml)
[![License: GPL (>= 3)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

castSDM (Causal Structure-Informed Species Distribution Modeling) is an R package for causal-aware species distribution modeling. It does not claim to turn observational SDMs into definitive causal identification. Instead, it brings causal structure learning, treatment-effect estimation, screening, feature construction, model fitting, spatial prediction, and interpretation into one workflow so that causal evidence can be used explicitly upstream rather than only discussed after the fact.

The package supports both end-to-end and modular use. You can run a full one-species workflow with `cast()`, process many species with `cast_batch()`, or call each stage separately depending on how much control you need.

## What castSDM does

- Learns putative dependency structure among environmental predictors with bootstrapped DAG learning or alternative structure learners.
- Estimates average treatment-effect summaries with double machine learning and optional sensitivity / identifiability checks.
- Uses DAG and ATE outputs to guide variable screening, structural role assignment, and causal feature engineering.
- Fits a causally informed neural model alongside standard SDM baselines including RF, MaxEnt, and BRT.
- Produces habitat suitability maps, spatial cross-validation summaries, inter-model consistency diagnostics, SHAP explanations, MC Dropout uncertainty, and optional CATE surfaces.
- Scales from focal species analysis to multi-species batch workflows.

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
  "pROC", "ggplot2", "patchwork", "sf", "terra", "future",
  "future.apply", "fastshap", "xgboost", "grf", "viridisLite"
))

# torch is installed separately if you want neural models / NOTEARS
# install.packages("torch")
# torch::install_torch()
```

## Workflow at a glance

| Stage | Main functions | Purpose |
|-------|----------------|---------|
| Utility | `get_env_vars()`, `cast_vif()` | Detect environmental predictors and optionally remove severe collinearity |
| Preparation | `cast_prepare()` | Validate data and create train/test splits |
| Causal structure and effects | `cast_dag()`, `cast_ate()`, `cast_evalue()`, `cast_backdoor()` | Learn putative predictor structure, estimate effect summaries, and add sensitivity / identifiability checks |
| Screening and representation | `cast_screen()`, `cast_roles()`, `cast_features()` | Select variables and construct feature spaces informed by DAG + ATE outputs |
| Fitting and evaluation | `cast_fit()`, `cast_evaluate()`, `cast_cv()` | Fit models and evaluate them under hold-out or spatial blocking |
| Spatial outputs | `cast_predict()`, `cast_consistency()`, `cast_uncertainty()`, `cast_cate()` | Map habitat suitability, compare models, quantify neural uncertainty, and estimate heterogeneous effects |
| Interpretation and reporting | `cast_shap_xgb()`, `cast_shap_fit()`, `cast_shap_write_csv()`, `cast_report()` | Explain fitted models and export a reproducible HTML summary |
| Wrappers | `cast()`, `cast_batch()` | Run one-species or multi-species end-to-end workflows |

## Model and interpretation options

### Fitted models

- `cast`: causally informed neural network using ATE-weighted inputs and DAG-guided interactions.
- `mlp_ate`: neural network using ATE-weighted predictors without DAG interaction features.
- `mlp`: neural baseline on standardized raw predictors.
- `rf`: random forest baseline.
- `maxent`: MaxEnt baseline.
- `brt`: boosted regression tree baseline.

### DAG structure learners

- `bootstrap_hc`: bootstrap-aggregated score-based DAG learning.
- `pc`: constraint-based PC learning.
- `bidag_bge`: Bayesian MAP search via BiDAG.
- `notears_linear`: linear NOTEARS backend.

### SHAP routes

- `cast_shap_xgb()`: XGBoost surrogate on raw environmental predictors for a standardized explanation layer.
- `cast_shap_fit()`: direct SHAP computation for fitted `rf` or `cast` models. For `cast`, explanations are in the engineered feature space, not only the raw environmental space.

## Dependency notes

Hard imports are currently limited to `cli`, `grid`, and `scales`.

Extended functionality is activated through suggested packages:

| Purpose | Packages |
|---------|----------|
| Collinearity | `car` |
| DAG / graph tools | `bnlearn`, `BiDAG`, `dagitty`, `igraph`, `ggraph`, `torch` |
| Model fitting | `ranger`, `gbm`, `maxnet`, `torch` |
| Evaluation | `pROC` |
| CATE | `grf` |
| SHAP | `xgboost`, `fastshap` |
| Spatial / plotting | `ggplot2`, `patchwork`, `sf`, `terra`, `viridisLite` |
| Parallel execution | `future`, `future.apply`, `pkgload` |
| Reporting / vignettes | `rmarkdown`, `knitr` |

Functions that rely on optional packages will emit informative errors if the required backend is not installed.

## Interpretation guidance

castSDM is best understood as a causal-aware SDM toolkit, not as a guarantee of causal identification from observational ecological data. In practice:

- DAG outputs should be treated as data-informed structural hypotheses, not confirmed mechanism.
- ATE and CATE estimates remain conditional on measured covariates and modelling assumptions.
- The main strength of the package is operational: it makes causal structure, effect summaries, screening, model comparison, and interpretation available inside one reproducible SDM workflow.

## Citation

```bibtex
@software{castSDM,
  author  = {Liqiang Q},
  title   = {castSDM: Causal Structure-Informed Species Distribution Modeling},
  year    = {2026},
  url     = {https://github.com/EldonQ/castSDM},
  version = {0.1.0}
}
```

## License

GPL (>= 3)