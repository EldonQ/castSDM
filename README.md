# castSDM <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EldonQ/castSDM/actions/workflows/R-CMD-check.yaml)
[![License: GPL (>=
3)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%203%29-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**castSDM** (Causal Structure-informed Species Distribution Modeling)
integrates causal inference methods with species distribution modeling
(SDM), enabling researchers to learn species-specific causal structures,
estimate variable-level causal effects, and build causally-informed
predictive models for habitat suitability mapping.

## Motivation

Traditional SDMs rely on correlative associations between species
occurrences and environmental predictors, making them vulnerable to
spurious correlations and poor transferability under novel conditions.
**castSDM** addresses these limitations by:

- Discovering causal structure among predictors via **Directed Acyclic
  Graphs (DAGs)**
- Quantifying causal effect sizes through **Double Machine Learning
  (DML)** based Average Treatment Effect (ATE) estimation
- Engineering ecologically interpretable features guided by causal
  theory
- Training **Causally-Informed Multi-Layer Perceptrons (CI-MLP)** that
  embed causal knowledge into model architecture
- Mapping **spatially heterogeneous treatment effects (CATE)** via
  causal forests

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("EldonQ/castSDM")
```

Or using `devtools`:

``` r
# install.packages("devtools")
devtools::install_github("EldonQ/castSDM")
```

## Quick Start

``` r
library(castSDM)

# Load example data
data(ovis_ammon)
data(china_env_grid)

# One-step pipeline
result <- cast(
  species_data = ovis_ammon,
  env_data     = china_env_grid,
  models       = c("cast", "rf", "maxent", "brt")
)

summary(result)
plot(result)
```

## Pipeline Overview

The **castSDM** workflow consists of 10 modular steps, each accessible
as a standalone function:

| Step | Function           | Description                           |
|-----:|--------------------|---------------------------------------|
|    0 | `cast_vif()`       | VIF-based collinearity screening      |
|    1 | `cast_prepare()`   | Data validation and train/test split  |
|    2 | `cast_dag()`       | Causal DAG learning (bootstrap HC)    |
|    3 | `cast_ate()`       | DML-based ATE estimation              |
|    4 | `cast_screen()`    | Adaptive variable screening           |
|    5 | `cast_roles()`     | Causal role assignment                |
|    6 | `cast_features()`  | Causal feature engineering            |
|    7 | `cast_fit()`       | Model fitting (CI-MLP / RF / MaxEnt / BRT) |
|    8 | `cast_evaluate()`  | Hold-out evaluation (AUC, TSS, CBI, SEDI, Kappa, PRAUC) |
|    8b| `cast_cv()`        | Spatial K-Fold cross-validation       |
|    9 | `cast_predict()`   | Spatial habitat suitability prediction|
|   10 | `cast_cate()`      | Spatially heterogeneous CATE mapping  |

All steps are also available through the unified `cast()` function.

## Step-by-Step Usage

``` r
library(castSDM)
data(ovis_ammon)

# 1. Data preparation
split <- cast_prepare(ovis_ammon, train_fraction = 0.7, seed = 42)

# 2. DAG learning
dag <- cast_dag(split$train, R = 100, seed = 42)

# 3. ATE estimation
ate <- cast_ate(split$train, K = 2, seed = 42)

# 4. Variable screening
screen <- cast_screen(dag, ate, split$train)

# 5. Causal role assignment
roles <- cast_roles(screen, dag)

# 6. Model fitting
fit <- cast_fit(split$train, screen = screen, dag = dag, ate = ate,
                models = c("cast", "rf"))

# 7. Evaluation
eval <- cast_evaluate(fit, split$test)
print(eval)

# 8. Spatial CV (optional, recommended)
cv <- cast_cv(ovis_ammon, screen = screen, dag = dag, ate = ate,
              k = 5, models = c("rf"))

# 9. Prediction
data(china_env_grid)
pred <- cast_predict(fit, china_env_grid)
plot(pred, model = "cast", basemap = "china")

# 10. CATE estimation
cate <- cast_cate(split$train, ate = ate, screen = screen, top_n = 3,
                  predict_data = china_env_grid, seed = 42)
plot(cate)
```

## Key Features

- **Causal inference integration**: DAG structure learning and DML-ATE
  estimation provide an ecologically grounded variable selection and
  feature engineering framework.
- **Multi-model comparison**: Fit CI-MLP alongside Random Forest,
  MaxEnt, and BRT within the same pipeline for fair benchmarking.
- **Spatial cross-validation**: Spatially blocked K-Fold CV with grid or
  cluster blocking strategies for honest evaluation under spatial
  autocorrelation.
- **Rich S3 methods**: `print()`, `summary()`, and `plot()` methods for
  every pipeline output, enabling rapid visual exploration.
- **Modular design**: Use individual functions for fine-grained control,
  or `cast()` for an end-to-end workflow.

## Example Data

The package includes two built-in datasets:

- **`ovis_ammon`**: Presence/background data for Argali (*Ovis ammon*)
  on the ISEA3H Res-9 hexagonal grid across China, with 19 bioclimatic
  variables.
- **`china_env_grid`**: Full environmental grid covering China at the
  same resolution, for spatial prediction.

## Dependencies

**castSDM** has minimal hard dependencies (`cli`, `coro`, `data.table`,
`rlang`). Extended functionality uses suggested packages:

- **DAG**: `bnlearn`, `igraph`, `ggraph`
- **Models**: `ranger` (RF), `maxnet` (MaxEnt), `gbm` (BRT), `torch`
  (CI-MLP)
- **CATE**: `grf` (causal forests)
- **Plotting**: `ggplot2`, `patchwork`, `sf`, `terra`, `viridisLite`
- **Evaluation**: `pROC`

## Citation

If you use **castSDM** in your research, please cite:

``` bibtex
@software{castSDM,
  author = {Liqiang},
  title = {castSDM: Causal Structure-Informed Species Distribution Modeling},
  year = {2026},
  url = {https://github.com/EldonQ/castSDM},
  version = {0.1.0}
}
```

## License

GPL (\>= 3)
