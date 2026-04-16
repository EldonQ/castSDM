#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' The \pkg{castSDM} package provides a complete pipeline for causal
#' structure-informed species distribution modeling:
#'
#' 1. **Data Preparation**: VIF-based collinearity screening ([cast_vif()])
#' 2. **Causal Structure Learning**: DAG discovery via bootstrap Hill-Climbing
#'    ([cast_dag()])
#' 3. **Causal Effect Estimation**: Double Machine Learning ATE ([cast_ate()])
#' 4. **Variable Screening**: Adaptive multi-criteria selection ([cast_screen()])
#' 5. **Causal Role Assignment**: Root/Mediator/Terminal classification
#'    ([cast_roles()])
#' 6. **Feature Engineering**: ATE-weighting + DAG interactions
#'    ([cast_features()])
#' 7. **Model Fitting**: CI-MLP and traditional SDMs ([cast_fit()])
#' 8. **Evaluation**: AUC, TSS metrics ([cast_evaluate()])
#' 9. **Prediction**: Spatial habitat suitability scores ([cast_predict()])
#' 10. **CATE Estimation**: Spatially heterogeneous treatment effects
#'     ([cast_cate()])
#'
#' @section Quick Start:
#' ```
#' result <- cast(species_data, env_data)
#' summary(result)
#' plot(result)
#' ```
#'
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom stats predict cor sd median var na.omit lm residuals
#' @importFrom utils head tail
NULL

# Suppress R CMD check NOTEs for non-standard evaluation variables
utils::globalVariables(c("self", ".data", "index", "edge_strength",
                         "from", "to", "evalue_point", "evalue_ci",
                         "significant"))
