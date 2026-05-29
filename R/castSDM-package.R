#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' The \pkg{castSDM} package provides a complete pipeline for causal-aware
#' species distribution modeling. The core innovation is **DAG-guided
#' variable selection via Markov Blanket** — a theoretically grounded
#' approach to selecting ecologically relevant predictors.
#'
#' **Pipeline steps:**
#'
#' 1. **Data Preparation**: Train/test splitting, VIF-based collinearity
#'    screening ([cast_prepare()], [cast_vif()])
#' 2. **Causal Structure Learning**: DAG discovery via PC algorithm,
#'    two-stage MB-First (IAMB + local PC), or bootstrap Hill-Climbing,
#'    with presence as a node ([cast_dag()])
#' 3. **Variable Selection**: Markov Blanket extraction + RF importance
#'    fusion ([cast_select()])
#' 4. **Model Fitting**: RF, BRT, MaxEnt, GAM ([cast_fit()]), plus ESM for
#'    rare species ([cast_esm()])
#' 5. **Evaluation**: AUC, TSS, CBI metrics ([cast_evaluate()]),
#'    spatial block cross-validation ([cast_cv()])
#' 6. **Prediction**: In-memory ([cast_predict()]) or tile-based
#'    ([cast_predict_tiled()]) spatial habitat suitability mapping
#' 7. **Ensemble**: Performance-weighted, best-model, or equal-weight
#'    ensemble prediction ([cast_ensemble()])
#' 8. **Future Projection**: Range change analysis under climate scenarios
#'    ([cast_project()])
#' 9. **Batch Workflows**: Multi-species runs with checkpoint/resume
#'    ([cast_batch()], [cast_batch_resume()]), YAML-driven config
#'    ([cast_run_from_config()]), and worker budget allocation
#'    ([cast_worker_budget()])
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
                         "from", "to"))
