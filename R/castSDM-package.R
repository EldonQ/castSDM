#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' The \pkg{castSDM} package provides a complete pipeline for causal-aware
#' species distribution modeling. The core innovation is **role-constrained
#' causal-core screening**: ecological roles are declared before selection,
#' and only direct candidates enter a shallow RF, invariance, conditional-
#' evidence, and redundancy ranking. The complete role and exclusion audit is
#' retained in the returned object.
#'
#' **Pipeline steps:**
#'
#' 1. **Data Preparation**: Train/test splitting, VIF-based collinearity
#'    screening ([cast_prepare()], [cast_vif()])
#' 2. **Causal Role Specification**: reviewed default or user overrides
#'    ([cast_causal_spec()])
#' 3. **Variable Selection**: role-constrained shallow RF, invariance,
#'    conditional evidence, and correlation redundancy control
#'    ([cast_select()]); DAG/MB learning remains optional ([cast_dag()])
#' 4. **Model Fitting**: RF, BRT, MaxEnt, GAM ([cast_fit()]), plus ESM for
#'    rare species ([cast_esm()])
#' 5. **Evaluation**: AUC, TSS, CBI metrics ([cast_evaluate()]),
#'    spatial block cross-validation ([cast_cv()])
#' 6. **Prediction**: In-memory ([cast_predict()]) or tile-based
#'    ([cast_predict_tiled()]) spatial habitat suitability mapping
#' 7. **Ensemble**: Performance-weighted, best-model, or equal-weight
#'    ensemble prediction ([cast_ensemble()])
#' 8. **CATE** (optional): Spatially heterogeneous treatment effects via
#'    causal forests ([cast_cate()], requires `grf`)
#' 9. **Future Projection**: Range change analysis under climate scenarios
#'    ([cast_project()])
#' 10. **Batch Workflows**: Multi-species runs with checkpoint/resume
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
