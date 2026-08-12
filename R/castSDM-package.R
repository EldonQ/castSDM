#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' \pkg{castSDM} is a focused species distribution modelling toolkit. It pairs
#' a conditional variable selector (conditional predictive impact with FDR
#' control by default; an optional double-machine-learning effect reporter)
#' with standard RF / BRT / MaxEnt / GAM learners, nested spatial
#' cross-validation, habitat-suitability prediction, performance-weighted
#' ensembles and future climate projection. The selector retains predictors
#' whose conditional contribution on occurrence survives FDR control; the
#' package supports interpretable driver analysis but does not claim
#' automatic causal discovery.
#'
#' **Pipeline steps:**
#'
#' 1. **Data Preparation**: train/test splitting and optional VIF collinearity
#'    screening ([cast_prepare()], [cast_vif()])
#' 2. **Variable Selection**: conditional predictive impact screening with
#'    FDR control ([cast_select()]); effect reporting via [cast_effect()]
#'    and counterfactual what-if maps ([cast_counterfactual()])
#' 3. **Model Fitting**: RF, BRT, MaxEnt, GAM ([cast_fit()]), plus ESM for
#'    rare species ([cast_esm()])
#' 4. **Evaluation**: AUC, TSS, CBI metrics ([cast_evaluate()]), nested spatial
#'    cross-validation with fold-specific selection ([cast_cv()])
#' 5. **Prediction**: In-memory ([cast_predict()]) or tile-based
#'    ([cast_predict_tiled()]) spatial habitat suitability mapping
#' 6. **Ensemble**: performance-weighted, best-model, or equal-weight
#'    ensemble prediction ([cast_ensemble()])
#' 7. **Future Projection**: range-change analysis under climate scenarios
#'    ([cast_project()])
#' 8. **Batch Workflows**: multi-species runs with checkpoint/resume
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
utils::globalVariables(c("self", ".data", "index"))
