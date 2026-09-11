#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' \pkg{castSDM} is a species-distribution-modelling workflow with an explicit
#' separation between two questions that are routinely conflated. Variable
#' selection ([cast_select()]) answers *which predictors to carry*: a
#' conventional two-stage screen — collinearity thinning followed by an
#' importance filter calibrated against a permutation null — serving parsimony
#' and projection robustness. Attribution answers *which predictors matter and
#' how*, and is handled separately by an audit that pairs interventional
#' effects ([cast_effect_table()], [cast_effect_map()], [cast_sensitivity()])
#' with a knockout necessity check ([cast_necessity()]). Neither product is a
#' claim of automatic causal discovery; the pairing is a reliability audit of
#' correlational attribution, and a driver is only reported as credible when
#' both agree.
#'
#' **Pipeline steps:**
#'
#' 1. **Data Preparation**: train/test splitting and optional VIF collinearity
#'    screening ([cast_prepare()], [cast_vif()])
#' 2. **Variable Selection**: two-stage collinearity + permutation-null
#'    importance screening ([cast_select()]), reported via [cast_importance()]
#' 3. **Model Fitting**: RF, BRT, MaxEnt, GAM ([cast_fit()])
#' 4. **Evaluation**: AUC, TSS, CBI metrics ([cast_evaluate()]), nested spatial
#'    cross-validation with fold-specific selection ([cast_cv()])
#' 5. **Prediction**: In-memory ([cast_predict()]) or tile-based
#'    ([cast_predict_tiled()]) spatial habitat suitability mapping
#' 6. **Ensemble**: performance-weighted, best-model, or equal-weight
#'    ensemble prediction ([cast_ensemble()])
#' 7. **Future Projection**: range-change analysis under climate scenarios
#'    ([cast_project()])
#' 8. **Attribution audit**: interventional effect sizes and maps
#'    ([cast_effect_table()], [cast_effect_map()]) paired with knockout
#'    necessity ([cast_necessity()])
#'
#' @section Quick Start:
#' ```
#' result <- cast(species_data, env_data)
#' summary(result)
#' plot(result)
#' ```
#'
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom stats predict cor sd median var na.omit lm
#' @importFrom utils head
NULL

# Suppress R CMD check NOTEs for non-standard evaluation variables
utils::globalVariables(c("self", ".data", "index"))
