#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' \pkg{castSDM} separates true ecological drivers from collinear bystanders
#' in species distribution modelling. Its conditional variable selector
#' (conditional predictive impact with Benjamini-Hochberg FDR control by
#' default) tests each predictor's contribution given every other predictor,
#' so collinear proxies are rejected where marginal screens retain them. The
#' selector is embedded in a complete workflow — standard RF / BRT / MaxEnt /
#' GAM learners, nested spatial cross-validation, habitat-suitability
#' prediction, performance-weighted ensembles and future climate projection —
#' with a predictor-level audit of every retention decision. The package
#' supports interpretable conditional selection but does not claim automatic
#' discovery of ecological causes.
#'
#' **Pipeline steps:**
#'
#' 1. **Data Preparation**: train/test splitting and optional VIF collinearity
#'    screening ([cast_prepare()], [cast_vif()])
#' 2. **Variable Selection**: conditional predictive impact screening with
#'    FDR control ([cast_select()]); importance reporting via [cast_importance()]
#'    and sensitivity what-if maps ([cast_sensitivity()])
#' 3. **Model Fitting**: RF, BRT, MaxEnt, GAM ([cast_fit()])
#' 4. **Evaluation**: AUC, TSS, CBI metrics ([cast_evaluate()]), nested spatial
#'    cross-validation with fold-specific selection ([cast_cv()])
#' 5. **Prediction**: In-memory ([cast_predict()]) or tile-based
#'    ([cast_predict_tiled()]) spatial habitat suitability mapping
#' 6. **Ensemble**: performance-weighted, best-model, or equal-weight
#'    ensemble prediction ([cast_ensemble()])
#' 7. **Future Projection**: range-change analysis under climate scenarios
#'    ([cast_project()])
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
