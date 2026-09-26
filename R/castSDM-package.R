#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' \pkg{castSDM} is a species-distribution-modelling workflow built around one
#' question: if a driver moved, what would happen to predicted suitability, and
#' where is that answer supported by data?
#'
#' Variable selection ([cast_select()]) answers *which predictors to carry*.
#' Stage 1 thins collinear predictors. Stage 2 runs a spatial forward search:
#' starting from prespecified predictors (if any), it repeatedly admits the
#' survivor that most improves the inner cross-validated loss of a probability
#' random forest — spatial folds when coordinates exist — and stops when
#' nothing improves it further (Meyer et al. 2018, 2019). There is no
#' predictor-count cap and no fallback set: an empty selection honestly means
#' nothing passed admission (admission needs at least two finite paired folds
#' and sufficient gain). Selection serves parsimony
#' for interpretation and projection, not causal discovery.
#'
#' The effect products ([cast_effect_table()], [cast_effect_map()])
#' report the interventional response to a single raw-unit shift, with
#' hard range-masking: rows or cells outside the training quantile box are
#' masked, not ranked.
#'
#' None of these is a claim of automatic causal discovery. Causal interpretation
#' requires consistency, a justified adjustment set, no uncontrolled confounding,
#' joint positivity, and adequate response and observation models. Quantile-box
#' coverage flags range extrapolation but cannot establish joint positivity.
#' Even agreement between products does not establish a causal driver.
#'
#' **Pipeline steps:**
#'
#' 1. **Data Preparation**: train/test splitting and optional VIF collinearity
#'    screening ([cast_prepare()], [cast_vif()])
#' 2. **Variable Selection**: collinearity thinning, then spatial forward
#'    selection on inner-CV loss ([cast_select()]), reported as an admission
#'    path via [cast_importance()]
#' 3. **Model Fitting**: RF, BRT, MaxEnt, GAM ([cast_fit()])
#' 4. **Evaluation**: AUC, TSS, CBI metrics ([cast_evaluate()]), nested spatial
#'    cross-validation with fold-specific selection ([cast_cv()])
#' 5. **Prediction**: In-memory ([cast_predict()]) or tile-based
#'    ([cast_predict_tiled()]) spatial habitat suitability mapping
#' 6. **Ensemble**: performance-weighted, best-model, or equal-weight
#'    ensemble prediction ([cast_ensemble()])
#' 7. **Future Projection**: range-change analysis under climate scenarios
#'    ([cast_project()])
#' 8. **Interventional attribution**: single-shift effect sizes and maps
#'    with hard range-masking
#'    ([cast_effect_table()], [cast_effect_map()])
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
