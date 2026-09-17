#' @keywords internal
"_PACKAGE"

#' @section Overview:
#' \pkg{castSDM} is a species-distribution-modelling workflow built around one
#' question: if a driver moved, what would happen to predicted suitability, and
#' where is that answer supported by data?
#'
#' Variable selection ([cast_select()]) answers *which predictors to carry*.
#' Stage 1 thins collinear predictors. Stage 2 ranks the survivors by the
#' **interventional effect**: the change each one causes in the fitted
#' probability when it is shifted while every other predictor stays at its
#' observed value, calibrated against a permuted-response null. The same
#' forest's permutation importance is reported alongside as a diagnostic,
#' because permuting a predictor breaks its correlation with the others and is
#' therefore governed by the model's extrapolation behaviour rather than by the
#' predictor's influence (Hooker, Mentch & Zhou 2021).
#'
#' The effect products ([cast_effect_table()], [cast_effect_map()],
#' [cast_dose_response()], [cast_sensitivity()],
#' [cast_effect_support()]) report the interventional response, its shape over
#' the size of the shift, and the fraction of observed predictor vectors still
#' inside the training support after the shift. The knockout audit
#' ([cast_necessity()]) asks the complementary question: does the model still
#' discriminate without this predictor?
#'
#' None of these is a claim of automatic causal discovery. The effect products
#' are model-based interventional contrasts under stated assumptions, and the
#' support column reports only the positivity assumption; the others are
#' assumptions, not results. Even agreement between products does not establish
#' a causal driver.
#'
#' **Pipeline steps:**
#'
#' 1. **Data Preparation**: train/test splitting and optional VIF collinearity
#'    screening ([cast_prepare()], [cast_vif()])
#' 2. **Variable Selection**: collinearity thinning, then an interventional
#'    effect filter against a permuted-response null ([cast_select()]),
#'    reported via [cast_importance()] with both attribution columns
#' 3. **Model Fitting**: RF, BRT, MaxEnt, GAM ([cast_fit()])
#' 4. **Evaluation**: AUC, TSS, CBI metrics ([cast_evaluate()]), nested spatial
#'    cross-validation with fold-specific selection ([cast_cv()])
#' 5. **Prediction**: In-memory ([cast_predict()]) or tile-based
#'    ([cast_predict_tiled()]) spatial habitat suitability mapping
#' 6. **Ensemble**: performance-weighted, best-model, or equal-weight
#'    ensemble prediction ([cast_ensemble()])
#' 7. **Future Projection**: range-change analysis under climate scenarios
#'    ([cast_project()])
#' 8. **Interventional attribution**: effect sizes and maps
#'    ([cast_effect_table()], [cast_effect_map()]), response shape
#'    ([cast_dose_response()]), positivity
#'    ([cast_effect_support()]), paired with knockout necessity
#'    ([cast_necessity()])
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
