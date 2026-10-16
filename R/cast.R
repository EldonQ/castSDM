#' Run the Full castSDM Pipeline
#'
#' One-step pipeline that executes the entire workflow: data splitting,
#' variable screening, model fitting, spatial cross-validation,
#' evaluation and ensemble prediction. Use [cast_attribution()] and
#' [cast_project()] to compare current and future environments afterwards.
#'
#' @param species_data A `data.frame` with columns: `lon`, `lat`, `presence`
#'   (0/1), and environmental variables.
#' @param env_data Optional `data.frame` of environmental variables for the
#'   full spatial grid (used for prediction). Must contain `lon`, `lat`, and
#'   the same environmental columns as `species_data`.
#' @param models Character vector of models to fit. Options: `"rf"`,
#'   `"maxent"`, `"brt"`, `"gam"`. Default `c("rf", "brt", "maxent", "gam")`.
#' @param train_fraction Numeric in (0, 1), strictly. Fraction of data for
#'   training. Values outside the open interval are rejected. Default `0.7`.
#' @param select_method Character. Variable screening method passed to
#'   [cast_select()]. Default `"two_stage"`.
#' @param select_num_trees Integer. Trees per forest in the stage-2 forward
#'   search. Default `300`.
#' @param select_metric Character. Inner-CV loss minimised by the forward
#'   search: `"brier"` (default) or `"auc"`.
#' @param select_tolerance Numeric. Minimum inner-CV loss improvement required
#'   to admit a predictor, on top of the 2-SE guard. Default `0`.
#' @param select_n_folds Integer. Inner folds for the stage-2 forward search.
#'   Default `3`.
#' @param select_corcut Correlation threshold for the `"covsel"` filter
#'   levels. Default `0.7` (Dormann et al. 2013).
#' @param select_var_id Character vector of variable-level identities aligned
#'   with the predictors (for `"covsel"`; `NULL` skips the variable level).
#' @param select_category Character vector of category-level identities
#'   aligned with the predictors (for `"covsel"`; `NULL` skips it).
#' @param select_algorithms Embedding learners for `"covsel"`: any of
#'   `"glm"`, `"gam"`, `"rf"`. Default all three.
#' @param select_ncov Target covariate count for `"covsel"`. Default `NULL` =
#'   `ceiling(log2(n_presences)) - 1` (Adde et al. 2023).
#' @param select_maxncov Maximum covariate count for `"covsel"`. Default `12`.
#' @param threshold_method Binary threshold rule, selected on training data and
#'   then applied to held-out predictions. See [cast_threshold()]. Default
#'   `"max_tss"`.
#' @param tune Logical. Per-engine hyperparameter grid search inside the
#'   nested spatial CV folds ([cast_cv()], argument `tune`). Default `FALSE`.
#' @param tune_folds Integer. Folds for the inner grid-search scoring.
#'   Default `3`.
#' @param prevalence_target Numeric in (0, 1), or `NULL`. Passed to
#'   [cast_fit()] (training fit, full-data refit) and [cast_cv()] (every
#'   fold): presence/background rows are weighted so the fitted prevalence
#'   equals the target (weight = target/n per class). Evaluation metrics
#'   remain unweighted. Default `NULL`.
#' @param num_threads Integer. Threads for the ranger learners. Default `1`.
#' @param do_cv Logical. Run spatial cross-validation. Default `TRUE`.
#' @param cv_k Integer. Number of spatial folds. Default `5`.
#' @param n_repeat Integer >= 1. Independent spatial-blocking repeats
#'   ([cast_cv()], argument `n_repeat`). Each repeat redraws the fold
#'   assignment; metrics are aggregated across repeats. Default `1`.
#' @param cv_block_method Character. Spatial blocking strategy. Default
#'   `"grid"` - grid cells grouped into spatially contiguous folds (the
#'   legacy interleaved packing is still available as `"grid_random"`; see
#'   [cast_cv()], which also offers a `buffer` exclusion band).
#' @param do_predict Logical. Generate spatial predictions. Default `TRUE`
#'   if `env_data` is provided.
#' @param do_ensemble Logical. Generate ensemble prediction. Default `TRUE`.
#' @param ensemble_method Character. Ensemble method: `"weighted"`,
#'   `"best"`, `"equal"`. Default `"weighted"`.
#' @param refit_full Logical. Refit the final models on the full data set
#'   (instead of the training split only) before spatial prediction and
#'   ensemble, standard SDM practice. Variable selection is still decided on
#'   the training split. Default `TRUE`.
#' @param seed Integer or `NULL`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#' @param select_keep Character vector passed to `keep` in [cast_select()],
#'   retained in the training screen and every nested CV fold. Prespecify these
#'   variables before evaluating outcomes; retention does not verify adjustment.
#'
#' @return A `cast_result` object (S3 class) containing pipeline outputs.
#'   Use [print()], [summary()], [plot()].
#'
#' @seealso [cast_select()], [cast_fit()], [cast_ensemble()],
#'   [cast_predict()]
#'
#' @export
cast <- function(species_data,
                 env_data = NULL,
                 models = c("rf", "brt", "maxent", "gam"),
                 train_fraction = 0.7,
                  select_method = "two_stage",
                  select_num_trees = 300L,
                   select_metric = "brier",
                   select_tolerance = 0,
                   select_n_folds = 3L,
                  select_corcut = 0.7,
                  select_var_id = NULL,
                  select_category = NULL,
                  select_algorithms = c("glm", "gam", "rf"),
                  select_ncov = NULL,
                  select_maxncov = 12L,
                  threshold_method = "max_tss",
                  tune = FALSE,
                  tune_folds = 3L,
                  prevalence_target = NULL,
                  num_threads = 1L,
                 do_cv = TRUE,
                 cv_k = 5L,
                 n_repeat = 1L,
                 cv_block_method = "grid",
                 do_predict = NULL,
                 do_ensemble = TRUE,
                 ensemble_method = "weighted",
                 refit_full = TRUE,
                 seed = NULL,
                 verbose = TRUE,
                 select_keep = character(0)) {
  cl <- match.call()
  do_predict <- do_predict %||% !is.null(env_data)

  if (verbose) cli::cli_h1("castSDM Pipeline")

  # === Step 1: Data Preparation ===
  if (verbose) cli::cli_h2("Step 1: Data Preparation")
  split <- cast_prepare(
    species_data, train_fraction = train_fraction, seed = seed,
    verbose = verbose
  )
  train_data <- split$train
  test_data  <- split$test
  if (verbose) {
    cli::cli_inform(
      "Train: {nrow(train_data)} | Test: {nrow(test_data)} | Vars: {length(split$env_vars)}"
    )
  }

  # === Step 2: Variable Selection ===
  if (verbose) cli::cli_h2("Step 2: Variable Selection")
  screen <- cast_select(
    train_data,
    method = select_method,
    num_trees = select_num_trees,
    metric = select_metric,
    tolerance = select_tolerance,
    n_folds = select_n_folds,
    corcut = select_corcut,
    var_id = select_var_id,
    category = select_category,
    algorithms = select_algorithms,
    ncov = select_ncov,
    maxncov = select_maxncov,
    keep = select_keep,
    seed = seed, verbose = verbose
  )

  # === Step 3: Model Fitting ===
  if (verbose) cli::cli_h2("Step 3: Model Fitting")
  fit <- cast_fit(
    train_data,
    screen = screen,
    models = models,
    tune = tune, tune_folds = tune_folds,
    prevalence_target = prevalence_target,
    num_threads = num_threads,
    seed = seed, verbose = verbose
  )

  # === Step 4: Cross-Validation ===
  cv_result <- NULL
  if (do_cv) {
    if (verbose) cli::cli_h2("Step 4: Spatial Cross-Validation")
    cv_select_args <- switch(select_method,
      covsel = list(corcut = select_corcut, var_id = select_var_id,
                    category = select_category, algorithms = select_algorithms,
                    ncov = select_ncov, maxncov = select_maxncov,
                    keep = select_keep),
      list(num_trees = select_num_trees, metric = select_metric,
           tolerance = select_tolerance, n_folds = select_n_folds,
           keep = select_keep))
    cv_result <- cast_cv(
        species_data,
        screen = screen,
        select_method = select_method,
         select_args = cv_select_args,
         k = cv_k, models = models,
         block_method = cv_block_method,
         n_repeat = n_repeat,
         threshold_method = threshold_method,
         tune = tune, tune_folds = tune_folds,
         prevalence_target = prevalence_target,
         num_threads = num_threads,
         seed = seed, verbose = verbose
    )
  }

  # === Step 5: Model Evaluation ===
  if (verbose) cli::cli_h2("Step 5: Model Evaluation")
  eval_result <- cast_evaluate(fit, test_data,
                               threshold_method = threshold_method)
  if (verbose) {
    if (!is.null(cv_result)) print(cv_result) else print(eval_result)
  }

  # === Step 6: Spatial Prediction & Ensemble ===
  pred_result <- NULL
  ensemble_result <- NULL
  fit_full <- NULL
  if (do_predict && !is.null(env_data)) {
    # Refit on the full data set for final maps: selection stays decided on
    # the training split, but published predictions reuse every record.
    if (isTRUE(refit_full)) {
      if (verbose) cli::cli_inform("Refitting final models on the full data set.")
      # Reuse the hyperparameters tuned on the training split: re-running the
      # grid on the full data would double grid compute and could silently
      # change the configuration of the published maps.
      tuned <- .cast_tuned_args(fit)
      fit_full <- cast_fit(
          species_data, screen = screen, models = models,
          rf_mtry = tuned$rf_mtry,
          brt_depth = tuned$brt_depth %||% 5L,
          brt_shrinkage = tuned$brt_shrinkage %||% 0.005,
          maxent_classes = tuned$maxent_classes,
          maxent_regmult = tuned$maxent_regmult,
          prevalence_target = prevalence_target,
          num_threads = num_threads, seed = seed, verbose = verbose
      )
    }
    fit_use <- if (!is.null(fit_full)) fit_full else fit

    if (verbose) cli::cli_h2("Step 6: Spatial Prediction")
    pred_result <- cast_predict(fit_use, env_data)
    if (verbose) {
      cli::cli_inform("Predicted {nrow(pred_result$predictions)} sites.")
    }

    if (do_ensemble && !is.null(cv_result)) {
      if (verbose) cli::cli_inform("Building ensemble ({ensemble_method})...")
      ensemble_result <- cast_ensemble(fit_use, cv_result, env_data,
                                       method = ensemble_method)
    } else if (do_ensemble) {
      if (verbose) {
        cli::cli_inform("Ensemble skipped: spatial CV unavailable (do_cv = FALSE or CV failed).")
      }
    }
  }

  if (verbose) cli::cli_h1("Pipeline Complete")

  new_cast_result(
    screen = screen,
    fit = fit, eval = eval_result,
    cv = cv_result,
    predict = pred_result,
    ensemble = ensemble_result,
    fit_full = fit_full,
    call = cl
  )
}

#' Extract Tuned Hyperparameters from a cast_fit
#'
#' Returns the fitted configuration, including explicit settings, so full-data
#' refitting preserves the configuration used for evaluation.
#'
#' @keywords internal
#' @noRd
.cast_tuned_args <- function(fit) {
  # Ask the registry which cast_fit() formals each engine contributes, and
  # in which field of `tune$best` the tuned value lives. This used to be a
  # chain of `if (mdl == "rf")` branches, so a new tunable engine silently
  # refit at its default unless this function was also edited.
  eng <- .cast_engines()
  out <- list()
  for (spec in eng) {
    if (length(spec$tune_args)) out[unname(spec$tune_args)] <- list(NULL)
  }
  for (mdl in names(fit$models %||% list())) {
    info <- fit$models[[mdl]]
    if (mdl == "rf" && !is.null(info$model$mtry)) out$rf_mtry <- info$model$mtry
    if (mdl == "brt" && !is.null(info$model$interaction.depth)) {
      out$brt_depth <- info$model$interaction.depth
      out$brt_shrinkage <- info$model$shrinkage
    }
    t <- fit$models[[mdl]]$tune
    spec <- eng[[mdl]]$tune_args
    if (is.null(t) || is.null(t$best) || is.null(spec)) next
    for (src in names(spec)) {
      val <- t$best[[src]]
      if (!is.null(val)) out[[unname(spec[[src]])]] <- val
    }
  }
  out
}
