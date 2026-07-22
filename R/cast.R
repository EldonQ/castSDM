#' Run the Full castSDM Pipeline
#'
#' One-step pipeline that executes the entire workflow: data splitting,
#' variable screening, model fitting, spatial cross-validation,
#' evaluation, ensemble prediction, and future projection.
#'
#' @param species_data A `data.frame` with columns: `lon`, `lat`, `presence`
#'   (0/1), and environmental variables.
#' @param env_data Optional `data.frame` of environmental variables for the
#'   full spatial grid (used for prediction). Must contain `lon`, `lat`, and
#'   the same environmental columns as `species_data`.
#' @param models Character vector of models to fit. Options: `"rf"`,
#'   `"maxent"`, `"brt"`, `"gam"`. Default `c("rf", "brt", "maxent", "gam")`.
#' @param train_fraction Numeric. Fraction of data for training. Default `0.7`.
#' @param select_min_vars Integer. Minimum retained variables. Default `3`.
#' @param select_method Character. Variable screening method passed to
#'   [cast_select()]. Default `"cpi"`.
#' @param select_num_trees Integer. Trees in the RF nuisance/benchmark forests.
#'   Default `300`.
#' @param select_max_vars Candidate ceiling for DML testing / RF output.
#'   Default `30`.
#' @param select_alpha Numeric. FDR level for the DML selector. Default `0.05`.
#' @param select_cor_threshold Numeric. Absolute correlation threshold for the
#'   RF benchmark. Default `0.8`.
#' @param do_cv Logical. Run spatial cross-validation. Default `TRUE`.
#' @param cv_k Integer. Number of spatial folds. Default `5`.
#' @param cv_block_method Character. Spatial blocking strategy. Default
#'   `"grid"`.
#' @param do_predict Logical. Generate spatial predictions. Default `TRUE`
#'   if `env_data` is provided.
#' @param do_ensemble Logical. Generate ensemble prediction. Default `TRUE`.
#' @param ensemble_method Character. Ensemble method: `"weighted"`,
#'   `"best"`, `"equal"`. Default `"weighted"`.
#' @param seed Integer or `NULL`.
#' @param verbose Logical. Print progress. Default `TRUE`.
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
                 select_min_vars = 3L,
                 select_method = "cpi",
                 select_num_trees = 300L,
                 select_max_vars = 30L,
                 select_alpha = 0.05,
                 select_cor_threshold = 0.8,
                 do_cv = TRUE,
                 cv_k = 5L,
                 cv_block_method = "grid",
                 do_predict = NULL,
                 do_ensemble = TRUE,
                 ensemble_method = "weighted",
                 seed = NULL,
                 verbose = TRUE) {
  do_predict <- do_predict %||% !is.null(env_data)

  if (verbose) cli::cli_h1("castSDM Pipeline")

  # === Step 1: Data Preparation ===
  if (verbose) cli::cli_h2("Step 1: Data Preparation")
  split <- cast_prepare(
    species_data, train_fraction = train_fraction, seed = seed
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
    alpha = select_alpha,
    min_vars = select_min_vars,
    num_trees = select_num_trees,
    max_candidates = select_max_vars,
    cor_threshold = select_cor_threshold,
    seed = seed, verbose = verbose
  )

  # === Step 3: Model Fitting ===
  if (verbose) cli::cli_h2("Step 3: Model Fitting")
  fit <- cast_fit(
    train_data,
    screen = screen,
    models = models,
    seed = seed, verbose = verbose
  )

  # === Step 4: Cross-Validation ===
  cv_result <- NULL
  if (do_cv) {
    if (verbose) cli::cli_h2("Step 4: Spatial Cross-Validation")
    cv_result <- tryCatch(
      cast_cv(
        species_data,
        screen = screen,
        select_method = select_method,
        select_args = list(
          alpha = select_alpha,
          min_vars = select_min_vars,
          max_candidates = select_max_vars,
          num_trees = select_num_trees,
          cor_threshold = select_cor_threshold
        ),
        k = cv_k, models = models,
        block_method = cv_block_method,
        seed = seed, verbose = verbose
      ),
      error = function(e) {
        cli::cli_warn(
          "Spatial CV failed ({e$message}). Falling back to hold-out eval."
        )
        NULL
      }
    )
  }

  # === Step 5: Model Evaluation ===
  if (verbose) cli::cli_h2("Step 5: Model Evaluation")
  eval_result <- cast_evaluate(fit, test_data)
  if (verbose) {
    if (!is.null(cv_result)) print(cv_result) else print(eval_result)
  }

  # === Step 6: Spatial Prediction & Ensemble ===
  pred_result <- NULL
  ensemble_result <- NULL
  if (do_predict && !is.null(env_data)) {
    if (verbose) cli::cli_h2("Step 6: Spatial Prediction")
    pred_result <- cast_predict(fit, env_data)
    if (verbose) {
      cli::cli_inform("Predicted {nrow(pred_result$predictions)} sites.")
    }

    if (do_ensemble && !is.null(cv_result)) {
      if (verbose) cli::cli_inform("Building ensemble ({ensemble_method})...")
      ensemble_result <- tryCatch(
        cast_ensemble(fit, cv_result, env_data, method = ensemble_method),
        error = function(e) {
          cli::cli_warn("Ensemble failed: {e$message}")
          NULL
        }
      )
    }
  }

  if (verbose) cli::cli_h1("Pipeline Complete")

  new_cast_result(
    screen = screen,
    fit = fit, eval = eval_result,
    cv = cv_result,
    predict = pred_result,
    ensemble = ensemble_result
  )
}
