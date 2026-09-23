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
#' @param select_method Character. Variable screening method passed to
#'   [cast_select()]. Default `"two_stage"`.
#' @param select_num_trees Integer. Trees per forest in the stage-2 importance
#'   filter. Default `300`.
#' @param select_n_perm Integer. Response permutations used to build the
#'   stage-2 null distribution. Default `49`.
#' @param select_ncov Integer or `NULL`. Maximum predictors retained by
#'   [cast_select()], ordered by the interventional effect. Default `NULL`
#'   selects `ceiling(log2(n_presence))`. Use `Inf` for no cap.
#' @param select_maxncov Integer. Upper bound on the automatic `select_ncov`.
#'   Default `12`.
#' @param num_threads Integer. Threads for the ranger learners. Default `1`.
#' @param do_cv Logical. Run spatial cross-validation. Default `TRUE`.
#' @param cv_k Integer. Number of spatial folds. Default `5`.
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
                  select_n_perm = 49L,
                  select_ncov = NULL,
                  select_maxncov = 12L,
                 num_threads = 1L,
                 do_cv = TRUE,
                 cv_k = 5L,
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
    n_perm = select_n_perm,
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
    num_threads = num_threads,
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
          num_trees = select_num_trees,
          n_perm = select_n_perm,
          ncov = select_ncov,
          maxncov = select_maxncov,
          keep = select_keep
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
  fit_full <- NULL
  if (do_predict && !is.null(env_data)) {
    # Refit on the full data set for final maps: selection stays decided on
    # the training split, but published predictions reuse every record.
    if (isTRUE(refit_full)) {
      if (verbose) cli::cli_inform("Refitting final models on the full data set.")
      fit_full <- tryCatch(
        cast_fit(
          species_data, screen = screen, models = models,
          num_threads = num_threads, seed = seed, verbose = verbose
        ),
        error = function(e) {
          cli::cli_warn("Full-data refit failed ({e$message}); using the training fit.")
          NULL
        }
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
      ensemble_result <- tryCatch(
        cast_ensemble(fit_use, cv_result, env_data, method = ensemble_method),
        error = function(e) {
          cli::cli_warn("Ensemble failed: {e$message}")
          NULL
        }
      )
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
