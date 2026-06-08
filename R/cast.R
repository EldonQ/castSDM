#' Run the Full castSDM Pipeline
#'
#' One-step pipeline that executes the entire workflow: data splitting,
#' DAG learning (with response node), DAG-guided variable selection via
#' Markov Blanket + RF importance, model fitting, spatial cross-validation,
#' evaluation, ensemble prediction, and optionally CATE estimation.
#'
#' @param species_data A `data.frame` with columns: `lon`, `lat`, `presence`
#'   (0/1), and environmental variables.
#' @param env_data Optional `data.frame` of environmental variables for the
#'   full spatial grid (used for prediction). Must contain `lon`, `lat`, and
#'   the same environmental columns as `species_data`.
#' @param models Character vector of models to fit. Options: `"rf"`,
#'   `"maxent"`, `"brt"`, `"gam"`. Default `c("rf", "brt", "maxent", "gam")`.
#' @param train_fraction Numeric. Fraction of data for training. Default `0.7`.
#' @param n_bootstrap Integer. Number of bootstrap replicates for DAG when
#'   `dag_structure_method = "bootstrap_hc"`. Default `100`.
#' @param dag_structure_method Character passed to [cast_dag()] as
#'   `structure_method`. Default `"mb_first"`.
#' @param dag_include_response Logical. Include response in DAG learning for
#'   Markov Blanket extraction. Default `TRUE`.
#' @param dag_response_as_sink Logical. If `TRUE` (default), forbids
#'   `presence -> environmental predictor` edges during DAG learning.
#' @param dag_pc_alpha Significance level for PC.
#' @param dag_pc_test Conditional-independence test for [cast_dag()].
#'   Default `NULL` (auto-selects `"mi-cg"` for mixed data).
#' @param dag_mb_method Character. MB discovery algorithm for
#'   `structure_method = "mb_first"`. Default `"fast.iamb"`.
#' @param dag_mb_alpha Numeric. Significance level for MB discovery stage.
#'   Default `0.05`.
#' @param dag_bidag_algorithm,dag_bidag_iterations BiDAG options.
#' @param strength_threshold Numeric. Minimum edge strength. Default `0.7`.
#' @param direction_threshold Numeric. Minimum direction consistency. Default
#'   `0.6`.
#' @param select_min_vars Integer. Minimum retained variables. Default `5`.
#' @param select_min_fraction Numeric. Minimum fraction of variables. Default
#'   `0.3`.
#' @param select_stability_reps Integer. Bootstrap repetitions for selection
#'   stability. Default `0`.
#' @param select_stability_threshold Numeric. Stability frequency threshold.
#'   Default `0.6`.
#' @param do_refute Logical. Run lightweight screen refutation diagnostics.
#'   Default `TRUE`.
#' @param refute_reps Integer. Number of refutation repetitions.
#' @param do_cv Logical. Run spatial cross-validation. Default `TRUE`.
#' @param cv_k Integer. Number of spatial folds. Default `5`.
#' @param cv_block_method Character. Spatial blocking strategy. Default
#'   `"grid"`.
#' @param do_predict Logical. Generate spatial predictions. Default `TRUE`
#'   if `env_data` is provided.
#' @param do_ensemble Logical. Generate ensemble prediction. Default `TRUE`.
#' @param ensemble_method Character. Ensemble method: `"weighted"`,
#'   `"best"`, `"equal"`. Default `"weighted"`.
#' @param do_cate Logical. Estimate spatial CATE via causal forests.
#'   Default `FALSE`. Requires package `grf`.
#' @param cate_top_n Integer. Top variables for CATE. Default `3`.
#' @param cate_n_trees Integer. Number of causal forest trees. Default
#'   `1000`.
#' @param cate_variables Character vector or `NULL`. Specific variables
#'   for CATE estimation. Default `NULL` (auto-select top variables).
#' @param blacklist,whitelist DAG edge constraints. Default `NULL`.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A `cast_result` object (S3 class) containing all pipeline outputs.
#'   Use [print()], [summary()], and [plot()] for inspection.
#'
#' @seealso [cast_dag()], [cast_select()], [cast_fit()], [cast_ensemble()],
#'   [cast_predict()], [cast_cate()]
#'
#' @export
cast <- function(species_data,
                 env_data = NULL,
                 models = c("rf", "brt", "maxent", "gam"),
                 train_fraction = 0.7,
                 n_bootstrap = 100L,
                 dag_structure_method = "mb_first",
                 dag_include_response = TRUE,
                 dag_response_as_sink = TRUE,
                 dag_pc_alpha = 0.05,
                 dag_pc_test = NULL,
                 dag_mb_method = "fast.iamb",
                 dag_mb_alpha = 0.05,
                 dag_bidag_algorithm = "order",
                 dag_bidag_iterations = NULL,
                 strength_threshold = 0.7,
                 direction_threshold = 0.6,
                 select_min_vars = 5L,
                 select_min_fraction = 0.3,
                 select_stability_reps = 0L,
                 select_stability_threshold = 0.6,
                 do_refute = TRUE,
                 refute_reps = 20L,
                 do_cv = TRUE,
                 cv_k = 5L,
                 cv_block_method = "grid",
                 do_predict = NULL,
                 do_ensemble = TRUE,
                 ensemble_method = "weighted",
                 do_cate = FALSE,
                 cate_top_n = 3L,
                 cate_n_trees = 1000L,
                 cate_variables = NULL,
                 blacklist = NULL,
                 whitelist = NULL,
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
  test_data <- split$test
  if (verbose) {
    cli::cli_inform(
      "Train: {nrow(train_data)} | Test: {nrow(test_data)} | Vars: {length(split$env_vars)}"
    )
  }

  # === Step 2: DAG Learning ===
  if (verbose) cli::cli_h2("Step 2: DAG Learning")
  dag <- cast_dag(
    train_data,
    include_response = dag_include_response,
    response_as_sink = dag_response_as_sink,
    R = n_bootstrap,
    strength_threshold = strength_threshold,
    direction_threshold = direction_threshold,
    seed = seed,
    verbose = verbose,
    structure_method = dag_structure_method,
    pc_alpha = dag_pc_alpha,
    pc_test = dag_pc_test,
    mb_method = dag_mb_method,
    mb_alpha = dag_mb_alpha,
    bidag_algorithm = dag_bidag_algorithm,
    bidag_iterations = dag_bidag_iterations,
    blacklist = blacklist,
    whitelist = whitelist
  )

  # === Step 3: Variable Selection (DAG MB + RF Importance) ===
  if (verbose) cli::cli_h2("Step 3: Variable Selection")
  screen <- cast_select(
    dag, train_data,
    min_vars = select_min_vars,
    min_fraction = select_min_fraction,
    stability_reps = select_stability_reps,
    stability_threshold = select_stability_threshold,
    seed = seed, verbose = verbose
  )

  refute_result <- NULL
  if (isTRUE(do_refute)) {
    if (verbose) cli::cli_h2("Step 3b: Screen Refutation")
    refute_result <- tryCatch(
      cast_refute_screen(
        dag, screen, train_data,
        reps = refute_reps,
        seed = seed,
        verbose = verbose
      ),
      error = function(e) {
        cli::cli_warn("Screen refutation failed: {e$message}")
        NULL
      }
    )
  }

  # === Step 4: Model Fitting ===
  if (verbose) cli::cli_h2("Step 4: Model Fitting")
  fit <- cast_fit(
    train_data,
    screen = screen, dag = dag,
    models = models,
    seed = seed, verbose = verbose
  )

  # === Step 5: Cross-Validation ===
  cv_result <- NULL
  if (do_cv) {
    if (verbose) cli::cli_h2("Step 5: Spatial Cross-Validation")
    cv_result <- tryCatch(
      cast_cv(
        species_data,
        screen = screen, dag = dag,
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

  # === Step 6: Model Evaluation ===
  if (verbose) cli::cli_h2("Step 6: Model Evaluation")
  eval_result <- cast_evaluate(fit, test_data)
  if (verbose) {
    if (!is.null(cv_result)) print(cv_result) else print(eval_result)
  }

  # === Step 7: Spatial Prediction & Ensemble ===
  pred_result <- NULL
  ensemble_result <- NULL
  if (do_predict && !is.null(env_data)) {
    if (verbose) cli::cli_h2("Step 7: Spatial Prediction")
    pred_result <- cast_predict(fit, env_data)
    if (verbose) {
      cli::cli_inform("Predicted {nrow(pred_result$predictions)} sites.")
    }

    # Ensemble prediction (requires CV for weights)
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

  # === Step 8: CATE (optional) ===
  cate_result <- NULL
  if (do_cate) {
    if (verbose) cli::cli_h2("Step 8: CATE Estimation")
    pred_for_cate <- if (!is.null(env_data)) env_data else NULL
    cate_result <- tryCatch(
      cast_cate(
        train_data,
        dag = dag, screen = screen,
        variables = cate_variables,
        top_n = cate_top_n,
        n_trees = cate_n_trees,
        predict_data = pred_for_cate,
        seed = seed, verbose = verbose
      ),
      error = function(e) {
        cli::cli_warn("CATE failed: {e$message}")
        NULL
      }
    )
  }

  if (verbose) cli::cli_h1("Pipeline Complete")

  new_cast_result(
    dag = dag, screen = screen,
    fit = fit, eval = eval_result,
    cv = cv_result,
    predict = pred_result,
    ensemble = ensemble_result,
    cate = cate_result,
    refute = refute_result
  )
}
