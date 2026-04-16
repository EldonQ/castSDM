#' Run the Full CAST Pipeline
#'
#' One-step pipeline that executes the entire CAST workflow: data splitting,
#' DAG learning, ATE estimation, variable screening, causal role assignment,
#' feature engineering, model fitting, evaluation, and optionally spatial
#' prediction and CATE estimation.
#'
#' @param species_data A `data.frame` with columns: `lon`, `lat`, `presence`
#'   (0/1), and environmental variables.
#' @param env_data Optional `data.frame` of environmental variables for the
#'   full spatial grid (used for prediction). Must contain `lon`, `lat`, and
#'   the same environmental columns as `species_data`.
#' @param models Character vector of models to fit. Options: `"cast"`,
#'   `"mlp_ate"`, `"mlp"`, `"rf"`, `"maxent"`, `"brt"`. Default is
#'   `c("cast", "rf", "maxent", "brt")`.
#' @param train_fraction Numeric. Fraction of data for training. Default `0.7`.
#' @param n_bootstrap Integer. Number of bootstrap replicates for DAG when
#'   `dag_structure_method = "bootstrap_hc"`. Default `100`.
#' @param dag_structure_method Character passed to [cast_dag()] as
#'   `structure_method`. Default `"bootstrap_hc"`. Alternatives: `"pc"`,
#'   `"bidag_bge"`, `"notears_linear"`.
#' @param dag_pc_alpha Significance level for PC.
#' @param dag_pc_test Conditional-independence test passed to
#'   [cast_dag()] as `pc_test` (default `"zf"`).
#' @param dag_bidag_algorithm,dag_bidag_iterations BiDAG options.
#' @param dag_notears_lambda NOTEARS L1 penalty; passed to \code{cast_dag()}.
#' @param dag_notears_max_iter Maximum NOTEARS optimization steps.
#' @param dag_notears_lr Adam learning rate for NOTEARS.
#' @param dag_notears_tol Acyclicity tolerance for NOTEARS.
#' @param dag_notears_rho_init Initial augmented-Lagrangian rho for NOTEARS.
#' @param dag_notears_alpha_mult Rho multiplier (every 100 steps) for NOTEARS.
#' @param strength_threshold Numeric. Minimum edge strength. Default `0.7`.
#' @param direction_threshold Numeric. Minimum direction consistency. Default
#'   `0.6`.
#' @param ate_folds Integer. DML cross-fitting folds. Default `2`.
#' @param ate_alpha Numeric. ATE significance level. Default `0.05`.
#' @param screen_min_vars Integer. Minimum retained variables. Default `5`.
#' @param do_cv Logical. Whether to run spatial k-fold cross-validation for
#'   honest model evaluation. Default `TRUE`. Results are stored in
#'   `$cv` of the returned object. Setting `FALSE` falls back to a single
#'   random hold-out evaluation via [cast_evaluate()].
#' @param cv_k Integer. Number of spatial folds. Default `5`. Pass `3` for
#'   small datasets (< 100 presences) or `10` for large ones.
#' @param cv_block_method Character. Spatial blocking strategy: `"grid"`
#'   (default) or `"cluster"`. See [cast_cv()] for details.
#' @param do_predict Logical. Whether to generate spatial predictions using
#'   `env_data`. Default `TRUE` if `env_data` is provided.
#' @param do_cate Logical. Whether to estimate spatial CATE. Default `FALSE`.
#' @param cate_top_n Integer. Top variables for CATE. Default `3`.
#' @param do_shap Logical. Whether to compute SHAP explanations after model
#'   fitting. Produces up to three `cast_shap` objects (XGBoost surrogate, RF,
#'   and CAST) stored in the returned result. Default `FALSE`.
#' @param shap_nrounds Integer. XGBoost boosting rounds for
#'   [cast_shap_xgb()]. Default `400`.
#' @param shap_test_fraction Numeric. Fraction held out for SHAP
#'   visualisation. Default `0.2`.
#' @param shap_fastshap_nsim Integer. Monte Carlo reps per feature in
#'   [cast_shap_fit()]. Default `60`.
#' @param shap_max_explain_rows Integer. Cap on explained rows in
#'   [cast_shap_fit()]. Default `80`.
#' @param blacklist A `data.frame` with columns `from` and `to` specifying
#'   forbidden directed edges in DAG learning. Default `NULL`.
#' @param whitelist A `data.frame` with columns `from` and `to` specifying
#'   required directed edges in DAG learning. Default `NULL`.
#' @param do_evalue Logical. Compute E-value sensitivity analysis on ATE
#'   estimates. Default `FALSE`.
#' @param do_backdoor Logical. Run backdoor criterion check via dagitty.
#'   Default `FALSE`.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A `cast_result` object (S3 class) containing all pipeline outputs.
#'   Use [print()], [summary()], and [plot()] for inspection.
#'
#' @seealso [cast_dag()], [cast_ate()], [cast_screen()], [cast_fit()],
#'   [cast_predict()], [cast_cate()], [cast_shap_xgb()], [cast_shap_fit()]
#'
#' @export
cast <- function(species_data,
                 env_data = NULL,
                 models = c("cast", "rf", "maxent", "brt"),
                 train_fraction = 0.7,
                 n_bootstrap = 100L,
                 dag_structure_method = "bootstrap_hc",
                 dag_pc_alpha = 0.05,
                 dag_pc_test = "zf",
                 dag_bidag_algorithm = "order",
                 dag_bidag_iterations = NULL,
                 dag_notears_lambda = 0.03,
                 dag_notears_max_iter = 2000L,
                 dag_notears_lr = 0.02,
                 dag_notears_tol = 1e-3,
                 dag_notears_rho_init = 0.1,
                 dag_notears_alpha_mult = 1.01,
                 strength_threshold = 0.7,
                 direction_threshold = 0.6,
                 ate_folds = 2L,
                 ate_alpha = 0.05,
                 screen_min_vars = 5L,
                 do_cv = TRUE,
                 cv_k = 5L,
                 cv_block_method = "grid",
                 do_predict = NULL,
                 do_cate = FALSE,
                 cate_top_n = 3L,
                 do_shap = FALSE,
                 shap_nrounds = 400L,
                 shap_test_fraction = 0.2,
                 shap_fastshap_nsim = 60L,
                 shap_max_explain_rows = 80L,
                 blacklist = NULL,
                 whitelist = NULL,
                 do_evalue = FALSE,
                 do_backdoor = FALSE,
                 seed = NULL,
                 verbose = TRUE) {
  do_predict <- do_predict %||% !is.null(env_data)

  if (verbose) cli::cli_h1("CAST Pipeline")

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
    R = n_bootstrap,
    strength_threshold = strength_threshold,
    direction_threshold = direction_threshold,
    seed = seed,
    verbose = verbose,
    structure_method = dag_structure_method,
    pc_alpha = dag_pc_alpha,
    pc_test = dag_pc_test,
    bidag_algorithm = dag_bidag_algorithm,
    bidag_iterations = dag_bidag_iterations,
    notears_lambda = dag_notears_lambda,
    notears_max_iter = dag_notears_max_iter,
    notears_lr = dag_notears_lr,
    notears_tol = dag_notears_tol,
    notears_rho_init = dag_notears_rho_init,
    notears_alpha_mult = dag_notears_alpha_mult,
    blacklist = blacklist,
    whitelist = whitelist
  )

  # === Step 3: ATE Estimation ===
  if (verbose) cli::cli_h2("Step 3: ATE Estimation")
  ate <- cast_ate(
    train_data,
    K = ate_folds, alpha = ate_alpha,
    seed = seed, verbose = verbose
  )

  # === Step 3b: E-value Sensitivity (optional) ===
  if (do_evalue) {
    if (verbose) cli::cli_h2("Step 3b: E-value Sensitivity Analysis")
    ate <- cast_evalue(ate, verbose = verbose)
  }

  # === Step 3c: Backdoor Criterion (optional) ===
  backdoor <- NULL
  if (do_backdoor) {
    if (verbose) cli::cli_h2("Step 3c: Backdoor Criterion Check")
    backdoor <- cast_backdoor(dag, verbose = verbose)
  }

  # === Step 4: Variable Screening ===
  if (verbose) cli::cli_h2("Step 4: Variable Screening")
  screen <- cast_screen(
    dag, ate, train_data,
    min_vars = screen_min_vars,
    seed = seed, verbose = verbose
  )

  # === Step 5: Causal Roles ===
  if (verbose) cli::cli_h2("Step 5: Causal Role Assignment")
  roles <- cast_roles(screen, dag)
  if (verbose) {
    role_tbl <- table(roles$roles$role)
    cli::cli_inform(
      "Roles: {paste(names(role_tbl), role_tbl, sep='=', collapse=', ')}"
    )
  }

  # === Step 6: Model Fitting ===
  if (verbose) cli::cli_h2("Step 6: Model Fitting")
  fit <- cast_fit(
    train_data,
    screen = screen, dag = dag, ate = ate,
    models = models,
    seed = seed, verbose = verbose
  )

  # === Step 7: Model Evaluation ===
  if (verbose) cli::cli_h2("Step 7: Model Evaluation")

  cv_result <- NULL
  if (do_cv) {
    # Spatial k-fold CV on the full dataset (preferred: more data, honest geo)
    if (verbose) {
      cli::cli_inform(
        "Running spatial {cv_k}-fold CV ({cv_block_method} blocks)..."
      )
    }
    cv_result <- tryCatch(
      cast_cv(
        species_data,
        screen       = screen,
        dag          = dag,
        ate          = ate,
        k            = cv_k,
        models       = models,
        block_method = cv_block_method,
        seed         = seed,
        verbose      = FALSE
      ),
      error = function(e) {
        cli::cli_warn(
          "Spatial CV failed ({e$message}). Falling back to hold-out eval."
        )
        NULL
      }
    )
  }

  # Always compute hold-out eval on the random test split (complementary)
  eval_result <- cast_evaluate(fit, test_data)
  if (verbose) {
    if (!is.null(cv_result)) print(cv_result) else print(eval_result)
  }

  # === Step 8: Spatial Prediction (optional) ===
  pred_result <- NULL
  if (do_predict && !is.null(env_data)) {
    if (verbose) cli::cli_h2("Step 8: Spatial Prediction")
    pred_result <- cast_predict(fit, env_data)
    if (verbose) {
      cli::cli_inform(
        "Predicted {nrow(pred_result$predictions)} sites."
      )
    }
  }

  # === Step 9: CATE (optional) ===
  cate_result <- NULL
  if (do_cate) {
    if (verbose) cli::cli_h2("Step 9: CATE Estimation")
    pred_for_cate <- if (!is.null(env_data)) env_data else NULL
    cate_result <- cast_cate(
      train_data,
      ate = ate, screen = screen,
      top_n = cate_top_n,
      predict_data = pred_for_cate,
      seed = seed, verbose = verbose
    )
  }

  # === Step 10: SHAP Explanations (optional) ===
  shap_result <- NULL
  if (do_shap) {
    if (verbose) cli::cli_h2("Step 10: SHAP Explanations")
    shap_result <- list()

    # XGBoost surrogate SHAP
    if (requireNamespace("xgboost", quietly = TRUE)) {
      if (verbose) cli::cli_inform("  Computing XGBoost TreeSHAP...")
      shap_result$xgb <- tryCatch(
        cast_shap_xgb(
          train_data,
          response = "presence",
          dag = dag,
          nrounds = shap_nrounds,
          test_fraction = shap_test_fraction,
          seed = seed,
          verbose = FALSE
        ),
        error = function(e) {
          if (verbose) cli::cli_warn("XGBoost SHAP failed: {e$message}")
          NULL
        }
      )
    }

    # RF SHAP (fastshap)
    if (requireNamespace("fastshap", quietly = TRUE) &&
        "rf" %in% names(fit$models) &&
        !is.null(fit$models$rf$model)) {
      if (verbose) cli::cli_inform("  Computing RF SHAP (fastshap)...")
      shap_result$rf <- tryCatch(
        cast_shap_fit(
          fit = fit,
          which = "rf",
          data = train_data,
          test_fraction = shap_test_fraction,
          seed = seed,
          fastshap_nsim = shap_fastshap_nsim,
          max_explain_rows = shap_max_explain_rows,
          verbose = FALSE
        ),
        error = function(e) {
          if (verbose) cli::cli_warn("RF SHAP failed: {e$message}")
          NULL
        }
      )
    }

    # CAST SHAP (fastshap on CI-MLP)
    torch_ok <- requireNamespace("torch", quietly = TRUE) &&
      tryCatch(torch::torch_is_installed(), error = function(e) FALSE)
    if (requireNamespace("fastshap", quietly = TRUE) &&
        isTRUE(torch_ok) &&
        "cast" %in% names(fit$models) &&
        !is.null(fit$models$cast$model)) {
      if (verbose) cli::cli_inform("  Computing CAST SHAP (fastshap)...")
      shap_result$cast <- tryCatch(
        cast_shap_fit(
          fit = fit,
          which = "cast",
          data = train_data,
          test_fraction = shap_test_fraction,
          seed = seed,
          fastshap_nsim = shap_fastshap_nsim,
          max_explain_rows = shap_max_explain_rows,
          verbose = FALSE
        ),
        error = function(e) {
          if (verbose) cli::cli_warn("CAST SHAP failed: {e$message}")
          NULL
        }
      )
    }

    n_ok <- sum(!vapply(shap_result, is.null, logical(1)))
    if (n_ok == 0L) shap_result <- NULL
    if (verbose && !is.null(shap_result)) {
      cli::cli_inform(
        "SHAP: {n_ok} model(s) explained ({.val {names(Filter(Negate(is.null), shap_result))}})"
      )
    }
  }

  if (verbose) cli::cli_h1("Pipeline Complete")

  new_cast_result(
    dag = dag, ate = ate, screen = screen, roles = roles,
    fit = fit, eval = eval_result,
    cv = cv_result,
    predict = pred_result, cate = cate_result,
    shap = shap_result,
    backdoor = backdoor
  )
}
