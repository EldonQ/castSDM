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
#'   `"fci"`, `"bidag_bge"`, `"notears_linear"`.
#' @param dag_pc_alpha,dag_fci_alpha Significance levels for PC / FCI.
#' @param dag_pc_test,dag_fci_test Conditional-independence tests passed to
#'   [cast_dag()] as `pc_test` / `fci_test` (default `"zf"`).
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
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A `cast_result` object (S3 class) containing all pipeline outputs.
#'   Use [print()], [summary()], and [plot()] for inspection.
#'
#' @seealso [cast_dag()], [cast_ate()], [cast_screen()], [cast_fit()],
#'   [cast_predict()], [cast_cate()]
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
                 dag_fci_alpha = 0.05,
                 dag_fci_test = "zf",
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
    fci_alpha = dag_fci_alpha,
    fci_test = dag_fci_test,
    bidag_algorithm = dag_bidag_algorithm,
    bidag_iterations = dag_bidag_iterations,
    notears_lambda = dag_notears_lambda,
    notears_max_iter = dag_notears_max_iter,
    notears_lr = dag_notears_lr,
    notears_tol = dag_notears_tol,
    notears_rho_init = dag_notears_rho_init,
    notears_alpha_mult = dag_notears_alpha_mult
  )

  # === Step 3: ATE Estimation ===
  if (verbose) cli::cli_h2("Step 3: ATE Estimation")
  ate <- cast_ate(
    train_data,
    K = ate_folds, alpha = ate_alpha,
    seed = seed, verbose = verbose
  )

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

  if (verbose) cli::cli_h1("Pipeline Complete")

  new_cast_result(
    dag = dag, ate = ate, screen = screen, roles = roles,
    fit = fit, eval = eval_result,
    cv = cv_result,
    predict = pred_result, cate = cate_result
  )
}
