#' Batch Multi-Species Modeling
#'
#' One-stop interface that runs the full castSDM pipeline on multiple
#' species, optionally in parallel. Every pipeline parameter can be
#' configured from this single function call. Each species is processed
#' independently and results (RDS + all diagnostic figures) are saved
#' to separate subdirectories. A comparative summary of model performance
#' across species is returned.
#'
#' @section Pipeline parameters:
#' All major stages of the castSDM pipeline are configurable through
#' explicit arguments. Parameters are grouped by prefix:
#' \itemize{
#'   \item **General**: `models`, `train_fraction`, `output_dir`,
#'     `fig_dpi`, `parallel`, `seed`, `verbose`.
#'     When `dev_package_root` is set and `parallel = TRUE`, species run on a
#'     [parallel::makeCluster()] backend: each worker calls [pkgload::load_all()]
#'     first (see `inst/examples/run_multi_species.R`). Without it, `parallel`
#'     uses [future.apply::future_lapply()] and the installed package must match
#'     your DAG API.
#'   \item **DAG** (`dag_*`): Bootstrap replicates, edge thresholds,
#'     structure learning algorithm and scoring criterion.
#'   \item **ATE** (`ate_*`): Cross-fitting folds, significance level,
#'     nuisance model trees, quantile discretization, Bonferroni,
#'     parallel ATE estimation.
#'   \item **Screen** (`screen_*`): Minimum retained variables, fraction,
#'     importance trees.
#'   \item **CV** (`cv_*`): Whether to run, folds, blocking method.
#'   \item **CATE** (`cate_*`): Whether to run, number of variables,
#'     causal forest trees.
#'   \item **Fit** (via `...`): Model hyper-parameters such as
#'     `n_epochs`, `n_runs`, `patience`, `rf_ntree`, `brt_n_trees`,
#'     `brt_depth`, `hidden_size`, `dropout`, `lr`, `batch_size`,
#'     `max_interactions`, `tune_grid` are forwarded to [cast_fit()].
#'     Use `fit_verbose` for [cast_fit()] console output (not `verbose` in `...`).
#'   \item **Spatial plots**: `plot_basemap` for CV / HSS / CATE maps.
#'   \item **SHAP** (`do_shap`, `shap_*`): optional XGBoost + fastshap figures.
#' }
#'
#' @param species_list A named list of `data.frame`s. Each element is a
#'   species dataset with `lon`, `lat`, `presence`, and environmental
#'   variables. Names are used as species identifiers.
#' @param env_data Optional shared `data.frame` for spatial prediction.
#'   If `NULL`, spatial prediction is skipped.
#' @param models Character vector. Models to fit per species. Default
#'   `c("cast", "rf", "maxent", "brt")`.
#' @param train_fraction Numeric. Train/test split ratio. Default `0.7`.
#' @param output_dir Character. Top-level directory for per-species
#'   outputs. Default `"castSDM_batch_output"`.
#' @param fig_dpi Integer. DPI for all saved figures. Default `300`.
#' @param parallel Logical. If `TRUE`, run species in parallel using
#'   \pkg{future}. Set a [future::plan()] before calling. Default `TRUE`.
#' @param seed Integer or `NULL`. Base seed. Each species gets
#'   `seed + species_index` for reproducibility.
#' @param verbose Logical. Progress messages for the batch driver itself.
#'   Default `TRUE`.
#' @param dev_package_root Character or `NULL`. Absolute path to the package
#'   source root (folder containing `DESCRIPTION`). When `parallel = TRUE`,
#'   each PSOCK worker runs [pkgload::load_all()] on this path before fitting
#'   (so code matches [devtools::load_all()] on the host). Example: pass the
#'   same path as in `inst/examples/run_multi_species.R` after locating the
#'   project root.
#' @param fit_verbose Logical. Passed to [cast_fit()] as `verbose` for each
#'   species. Default `FALSE`. (Do not pass `verbose` in `...` to `cast_batch`;
#'   it would clash with `verbose` here; use `fit_verbose` instead.)
#' @param response Character. Response column name in each species
#'   `data.frame`. Default `"presence"`.
#' @param prepare_env_vars,prepare_verbose Passed to [cast_prepare()] as
#'   `env_vars` and `verbose`. Defaults `NULL` and `FALSE`.
#' @param dag_env_vars,dag_verbose Passed to [cast_dag()] as `env_vars` and
#'   `verbose`. Defaults `NULL` and `FALSE`.
#' @param ate_variables,ate_verbose Passed to [cast_ate()] as `variables` and
#'   `verbose`. Defaults `NULL` and `FALSE`.
#' @param screen_verbose Passed to [cast_screen()] as `verbose`. Default `FALSE`.
#' @param cv_models Character vector or `NULL`. Models for [cast_cv()]. If
#'   `NULL`, uses the same vector as `models`.
#' @param cv_n_epochs,cv_n_runs,cv_rf_ntree,cv_brt_n_trees,cv_parallel,cv_verbose
#'   Passed to [cast_cv()]. Defaults match that function.
#' @param predict_models Passed to [cast_predict()] as `models`. Default `NULL`
#'   (all fitted models).
#' @param plot_basemap Character. `basemap` for spatial plots (CV, HSS, CATE).
#'   Default `"world"`.
#' @param cate_variables,cate_verbose Passed to [cast_cate()]. Defaults `NULL`
#'   and `FALSE`.
#' @param cate_point_size Numeric. Point size for saved [plot.cast_cate()] maps.
#'   Default `0.45` (matches `run_ovis_ammon.R`).
#' @param eval_response Character or `NULL`. Passed to [cast_evaluate()] as
#'   `response`. If `NULL`, uses `response`.
#' @param var_labels Named character vector or `NULL`. Optional labels for
#'   [plot.cast_dag()], [plot.cast_ate()], [plot.cast_screen()], [plot.cast_cate()].
#' @param do_shap Logical. If `TRUE`, saves SHAP figures (XGBoost surrogate,
#'   RF and CAST via [cast_shap_xgb()] / [cast_shap_fit()]) when dependencies
#'   are available. Default `FALSE` (can be slow in parallel batches).
#' @param shap_nrounds,shap_max_depth,shap_eta,shap_subsample,shap_colsample_bytree,shap_test_fraction,shap_verbose
#'   Passed to [cast_shap_xgb()].
#' @param shap_plot_top_n Integer. Passed to SHAP plot methods.
#' @param shap_fastshap_nsim,shap_max_explain_rows Passed to [cast_shap_fit()].
#'
#' @param dag_R Integer. Number of DAG bootstrap replicates.
#'   Default `100`.
#' @param dag_structure_method Character. Passed to [cast_dag()] as
#'   `structure_method`. Default `"bootstrap_hc"`.
#' @param dag_pc_alpha PC alpha level. Default `0.05`.
#' @param dag_pc_test Passed to [cast_dag()] as `pc_test`. Default `"zf"`.
#' @param dag_bidag_algorithm,dag_bidag_iterations BiDAG controls.
#' @param dag_notears_lambda NOTEARS L1 penalty (see \code{cast_dag()}).
#' @param dag_notears_max_iter Maximum NOTEARS optimization steps.
#' @param dag_notears_lr Adam learning rate for NOTEARS.
#' @param dag_notears_tol Acyclicity tolerance for NOTEARS.
#' @param dag_notears_rho_init Initial augmented-Lagrangian rho for NOTEARS.
#' @param dag_notears_alpha_mult Rho multiplier for NOTEARS.
#' @param dag_algorithm Character. Passed to [cast_dag()] as `algorithm`.
#'   **Only** used when `dag_structure_method = "bootstrap_hc"` (score-based
#'   learner inside each bootstrap, e.g. `"hc"`, `"tabu"`). Default `"hc"`.
#' @param dag_score Character. Passed to [cast_dag()] as `score` for
#'   bootstrap HC only. Default `"bic-g"`.
#' @param dag_strength_threshold Numeric. Minimum bootstrap edge
#'   strength to retain. Default `0.7`.
#' @param dag_direction_threshold Numeric. Minimum direction
#'   consistency to retain. Default `0.6`.
#' @param dag_max_rows Integer. Subsample size for DAG learning.
#'   Default `8000`.
#' @param cate_hss_model,cate_hss_threshold When saving CATE maps, pass to
#'   [plot.cast_cate()] to mask by `HSS_<model> >= threshold`. Defaults
#'   `cast` and `0.1`. Set `cate_hss_model = NULL` to disable masking.
#'
#' @param ate_K Integer. DML cross-fitting folds. Default `5`.
#' @param ate_alpha Numeric. Significance level for ATE tests.
#'   Default `0.05`.
#' @param ate_num_trees Integer. Trees per nuisance RF model.
#'   Default `300`.
#' @param ate_quantile_cuts Numeric vector. Quantile thresholds for
#'   binarization. Default `c(0.25, 0.50, 0.75)`.
#' @param ate_bonferroni Logical. Apply Bonferroni correction.
#'   Default `TRUE`.
#' @param ate_parallel Logical. Parallel ATE estimation within each
#'   species (via \pkg{future}). Default `FALSE`.
#'
#' @param screen_min_vars Integer. Minimum variables to retain.
#'   Default `5`.
#' @param screen_min_fraction Numeric. Minimum fraction of variables.
#'   Default `0.5`.
#' @param screen_num_trees Integer. RF trees for importance ranking.
#'   Default `300`.
#'
#' @param do_cv Logical. Run spatial block CV. Default `TRUE`.
#' @param cv_k Integer. Number of spatial CV folds. Default `5`.
#' @param cv_block_method Character. Blocking strategy: `"grid"` or
#'   `"cluster"`. Default `"grid"`.
#'
#' @param do_cate Logical. Estimate CATE via causal forests.
#'   Default `TRUE`.
#' @param cate_top_n Integer. Number of top variables for CATE.
#'   Default `3`.
#' @param cate_n_trees Integer. Trees in each causal forest.
#'   Default `1000`.
#'
#' @param ... Additional arguments forwarded to [cast_fit()] for
#'   model hyper-parameter tuning (e.g., `n_epochs`, `n_runs`,
#'   `rf_ntree`, `brt_n_trees`, `tune_grid`).
#'
#' @return A `cast_batch` object with components:
#' \describe{
#'   \item{`species_metrics`}{`data.frame` with per-species per-model
#'     evaluation metrics (from spatial CV when available).}
#'   \item{`species`}{Character vector of species names.}
#'   \item{`models`}{Character vector of model names.}
#'   \item{`results`}{Named list of per-species result lists.}
#'   \item{`output_dir`}{Output directory path.}
#' }
#'
#' @seealso [cast()], [plot.cast_batch()], [cast_consistency()],
#'   [cast_fit()], [cast_dag()], [cast_ate()], [cast_screen()],
#'   [cast_cv()], [cast_cate()], [cast_shap_xgb()], [cast_shap_fit()]
#'
#' @export
cast_batch <- function(species_list,
                       env_data    = NULL,
                       models      = c("cast", "rf", "maxent", "brt"),
                       train_fraction = 0.7,
                       output_dir  = "castSDM_batch_output",
                       fig_dpi     = 300L,
                       parallel    = TRUE,
                       seed        = NULL,
                       verbose     = TRUE,
                       fit_verbose = FALSE,
                       # ── DAG ──
                       dag_R                  = 100L,
                       dag_structure_method   = "bootstrap_hc",
                       dag_pc_alpha           = 0.05,
                       dag_pc_test            = "zf",
                       dag_bidag_algorithm    = "order",
                       dag_bidag_iterations   = NULL,
                       dag_notears_lambda     = 0.03,
                       dag_notears_max_iter   = 2000L,
                       dag_notears_lr         = 0.02,
                       dag_notears_tol        = 1e-3,
                       dag_notears_rho_init   = 0.1,
                       dag_notears_alpha_mult = 1.01,
                       dag_algorithm          = "hc",
                       dag_score              = "bic-g",
                       dag_strength_threshold = 0.7,
                       dag_direction_threshold = 0.6,
                       dag_max_rows           = 8000L,
                       # ── ATE ──
                       ate_K              = 5L,
                       ate_alpha          = 0.05,
                       ate_num_trees      = 300L,
                       ate_quantile_cuts  = c(0.25, 0.50, 0.75),
                       ate_bonferroni     = TRUE,
                       ate_parallel       = FALSE,
                       # ── Screen ──
                       screen_min_vars     = 5L,
                       screen_min_fraction = 0.5,
                       screen_num_trees    = 300L,
                       # ── CV ──
                       do_cv           = TRUE,
                       cv_k            = 5L,
                       cv_block_method = "grid",
                       # ── CATE ──
                       do_cate      = TRUE,
                       cate_top_n   = 3L,
                       cate_n_trees = 1000L,
                       cate_hss_model = "cast",
                       cate_hss_threshold = 0.1,
                       response = "presence",
                       prepare_env_vars = NULL,
                       prepare_verbose = FALSE,
                       dag_env_vars = NULL,
                       dag_verbose = FALSE,
                       ate_variables = NULL,
                       ate_verbose = FALSE,
                       screen_verbose = FALSE,
                       cv_models = NULL,
                       cv_n_epochs = 200L,
                       cv_n_runs = 3L,
                       cv_rf_ntree = 300L,
                       cv_brt_n_trees = 500L,
                       cv_parallel = FALSE,
                       cv_verbose = FALSE,
                       predict_models = NULL,
                       plot_basemap = "world",
                       cate_variables = NULL,
                       cate_verbose = FALSE,
                       cate_point_size = 0.45,
                       eval_response = NULL,
                       var_labels = NULL,
                       do_shap = FALSE,
                       shap_nrounds = 200L,
                       shap_max_depth = 6L,
                       shap_eta = 0.05,
                       shap_subsample = 0.8,
                       shap_colsample_bytree = 0.8,
                       shap_test_fraction = 0.2,
                       shap_verbose = FALSE,
                       shap_plot_top_n = 15L,
                       shap_fastshap_nsim = 40L,
                       shap_max_explain_rows = 50L,
                       dev_package_root = NULL,
                       ...) {

  if (!is.list(species_list) || is.null(names(species_list))) {
    cli::cli_abort("{.arg species_list} must be a named list of data.frames.")
  }
  sp_names <- names(species_list)
  n_sp <- length(sp_names)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (verbose) {
    cli::cli_h1("CAST Batch: {n_sp} species")
    cli::cli_inform("Models: {.val {models}}")
    cli::cli_inform("Output: {output_dir}")
  }

  dev_root_workers <- NULL
  if (!is.null(dev_package_root) && nzchar(as.character(dev_package_root)[1])) {
    dev_root_workers <- tryCatch(
      normalizePath(
        as.character(dev_package_root)[1],
        winslash = "/", mustWork = FALSE
      ),
      error = function(e) NA_character_
    )
    if (!is.na(dev_root_workers) && nzchar(dev_root_workers) &&
        file.exists(file.path(dev_root_workers, "DESCRIPTION"))) {
      Sys.setenv(CASTSDM_ROOT = dev_root_workers)
    } else {
      dev_root_workers <- NULL
    }
  }

  if (isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE) &&
      is.null(dev_root_workers) && !nzchar(Sys.getenv("CASTSDM_ROOT", ""))) {
    root <- tryCatch(
      normalizePath(
        as.character(getNamespaceInfo(asNamespace("castSDM"), "path")),
        winslash = "/", mustWork = FALSE
      ),
      error = function(e) NA_character_
    )
    if (!is.na(root) && nzchar(root) &&
        file.exists(file.path(root, "DESCRIPTION"))) {
      np <- root
      in_lib <- any(vapply(.libPaths(), function(lib) {
        nl <- tryCatch(
          normalizePath(lib, winslash = "/", mustWork = FALSE),
          error = function(e) ""
        )
        nzchar(nl) && (identical(np, nl) || startsWith(np, paste0(nl, "/")))
      }, logical(1L)))
      if (!isTRUE(in_lib)) Sys.setenv(CASTSDM_ROOT = root)
    }
  }
  if (verbose && isTRUE(parallel) && nzchar(Sys.getenv("CASTSDM_ROOT", ""))) {
    cli::cli_inform(
      "Parallel: {.envvar CASTSDM_ROOT} = {.path {Sys.getenv('CASTSDM_ROOT')}}"
    )
  }

  fit_args <- list(...)
  if ("verbose" %in% names(fit_args)) {
    cli::cli_abort(
      "{.arg verbose} in {.arg ...} is ambiguous: use {.arg fit_verbose} for {.fn cast_fit} logging."
    )
  }
  cast_fit_names <- names(formals(cast_fit))
  if (length(fit_args)) {
    bad <- !names(fit_args) %in% cast_fit_names
    if (any(bad)) {
      cli::cli_warn(
        "Dropping {.arg ...} names not in {.fn cast_fit}: {.val {names(fit_args)[bad]}} (use {.code devtools::load_all()} on current source if you meant batch options)."
      )
      fit_args <- fit_args[!bad]
    }
  }

  eval_resp <- if (is.null(eval_response)) response else eval_response
  cv_models_use <- if (is.null(cv_models)) models else cv_models

  # Pack all pipeline config so workers receive a single list
  cfg <- list(
    response = response,
    eval_response = eval_resp,
    prepare_env_vars = prepare_env_vars,
    prepare_verbose = prepare_verbose,
    train_fraction = train_fraction,
    dag_env_vars = dag_env_vars,
    dag_verbose = dag_verbose,
    dag_R = dag_R,
    dag_structure_method = dag_structure_method,
    dag_pc_alpha = dag_pc_alpha,
    dag_pc_test = dag_pc_test,
    dag_bidag_algorithm = dag_bidag_algorithm,
    dag_bidag_iterations = dag_bidag_iterations,
    dag_notears_lambda = dag_notears_lambda,
    dag_notears_max_iter = dag_notears_max_iter,
    dag_notears_lr = dag_notears_lr,
    dag_notears_tol = dag_notears_tol,
    dag_notears_rho_init = dag_notears_rho_init,
    dag_notears_alpha_mult = dag_notears_alpha_mult,
    dag_algorithm = dag_algorithm,
    dag_score = dag_score,
    dag_strength_threshold = dag_strength_threshold,
    dag_direction_threshold = dag_direction_threshold,
    dag_max_rows = dag_max_rows,
    ate_variables = ate_variables,
    ate_verbose = ate_verbose,
    ate_K = ate_K, ate_alpha = ate_alpha,
    ate_num_trees = ate_num_trees,
    ate_quantile_cuts = ate_quantile_cuts,
    ate_bonferroni = ate_bonferroni,
    ate_parallel = ate_parallel,
    screen_verbose = screen_verbose,
    screen_min_vars = screen_min_vars,
    screen_min_fraction = screen_min_fraction,
    screen_num_trees = screen_num_trees,
    do_cv = do_cv, cv_k = cv_k, cv_block_method = cv_block_method,
    cv_models = cv_models_use,
    cv_n_epochs = cv_n_epochs,
    cv_n_runs = cv_n_runs,
    cv_rf_ntree = cv_rf_ntree,
    cv_brt_n_trees = cv_brt_n_trees,
    cv_parallel = cv_parallel,
    cv_verbose = cv_verbose,
    predict_models = predict_models,
    plot_basemap = plot_basemap,
    do_cate = do_cate, cate_top_n = cate_top_n,
    cate_n_trees = cate_n_trees,
    cate_variables = cate_variables,
    cate_verbose = cate_verbose,
    cate_point_size = cate_point_size,
    cate_hss_model = cate_hss_model,
    cate_hss_threshold = cate_hss_threshold,
    var_labels = var_labels,
    do_shap = do_shap,
    shap_nrounds = shap_nrounds,
    shap_max_depth = shap_max_depth,
    shap_eta = shap_eta,
    shap_subsample = shap_subsample,
    shap_colsample_bytree = shap_colsample_bytree,
    shap_test_fraction = shap_test_fraction,
    shap_verbose = shap_verbose,
    shap_plot_top_n = shap_plot_top_n,
    shap_fastshap_nsim = shap_fastshap_nsim,
    shap_max_explain_rows = shap_max_explain_rows,
    fit_verbose = fit_verbose
  )

  if (parallel && !is.null(dev_root_workers)) {
    nwrk <- tryCatch(
      if (requireNamespace("future", quietly = TRUE)) {
        future::nbrOfWorkers()
      } else {
        NA_integer_
      },
      error = function(e) NA_integer_
    )
    if (is.na(nwrk) || nwrk < 1L) {
      nwrk <- max(1L, parallel::detectCores() - 1L)
    }
    nwrk <- min(as.integer(nwrk), n_sp)
    if (verbose) {
      cli::cli_inform(
        "Parallel (dev): PSOCK with {nwrk} workers; pkgload::load_all() on each."
      )
    }
    root_value <- as.character(dev_root_workers)[1L]
    cl <- NULL
    cl <- parallel::makeCluster(nwrk)
    results <- tryCatch(
      {
        parallel::clusterExport(cl, "root_value", envir = environment())
        ok <- parallel::clusterEvalQ(cl, {
          if (!requireNamespace("pkgload", quietly = TRUE)) {
            FALSE
          } else {
            suppressPackageStartupMessages(pkgload::load_all(root_value, quiet = TRUE))
            TRUE
          }
        })
        if (!all(vapply(ok, isTRUE, logical(1L)))) {
          cli::cli_abort(
            "Parallel with {.arg dev_package_root} requires {.pkg pkgload} on workers."
          )
        }
        eb <- environment()
        parallel::clusterExport(
          cl,
          c(
            "species_list", "sp_names", "env_data", "models", "output_dir",
            "fig_dpi", "seed", "cfg", "fit_args"
          ),
          envir = eb
        )
        parallel::parLapply(cl, seq_along(sp_names), function(ii) {
          sp <- sp_names[[ii]]
          sd <- species_list[[sp]]
          seed_i <- if (!is.null(seed)) seed + ii else NULL
          worker_run <- utils::getFromNamespace(
            ".cast_batch_run_one_species", "castSDM"
          )
          worker_run(
            sp, sd, env_data, models,
            output_dir, fig_dpi, seed_i,
            cfg, fit_args, FALSE, NULL
          )
        })
      },
      finally = if (!is.null(cl)) parallel::stopCluster(cl)
    )
  } else if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    if (verbose) cli::cli_inform("Running in parallel (future)...")
    results <- future.apply::future_lapply(
      seq_along(sp_names),
      function(ii) {
        sp <- sp_names[[ii]]
        sd <- species_list[[sp]]
        seed_i <- if (!is.null(seed)) seed + ii else NULL
        .cast_batch_run_one_species(
          sp, sd, env_data, models,
          output_dir, fig_dpi, seed_i,
          cfg, fit_args, TRUE, dev_root_workers
        )
      },
      future.seed = TRUE
    )
  } else {
    results <- vector("list", n_sp)
    names(results) <- sp_names
    for (i in seq_along(sp_names)) {
      sp <- sp_names[i]
      if (verbose) cli::cli_inform("[{i}/{n_sp}] Processing {.val {sp}}...")
      seed_i <- if (!is.null(seed)) seed + i else NULL
      results[[i]] <- .cast_batch_run_one_species(
        sp, species_list[[sp]], env_data, models, output_dir,
        fig_dpi, seed_i, cfg, fit_args, FALSE, dev_root_workers
      )
    }
  }

  names(results) <- sp_names

  # Collect metrics across species
  metrics_rows <- list()
  for (sp in sp_names) {
    r <- results[[sp]]
    if (is.null(r)) next

    if (!is.null(r$cv) && !is.null(r$cv$fold_metrics) &&
        nrow(r$cv$fold_metrics) > 0) {
      fm <- r$cv$fold_metrics
      fm$species <- sp
      metrics_rows[[sp]] <- fm
    } else if (!is.null(r$eval) && !is.null(r$eval$metrics)) {
      em <- r$eval$metrics
      em$species <- sp
      for (mcol in c("auc", "tss", "cbi", "sedi", "kappa", "prauc")) {
        mean_col <- paste0(mcol, "_mean")
        if (mean_col %in% names(em)) em[[mcol]] <- em[[mean_col]]
      }
      em$fold <- 0L
      metrics_rows[[sp]] <- em[, intersect(
        c("fold", "model", "auc", "tss", "cbi", "sedi", "kappa",
          "prauc", "species"),
        names(em)
      ), drop = FALSE]
    }
  }

  species_metrics <- if (length(metrics_rows) > 0) {
    do.call(rbind, metrics_rows)
  } else {
    data.frame()
  }
  rownames(species_metrics) <- NULL

  n_ok <- sum(!vapply(results, is.null, logical(1)))
  if (verbose) {
    cli::cli_inform("Batch complete: {n_ok}/{n_sp} species succeeded.")
  }

  new_cast_batch(
    species_metrics = species_metrics,
    species         = sp_names,
    models          = models,
    results         = results,
    output_dir      = output_dir
  )
}


#' Plot Multi-Species Model Performance Comparison
#'
#' Generates boxplots with jittered points comparing model performance
#' metrics across species. Uses a clean black-and-white style with
#' grayscale fill.
#'
#' @param x A `cast_batch` object.
#' @param metrics Character vector. Metrics to show. Default
#'   `c("auc", "tss", "cbi")`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (faceted by metric).
#' @export
plot.cast_batch <- function(x, metrics = c("auc", "tss", "cbi"), ...) {
  check_suggested("ggplot2", "for plotting")

  sm <- x$species_metrics
  if (is.null(sm) || nrow(sm) == 0) {
    cli::cli_abort("No species metrics available for plotting.")
  }

  present_metrics <- intersect(metrics, names(sm))
  if (length(present_metrics) == 0) {
    cli::cli_abort("None of the requested metrics found in species_metrics.")
  }

  long_rows <- list()
  for (mc in present_metrics) {
    long_rows[[mc]] <- data.frame(
      species = sm$species,
      model   = sm$model,
      metric  = toupper(mc),
      value   = sm[[mc]],
      stringsAsFactors = FALSE
    )
  }
  long <- do.call(rbind, long_rows)
  long <- long[!is.na(long$value), ]

  long$metric <- factor(long$metric, levels = toupper(present_metrics))
  long$model  <- factor(long$model)

  n_models <- length(levels(long$model))
  gray_fills <- grDevices::gray.colors(n_models, start = 0.4, end = 0.9)

  ggplot2::ggplot(long, ggplot2::aes(
    x = .data$model, y = .data$value, fill = .data$model
  )) +
    ggplot2::geom_boxplot(
      width = 0.6, outlier.shape = NA, alpha = 0.8,
      color = "black", linewidth = 0.4
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(shape = .data$species),
      width = 0.15, size = 1.8, alpha = 0.7, color = "black"
    ) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = gray_fills, guide = "none") +
    ggplot2::scale_shape_manual(
      values = c(16, 17, 15, 18, 8, 3, 4, 1, 2, 0),
      name = "Species"
    ) +
    ggplot2::labs(
      title    = "Multi-species Model Performance Comparison",
      subtitle = sprintf("%d species | Spatial CV", length(x$species)),
      x = "", y = "Score"
    ) +
    ggplot2::theme_minimal(base_size = 11, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(
        fill = NA, color = "black", linewidth = 0.5
      ),
      strip.text         = ggplot2::element_text(face = "bold", size = 10),
      axis.title         = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle      = ggplot2::element_text(
        hjust = 0.5, color = "grey40", size = 9
      ),
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position    = "bottom",
      legend.text        = ggplot2::element_text(size = 8)
    )
}
