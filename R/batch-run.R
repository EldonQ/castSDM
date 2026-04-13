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
#' @param verbose Logical. Default `TRUE`.
#'
#' @param dag_R Integer. Number of DAG bootstrap replicates.
#'   Default `100`.
#' @param dag_structure_method Character. Passed to [cast_dag()] as
#'   `structure_method`. Default `"bootstrap_hc"`.
#' @param dag_pc_alpha,dag_fci_alpha PC/FCI alpha levels. Defaults `0.05`.
#' @param dag_bidag_algorithm,dag_bidag_iterations BiDAG controls.
#' @param dag_notears_lambda NOTEARS L1 penalty (see \code{cast_dag()}).
#' @param dag_notears_max_iter Maximum NOTEARS optimization steps.
#' @param dag_notears_lr Adam learning rate for NOTEARS.
#' @param dag_notears_tol Acyclicity tolerance for NOTEARS.
#' @param dag_notears_rho_init Initial augmented-Lagrangian rho for NOTEARS.
#' @param dag_notears_alpha_mult Rho multiplier for NOTEARS.
#' @param dag_algorithm Character. Structure learning algorithm
#'   (`"hc"`). Default `"hc"`.
#' @param dag_score Character. Scoring criterion for structure
#'   learning. Default `"bic-g"`.
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
#'   [cast_cv()], [cast_cate()]
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
                       # ── DAG ──
                       dag_R                  = 100L,
                       dag_structure_method   = "bootstrap_hc",
                       dag_pc_alpha           = 0.05,
                       dag_fci_alpha          = 0.05,
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

  fit_args <- list(...)

  # Pack all pipeline config so workers receive a single list
  cfg <- list(
    train_fraction = train_fraction,
    dag_R = dag_R,
    dag_structure_method = dag_structure_method,
    dag_pc_alpha = dag_pc_alpha,
    dag_fci_alpha = dag_fci_alpha,
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
    ate_K = ate_K, ate_alpha = ate_alpha,
    ate_num_trees = ate_num_trees,
    ate_quantile_cuts = ate_quantile_cuts,
    ate_bonferroni = ate_bonferroni,
    ate_parallel = ate_parallel,
    screen_min_vars = screen_min_vars,
    screen_min_fraction = screen_min_fraction,
    screen_num_trees = screen_num_trees,
    do_cv = do_cv, cv_k = cv_k, cv_block_method = cv_block_method,
    do_cate = do_cate, cate_top_n = cate_top_n,
    cate_n_trees = cate_n_trees,
    cate_hss_model = cate_hss_model,
    cate_hss_threshold = cate_hss_threshold
  )

  run_one_species <- function(sp_name, sp_data, env_data, models,
                              output_dir, fig_dpi, seed_i,
                              cfg, fit_args) {
    sp_dir  <- file.path(output_dir, sp_name)
    fig_dir <- file.path(sp_dir, "figures")
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

    save_fig <- function(p, fname, w = 10, h = 7) {
      if (is.null(p)) return(invisible(NULL))
      tryCatch(
        ggplot2::ggsave(file.path(fig_dir, fname), p,
                        width = w, height = h, dpi = fig_dpi),
        error = function(e) NULL
      )
    }

    result <- tryCatch({
      split <- cast_prepare(sp_data, train_fraction = cfg$train_fraction,
                            seed = seed_i)

      dag <- cast_dag(
        split$train,
        R = cfg$dag_R,
        algorithm = cfg$dag_algorithm,
        score = cfg$dag_score,
        strength_threshold = cfg$dag_strength_threshold,
        direction_threshold = cfg$dag_direction_threshold,
        max_rows = cfg$dag_max_rows,
        seed = seed_i,
        verbose = FALSE,
        structure_method = cfg$dag_structure_method,
        pc_alpha = cfg$dag_pc_alpha,
        fci_alpha = cfg$dag_fci_alpha,
        bidag_algorithm = cfg$dag_bidag_algorithm,
        bidag_iterations = cfg$dag_bidag_iterations,
        notears_lambda = cfg$dag_notears_lambda,
        notears_max_iter = cfg$dag_notears_max_iter,
        notears_lr = cfg$dag_notears_lr,
        notears_tol = cfg$dag_notears_tol,
        notears_rho_init = cfg$dag_notears_rho_init,
        notears_alpha_mult = cfg$dag_notears_alpha_mult
      )

      ate <- cast_ate(split$train, K = cfg$ate_K,
                       alpha = cfg$ate_alpha,
                       num_trees = cfg$ate_num_trees,
                       quantile_cuts = cfg$ate_quantile_cuts,
                       bonferroni = cfg$ate_bonferroni,
                       parallel = cfg$ate_parallel,
                       seed = seed_i, verbose = FALSE)

      screen <- cast_screen(dag, ate, split$train,
                             min_vars = cfg$screen_min_vars,
                             min_fraction = cfg$screen_min_fraction,
                             num_trees = cfg$screen_num_trees,
                             seed = seed_i, verbose = FALSE)

      roles <- cast_roles(screen, dag)

      fit_call_args <- c(
        list(data = split$train, screen = screen, dag = dag, ate = ate,
             models = models, seed = seed_i, verbose = FALSE),
        fit_args
      )
      fit <- do.call(cast_fit, fit_call_args)

      eval_result <- cast_evaluate(fit, split$test)

      cv_result <- NULL
      if (cfg$do_cv) {
        cv_result <- tryCatch(
          cast_cv(sp_data, screen = screen, dag = dag, ate = ate,
                  k = cfg$cv_k, models = models,
                  block_method = cfg$cv_block_method,
                  seed = seed_i, verbose = FALSE),
          error = function(e) NULL
        )
      }

      pred_result <- NULL
      if (!is.null(env_data)) {
        pred_result <- tryCatch(
          cast_predict(fit, env_data),
          error = function(e) NULL
        )
      }

      cate_result <- NULL
      if (cfg$do_cate) {
        cate_result <- tryCatch({
          check_suggested("grf", "for CATE")
          pred_for_cate <- if (!is.null(env_data)) env_data else NULL
          cast_cate(split$train, ate = ate, screen = screen,
                    predict_data = pred_for_cate, top_n = cfg$cate_top_n,
                    n_trees = cfg$cate_n_trees,
                    seed = seed_i, verbose = FALSE)
        }, error = function(e) NULL)
      }

      cons_result <- NULL
      if (!is.null(pred_result) && length(pred_result$models) >= 2) {
        cons_result <- tryCatch(cast_consistency(pred_result),
                                error = function(e) NULL)
      }

      # ── Save ALL plots ────────────────────────────────────────────────
      check_suggested("ggplot2", "for plotting")

      if (requireNamespace("ggraph", quietly = TRUE) &&
          requireNamespace("igraph", quietly = TRUE)) {
        p <- tryCatch(plot(dag, roles = roles, screen = screen,
                           species = sp_name), error = function(e) NULL)
        save_fig(p, "causal_dag.png", w = 14, h = 10)
      }

      p <- tryCatch(plot(ate), error = function(e) NULL)
      save_fig(p, "ate_forest_plot.png", w = 10, h = 7)

      p <- tryCatch(plot(screen), error = function(e) NULL)
      save_fig(p, "variable_screening.png", w = 10, h = 7)

      p <- tryCatch(plot(eval_result), error = function(e) NULL)
      save_fig(p, "model_evaluation.png", w = 10, h = 6)

      if (!is.null(cv_result)) {
        p <- tryCatch(
          plot(cv_result, lon = sp_data$lon, lat = sp_data$lat,
               metric = "auc", basemap = "world"),
          error = function(e) NULL
        )
        save_fig(p, "spatial_cv_map.png", w = 14, h = 8)
      }

      if (!is.null(pred_result) && requireNamespace("sf", quietly = TRUE)) {
        for (mdl in pred_result$models) {
          p <- tryCatch(
            plot(pred_result, model = mdl, basemap = "world",
                 title = sprintf("%s HSS (%s)",
                                 gsub("_", " ", sp_name), mdl)),
            error = function(e) NULL
          )
          save_fig(p, sprintf("HSS_%s.png", mdl), w = 14, h = 8)
        }
      }

      if (!is.null(cate_result) && requireNamespace("sf", quietly = TRUE)) {
        for (cv in cate_result$variables) {
          p <- tryCatch({
            pm <- cfg$cate_hss_model
            if (!is.null(pm) && !is.null(pred_result) &&
                pm %in% pred_result$models) {
              plot(
                cate_result,
                variable = cv,
                species = sp_name,
                basemap = "world",
                legend_position = "bottom",
                hss_predict = pred_result,
                hss_model = pm,
                hss_threshold = cfg$cate_hss_threshold
              )
            } else {
              plot(
                cate_result,
                variable = cv,
                species = sp_name,
                basemap = "world",
                legend_position = "bottom"
              )
            }
          }, error = function(e) NULL)
          save_fig(p, sprintf("CATE_%s.png", cv), w = 14, h = 8)
        }
      }

      if (!is.null(cons_result) &&
          requireNamespace("patchwork", quietly = TRUE)) {
        p <- tryCatch(plot(cons_result, species = sp_name),
                       error = function(e) NULL)
        save_fig(p, "model_consistency.png", w = 16, h = 5.5)
      }

      # ── Save RDS ──────────────────────────────────────────────────────
      sp_result <- list(
        dag = dag, ate = ate, screen = screen, roles = roles,
        fit = fit, eval = eval_result, cv = cv_result,
        predict = pred_result, cate = cate_result,
        consistency = cons_result
      )
      saveRDS(sp_result, file.path(sp_dir, "cast_result.rds"))

      sp_result
    }, error = function(e) {
      warning(sprintf("Species '%s' failed: %s", sp_name, e$message))
      NULL
    })

    result
  }

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    if (verbose) cli::cli_inform("Running in parallel...")
    results <- future.apply::future_mapply(
      run_one_species,
      sp_name  = sp_names,
      sp_data  = species_list,
      seed_i   = if (!is.null(seed)) seed + seq_along(sp_names) else
                   rep(list(NULL), n_sp),
      MoreArgs = list(
        env_data   = env_data,
        models     = models,
        output_dir = output_dir,
        fig_dpi    = fig_dpi,
        cfg        = cfg,
        fit_args   = fit_args
      ),
      SIMPLIFY = FALSE,
      future.seed = TRUE
    )
  } else {
    results <- vector("list", n_sp)
    names(results) <- sp_names
    for (i in seq_along(sp_names)) {
      sp <- sp_names[i]
      if (verbose) cli::cli_inform("[{i}/{n_sp}] Processing {.val {sp}}...")
      seed_i <- if (!is.null(seed)) seed + i else NULL
      results[[i]] <- run_one_species(
        sp, species_list[[sp]], env_data, models, output_dir,
        fig_dpi, seed_i, cfg, fit_args
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
