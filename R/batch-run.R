#' Batch Multi-Species Modeling
#'
#' Runs the full castSDM pipeline on multiple species, optionally in
#' parallel. Each species is processed independently and results are saved
#' to separate subdirectories with all diagnostic plots. A comparative
#' summary of model performance across species is returned.
#'
#' @param species_list A named list of `data.frame`s. Each element is a
#'   species dataset with `lon`, `lat`, `presence`, and environmental
#'   variables. Names are used as species identifiers.
#' @param env_data Optional shared `data.frame` for spatial prediction.
#'   If `NULL`, spatial prediction is skipped.
#' @param models Character vector. Models to fit per species. Default
#'   `c("cast", "rf", "maxent", "brt")`.
#' @param output_dir Character. Directory for per-species outputs.
#'   Default `"castSDM_batch_output"`.
#' @param do_cv Logical. Run spatial CV per species. Default `TRUE`.
#' @param cv_k Integer. Number of spatial CV folds. Default `5`.
#' @param n_bootstrap Integer. DAG bootstrap replicates. Default `100`.
#' @param fig_dpi Integer. DPI for all saved figures. Default `300`.
#' @param parallel Logical. If `TRUE`, run species in parallel using
#'   \pkg{future}. Set a [future::plan()] before calling. Default `TRUE`.
#' @param seed Integer or `NULL`. Base seed. Each species gets
#'   `seed + species_index` for reproducibility.
#' @param verbose Logical. Default `TRUE`.
#' @param ... Additional arguments passed to [cast_fit()] (e.g.,
#'   `n_epochs`, `n_runs`, `tune_grid`).
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
#' @seealso [cast()], [plot.cast_batch()], [cast_consistency()]
#'
#' @export
cast_batch <- function(species_list,
                       env_data    = NULL,
                       models      = c("cast", "rf", "maxent", "brt"),
                       output_dir  = "castSDM_batch_output",
                       do_cv       = TRUE,
                       cv_k        = 5L,
                       n_bootstrap = 100L,
                       fig_dpi     = 300L,
                       parallel    = TRUE,
                       seed        = NULL,
                       verbose     = TRUE,
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

  run_one_species <- function(sp_name, sp_data, env_data, models, output_dir,
                              do_cv, cv_k, n_bootstrap, fig_dpi, seed_i,
                              fit_args) {
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
      # Step 1: Prepare
      split <- cast_prepare(sp_data, train_fraction = 0.7, seed = seed_i)

      # Step 2: DAG
      dag <- cast_dag(split$train, R = n_bootstrap, seed = seed_i,
                       verbose = FALSE)

      # Step 3: ATE
      ate <- cast_ate(split$train, K = 5L, seed = seed_i, verbose = FALSE)

      # Step 4: Screen
      screen <- cast_screen(dag, ate, split$train, seed = seed_i,
                             verbose = FALSE)

      # Step 5: Roles
      roles <- cast_roles(screen, dag)

      # Step 6: Fit
      fit_call_args <- c(
        list(data = split$train, screen = screen, dag = dag, ate = ate,
             models = models, seed = seed_i, verbose = FALSE),
        fit_args
      )
      fit <- do.call(cast_fit, fit_call_args)

      # Step 7: Evaluate (hold-out)
      eval_result <- cast_evaluate(fit, split$test)

      # Step 8: Spatial CV
      cv_result <- NULL
      if (do_cv) {
        cv_result <- tryCatch(
          cast_cv(sp_data, screen = screen, dag = dag, ate = ate,
                  k = cv_k, models = models, block_method = "grid",
                  seed = seed_i, verbose = FALSE),
          error = function(e) NULL
        )
      }

      # Step 9: Predict
      pred_result <- NULL
      if (!is.null(env_data)) {
        pred_result <- tryCatch(
          cast_predict(fit, env_data),
          error = function(e) NULL
        )
      }

      # Step 10: CATE
      cate_result <- tryCatch({
        check_suggested("grf", "for CATE")
        pred_for_cate <- if (!is.null(env_data)) env_data else NULL
        cast_cate(split$train, ate = ate, screen = screen,
                  predict_data = pred_for_cate, top_n = 3L,
                  seed = seed_i, verbose = FALSE)
      }, error = function(e) NULL)

      # Step 11: Consistency
      cons_result <- NULL
      if (!is.null(pred_result) && length(pred_result$models) >= 2) {
        cons_result <- tryCatch(cast_consistency(pred_result),
                                error = function(e) NULL)
      }

      # ── Save ALL plots ────────────────────────────────────────────────
      check_suggested("ggplot2", "for plotting")

      # DAG
      if (requireNamespace("ggraph", quietly = TRUE) &&
          requireNamespace("igraph", quietly = TRUE)) {
        p <- tryCatch(plot(dag, roles = roles, screen = screen,
                           species = sp_name), error = function(e) NULL)
        save_fig(p, "causal_dag.png", w = 14, h = 10)
      }

      # ATE
      p <- tryCatch(plot(ate), error = function(e) NULL)
      save_fig(p, "ate_forest_plot.png", w = 10, h = 7)

      # Screen
      p <- tryCatch(plot(screen), error = function(e) NULL)
      save_fig(p, "variable_screening.png", w = 10, h = 7)

      # Eval
      p <- tryCatch(plot(eval_result), error = function(e) NULL)
      save_fig(p, "model_evaluation.png", w = 10, h = 6)

      # CV
      if (!is.null(cv_result)) {
        p <- tryCatch(
          plot(cv_result, lon = sp_data$lon, lat = sp_data$lat,
               metric = "auc", basemap = "world"),
          error = function(e) NULL
        )
        save_fig(p, "spatial_cv_map.png", w = 14, h = 8)
      }

      # HSS maps (one per model)
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

      # CATE maps (one per variable)
      if (!is.null(cate_result) && requireNamespace("sf", quietly = TRUE)) {
        for (cv in cate_result$variables) {
          p <- tryCatch(
            plot(cate_result, variable = cv, species = sp_name,
                 basemap = "world", legend_position = "bottom"),
            error = function(e) NULL
          )
          save_fig(p, sprintf("CATE_%s.png", cv), w = 14, h = 8)
        }
      }

      # Consistency heatmap
      if (!is.null(cons_result) &&
          requireNamespace("patchwork", quietly = TRUE)) {
        p <- tryCatch(plot(cons_result, species = sp_name),
                       error = function(e) NULL)
        save_fig(p, "model_consistency.png", w = 16, h = 5.5)
      }

      # ── Save RDS results ──────────────────────────────────────────────
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
        env_data    = env_data,
        models      = models,
        output_dir  = output_dir,
        do_cv       = do_cv,
        cv_k        = cv_k,
        n_bootstrap = n_bootstrap,
        fig_dpi     = fig_dpi,
        fit_args    = fit_args
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
        do_cv, cv_k, n_bootstrap, fig_dpi, seed_i, fit_args
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
