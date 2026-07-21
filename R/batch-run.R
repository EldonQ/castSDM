#' Batch Multi-Species Modeling
#'
#' One-stop interface that runs the full castSDM pipeline on multiple
#' species, optionally in parallel. Results (RDS + diagnostic figures) are
#' saved to separate subdirectories per species.
#'
#' @param species_list A named list of `data.frame`s.
#' @param env_data Optional shared `data.frame` for spatial prediction.
#' @param models Character vector. Default `c("rf", "brt", "maxent", "gam")`.
#' @param train_fraction Numeric. Default `0.7`.
#' @param output_dir Character. Default `"castSDM_batch_output"`.
#' @param fig_dpi Integer. Default `600`.
#' @param parallel Logical. Default `TRUE`.
#' @param seed Integer or `NULL`.
#' @param verbose Logical. Default `TRUE`.
#' @param fit_verbose Logical. Default `FALSE`.
#' @param select_min_vars Integer. Default `3`.
#' @param select_num_trees Integer. Default `300`.
#' @param select_method Character. Variable screening method passed to
#'   [cast_select()]. Default `"dml"`.
#' @param select_max_vars Candidate ceiling for DML testing / RF output.
#'   Default `30`.
#' @param select_alpha Numeric. FDR level for the DML selector. Default `0.05`.
#' @param select_cor_threshold Numeric. Redundant-proxy correlation threshold.
#'   Default `0.8`.
#' @param do_cv Logical. Default `TRUE`.
#' @param cv_k Integer. Default `5`.
#' @param cv_block_method Character. Default `"grid"`.
#' @param do_predict Logical. Default `TRUE`.
#' @param do_ensemble Logical. Default `TRUE`.
#' @param ensemble_method Character. `"weighted"`, `"best"`, or `"equal"`.
#'   Default `"weighted"`.
#' @param response Character. Default `"presence"`.
#' @param prepare_env_vars,prepare_verbose Passed to [cast_prepare()].
#' @param select_verbose Passed to [cast_select()]. Default `FALSE`.
#' @param cv_models Character or `NULL`.
#' @param cv_rf_ntree,cv_brt_n_trees,cv_parallel,cv_verbose CV controls.
#' @param predict_models Passed to [cast_predict()].
#' @param plot_basemap Character. Default `"world"`.
#' @param eval_response Character or `NULL`.
#' @param var_labels Named character vector or `NULL`.
#' @param dev_package_root Character or `NULL`.
#' @param raster_stack Optional `terra::SpatRaster` for raster-based
#'   ensemble prediction via [cast_ensemble_raster()].
#' @param future_rasters Optional named list of `terra::SpatRaster` stacks
#'   for future projection via [cast_project_raster()].
#' @param raster_mask Optional `terra::SpatRaster` prediction mask.
#' @param raster_compression Character. Default `"LZW"`.
#' @param overwrite_rasters Logical. Default `FALSE`.
#' @param ... Additional arguments forwarded to [cast_fit()].
#'
#' @return A `cast_batch` object.
#' @seealso [cast()], [cast_fit()], [cast_select()], [cast_cv()]
#' @export
cast_batch <- function(species_list,
                      env_data    = NULL,
                      models      = c("rf", "brt", "maxent", "gam"),
                      train_fraction = 0.7,
                      output_dir  = "castSDM_batch_output",
                      fig_dpi     = 600L,
                      parallel    = TRUE,
                      seed        = NULL,
                      verbose     = TRUE,
                      fit_verbose = FALSE,
                      # -- Selection --
                      select_min_vars     = 3L,
                      select_num_trees    = 300L,
                      select_method = "dml",
                      select_max_vars = 30L,
                      select_alpha = 0.05,
                      select_cor_threshold = 0.8,
                      # -- CV --
                      do_cv           = TRUE,
                      cv_k            = 5L,
                      cv_block_method = "grid",
                      do_predict = TRUE,
                      do_ensemble = TRUE,
                      ensemble_method = "weighted",
                      response = "presence",
                      prepare_env_vars = NULL,
                      prepare_verbose = FALSE,
                      select_verbose = FALSE,
                      cv_models = NULL,
                      cv_rf_ntree = 300L,
                      cv_brt_n_trees = 500L,
                      cv_parallel = FALSE,
                      cv_verbose = FALSE,
                      predict_models = NULL,
                      plot_basemap = "world",
                      eval_response = NULL,
                      var_labels = NULL,
                      dev_package_root = NULL,
                      raster_stack = NULL,
                      future_rasters = NULL,
                      raster_mask = NULL,
                      raster_compression = "LZW",
                      overwrite_rasters = FALSE,
                      ...) {

  if (!is.list(species_list) || is.null(names(species_list))) {
    cli::cli_abort("{.arg species_list} must be a named list of data.frames.")
  }
  sp_names <- names(species_list)
  n_sp <- length(sp_names)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (verbose) {
    cli::cli_h1("castSDM Batch: {n_sp} species")
    cli::cli_inform("Models: {.val {models}}")
    cli::cli_inform("Output: {output_dir}")
  }

  dev_root_workers <- NULL
  if (!is.null(dev_package_root) && nzchar(as.character(dev_package_root)[1])) {
    dev_root_workers <- tryCatch(
      normalizePath(as.character(dev_package_root)[1],
                    winslash = "/", mustWork = FALSE),
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
        nl <- tryCatch(normalizePath(lib, winslash = "/", mustWork = FALSE),
                       error = function(e) "")
        nzchar(nl) && (identical(np, nl) || startsWith(np, paste0(nl, "/")))
      }, logical(1L)))
      if (!isTRUE(in_lib)) Sys.setenv(CASTSDM_ROOT = root)
    }
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
        "Dropping {.arg ...} names not in {.fn cast_fit}: {.val {names(fit_args)[bad]}}."
      )
      fit_args <- fit_args[!bad]
    }
  }

  eval_resp <- if (is.null(eval_response)) response else eval_response
  cv_models_use <- if (is.null(cv_models)) models else cv_models

  cfg <- list(
    response = response,
    eval_response = eval_resp,
    prepare_env_vars = prepare_env_vars,
    prepare_verbose = prepare_verbose,
    train_fraction = train_fraction,
    select_min_vars = select_min_vars,
    select_num_trees = select_num_trees,
    select_method = select_method,
    select_max_vars = select_max_vars,
    select_alpha = select_alpha,
    select_cor_threshold = select_cor_threshold,
    select_verbose = select_verbose,
    do_cv = do_cv, cv_k = cv_k, cv_block_method = cv_block_method,
    cv_models = cv_models_use,
    cv_rf_ntree = cv_rf_ntree,
    cv_brt_n_trees = cv_brt_n_trees,
    cv_parallel = cv_parallel,
    cv_verbose = cv_verbose,
    predict_models = predict_models,
    do_predict = do_predict,
    do_ensemble = do_ensemble,
    ensemble_method = ensemble_method,
    plot_basemap = plot_basemap,
    var_labels = var_labels,
    fit_verbose = fit_verbose,
    raster_stack = raster_stack,
    future_rasters = future_rasters,
    raster_mask = raster_mask,
    raster_compression = raster_compression,
    overwrite_rasters = overwrite_rasters
  )

  if (parallel && !is.null(dev_root_workers)) {
    nwrk <- tryCatch(
      if (requireNamespace("future", quietly = TRUE)) future::nbrOfWorkers()
      else NA_integer_,
      error = function(e) NA_integer_
    )
    if (is.na(nwrk) || nwrk < 1L) nwrk <- max(1L, parallel::detectCores() - 1L)
    nwrk <- min(as.integer(nwrk), n_sp)
    if (verbose) {
      cli::cli_inform(
        "Parallel (dev): PSOCK with {nwrk} workers; pkgload::load_all() on each."
      )
    }
    root_value <- as.character(dev_root_workers)[1L]
    cl <- parallel::makeCluster(nwrk)
    results <- tryCatch(
      {
        parallel::clusterExport(cl, "root_value", envir = environment())
        ok <- parallel::clusterEvalQ(cl, {
          if (!requireNamespace("pkgload", quietly = TRUE)) FALSE
          else {
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
          c("species_list", "sp_names", "env_data", "models", "output_dir",
            "fig_dpi", "seed", "cfg", "fit_args"),
          envir = eb
        )
        parallel::parLapply(cl, seq_along(sp_names), function(ii) {
          sp <- sp_names[[ii]]
          sd <- species_list[[sp]]
          seed_i <- if (!is.null(seed)) seed + ii else NULL
          worker_run <- utils::getFromNamespace(".cast_batch_run_one_species", "castSDM")
          worker_run(sp, sd, env_data, models,
                     output_dir, fig_dpi, seed_i,
                     cfg, fit_args, FALSE, NULL)
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
        .cast_batch_run_one_species(sp, sd, env_data, models,
                                   output_dir, fig_dpi, seed_i,
                                   cfg, fit_args, TRUE, dev_root_workers)
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
      results[i] <- list(.cast_batch_run_one_species(
        sp, species_list[[sp]], env_data, models, output_dir,
        fig_dpi, seed_i, cfg, fit_args, FALSE, dev_root_workers
      ))
    }
  }

  names(results) <- sp_names

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
      for (mcol in c("auc", "tss", "cbi")) {
        mean_col <- paste0(mcol, "_mean")
        if (mean_col %in% names(em)) em[[mcol]] <- em[[mean_col]]
      }
      em$fold <- 0L
      metrics_rows[[sp]] <- em[, intersect(
        c("fold", "model", "auc", "tss", "cbi", "species"), names(em)
      ), drop = FALSE]
    }
  }

  species_metrics <- if (length(metrics_rows) > 0) do.call(rbind, metrics_rows)
                     else data.frame()
  rownames(species_metrics) <- NULL

  n_ok <- sum(!vapply(results, is.null, logical(1)))
  if (verbose) cli::cli_inform("Batch complete: {n_ok}/{n_sp} species succeeded.")

  new_cast_batch(
    species_metrics = species_metrics,
    species         = sp_names,
    models          = models,
    results         = results,
    output_dir      = output_dir
  )
}


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
      species = sm$species, model = sm$model,
      metric = toupper(mc), value = sm[[mc]],
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
      ggplot2::aes(group = .data$species),
      width = 0.15, size = 1.5, alpha = 0.55, color = "black", shape = 16
    ) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_manual(values = gray_fills, guide = "none") +
    ggplot2::labs(
      title    = "Multi-species Model Performance Comparison",
      subtitle = sprintf("%d species | Spatial CV", length(x$species)),
      x = "", y = "Score"
    ) +
    ggplot2::theme_minimal(
      base_size = 11,
      base_family = getOption("castSDM.font_family", "Arial")
    ) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border       = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5),
      strip.text         = ggplot2::element_text(face = "bold", size = 10),
      axis.title         = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle      = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position    = "bottom",
      legend.text        = ggplot2::element_text(size = 8)
    )
}
