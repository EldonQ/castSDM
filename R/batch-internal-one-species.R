# One species pipeline for cast_batch(); kept at package level so PSOCK
# workers can utils::getFromNamespace() after pkgload::load_all (future::
# multisession cannot).
.cast_batch_run_one_species <- function(sp_name, sp_data, env_data, models,
                            output_dir, fig_dpi, seed_i,
                            cfg, fit_args, parallel_batch, dev_root) {
  if (isTRUE(parallel_batch)) {
    root <- dev_root
    if (is.null(root) || !nzchar(as.character(root)[1])) {
      root <- Sys.getenv("CASTSDM_ROOT", "")
    }
    root <- tryCatch(
      normalizePath(as.character(root)[1], winslash = "/", mustWork = FALSE),
      error = function(e) ""
    )
    if (nzchar(root) && file.exists(file.path(root, "DESCRIPTION"))) {
      if (requireNamespace("pkgload", quietly = TRUE)) {
        suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE))
      }
    }
  }
  sp_dir  <- file.path(output_dir, sp_name)
  fig_dir <- file.path(sp_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  save_fig <- function(p, fname, w = 10, h = 7) {
    if (is.null(p)) return(invisible(NULL))
    tryCatch(
      cast_safe_ggsave(
        file.path(fig_dir, fname), p,
        width = w, height = h, dpi = fig_dpi,
        bg = "transparent", limitsize = FALSE
      ),
      error = function(e) NULL
    )
  }

  result <- tryCatch({
    # Run signature: every step cache is keyed on this, so changing any
    # input (data, config, seed, models) invalidates all stale caches.
    sig <- .cast_digest(list(
      species = sp_name, data = sp_data, models = models,
      seed = seed_i, cfg = cfg, fit_args = fit_args
    ))

    split <- cast_run_step("prepare", output_dir, sp_name,
      cast_prepare(
        sp_data,
        train_fraction = cfg$train_fraction,
        seed = seed_i,
        env_vars = cfg$prepare_env_vars,
        verbose = cfg$prepare_verbose
      ),
      params = list(signature = sig)
    )

    select_suffix <- paste0("_", cfg$select_method %||% "cpi")
    screen <- cast_run_step(paste0("select", select_suffix), output_dir, sp_name,
      cast_select(
        split$train,
        response = cfg$response,
        method = cfg$select_method %||% "cpi",
        alpha = cfg$select_alpha %||% 0.05,
        min_vars = cfg$select_min_vars %||% 3L,
        num_trees = cfg$select_num_trees %||% 300L,
        max_candidates = cfg$select_max_vars %||% 30L,
        dml_folds = cfg$select_dml_folds %||% 5L,
        cor_threshold = cfg$select_cor_threshold %||% 0.8,
        num_threads = cfg$num_threads %||% 1L,
        seed = seed_i,
        verbose = cfg$select_verbose %||% FALSE
      ),
      params = list(signature = sig)
    )

    fit_call_args <- utils::modifyList(
      list(
        data = split$train,
        screen = screen,
        models = models,
        response = cfg$response,
        num_threads = cfg$num_threads %||% 1L,
        seed = seed_i,
        verbose = cfg$fit_verbose
      ),
      fit_args
    )
    fit <- cast_run_step(paste0("fit", select_suffix), output_dir, sp_name,
      do.call(cast_fit, fit_call_args),
      params = list(signature = sig)
    )

    eval_result <- cast_run_step(paste0("eval", select_suffix), output_dir, sp_name,
      cast_evaluate(fit, split$test, response = cfg$eval_response),
      params = list(signature = sig)
    )

    cv_result <- NULL
    if (cfg$do_cv) {
      cv_result <- cast_run_step(paste0("cv", select_suffix), output_dir, sp_name,
        tryCatch(
          cast_cv(
            sp_data,
            screen = screen,
            select_method = cfg$select_method %||% "cpi",
            select_args = list(
              alpha = cfg$select_alpha %||% 0.05,
              min_vars = cfg$select_min_vars %||% 3L,
              max_candidates = cfg$select_max_vars %||% 30L,
              num_trees = cfg$select_num_trees %||% 300L,
              dml_folds = cfg$select_dml_folds %||% 5L,
              cor_threshold = cfg$select_cor_threshold %||% 0.8,
              num_threads = cfg$num_threads %||% 1L
            ),
            k = cfg$cv_k, models = cfg$cv_models,
            block_method = cfg$cv_block_method,
            response = cfg$response,
            rf_ntree = cfg$cv_rf_ntree,
            brt_n_trees = cfg$cv_brt_n_trees,
            parallel = cfg$cv_parallel,
            seed = seed_i,
            verbose = cfg$cv_verbose
          ),
          error = function(e) NULL
        ),
        params = list(signature = sig)
      )
    }

    # Refit on the full data set for final maps (standard SDM practice:
    # selection is decided on the training split, but the published
    # predictions reuse every record).
    fit_full <- NULL
    if (isTRUE(cfg$refit_full %||% TRUE)) {
      fit_full_args <- utils::modifyList(
        list(
          data = sp_data, screen = screen, models = models,
          response = cfg$response, num_threads = cfg$num_threads %||% 1L,
          seed = seed_i, verbose = cfg$fit_verbose
        ),
        fit_args
      )
      fit_full <- cast_run_step(paste0("fitfull", select_suffix),
        output_dir, sp_name,
        do.call(cast_fit, fit_full_args),
        params = list(signature = sig)
      )
    }
    fit_use <- if (!is.null(fit_full)) fit_full else fit

    pred_result <- NULL
    ensemble_result <- NULL
    if (isTRUE(cfg$do_predict) && !is.null(env_data)) {
      pred_result <- cast_run_step("predict", output_dir, sp_name,
        tryCatch(
          cast_predict(fit_use, env_data, models = cfg$predict_models),
          error = function(e) NULL
        ),
        params = list(signature = sig)
      )

      if (isTRUE(cfg$do_ensemble) && !is.null(pred_result) && !is.null(cv_result)) {
        ensemble_result <- cast_run_step("ensemble", output_dir, sp_name,
          tryCatch(
            cast_ensemble(fit_use, cv_result, env_data,
                          method = cfg$ensemble_method %||% "weighted"),
            error = function(e) NULL
          ),
          params = list(signature = sig)
        )
      } else if (isTRUE(cfg$do_ensemble) && is.null(cv_result)) {
        cli::cli_inform(
          "[{sp_name}] Ensemble skipped: spatial CV unavailable."
        )
      }
    }

    # -- Raster Ensemble & Future Projection (optional) --
    raster_result <- NULL
    if (!is.null(cfg$raster_stack) && !is.null(cv_result) &&
        requireNamespace("terra", quietly = TRUE)) {
      raster_dir <- file.path(sp_dir, "rasters")
      raster_result <- tryCatch(
        cast_ensemble_raster(
          fit_use, cv_result, cfg$raster_stack,
          output_dir = raster_dir,
          method = "weighted",
          models = cfg$predict_models,
          mask = cfg$raster_mask,
          prefix = "current",
          overwrite = cfg$overwrite_rasters %||% FALSE,
          compression = cfg$raster_compression %||% "LZW",
          verbose = FALSE
        ),
        error = function(e) {
          warning(sprintf("Raster ensemble for '%s' failed: %s", sp_name, e$message))
          NULL
        }
      )

      if (!is.null(raster_result) && !is.null(cfg$future_rasters)) {
        tryCatch(
          cast_project_raster(
            fit_use, cv_result,
            current_raster = cfg$raster_stack,
            future_rasters = cfg$future_rasters,
            output_dir = sp_dir,
            method = "weighted",
            models = cfg$predict_models,
            mask = cfg$raster_mask,
            overwrite = cfg$overwrite_rasters %||% FALSE,
            compression = cfg$raster_compression %||% "LZW",
            verbose = FALSE
          ),
          error = function(e) {
            warning(sprintf("Future projection for '%s' failed: %s", sp_name, e$message))
          }
        )
      }
    }

    # -- Save diagnostic plots --
    check_suggested("ggplot2", "for plotting")
    vl <- cfg$var_labels
    bm <- cfg$plot_basemap

    p <- tryCatch(plot(screen, var_labels = vl), error = function(e) NULL)
    save_fig(p, "variable_selection.png", w = 10, h = 7)

    p <- tryCatch(plot(eval_result), error = function(e) NULL)
    save_fig(p, "model_evaluation.png", w = 10, h = 6)

    if (!is.null(cv_result)) {
      p <- tryCatch(
        plot(cv_result, lon = sp_data$lon, lat = sp_data$lat,
             metric = "auc", basemap = bm),
        error = function(e) NULL
      )
      save_fig(p, "spatial_cv_map.png", w = 14, h = 8)
    }

    if (!is.null(ensemble_result) && requireNamespace("sf", quietly = TRUE)) {
      p <- tryCatch(
        plot(ensemble_result, basemap = bm,
             title = sprintf("%s Ensemble HSS", gsub("_", " ", sp_name))),
        error = function(e) NULL
      )
      save_fig(p, "ensemble_HSS.png", w = 14, h = 8)
    } else if (!is.null(pred_result) && requireNamespace("sf", quietly = TRUE)) {
      best_model <- pred_result$models[1]
      p <- tryCatch(
        plot(pred_result, model = best_model, basemap = bm,
             title = sprintf("%s HSS (%s)", gsub("_", " ", sp_name), best_model)),
        error = function(e) NULL
      )
      save_fig(p, sprintf("HSS_%s.png", best_model), w = 14, h = 8)
    }

    # -- Save RDS --
    sp_result <- new_cast_result(
      screen = screen,
      fit = fit,
      eval = eval_result,
      cv = cv_result,
      predict = pred_result,
      ensemble = ensemble_result
    )
    sp_result$fit_full <- fit_full
    sp_result$split <- split
    saveRDS(sp_result, file.path(sp_dir, "cast_result.rds"))

    sp_result
  }, error = function(e) {
    warning(sprintf("Species '%s' failed: %s", sp_name, e$message))
    NULL
  })

  result
}
