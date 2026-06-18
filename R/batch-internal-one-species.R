# One species pipeline for cast_batch(); must live at package level so PSOCK workers
# can utils::getFromNamespace() after pkgload::load_all (future::multisession cannot).
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
    shared_suffix <- if (!is.null(cfg$shared_dag)) "_shared_dag" else ""

    split <- cast_run_step("prepare", output_dir, sp_name,
      cast_prepare(
        sp_data,
        train_fraction = cfg$train_fraction,
        seed = seed_i,
        env_vars = cfg$prepare_env_vars,
        verbose = cfg$prepare_verbose
      )
    )

    dag <- cfg$shared_dag
    if (is.null(dag)) {
      dag <- cast_run_step("dag", output_dir, sp_name,
        cast_dag(
          split$train,
          response = cfg$response,
          include_response = cfg$dag_include_response %||% TRUE,
          response_as_sink = cfg$dag_response_as_sink %||% TRUE,
          env_vars = cfg$dag_env_vars,
          R = cfg$dag_R,
          algorithm = cfg$dag_algorithm,
          score = cfg$dag_score,
          strength_threshold = cfg$dag_strength_threshold,
          direction_threshold = cfg$dag_direction_threshold,
          max_rows = cfg$dag_max_rows,
          seed = seed_i,
          verbose = cfg$dag_verbose,
          structure_method = cfg$dag_structure_method,
          pc_alpha = cfg$dag_pc_alpha,
          pc_test = cfg$dag_pc_test,
          mb_method = cfg$dag_mb_method,
          mb_alpha = cfg$dag_mb_alpha,
          bidag_algorithm = cfg$dag_bidag_algorithm,
          bidag_iterations = cfg$dag_bidag_iterations
        )
      )
    }

    select_suffix <- paste0("_", cfg$select_method %||% "invariant_screen", shared_suffix)
    screen <- cast_run_step(paste0("select", select_suffix), output_dir, sp_name,
      cast_select(
        dag, split$train,
        response = cfg$response,
        method = cfg$select_method %||% "invariant_screen",
        min_vars = cfg$select_min_vars %||% 5L,
        min_fraction = cfg$select_min_fraction %||% 0,
        num_trees = cfg$select_num_trees %||% 300L,
        stability_reps = cfg$select_stability_reps %||% 0L,
        stability_threshold = cfg$select_stability_threshold %||% 0.6,
        max_vars = cfg$select_max_vars,
        cor_threshold = cfg$select_cor_threshold %||% 0.8,
        seed = seed_i,
        verbose = cfg$select_verbose %||% FALSE
      )
    )

    refute_result <- NULL
    if (isTRUE(cfg$do_refute)) {
      refute_result <- cast_run_step(paste0("refute", select_suffix), output_dir, sp_name,
        tryCatch(
          cast_refute_screen(
            dag, screen, split$train,
            response = cfg$response,
            reps = cfg$refute_reps %||% 10L,
            num_trees = cfg$refute_num_trees %||% 80L,
            seed = seed_i,
            verbose = cfg$select_verbose %||% FALSE
          ),
          error = function(e) NULL
        )
      )
    }

    fit_call_args <- utils::modifyList(
      list(
        data = split$train,
        screen = screen,
        dag = dag,
        models = models,
        response = cfg$response,
        seed = seed_i,
        verbose = cfg$fit_verbose
      ),
      fit_args
    )
    fit <- cast_run_step(paste0("fit", select_suffix), output_dir, sp_name,
      do.call(cast_fit, fit_call_args)
    )

    eval_result <- cast_run_step(paste0("eval", select_suffix), output_dir, sp_name,
      cast_evaluate(fit, split$test, response = cfg$eval_response)
    )

    cv_result <- NULL
    if (cfg$do_cv) {
      cv_result <- cast_run_step(paste0("cv", select_suffix), output_dir, sp_name,
        tryCatch(
          cast_cv(
            sp_data,
            screen = screen, dag = dag,
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
        )
      )
    }

    pred_result <- NULL
    ensemble_result <- NULL
    if (isTRUE(cfg$do_predict) && !is.null(env_data)) {
      pred_result <- cast_run_step(paste0("predict", shared_suffix), output_dir, sp_name,
        tryCatch(
          cast_predict(fit, env_data, models = cfg$predict_models),
          error = function(e) NULL
        )
      )

      # Ensemble prediction
      if (isTRUE(cfg$do_ensemble) && !is.null(pred_result) && !is.null(cv_result)) {
        ensemble_result <- cast_run_step(paste0("ensemble", shared_suffix), output_dir, sp_name,
          tryCatch(
            cast_ensemble(fit, cv_result, env_data, method = cfg$ensemble_method %||% "weighted"),
            error = function(e) NULL
          )
        )
      }
    }

    # -- Raster Ensemble & Future Projection (optional) ----------
    raster_result <- NULL
    if (!is.null(cfg$raster_stack) && !is.null(cv_result) &&
        requireNamespace("terra", quietly = TRUE)) {
      raster_dir <- file.path(sp_dir, "rasters")
      raster_result <- tryCatch(
        cast_ensemble_raster(
          fit, cv_result, cfg$raster_stack,
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

      # Future projections
      if (!is.null(raster_result) && !is.null(cfg$future_rasters)) {
        tryCatch(
          cast_project_raster(
            fit, cv_result,
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

    # -- CATE (optional) ----------------------------------------------
    cate_result <- NULL
    if (isTRUE(cfg$do_cate)) {
      cate_result <- cast_run_step(paste0("cate", shared_suffix), output_dir, sp_name,
        tryCatch({
          check_suggested("grf", "for CATE")
          pred_for_cate <- if (!is.null(env_data)) env_data else NULL
          cast_cate(
            split$train,
            dag = dag, screen = screen,
            variables = cfg$cate_variables,
            response = cfg$response,
            predict_data = pred_for_cate,
            top_n = cfg$cate_top_n,
            n_trees = cfg$cate_n_trees,
            seed = seed_i,
            verbose = cfg$cate_verbose
          )
        }, error = function(e) NULL)
      )
    }

    # -- Save diagnostic plots ----------------------------------------
    check_suggested("ggplot2", "for plotting")

    vl <- cfg$var_labels
    bm <- cfg$plot_basemap

    # DAG structure plot
    if (requireNamespace("ggraph", quietly = TRUE) &&
        requireNamespace("igraph", quietly = TRUE)) {
      p <- tryCatch(
        plot(dag, screen = screen, species = sp_name, var_labels = vl),
        error = function(e) NULL
      )
      save_fig(p, "causal_dag.png", w = 14, h = 10)
    }

    # Variable selection plot
    p <- tryCatch(plot(screen, var_labels = vl), error = function(e) NULL)
    save_fig(p, "variable_selection.png", w = 10, h = 7)

    # Model evaluation plot
    p <- tryCatch(plot(eval_result), error = function(e) NULL)
    save_fig(p, "model_evaluation.png", w = 10, h = 6)

    # Spatial CV map
    if (!is.null(cv_result)) {
      p <- tryCatch(
        plot(
          cv_result, lon = sp_data$lon, lat = sp_data$lat,
          metric = "auc", basemap = bm
        ),
        error = function(e) NULL
      )
      save_fig(p, "spatial_cv_map.png", w = 14, h = 8)
    }

    # Ensemble HSS map (single final output, not per-algorithm)
    if (!is.null(ensemble_result) && requireNamespace("sf", quietly = TRUE)) {
      p <- tryCatch(
        plot(
          ensemble_result, basemap = bm,
          title = sprintf("%s Ensemble HSS", gsub("_", " ", sp_name))
        ),
        error = function(e) NULL
      )
      save_fig(p, "ensemble_HSS.png", w = 14, h = 8)
    } else if (!is.null(pred_result) && requireNamespace("sf", quietly = TRUE)) {
      # Fallback: best single model if ensemble not available
      best_model <- pred_result$models[1]
      p <- tryCatch(
        plot(
          pred_result, model = best_model, basemap = bm,
          title = sprintf("%s HSS (%s)", gsub("_", " ", sp_name), best_model)
        ),
        error = function(e) NULL
      )
      save_fig(p, sprintf("HSS_%s.png", best_model), w = 14, h = 8)
    }

    # CATE spatial maps (optional)
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
              basemap = bm,
              var_labels = vl,
              point_size = cfg$cate_point_size %||% 0.45,
              legend_position = "bottom",
              hss_predict = pred_result,
              hss_model = pm,
              hss_threshold = cfg$cate_hss_threshold %||% 0.1
            )
          } else {
            plot(
              cate_result,
              variable = cv,
              species = sp_name,
              basemap = bm,
              var_labels = vl,
              point_size = cfg$cate_point_size %||% 0.45,
              legend_position = "bottom"
            )
          }
        }, error = function(e) NULL)
        save_fig(p, sprintf("CATE_%s.png", cv), w = 14, h = 8)
      }
    }

    # -- Save RDS ------------------------------------------------------
    sp_result <- list(
      split = split,
      dag = dag, screen = screen,
      fit = fit, eval = eval_result, cv = cv_result,
      predict = pred_result, ensemble = ensemble_result,
      cate = cate_result,
      refute = refute_result
    )
    saveRDS(sp_result, file.path(sp_dir, "cast_result.rds"))

    sp_result
  }, error = function(e) {
    warning(sprintf("Species '%s' failed: %s", sp_name, e$message))
    NULL
  })

  result
}
