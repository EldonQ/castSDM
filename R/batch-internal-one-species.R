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
      ggplot2::ggsave(
        file.path(fig_dir, fname), p,
        width = w, height = h, dpi = fig_dpi,
        bg = "white", limitsize = FALSE
      ),
      error = function(e) NULL
    )
  }

  result <- tryCatch({
    split <- cast_prepare(
      sp_data,
      train_fraction = cfg$train_fraction,
      seed = seed_i,
      env_vars = cfg$prepare_env_vars,
      verbose = cfg$prepare_verbose
    )

    dag <- cast_dag(
      split$train,
      response = cfg$response,
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
      bidag_algorithm = cfg$dag_bidag_algorithm,
      bidag_iterations = cfg$dag_bidag_iterations,
      notears_lambda = cfg$dag_notears_lambda,
      notears_max_iter = cfg$dag_notears_max_iter,
      notears_lr = cfg$dag_notears_lr,
      notears_tol = cfg$dag_notears_tol,
      notears_rho_init = cfg$dag_notears_rho_init,
      notears_alpha_mult = cfg$dag_notears_alpha_mult
    )

    ate <- cast_ate(
      split$train,
      response = cfg$response,
      variables = cfg$ate_variables,
      K = cfg$ate_K,
      alpha = cfg$ate_alpha,
      num_trees = cfg$ate_num_trees,
      quantile_cuts = cfg$ate_quantile_cuts,
      bonferroni = cfg$ate_bonferroni,
      parallel = cfg$ate_parallel,
      seed = seed_i,
      verbose = cfg$ate_verbose
    )

    screen <- cast_screen(
      dag, ate, split$train,
      response = cfg$response,
      min_vars = cfg$screen_min_vars,
      min_fraction = cfg$screen_min_fraction,
      num_trees = cfg$screen_num_trees,
      seed = seed_i,
      verbose = cfg$screen_verbose
    )

    roles <- cast_roles(screen, dag)

    fit_call_args <- utils::modifyList(
      list(
        data = split$train,
        screen = screen,
        dag = dag,
        ate = ate,
        models = models,
        response = cfg$response,
        seed = seed_i,
        verbose = cfg$fit_verbose
      ),
      fit_args
    )
    fit <- do.call(cast_fit, fit_call_args)

    eval_result <- cast_evaluate(
      fit, split$test,
      response = cfg$eval_response
    )

    cv_result <- NULL
    if (cfg$do_cv) {
      cv_result <- tryCatch(
        cast_cv(
          sp_data,
          screen = screen, dag = dag, ate = ate,
          k = cfg$cv_k, models = cfg$cv_models,
          block_method = cfg$cv_block_method,
          response = cfg$response,
          n_epochs = cfg$cv_n_epochs,
          n_runs = cfg$cv_n_runs,
          rf_ntree = cfg$cv_rf_ntree,
          brt_n_trees = cfg$cv_brt_n_trees,
          parallel = cfg$cv_parallel,
          seed = seed_i,
          verbose = cfg$cv_verbose
        ),
        error = function(e) NULL
      )
    }

    pred_result <- NULL
    if (!is.null(env_data)) {
      pred_result <- tryCatch(
        cast_predict(fit, env_data, models = cfg$predict_models),
        error = function(e) NULL
      )
    }

    cate_result <- NULL
    if (cfg$do_cate) {
      cate_result <- tryCatch({
        check_suggested("grf", "for CATE")
        pred_for_cate <- if (!is.null(env_data)) env_data else NULL
        cast_cate(
          split$train,
          variables = cfg$cate_variables,
          ate = ate, screen = screen,
          response = cfg$response,
          predict_data = pred_for_cate,
          top_n = cfg$cate_top_n,
          n_trees = cfg$cate_n_trees,
          seed = seed_i,
          verbose = cfg$cate_verbose
        )
      }, error = function(e) NULL)
    }

    cons_result <- NULL
    if (!is.null(pred_result) && length(pred_result$models) >= 2) {
      cons_result <- tryCatch(cast_consistency(pred_result),
                              error = function(e) NULL)
    }

    # ── Save ALL plots ────────────────────────────────────────────────
    check_suggested("ggplot2", "for plotting")

    vl <- cfg$var_labels
    bm <- cfg$plot_basemap

    if (requireNamespace("ggraph", quietly = TRUE) &&
        requireNamespace("igraph", quietly = TRUE)) {
      p <- tryCatch(
        plot(
          dag, roles = roles, screen = screen,
          species = sp_name, var_labels = vl
        ),
        error = function(e) NULL
      )
      save_fig(p, "causal_dag.png", w = 14, h = 10)
    }

    p <- tryCatch(plot(ate, var_labels = vl), error = function(e) NULL)
    save_fig(p, "ate_forest_plot.png", w = 10, h = 7)

    p <- tryCatch(plot(screen, var_labels = vl), error = function(e) NULL)
    save_fig(p, "variable_screening.png", w = 10, h = 7)

    p <- tryCatch(plot(eval_result), error = function(e) NULL)
    save_fig(p, "model_evaluation.png", w = 10, h = 6)

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

    if (!is.null(pred_result) && requireNamespace("sf", quietly = TRUE)) {
      for (mdl in pred_result$models) {
        p <- tryCatch(
          plot(
            pred_result, model = mdl, basemap = bm,
            title = sprintf("%s HSS (%s)",
                            gsub("_", " ", sp_name), mdl)
          ),
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
              basemap = bm,
              var_labels = vl,
              point_size = cfg$cate_point_size,
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
              basemap = bm,
              var_labels = vl,
              point_size = cfg$cate_point_size,
              legend_position = "bottom"
            )
          }
        }, error = function(e) NULL)
        save_fig(p, sprintf("CATE_%s.png", cv), w = 14, h = 8)
      }
    }

    if (!is.null(cons_result) &&
        requireNamespace("patchwork", quietly = TRUE)) {
      p <- tryCatch(
        plot(
          cons_result,
          species = sp_name,
          font_family = "sans",
          use_bold = TRUE,
          font_base = 11L,
          font_main_title = 14L,
          font_panel_title = 11L,
          font_axis = 8L,
          font_cell_value = 3.5,
          value_decimals = 3L,
          tile_linewidth = 0.5,
          text_white_above = 0.7
        ),
        error = function(e) NULL
      )
      save_fig(p, "model_consistency.png", w = 16, h = 5.5)
    }

    # ── SHAP (optional; pairs + 2x3 panel) ────────────────────────────
    save_cast_batch_shap_outputs(
      train_df = split$train,
      fit = fit,
      dag = dag,
      screen = screen,
      cfg = cfg,
      fig_dir = fig_dir,
      fig_dpi = fig_dpi,
      seed_i = seed_i
    )

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


