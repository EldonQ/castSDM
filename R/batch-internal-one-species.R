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

    screen <- cast_run_step(paste0("select", shared_suffix), output_dir, sp_name,
      cast_select(
        dag, split$train,
        response = cfg$response,
        min_vars = cfg$select_min_vars %||% 5L,
        min_fraction = cfg$select_min_fraction %||% 0.3,
        num_trees = cfg$select_num_trees %||% 300L,
        seed = seed_i,
        verbose = cfg$select_verbose %||% FALSE
      )
    )

    roles <- cast_roles(screen, dag)

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
    fit <- cast_run_step(paste0("fit", shared_suffix), output_dir, sp_name,
      do.call(cast_fit, fit_call_args)
    )

    eval_result <- cast_run_step(paste0("eval", shared_suffix), output_dir, sp_name,
      cast_evaluate(fit, split$test, response = cfg$eval_response)
    )

    cv_result <- NULL
    if (cfg$do_cv) {
      cv_result <- cast_run_step(paste0("cv", shared_suffix), output_dir, sp_name,
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
    if (!is.null(env_data)) {
      pred_result <- cast_run_step(paste0("predict", shared_suffix), output_dir, sp_name,
        tryCatch(
          cast_predict(fit, env_data, models = cfg$predict_models),
          error = function(e) NULL
        )
      )
    }

    cate_result <- NULL
    if (cfg$do_cate) {
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

    cons_result <- NULL
    if (!is.null(pred_result) && length(pred_result$models) >= 2) {
      cons_result <- cast_run_step(paste0("consistency", shared_suffix), output_dir, sp_name,
        tryCatch(cast_consistency(pred_result), error = function(e) NULL)
      )
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

    p <- tryCatch(plot(screen, var_labels = vl), error = function(e) NULL)
    save_fig(p, "variable_selection.png", w = 10, h = 7)

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
      dag = dag, screen = screen, roles = roles,
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
