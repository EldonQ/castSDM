.cast_batch_shared_dag_data <- function(species_list,
                                        env_data,
                                        env_vars,
                                        response,
                                        max_rows,
                                        seed = NULL) {
  env_vars <- unique(env_vars)
  if (!is.null(env_data) && is.data.frame(env_data) &&
      nrow(env_data) > 0L && all(env_vars %in% names(env_data))) {
    out <- as.data.frame(env_data[, env_vars, drop = FALSE])
  } else {
    rows_per_species <- max(1L, ceiling(max_rows / max(1L, length(species_list))))
    chunks <- vector("list", length(species_list))
    if (!is.null(seed)) set.seed(seed)
    for (i in seq_along(species_list)) {
      d <- species_list[[i]]
      if (!is.data.frame(d) || !all(env_vars %in% names(d))) next
      take <- min(nrow(d), rows_per_species)
      idx <- if (nrow(d) > take) sample.int(nrow(d), take) else seq_len(nrow(d))
      chunks[[i]] <- as.data.frame(d[idx, env_vars, drop = FALSE])
    }
    chunks <- chunks[!vapply(chunks, is.null, logical(1L))]
    if (!length(chunks)) {
      cli::cli_abort(
        "Could not build shared DAG data: {.arg env_vars} are absent from both {.arg env_data} and {.arg species_list}."
      )
    }
    out <- do.call(rbind, chunks)
  }
  for (nm in names(out)) out[[nm]] <- suppressWarnings(as.numeric(out[[nm]]))
  out <- stats::na.omit(out)
  if (nrow(out) > max_rows) {
    if (!is.null(seed)) set.seed(seed + 7919L)
    out <- out[sample.int(nrow(out), max_rows), , drop = FALSE]
  }
  if (nrow(out) < 10L) {
    cli::cli_abort("Shared DAG data has fewer than 10 complete rows.")
  }
  out
}

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
#' @param fig_dpi Integer. Default `300`.
#' @param parallel Logical. Default `TRUE`.
#' @param seed Integer or `NULL`.
#' @param verbose Logical. Default `TRUE`.
#' @param fit_verbose Logical. Default `FALSE`.
#' @param dag_R Integer. Bootstrap replicates. Default `100`.
#' @param dag_structure_method Character. Default `"mb_first"`.
#' @param dag_include_response Logical. Default `TRUE`.
#' @param dag_pc_alpha Numeric. Default `0.05`.
#' @param dag_pc_test Character or `NULL`. Default `NULL`.
#' @param dag_mb_method Character. MB discovery algorithm. Default `"fast.iamb"`.
#' @param dag_mb_alpha Numeric. MB discovery significance. Default `0.05`.
#' @param dag_bidag_algorithm,dag_bidag_iterations BiDAG controls.
#' @param dag_algorithm Character. Default `"hc"`.
#' @param dag_score Character or `NULL`. Default `NULL`.
#' @param dag_strength_threshold Numeric. Default `0.7`.
#' @param dag_direction_threshold Numeric. Default `0.6`.
#' @param dag_max_rows Integer. Default `8000`.
#' @param select_min_vars Integer. Default `5`.
#' @param select_min_fraction Numeric. Default `0.3`.
#' @param select_num_trees Integer. Default `300`.
#' @param do_cv Logical. Default `TRUE`.
#' @param cv_k Integer. Default `5`.
#' @param cv_block_method Character. Default `"grid"`.
#' @param response Character. Default `"presence"`.
#' @param prepare_env_vars,prepare_verbose Passed to [cast_prepare()].
#' @param dag_env_vars,dag_verbose Passed to [cast_dag()].
#' @param select_verbose Passed to [cast_select()]. Default `FALSE`.
#' @param cv_models Character or `NULL`. Models for CV.
#' @param cv_rf_ntree,cv_brt_n_trees,cv_parallel,cv_verbose CV controls.
#' @param predict_models Passed to [cast_predict()].
#' @param plot_basemap Character. Default `"world"`.
#' @param eval_response Character or `NULL`.
#' @param var_labels Named character vector or `NULL`.
#' @param dev_package_root Character or `NULL`.
#' @param learn_shared_dag Logical. Default `FALSE`.
#' @param shared_dag Optional precomputed [cast_dag].
#' @param shared_dag_data Optional `data.frame`.
#' @param ... Additional arguments forwarded to [cast_fit()].
#'
#' @return A `cast_batch` object.
#' @seealso [cast()], [cast_fit()], [cast_dag()], [cast_select()],
#'   [cast_cv()]
#' @export
cast_batch <- function(species_list,
                       env_data    = NULL,
                       models      = c("rf", "brt", "maxent", "gam"),
                       train_fraction = 0.7,
                       output_dir  = "castSDM_batch_output",
                       fig_dpi     = 300L,
                       parallel    = TRUE,
                       seed        = NULL,
                       verbose     = TRUE,
                       fit_verbose = FALSE,
                       # ── DAG ──
                       dag_R                  = 100L,
                       dag_structure_method   = "mb_first",
                       dag_include_response   = TRUE,
                       dag_pc_alpha           = 0.05,
                       dag_pc_test            = NULL,
                       dag_mb_method          = "fast.iamb",
                       dag_mb_alpha           = 0.05,
                       dag_bidag_algorithm    = "order",
                       dag_bidag_iterations   = NULL,
                       dag_algorithm          = "hc",
                       dag_score              = NULL,
                       dag_strength_threshold = 0.7,
                       dag_direction_threshold = 0.6,
                       dag_max_rows           = 8000L,
                       # ── Selection ──
                       select_min_vars     = 5L,
                       select_min_fraction = 0.3,
                       select_num_trees    = 300L,
                       # ── CV ──
                       do_cv           = TRUE,
                       cv_k            = 5L,
                       cv_block_method = "grid",
                       response = "presence",
                       prepare_env_vars = NULL,
                       prepare_verbose = FALSE,
                       dag_env_vars = NULL,
                       dag_verbose = FALSE,
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
                       learn_shared_dag = FALSE,
                       shared_dag = NULL,
                       shared_dag_data = NULL,
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

  if (!is.null(shared_dag) && !inherits(shared_dag, "cast_dag")) {
    cli::cli_abort("{.arg shared_dag} must be a {.cls cast_dag} object or {.val NULL}.")
  }

  if (is.null(shared_dag) && isTRUE(learn_shared_dag)) {
    shared_dir <- file.path(output_dir, ".shared")
    dir.create(shared_dir, showWarnings = FALSE, recursive = TRUE)
    shared_path <- file.path(shared_dir, "shared_dag.rds")
    shared_dag <- tryCatch(
      if (file.exists(shared_path)) readRDS(shared_path) else NULL,
      error = function(e) NULL
    )
    if (!inherits(shared_dag, "cast_dag")) {
      dag_vars <- dag_env_vars %||% prepare_env_vars %||%
        get_env_vars(species_list[[sp_names[1L]]], response = response)
      dag_data <- shared_dag_data %||%
        .cast_batch_shared_dag_data(
          species_list = species_list,
          env_data = env_data,
          env_vars = dag_vars,
          response = response,
          max_rows = max(dag_max_rows * 4L, dag_max_rows),
          seed = seed
        )
      if (identical(dag_structure_method, "mb_first") &&
          isTRUE(dag_include_response) &&
          !response %in% names(dag_data)) {
        cli::cli_abort(c(
          "{.code learn_shared_dag = TRUE} is incompatible with response-focused {.code structure_method = \"mb_first\"} when shared data has no response column.",
          "i" = "For robust species-specific Markov Blanket selection, set {.code learn_shared_dag = FALSE}.",
          "i" = "Use shared DAGs only for predictor-only environmental structure, e.g. {.code dag_include_response = FALSE} with {.code structure_method = \"pc\"} or {.code \"bootstrap_hc\"}."
        ))
      }
      if (verbose) {
        cli::cli_inform(
          "Learning shared DAG once: {length(dag_vars)} env vars, {nrow(dag_data)} rows."
        )
      }
      shared_dag <- cast_run_step("shared_dag", output_dir, ".shared",
        cast_dag(
          dag_data,
          response = response,
          include_response = dag_include_response,
          env_vars = dag_vars,
          R = dag_R,
          algorithm = dag_algorithm,
          score = dag_score,
          strength_threshold = dag_strength_threshold,
          direction_threshold = dag_direction_threshold,
          max_rows = dag_max_rows,
          seed = seed,
          verbose = dag_verbose,
          structure_method = dag_structure_method,
          pc_alpha = dag_pc_alpha,
          pc_test = dag_pc_test,
          mb_method = dag_mb_method,
          mb_alpha = dag_mb_alpha,
          bidag_algorithm = dag_bidag_algorithm,
          bidag_iterations = dag_bidag_iterations
        )
      )
      saveRDS(shared_dag, shared_path)
    } else if (verbose) {
      cli::cli_inform("Shared DAG cache hit: {.path {shared_path}}")
    }
  }

  # Pack all pipeline config
  cfg <- list(
    response = response,
    eval_response = eval_resp,
    prepare_env_vars = prepare_env_vars,
    prepare_verbose = prepare_verbose,
    train_fraction = train_fraction,
    dag_include_response = dag_include_response,
    dag_env_vars = dag_env_vars,
    dag_verbose = dag_verbose,
    dag_R = dag_R,
    dag_structure_method = dag_structure_method,
    dag_pc_alpha = dag_pc_alpha,
    dag_pc_test = dag_pc_test,
    dag_mb_method = dag_mb_method,
    dag_mb_alpha = dag_mb_alpha,
    dag_bidag_algorithm = dag_bidag_algorithm,
    dag_bidag_iterations = dag_bidag_iterations,
    dag_algorithm = dag_algorithm,
    dag_score = dag_score,
    dag_strength_threshold = dag_strength_threshold,
    dag_direction_threshold = dag_direction_threshold,
    dag_max_rows = dag_max_rows,
    select_min_vars = select_min_vars,
    select_min_fraction = select_min_fraction,
    select_num_trees = select_num_trees,
    select_verbose = select_verbose,
    do_cv = do_cv, cv_k = cv_k, cv_block_method = cv_block_method,
    cv_models = cv_models_use,
    cv_rf_ntree = cv_rf_ntree,
    cv_brt_n_trees = cv_brt_n_trees,
    cv_parallel = cv_parallel,
    cv_verbose = cv_verbose,
    predict_models = predict_models,
    plot_basemap = plot_basemap,
    var_labels = var_labels,
    fit_verbose = fit_verbose,
    shared_dag = shared_dag
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
      results[i] <- list(.cast_batch_run_one_species(
        sp, species_list[[sp]], env_data, models, output_dir,
        fig_dpi, seed_i, cfg, fit_args, FALSE, dev_root_workers
      ))
    }
  }

  names(results) <- sp_names

  # Collect metrics
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
        c("fold", "model", "auc", "tss", "cbi", "species"),
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
    ggplot2::theme_minimal(
      base_size = 11,
      base_family = getOption("castSDM.font_family", "Arial")
    ) +
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
