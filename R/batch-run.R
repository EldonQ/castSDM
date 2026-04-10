#' Batch Multi-Species Modeling
#'
#' Runs the full castSDM pipeline on multiple species, optionally in
#' parallel. Each species is processed independently and results are saved
#' to separate subdirectories. A comparative summary of model performance
#' across species is returned.
#'
#' @param species_list A named list of `data.frame`s. Each element is a
#'   species dataset with `lon`, `lat`, `presence`, and environmental
#'   variables. Names are used as species identifiers.
#' @param env_data Optional shared `data.frame` for spatial prediction.
#'   If `NULL`, spatial prediction is skipped.
#' @param models Character vector. Models to fit per species. Default
#'   `c("cast", "rf", "maxent", "brt")`.
#' @param output_dir Character. Directory for per-species RDS outputs.
#'   Default `"castSDM_batch_output"`.
#' @param do_cv Logical. Run spatial CV per species. Default `TRUE`.
#' @param cv_k Integer. Number of spatial CV folds. Default `5`.
#' @param parallel Logical. If `TRUE`, run species in parallel using
#'   \pkg{future}. Set a [future::plan()] before calling. Default `TRUE`.
#' @param seed Integer or `NULL`. Base seed. Each species gets
#'   `seed + species_index` for reproducibility.
#' @param verbose Logical. Default `TRUE`.
#' @param ... Additional arguments passed to [cast()].
#'
#' @return A `cast_batch` object with components:
#' \describe{
#'   \item{`species_metrics`}{`data.frame` with per-species per-model
#'     evaluation metrics (from spatial CV when available).}
#'   \item{`species`}{Character vector of species names.}
#'   \item{`models`}{Character vector of model names.}
#'   \item{`results`}{Named list of `cast_result` objects (see [cast()]).}
#'   \item{`output_dir`}{Output directory path.}
#' }
#'
#' @seealso [cast()], [plot.cast_batch()], [cast_consistency()]
#'
#' @export
cast_batch <- function(species_list,
                       env_data   = NULL,
                       models     = c("cast", "rf", "maxent", "brt"),
                       output_dir = "castSDM_batch_output",
                       do_cv      = TRUE,
                       cv_k       = 5L,
                       parallel   = TRUE,
                       seed       = NULL,
                       verbose    = TRUE,
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

  # Worker: process one species
  run_one_species <- function(sp_name, sp_data, env_data, models, output_dir,
                              do_cv, cv_k, seed_i, ...) {
    sp_dir <- file.path(output_dir, sp_name)
    dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)

    result <- tryCatch(
      cast(
        species_data = sp_data,
        env_data     = env_data,
        models       = models,
        do_cv        = do_cv,
        cv_k         = cv_k,
        seed         = seed_i,
        verbose      = FALSE,
        ...
      ),
      error = function(e) {
        warning(sprintf("Species '%s' failed: %s", sp_name, e$message))
        NULL
      }
    )

    if (!is.null(result)) {
      saveRDS(result, file.path(sp_dir, "cast_result.rds"))
    }

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
        do_cv      = do_cv,
        cv_k       = cv_k,
        ...
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
        do_cv, cv_k, seed_i, ...
      )
    }
  }

  names(results) <- sp_names

  # Collect metrics across species
  metrics_rows <- list()
  for (sp in sp_names) {
    r <- results[[sp]]
    if (is.null(r)) next

    # Prefer spatial CV metrics; fall back to hold-out eval
    if (!is.null(r$cv) && !is.null(r$cv$fold_metrics) &&
        nrow(r$cv$fold_metrics) > 0) {
      fm <- r$cv$fold_metrics
      fm$species <- sp
      metrics_rows[[sp]] <- fm
    } else if (!is.null(r$eval) && !is.null(r$eval$metrics)) {
      em <- r$eval$metrics
      em$species <- sp
      # Rename _mean columns for consistency
      for (mcol in c("auc", "tss", "cbi", "sedi", "kappa", "prauc")) {
        mean_col <- paste0(mcol, "_mean")
        if (mean_col %in% names(em)) {
          em[[mcol]] <- em[[mean_col]]
        }
      }
      em$fold <- 0L
      metrics_rows[[sp]] <- em[, intersect(
        c("fold", "model", "auc", "tss", "cbi", "sedi", "kappa", "prauc", "species"),
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

  if (verbose) {
    cli::cli_inform("Batch complete: {sum(!vapply(results, is.null, logical(1)))}/{n_sp} species succeeded.")
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

  # Pivot to long format
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
