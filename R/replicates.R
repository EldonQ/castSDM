# Replicate design ----------------------------------------------------------
#
# cast_rep() repeats background (pseudo-absence) sampling n_reps times and
# runs the full castSDM pipeline on every replicate. The returned object
# aggregates evaluation metrics (mean +/- SD across replicates), the
# selection frequency of each predictor, and the per-cell prediction mean
# and standard deviation, giving the replicate-level variance that a single
# background draw cannot provide (cf. N-SDM's n_reps pseudoabsence
# replicates; Barbet-Massin et al. 2012).

#' Replicated Background-Sampling Pipeline
#'
#' Runs [cast_background()] + [cast()] for `n_reps` independent background
#' draws and aggregates the results. Replicates quantify the sensitivity of
#' evaluation metrics, variable selection, and predicted suitability to the
#' arbitrary choice of pseudo-absence locations - the dominant source of
#' run-to-run variance in presence-background SDMs.
#'
#' @param occurrences A `data.frame` with `lon`, `lat` columns (presence
#'   locations). May carry extra columns, which are preserved.
#' @param raster_stack A `terra::SpatRaster` of environmental variables
#'   used for background sampling and value extraction.
#' @param env_data Optional prediction grid (`lon`, `lat` + predictors). When
#'   supplied, per-replicate predictions are stored and aggregated into
#'   `prediction` (mean + SD per cell).
#' @param n_reps Number of background replicates. Default `10`.
#' @param models Models passed to [cast()]. Default `c("rf")`.
#' @param select_method Selection method passed to [cast()]. Default `"cpi"`.
#' @param seed Integer or `NULL`. Replicate `r` uses `seed + r`.
#' @param verbose Logical. Default `TRUE`.
#' @param ... Further arguments passed to [cast()] (e.g. `do_cv`, `cv_k`,
#'   `select_min_vars`), and `background` arguments passed to
#'   [cast_background()] (`ratio`, `min_bg`, `max_bg`, `strategy`,
#'   `cell_thin`, `exclude_presence`).
#'
#' @return A `cast_rep` object with components:
#' \describe{
#'   \item{reps}{List of per-replicate [cast()] results.}
#'   \item{metrics}{Aggregated evaluation metrics (mean and SD across
#'     replicates) per model.}
#'   \item{selection_freq}{Data.frame with each predictor's selection
#'     frequency across replicates.}
#'   \item{prediction}{`data.frame` with `lon`, `lat`, `hss_mean`,
#'     `hss_sd` aggregated from the replicate predictions (only when
#'     `env_data` is supplied).}
#'   \item{n_reps}{Number of completed replicates.}
#' }
#'
#' @references
#' Barbet-Massin, M. et al. (2012). Selecting pseudo-absences for species
#' distribution models: how, where and how many? *Methods in Ecology and
#' Evolution*, 3(2), 327-338.
#'
#' @seealso [cast_background()], [cast()], [plot.cast_rep()]
#' @export
cast_rep <- function(occurrences, raster_stack, env_data = NULL,
                     n_reps = 10, models = c("rf"),
                     select_method = "cpi",
                     seed = NULL, verbose = TRUE, ...) {
  check_suggested("terra", "for background raster sampling")
  dots <- list(...)

  bg_args <- dots[intersect(
    names(dots),
    c("n_bg", "ratio", "min_bg", "max_bg", "strategy", "cell_thin",
      "exclude_presence")
  )]
  cast_args <- dots[setdiff(names(dots), names(bg_args))]

  n_reps <- max(1L, as.integer(n_reps))
  reps <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    seed_r <- if (is.null(seed)) NULL else seed + r
    if (verbose) cli::cli_inform("Replicate {r}/{n_reps}...")

    bg <- do.call(cast_background, c(
      list(occurrences = occurrences, raster_stack = raster_stack,
           seed = seed_r, verbose = FALSE),
      bg_args
    ))
    reps[[r]] <- tryCatch(
      do.call(cast, c(
        list(species_data = bg, env_data = env_data,
             models = models, select_method = select_method,
             seed = seed_r, verbose = FALSE),
        cast_args
      )),
      error = function(e) {
        if (verbose) cli::cli_warn("Replicate {r} failed: {e$message}")
        NULL
      }
    )
  }
  ok <- !vapply(reps, is.null, logical(1))
  reps <- reps[ok]
  n_done <- length(reps)
  if (!n_done) cli::cli_abort("All replicates failed.")

  # -- Aggregated evaluation metrics (hold-out, per model) ----------------
  metric_rows <- lapply(seq_along(reps), function(i) {
    m <- reps[[i]]$eval$metrics
    m$rep <- i
    m
  })
  metrics_long <- do.call(rbind, metric_rows)
  agg <- lapply(unique(metrics_long$model), function(md) {
    z <- metrics_long[metrics_long$model == md, ]
    data.frame(
      model = md,
      auc_mean = mean(z$auc_mean, na.rm = TRUE),
      auc_sd = stats::sd(z$auc_mean, na.rm = TRUE),
      tss_mean = mean(z$tss_mean, na.rm = TRUE),
      tss_sd = stats::sd(z$tss_mean, na.rm = TRUE),
      cbi_mean = mean(z$cbi_mean, na.rm = TRUE),
      cbi_sd = stats::sd(z$cbi_mean, na.rm = TRUE),
      n_reps = nrow(z),
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, agg)
  rownames(metrics) <- NULL

  # -- Selection frequency ------------------------------------------------
  all_vars <- unique(unlist(lapply(reps, function(r) r$screen$scores$variable)))
  freq <- vapply(all_vars, function(v) {
    mean(vapply(reps, function(r) v %in% r$screen$selected, logical(1)))
  }, numeric(1))
  selection_freq <- data.frame(
    variable = all_vars,
    freq = unname(freq),
    stringsAsFactors = FALSE
  )
  selection_freq <- selection_freq[order(-selection_freq$freq), ]
  rownames(selection_freq) <- NULL

  # -- Aggregated prediction surface --------------------------------------
  prediction <- NULL
  if (!is.null(env_data) && all(vapply(reps, function(r) !is.null(r$predict),
                                       logical(1)))) {
    cols <- grep("^HSS_", names(reps[[1]]$predict$predictions), value = TRUE)
    if (length(cols)) {
      mat <- do.call(cbind, lapply(reps, function(r) {
        r$predict$predictions[[cols[1]]]
      }))
      prediction <- data.frame(
        lon = reps[[1]]$predict$predictions$lon,
        lat = reps[[1]]$predict$predictions$lat,
        hss_mean = rowMeans(mat, na.rm = TRUE),
        hss_sd = if (ncol(mat) > 1L) apply(mat, 1L, stats::sd, na.rm = TRUE)
                 else rep(NA_real_, nrow(mat)),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- list(
    reps = reps,
    metrics = metrics,
    selection_freq = selection_freq,
    prediction = prediction,
    n_reps = n_done
  )
  class(out) <- "cast_rep"
  out
}

#' @export
print.cast_rep <- function(x, ...) {
  cli::cli_h1("castSDM Replicate Design")
  cli::cli_ul(c(
    "Completed replicates: {x$n_reps}",
    "Models: {.val {x$metrics$model}}"
  ))
  cli::cli_h2("Evaluation metrics (mean +/- SD across replicates)")
  m <- x$metrics
  for (i in seq_len(nrow(m))) {
    cli::cli_li(paste0(
      m$model[i],
      " | AUC=", round(m$auc_mean[i], 3), " (", round(m$auc_sd[i], 3), ")",
      " TSS=", round(m$tss_mean[i], 3), " (", round(m$tss_sd[i], 3), ")",
      " CBI=", round(m$cbi_mean[i], 3)
    ))
  }
  invisible(x)
}

#' Plot Replicate Results
#'
#' Two-panel figure: (left) per-model evaluation metrics with replicate
#' spread; (right) prediction mean and uncertainty (SD) when a prediction
#' grid was supplied.
#'
#' @param x A `cast_rep` object.
#' @param ... Ignored.
#' @return A `ggplot` object (or `patchwork` combination).
#' @export
plot.cast_rep <- function(x, ...) {
  check_suggested("ggplot2", "for plotting")

  m <- x$metrics
  long <- do.call(rbind, lapply(seq_len(nrow(m)), function(i) {
    data.frame(
      model = m$model[i],
      metric = c("AUC", "TSS", "CBI"),
      mean = c(m$auc_mean[i], m$tss_mean[i], m$cbi_mean[i]),
      sd = c(m$auc_sd[i], m$tss_sd[i], m$cbi_sd[i])
    )
  }))

  p1 <- ggplot2::ggplot(long, ggplot2::aes(x = .data$metric, y = .data$mean,
                                           fill = .data$model)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.7), width = 0.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$mean - .data$sd, ymax = .data$mean + .data$sd),
      position = ggplot2::position_dodge(0.7), width = 0.25
    ) +
    ggplot2::labs(title = "Evaluation across replicates",
                  subtitle = sprintf("%d replicates", x$n_reps),
                  x = "", y = "Score") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_cast()

  if (is.null(x$prediction)) return(p1)

  p2 <- ggplot2::ggplot(x$prediction,
                        ggplot2::aes(x = .data$hss_mean, y = .data$hss_sd)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6, color = "#256E92") +
    ggplot2::labs(title = "Prediction uncertainty across replicates",
                  x = "Mean HSS", y = "HSS SD") +
    theme_cast()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p1 + p2 + patchwork::plot_layout(widths = c(1, 1))
  } else {
    p1
  }
}
