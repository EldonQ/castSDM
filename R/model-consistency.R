#' Inter-model Spatial Consistency Metrics
#'
#' Computes pairwise spatial consistency metrics between model predictions,
#' following the validation framework of Wu et al. (2025). Three metrics
#' are computed for every pair of models:
#' \describe{
#'   \item{Cosine similarity}{Directional agreement between suitability
#'     vectors. Range \[0, 1\], 1 = identical direction.}
#'   \item{Warren's I}{Niche overlap metric (Warren et al. 2008). Measures
#'     the shared "probability mass" between two normalized suitability
#'     surfaces. Range \[0, 1\], 1 = perfect overlap.}
#'   \item{Pearson's r}{Linear correlation between suitability scores.
#'     Range \[-1, 1\].}
#' }
#'
#' @param x A [cast_predict] object containing HSS columns for multiple
#'   models.
#' @param models Character vector. Subset of models to compare. Default
#'   `NULL` uses all models present in `x`.
#'
#' @return A `cast_consistency` object with components:
#' \describe{
#'   \item{`metrics`}{`data.frame` with columns: `model_a`, `model_b`,
#'     `cosine_sim`, `warrens_I`, `pearson_r`, `pearson_p`.}
#'   \item{`models`}{Character vector of models compared.}
#' }
#'
#' @references
#' Warren, D.L., Glor, R.E. & Turelli, M. (2008). Environmental niche
#' equivalency versus conservatism. *Evolution*, 62, 2868-2883.
#'
#' @seealso [cast_predict()], [plot.cast_consistency()]
#'
#' @export
cast_consistency <- function(x, models = NULL) {
  if (!inherits(x, "cast_predict")) {
    cli::cli_abort("{.arg x} must be a {.cls cast_predict} object.")
  }

  pred <- x$predictions
  avail <- if (!is.null(models)) {
    models[paste0("HSS_", models) %in% names(pred)]
  } else {
    x$models
  }

  if (length(avail) < 2) {
    cli::cli_abort("At least 2 models required for consistency metrics.")
  }

  n <- length(avail)
  results <- list()

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      p1 <- as.numeric(pred[[paste0("HSS_", avail[i])]])
      p2 <- as.numeric(pred[[paste0("HSS_", avail[j])]])

      cos_val <- cosine_sim(p1, p2)
      wi_val  <- warrens_i(p1, p2)
      pr      <- pearson_r(p1, p2)

      results[[length(results) + 1]] <- data.frame(
        model_a    = avail[i],
        model_b    = avail[j],
        cosine_sim = cos_val,
        warrens_I  = wi_val,
        pearson_r  = pr$r,
        pearson_p  = pr$p,
        stringsAsFactors = FALSE
      )
    }
  }

  metrics_df <- do.call(rbind, results)
  new_cast_consistency(metrics = metrics_df, models = avail)
}


# ========================================================================
# Internal metric functions
# ========================================================================

#' Cosine Similarity
#' @keywords internal
#' @noRd
cosine_sim <- function(p1, p2) {
  valid <- !is.na(p1) & !is.na(p2)
  a <- p1[valid]
  b <- p2[valid]
  na <- sqrt(sum(a^2))
  nb <- sqrt(sum(b^2))
  if (na < 1e-10 || nb < 1e-10) return(NA_real_)
  sum(a * b) / (na * nb)
}

#' Warren's I Niche Overlap
#' @keywords internal
#' @noRd
warrens_i <- function(p1, p2) {
  valid <- !is.na(p1) & !is.na(p2) & (p1 >= 0) & (p2 >= 0)
  a <- p1[valid]
  b <- p2[valid]
  sa <- sum(a)
  sb <- sum(b)
  if (sa < 1e-10 || sb < 1e-10) return(NA_real_)
  sum(sqrt((a / sa) * (b / sb)))
}

#' Pearson Correlation
#' @keywords internal
#' @noRd
pearson_r <- function(p1, p2) {
  valid <- !is.na(p1) & !is.na(p2)
  a <- p1[valid]
  b <- p2[valid]
  if (length(a) < 10) return(list(r = NA_real_, p = NA_real_))
  ct <- stats::cor.test(a, b, method = "pearson")
  list(r = as.numeric(ct$estimate), p = as.numeric(ct$p.value))
}


# ========================================================================
# Plot method
# ========================================================================

#' Plot Inter-model Consistency Heatmaps
#'
#' Generates a 1x3 panel of heatmap matrices showing pairwise Cosine
#' similarity, Warren's I, and Pearson's r between all model pairs.
#' Follows the fig8 Part C visualization style with annotated cell values
#' and a YlGnBu color scale.
#'
#' @param x A `cast_consistency` object.
#' @param species Optional character. Species name for the title.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (or `patchwork` composite).
#' @export
plot.cast_consistency <- function(x, species = NULL, ...) {
  check_suggested("ggplot2", "for plotting")
  check_suggested("patchwork", "for multi-panel layout")

  models <- x$models
  metrics_df <- x$metrics
  n <- length(models)

  metric_defs <- list(
    list(key = "cosine_sim", label = "Cosine Similarity"),
    list(key = "warrens_I",  label = "Warren's I"),
    list(key = "pearson_r",  label = "Pearson's r")
  )

  panels <- list()

  for (mi in seq_along(metric_defs)) {
    mkey   <- metric_defs[[mi]]$key
    mlabel <- metric_defs[[mi]]$label

    # Build symmetric matrix
    mat_df <- expand.grid(
      model_x = factor(models, levels = models),
      model_y = factor(models, levels = rev(models)),
      stringsAsFactors = FALSE
    )
    mat_df$value <- NA_real_

    for (ri in seq_len(nrow(mat_df))) {
      mx <- as.character(mat_df$model_x[ri])
      my <- as.character(mat_df$model_y[ri])
      if (mx == my) {
        mat_df$value[ri] <- 1.0
      } else {
        row <- metrics_df[
          (metrics_df$model_a == mx & metrics_df$model_b == my) |
          (metrics_df$model_a == my & metrics_df$model_b == mx), ,
          drop = FALSE
        ]
        if (nrow(row) > 0) {
          val <- row[[mkey]][1]
          mat_df$value[ri] <- if (is.na(val)) 0 else val
        }
      }
    }

    mat_df$text_color <- ifelse(
      !is.na(mat_df$value) & mat_df$value > 0.7, "white", "black"
    )
    mat_df$label_text <- ifelse(
      is.na(mat_df$value), "", sprintf("%.3f", mat_df$value)
    )

    p <- ggplot2::ggplot(mat_df, ggplot2::aes(
      x = .data$model_x, y = .data$model_y, fill = .data$value
    )) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label_text, color = .data$text_color),
        size = 3.5, fontface = "bold", show.legend = FALSE
      ) +
      ggplot2::scale_fill_gradientn(
        colors = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1",
                   "#6BAED6", "#4292C6", "#2171B5", "#084594"),
        limits = c(0, 1), na.value = "grey90", name = NULL
      ) +
      ggplot2::scale_color_identity() +
      ggplot2::labs(title = mlabel, x = "", y = "") +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", hjust = 0.5, size = 11
        ),
        axis.text.x = ggplot2::element_text(
          angle = 45, hjust = 1, size = 8
        ),
        axis.text.y = ggplot2::element_text(size = 8),
        panel.grid = ggplot2::element_blank(),
        legend.position = if (mi == length(metric_defs)) "right" else "none"
      )

    panels[[mi]] <- p
  }

  sp_title <- if (!is.null(species)) {
    sprintf("Inter-model Spatial Consistency \u2014 %s", gsub("_", " ", species))
  } else {
    "Inter-model Spatial Consistency"
  }

  combined <- panels[[1]] + panels[[2]] + panels[[3]] +
    patchwork::plot_layout(nrow = 1) +
    patchwork::plot_annotation(
      title = sp_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = 14, hjust = 0.5
        )
      )
    )

  combined
}
