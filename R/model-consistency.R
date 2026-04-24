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
#' Default typography and colour layout follow the high-resolution reference
#' style used in `inst/examples/replot_model_consistency_from_rds.R` (shared
#' colour scale via \pkg{patchwork}, framed vertical colour bar, larger axis
#' fonts). All parameters are exposed so you can reproduce publication exports
#' without sourcing the example script.
#'
#' @param x A `cast_consistency` object.
#' @param species Optional character. Species name for the title.
#' @param font_family Base font family (passed to ggplot2). Default `"Arial"`.
#' @param use_bold Logical. Use bold face for titles, axes, legend, and cell
#'   values. Default `TRUE`.
#' @param font_base Base theme size. Default `14`.
#' @param font_main_title Outer patchwork title size. Default `20`.
#' @param font_panel_title Panel title size. Default `28`.
#' @param font_axis Axis label size (model names). Default `25`.
#' @param font_cell_value Text size inside tiles. Default `7`.
#' @param value_decimals Decimals printed in each cell. Default `4L`.
#' @param tile_linewidth Grid line width between tiles. Default `0.8`.
#' @param text_white_above Values above this threshold use white cell text.
#'   Default `0.7`.
#' @param consistency_colors Colours for the shared gradient (low to high).
#' @param fill_limits Numeric vector of length 2, fill scale limits. Default
#'   `c(0, 1)`.
#' @param color_stops Positions in `[0, 1]` matching `consistency_colors`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object (or `patchwork` composite).
#' @export
plot.cast_consistency <- function(x,
                                  species = NULL,
                                  font_family = "Arial",
                                  use_bold = TRUE,
                                  font_base = 14L,
                                  font_main_title = 24L,
                                  font_panel_title = 24L,
                                  font_axis = 20L,
                                  font_cell_value = 5,
                                  value_decimals = 4L,
                                  tile_linewidth = 0.8,
                                  text_white_above = 0.7,
                                  consistency_colors = c(
                                    "#F7FBFF", "#DEEBF7", "#C6DBEF",
                                    "#9ECAE1", "#6BAED6", "#4292C6",
                                    "#2171B5", "#084594"
                                  ),
                                  fill_limits = c(0, 1),
                                  color_stops = seq(0, 1,
                                                    length.out = length(
                                                      consistency_colors
                                                    )),
                                  ...) {
  check_suggested("ggplot2", "for plotting")
  check_suggested("patchwork", "for multi-panel layout")

  models <- x$models
  metrics_df <- x$metrics

  metric_defs <- list(
    list(key = "cosine_sim", label = "Cosine Similarity"),
    list(key = "warrens_I",  label = "Warren's I"),
    list(key = "pearson_r",  label = "Pearson's r")
  )

  face_plain_or_bold <- if (isTRUE(use_bold)) "bold" else "plain"

  fmt_val <- function(v) {
    ifelse(
      is.na(v), "NA",
      format(round(v, value_decimals), nsmall = value_decimals)
    )
  }

  panels <- vector("list", length(metric_defs))

  for (mi in seq_along(metric_defs)) {
    mkey   <- metric_defs[[mi]]$key
    mlabel <- metric_defs[[mi]]$label

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
            (metrics_df$model_a == my & metrics_df$model_b == mx),
          ,
          drop = FALSE
        ]
        if (nrow(row) > 0) {
          val <- row[[mkey]][1]
          mat_df$value[ri] <- val
        }
      }
    }

    mat_df$text_color <- ifelse(
      !is.na(mat_df$value) & mat_df$value > text_white_above, "white", "black"
    )
    mat_df$label_text <- fmt_val(mat_df$value)

    p <- ggplot2::ggplot(
      mat_df,
      ggplot2::aes(x = .data$model_x, y = .data$model_y, fill = .data$value)
    ) +
      ggplot2::geom_tile(
        color = "white",
        linewidth = tile_linewidth
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label_text, color = .data$text_color),
        size = font_cell_value,
        family = font_family,
        fontface = face_plain_or_bold,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_identity() +
      ggplot2::labs(title = mlabel, x = NULL, y = NULL) +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::theme_minimal(
        base_size = font_base,
        base_family = font_family
      ) +
      ggplot2::theme(
        text = ggplot2::element_text(family = font_family),
        plot.title = ggplot2::element_text(
          family = font_family,
          face = face_plain_or_bold,
          hjust = 0.5,
          size = font_panel_title
        ),
        plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
        panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
        axis.text.x = ggplot2::element_text(
          family = font_family,
          face = face_plain_or_bold,
          angle = 45,
          hjust = 1,
          size = font_axis
        ),
        axis.text.y = ggplot2::element_text(
          family = font_family,
          face = face_plain_or_bold,
          size = font_axis
        ),
        legend.text = ggplot2::element_text(
          family = font_family,
          face = face_plain_or_bold
        ),
        legend.title = ggplot2::element_text(
          family = font_family,
          face = face_plain_or_bold
        ),
        legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
        legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
        panel.grid = ggplot2::element_blank(),
        legend.position = "right",
        plot.margin = ggplot2::margin(8, 8, 8, 8)
      )

    panels[[mi]] <- p
  }

  sp_title <- if (!is.null(species)) {
    sprintf("Inter-model Spatial Consistency \u2014 %s", gsub("_", " ", species))
  } else {
    "Inter-model Spatial Consistency"
  }

  combined <- panels[[1]] + panels[[2]] + panels[[3]] +
    patchwork::plot_layout(nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = sp_title,
      theme = ggplot2::theme(
        text = ggplot2::element_text(family = font_family),
        plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
        plot.title = ggplot2::element_text(
          family = font_family,
          face = face_plain_or_bold,
          size = font_main_title,
          hjust = 0.5
        )
      )
    ) &
    ggplot2::scale_fill_gradientn(
      colours   = consistency_colors,
      values    = color_stops,
      limits    = fill_limits,
      breaks    = seq(fill_limits[1], fill_limits[2], by = 0.1),
      na.value  = "grey90",
      name      = NULL,
      oob       = scales::squish,
      guide     = ggplot2::guide_colorbar(
        frame.colour    = "black",
        ticks.colour    = "black",
        ticks.linewidth = 0.5,
        frame.linewidth = 0.5,
        barheight       = grid::unit(10, "cm"),
        barwidth        = grid::unit(0.6, "cm")
      )
    )

  combined
}
