# ==============================================================================
# Re-plot model consistency heatmaps from a saved cast_result.rds
#
# Does NOT re-train models. Loads predictions from RDS, calls cast_consistency()
# (lightweight matrix operations), and draws a three-panel figure with sizes
# you control below.
#
# Usage:
#   1. Edit the CONFIG block (paths, fonts, dimensions).
#   2. Run:  Rscript inst/examples/replot_model_consistency_from_rds.R
#      or source() this file from RStudio.
#
# Requires: castSDM, ggplot2, patchwork
#
# Typography: set FONT_FAMILY = "Arial" and USE_BOLD = FALSE for Arial Regular.
# Set USE_BOLD <- TRUE to use bold for outer title, panel titles, axes, legend,
# and cell values. On Linux, install fonts or set FONT_FAMILY to a sans alias.
# ==============================================================================

# ──────────────────────────────────────────────────────────────────────────────
# CONFIG — edit only this block (manual font / layout tuning)
# ──────────────────────────────────────────────────────────────────────────────

RDS_PATH <- "E:/Package/castSDM_multi_species/Ovis_ammon/cast_result.rds"
OUT_PATH <- "E:/Package/castSDM_multi_species/Ovis_ammon/figures/model_consistency_replot.png"

# Species label in main title (underscores → spaces)
SPECIES_NAME <- "Ovis ammon"

# Output size (inches) and resolution
WIDTH_IN   <- 18
HEIGHT_IN  <- 7
DPI        <- 1200

# Font: all text uses this family (Windows: "Arial"; macOS often same; Linux may need "sans")
FONT_FAMILY <- "Arial"

# Bold switch: FALSE = Arial Regular everywhere; TRUE = bold for titles, axes, and cell numbers
USE_BOLD <- TRUE

# Font sizes (ggplot2 / grid units; increase for larger text)
FONT_BASE          <- 14L   # theme_minimal(base_size)
FONT_MAIN_TITLE    <- 20L   # patchwork outer title
FONT_PANEL_TITLE   <- 28L   # each panel title (Cosine / Warren / Pearson)
FONT_AXIS          <- 25L  # model names on axes
FONT_CELL_VALUE    <- 7  # numbers inside tiles (was ~3.5 in package plot)

# How many decimals to print in each cell
VALUE_DECIMALS <- 4L

# Tile border
TILE_LINEWIDTH <- 0.8

# Text colour switch inside cells (same logic as package: high values → white)
TEXT_WHITE_ABOVE <- 0.7

# Optional: only use these models (NULL = all models in $predict)
MODELS_SUBSET <- NULL

# Heatmap colour ramp — must match castSDM::plot.cast_consistency() (model-consistency.R)
CONSISTENCY_COLORS <- c(
  "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1",
  "#6BAED6", "#4292C6", "#2171B5", "#084594"
)
FILL_LIMITS <- c(0, 1)
# Anchor each colour to an exact position on [0, 1] so tiles and colorbar use one mapping
COLOR_STOPS <- seq(0, 1, length.out = length(CONSISTENCY_COLORS))

# ──────────────────────────────────────────────────────────────────────────────
# Code (no need to edit below unless you customize the plot further)
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(castSDM)
})

if (!file.exists(RDS_PATH)) {
  stop("RDS not found: ", RDS_PATH, call. = FALSE)
}

res <- readRDS(RDS_PATH)
pred <- res$predict
if (is.null(pred) || !inherits(pred, "cast_predict")) {
  stop(
    "RDS must contain a cast_predict object in $predict (spatial predictions).",
    call. = FALSE
  )
}

cons <- cast_consistency(pred, models = MODELS_SUBSET)
models    <- cons$models
metrics_df <- cons$metrics
n <- length(models)

metric_defs <- list(
  list(key = "cosine_sim", label = "Cosine Similarity"),
  list(key = "warrens_I",  label = "Warren's I"),
  list(key = "pearson_r",  label = "Pearson's r")
)

fmt_val <- function(v) {
  ifelse(is.na(v), "", format(round(v, VALUE_DECIMALS), nsmall = VALUE_DECIMALS))
}

face_plain_or_bold <- if (isTRUE(USE_BOLD)) "bold" else "plain"

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
        mat_df$value[ri] <- if (is.na(val)) 0 else val
      }
    }
  }

  mat_df$text_color <- ifelse(
    !is.na(mat_df$value) & mat_df$value > TEXT_WHITE_ABOVE, "white", "black"
  )
  mat_df$label_text <- fmt_val(mat_df$value)

  p <- ggplot(mat_df, aes(x = .data$model_x, y = .data$model_y, fill = .data$value)) +
    geom_tile(color = "white", linewidth = TILE_LINEWIDTH) +
    geom_text(
      aes(label = .data$label_text, color = .data$text_color),
      size = FONT_CELL_VALUE,
      family = FONT_FAMILY,
      fontface = face_plain_or_bold,
      show.legend = FALSE
    ) +
    scale_color_identity() +
    labs(title = mlabel, x = NULL, y = NULL) +
    coord_fixed(expand = FALSE) +
    theme_minimal(base_size = FONT_BASE, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY),
      plot.title = element_text(
        family = FONT_FAMILY,
        face = face_plain_or_bold,
        hjust = 0.5,
        size = FONT_PANEL_TITLE
      ),
      axis.text.x = element_text(
        family = FONT_FAMILY,
        face = face_plain_or_bold,
        angle = 45,
        hjust = 1,
        size = FONT_AXIS
      ),
      axis.text.y = element_text(
        family = FONT_FAMILY,
        face = face_plain_or_bold,
        size = FONT_AXIS
      ),
      legend.text = element_text(family = FONT_FAMILY, face = face_plain_or_bold),
      legend.title = element_text(family = FONT_FAMILY, face = face_plain_or_bold),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.margin = margin(8, 8, 8, 8)
    )

  panels[[mi]] <- p
}

sp_title <- sprintf("Inter-model Spatial Consistency — %s", SPECIES_NAME)

# One fill scale for all panels: same limits, stops, and guide as package plot → tiles match colorbar
combined <- panels[[1]] + panels[[2]] + panels[[3]] +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(
    title = sp_title,
    theme = theme(
      text = element_text(family = FONT_FAMILY),
      plot.title = element_text(
        family = FONT_FAMILY,
        face = face_plain_or_bold,
        size = FONT_MAIN_TITLE,
        hjust = 0.5
      )
    )
  ) &
  scale_fill_gradientn(
    colours   = CONSISTENCY_COLORS,
    values    = COLOR_STOPS,
    limits    = FILL_LIMITS,
    breaks    = seq(0, 1, by = 0.1),
    na.value  = "grey90",
    name      = NULL,
    oob       = squish,
    guide     = guide_colorbar(
      frame.colour    = "black",
      ticks.colour    = "black",
      ticks.linewidth = 0.5,
      frame.linewidth = 0.5,
      barheight       = grid::unit(10, "cm"),
      barwidth        = grid::unit(0.6, "cm")
    )
  )

dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_PATH, combined, width = WIDTH_IN, height = HEIGHT_IN, dpi = DPI)
message("Saved: ", OUT_PATH)
