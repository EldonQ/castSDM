#' XGBoost + SHAP explainability for species distribution inputs
#'
#' Trains a gradient-boosted tree on `presence` ~ environmental variables (the
#' same columns used elsewhere in castSDM). Computes SHAP contributions and
#' SHAP interaction values on a held-out subset, then [plot.cast_shap()] can
#' reproduce the **Impact intensity** network and **Impact direction** waterfall
#' layouts from the Python reference (see `inst/examples/../Ref/DAG图.md`).
#'
#' @param data A `data.frame` with `presence` (0/1) and numeric predictors.
#' @param response Name of the binary response column. Default `"presence"`.
#' @param env_vars Predictor column names; `NULL` uses [get_env_vars()].
#' @param screen Optional [cast_screen]; when set, only `screen$selected`
#'   variables are used (aligned with CAST screening).
#' @param nrounds Maximum boosting rounds. Default `400`.
#' @param max_depth,eta,subsample,colsample_bytree Passed to \pkg{xgboost}
#'   (see [xgboost::xgb.train()]).
#' @param test_fraction Fraction held out for SHAP visualization. Default `0.2`.
#' @param seed Optional integer seed.
#' @param verbose Passed to [xgboost::xgb.train()].
#'
#' @return A `cast_shap` object: `shap` (matrix, rows = test samples),
#'   `shap_interaction` (averaged absolute interaction matrix `p x p`),
#'   `feature_names`, `base_score`, `xgb_model`, `expected_value` (mean
#'   prediction on probability scale for reference), `train_matrix`, `label`.
#'
#' @seealso [plot.cast_shap()], [cast_fit()], [cast_screen()]
#'
#' @export
cast_shap_xgb <- function(data,
                          response = "presence",
                          env_vars = NULL,
                          screen = NULL,
                          nrounds = 400L,
                          max_depth = 6L,
                          eta = 0.05,
                          subsample = 0.8,
                          colsample_bytree = 0.8,
                          test_fraction = 0.2,
                          seed = NULL,
                          verbose = FALSE) {
  check_suggested("xgboost", "for SHAP explainability")

  env_vars <- env_vars %||% get_env_vars(data, response = response)
  if (!is.null(screen)) {
    env_vars <- intersect(env_vars, screen$selected)
  }
  if (length(env_vars) < 2L) {
    cli::cli_abort("Need at least 2 predictor columns for SHAP.")
  }

  df <- data[, unique(c(env_vars, response)), drop = FALSE]
  df <- stats::na.omit(df)
  y <- as.numeric(df[[response]])
  X <- as.matrix(df[, env_vars, drop = FALSE])
  storage.mode(X) <- "double"
  storage.mode(y) <- "double"

  n <- nrow(X)
  if (n < 30L) {
    cli::cli_abort("Too few complete rows ({n}) for stable SHAP.")
  }

  if (!is.null(seed)) set.seed(seed)
  n_te <- max(5L, floor(n * test_fraction))
  ii <- sample.int(n, size = n - n_te)
  te <- setdiff(seq_len(n), ii)
  X_tr <- X[ii, , drop = FALSE]
  y_tr <- y[ii]
  X_te <- X[te, , drop = FALSE]
  y_te <- y[te]

  dtrain <- xgboost::xgb.DMatrix(X_tr, label = y_tr)
  dtest <- xgboost::xgb.DMatrix(X_te, label = y_te)

  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = max_depth,
    eta = eta,
    subsample = subsample,
    colsample_bytree = colsample_bytree
  )

  fit <- xgboost::xgb.train(
    params,
    data = dtrain,
    nrounds = nrounds,
    verbose = verbose
  )

  shap_mat <- predict(fit, dtest, predcontrib = TRUE)
  if (!is.matrix(shap_mat)) {
    shap_mat <- matrix(shap_mat, nrow = nrow(X_te))
  }
  bias_col <- which(colnames(shap_mat) %in% c("BIAS", "bias"))
  if (length(bias_col) == 0L) bias_col <- ncol(shap_mat)
  shap_only <- shap_mat[, setdiff(seq_len(ncol(shap_mat)), bias_col), drop = FALSE]
  if (is.null(colnames(shap_only))) {
    colnames(shap_only) <- env_vars
  }

  p <- ncol(X_te)
  shap_int_raw <- tryCatch(
    predict(fit, dtest, predinteraction = TRUE),
    error = function(e) NULL
  )

  mean_inter <- matrix(0, nrow = p, ncol = p)
  dimnames(mean_inter) <- list(colnames(X_te), colnames(X_te))

  if (!is.null(shap_int_raw)) {
    tot <- length(shap_int_raw)
    per <- tot / nrow(X_te)
    sqr <- as.integer(round(sqrt(per)))
    if (sqr > 1L && abs(sqr * sqr - per) < 1e-4) {
      arr <- array(shap_int_raw, dim = c(nrow(X_te), sqr, sqr))
      arr <- arr[, -1L, -1L, drop = FALSE]
      p2 <- dim(arr)[2]
      for (i in seq_len(p2)) {
        for (j in seq_len(p2)) {
          mean_inter[i, j] <- mean(abs(arr[, i, j]))
        }
      }
    } else {
      Ms <- abs(shap_only)
      mean_inter <- (t(Ms) %*% Ms) / pmax(1L, nrow(Ms))
    }
  }

  base_score <- if (length(bias_col)) {
    mean(shap_mat[, bias_col[1]])
  } else {
    mean(predict(fit, dtest))
  }

  structure(
    list(
      shap = shap_only,
      shap_interaction = mean_inter,
      feature_names = colnames(X_te),
      base_score = base_score,
      xgb_model = fit,
      expected_value = mean(predict(fit, dtest)),
      train_matrix = X_te,
      label = y_te,
      response = response
    ),
    class = "cast_shap"
  )
}


#' Plot SHAP diagnostics (interaction network or waterfall)
#'
#' @param x Object returned by \code{\link{cast_shap_xgb}}.
#' @param type `"interaction_network"` (circular impact-intensity style) or
#'   `"waterfall"` (mean SHAP contributions, impact-direction style).
#' @param top_n Number of features in the figure. Default `20`.
#' @param pos_color,neg_color Waterfall bar colours. Defaults match the
#'   reference (`#f1696d`, `#44c3df`).
#' @param ... Reserved.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_shap <- function(x,
                           type = c("interaction_network", "waterfall"),
                           top_n = 20L,
                           pos_color = "#f1696d",
                           neg_color = "#44c3df",
                           ...) {
  type <- match.arg(type)
  check_suggested("ggplot2", "for plotting")

  if (type == "waterfall") {
    return(.plot_shap_waterfall(x, top_n, pos_color, neg_color))
  }
  .plot_shap_interaction_network(x, top_n)
}


#' @keywords internal
#' @noRd
.plot_shap_interaction_network <- function(x, top_n) {
  feats <- x$feature_names
  imp <- colMeans(abs(x$shap))
  inter <- x$shap_interaction
  if (is.null(inter) || !is.matrix(inter)) {
    cli::cli_abort("Interaction matrix missing; check xgboost predinteraction.")
  }

  if (top_n < length(feats)) {
    ord <- order(imp, decreasing = TRUE)[seq_len(top_n)]
    feats <- feats[ord]
    imp <- imp[ord]
    inter <- inter[ord, ord, drop = FALSE]
  }

  n <- length(feats)
  ang <- seq(0, 2 * pi, length.out = n + 1)[seq_len(n)]
  pos <- data.frame(
    x = cos(ang - pi / 2),
    y = sin(ang - pi / 2),
    feature = feats,
    importance = imp,
    stringsAsFactors = FALSE
  )
  pos$label_x <- pos$x * 1.18
  pos$label_y <- pos$y * 1.18

  max_imp <- max(imp, 1e-8)
  pair <- list()
  pk <- 0L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i >= j) next
      v <- (inter[i, j] + inter[j, i]) / 2
      if (v <= 0) next
      pk <- pk + 1L
      pair[[pk]] <- data.frame(
        x = pos$x[i], y = pos$y[i],
        xend = pos$x[j], yend = pos$y[j],
        w = v,
        stringsAsFactors = FALSE
      )
    }
  }
  seg <- if (pk > 0L) do.call(rbind, pair) else NULL
  wvals <- if (!is.null(seg)) seg$w else numeric(0)
  wmax <- if (length(wvals)) max(wvals) else 1
  wmin <- if (length(wvals)) min(wvals[wvals > 0], na.rm = TRUE) else 0
  if (!is.finite(wmin)) wmin <- 0

  greens <- grDevices::colorRampPalette(
    c("#edf8e9", "#31a354")
  )(64)

  # Single ggplot + one coord_fixed (avoids "coordinate system already present"
  # noise). scale_linewidth() does not support `oob` on all ggplot2 versions.
  p <- ggplot2::ggplot(
    pos,
    ggplot2::aes(x = .data$x, y = .data$y)
  )

  if (!is.null(seg) && nrow(seg) > 0) {
    seg$lw_map <- 0.4 + 4.6 * (seg$w - wmin) / (wmax - wmin + 1e-8)
    seg$w_n <- (seg$w - wmin) / (wmax - wmin + 1e-8)
    purples <- grDevices::colorRampPalette(
      c("#fcfbfd", "#dadaeb", "#9e9ac8", "#6a51a3", "#3f007d")
    )(64)
    p <- p +
      ggplot2::geom_segment(
        data = seg,
        ggplot2::aes(
          x = .data$x, y = .data$y,
          xend = .data$xend, yend = .data$yend,
          linewidth = .data$lw_map,
          colour = .data$w_n
        ),
        inherit.aes = FALSE,
        lineend = "round",
        alpha = 0.92
      ) +
      ggplot2::scale_linewidth(
        name = "Vint",
        range = c(0.35, 5.2),
        guide = "none"
      ) +
      ggplot2::scale_colour_gradientn(
        colours = purples,
        limits = c(0, 1),
        name = "Vint",
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0.25", "0.5", "0.75", "1"),
        oob = scales::squish
      )
  }

  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$importance,
        fill = .data$importance
      ),
      shape = 21,
      color = "#2f2f2f",
      stroke = 0.4
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$label_x,
        y = .data$label_y,
        label = .data$feature
      ),
      inherit.aes = TRUE,
      size = 3.2,
      fontface = "bold",
      hjust = 0.5
    ) +
    ggplot2::scale_fill_gradientn(
      colours = greens,
      limits = c(0, max_imp),
      name = "Vimp",
      oob = scales::squish
    ) +
    ggplot2::scale_size(range = c(2, 10), guide = "none") +
    ggplot2::coord_fixed(
      xlim = c(-1.55, 1.55),
      ylim = c(-1.55, 1.55),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = sprintf("Impact intensity (Top %d features)", n),
      x = NULL, y = NULL
    )

  if (!is.null(seg) && nrow(seg) > 0) {
    p <- p + ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title = "Vimp",
        barwidth = ggplot2::unit(3.2, "cm"),
        barheight = ggplot2::unit(0.35, "cm"),
        title.position = "top",
        order = 1L
      ),
      colour = ggplot2::guide_colorbar(
        title = "Vint",
        barwidth = ggplot2::unit(3.2, "cm"),
        barheight = ggplot2::unit(0.35, "cm"),
        title.position = "top",
        order = 2L
      )
    )
  } else {
    p <- p + ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title = "Vimp",
        barwidth = ggplot2::unit(3.2, "cm"),
        barheight = ggplot2::unit(0.35, "cm"),
        title.position = "top"
      )
    )
  }

  p <- p +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0, size = 12
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.justification = "center",
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  p
}


#' @keywords internal
#' @noRd
.plot_shap_waterfall <- function(x, top_n, pos_color, neg_color) {
  sh <- x$shap
  feats <- colnames(sh)
  ms <- colMeans(sh)
  ord <- order(abs(ms), decreasing = TRUE)[seq_len(min(top_n, length(ms)))]
  ms <- ms[ord]
  feats <- feats[ord]
  base <- as.numeric(x$base_score)[1]

  ends <- cumsum(ms) + base
  starts <- c(base, utils::head(ends, -1L))
  df <- data.frame(
    feature = feats,
    start = starts,
    end = ends,
    val = ms,
    ymin = pmin(starts, ends),
    ymax = pmax(starts, ends),
    x = seq_along(feats),
    direction = ifelse(ms >= 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )
  final_val <- utils::tail(ends, 1L)
  df$val_lab <- ifelse(
    df$val >= 0,
    sprintf("+%.4g", df$val),
    sprintf("%.4g", df$val)
  )

  y_rng <- range(c(df$ymin, df$ymax, base, final_val))
  y_pad <- max(diff(y_rng) * 0.06, 1e-6)
  y_top <- max(df$ymax) + y_pad * 2.2
  y_bot <- min(df$ymin) - y_pad * 0.5

  bridge <- if (nrow(df) > 1L) {
    data.frame(
      x = utils::head(df$x, -1L) + 0.22,
      xend = utils::tail(df$x, -1L) - 0.22,
      y = utils::head(df$end, -1L),
      yend = utils::head(df$end, -1L),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_vline(
      xintercept = df$x,
      color = "#e0e0e0",
      linewidth = 0.28
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = base),
      linetype = "dotted",
      color = "gray40",
      linewidth = 0.35
    )

  if (!is.null(bridge)) {
    p <- p + ggplot2::geom_segment(
      data = bridge,
      ggplot2::aes(
        x = .data$x, xend = .data$xend,
        y = .data$y, yend = .data$yend
      ),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.35,
      lineend = "round"
    )
  }

  p <- p +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$x - 0.18,
        xmax = .data$x + 0.18,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$direction
      ),
      color = "#2a2a2a",
      linewidth = 0.28
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = ifelse(.data$val >= 0, .data$ymax + y_pad * 0.55, .data$ymin - y_pad * 0.55),
        label = .data$val_lab
      ),
      size = 3.1,
      color = "#1a1a1a",
      inherit.aes = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c("Positive" = pos_color, "Negative" = neg_color),
      breaks = c("Positive", "Negative"),
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = df$x,
      labels = df$feature,
      expand = ggplot2::expansion(mult = 0.05)
    ) +
    ggplot2::expand_limits(
      y = c(y_bot, y_top),
      x = c(0.55, max(df$x) + 0.45)
    ) +
    ggplot2::annotate(
      "text",
      x = 0.65,
      y = y_top,
      label = sprintf("E[f(x)] = %.4f", base),
      hjust = 0,
      size = 3.4,
      color = "gray20",
      fontface = "plain"
    ) +
    ggplot2::annotate(
      "text",
      x = max(df$x) + 0.35,
      y = y_top,
      label = sprintf("f(x) = %.4f", final_val),
      hjust = 1,
      size = 3.4,
      color = "gray20",
      fontface = "plain"
    ) +
    ggplot2::labs(
      title = "Impact direction",
      x = NULL,
      y = "SHAP value"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 8.5,
        color = "gray15"
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0,
        size = 13.5,
        margin = ggplot2::margin(b = 4)
      ),
      plot.margin = ggplot2::margin(10, 14, 10, 10),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        color = "gray35",
        size = 9.5,
        angle = 90
      )
    )

  p
}


#' Write SHAP summary tables to CSV
#'
#' @param x Object returned by \code{\link{cast_shap_xgb}}.
#' @param dir Output directory (created if needed).
#' @param prefix File name prefix. Default `"shap"`.
#'
#' @return Invisibly, paths written.
#' @export
cast_shap_write_csv <- function(x, dir, prefix = "shap") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  imp <- data.frame(
    feature = x$feature_names,
    importance = colMeans(abs(x$shap)),
    stringsAsFactors = FALSE
  )
  imp <- imp[order(imp$importance, decreasing = TRUE), ]
  p1 <- file.path(dir, paste0(prefix, "_feature_importance.csv"))
  utils::write.csv(imp, p1, row.names = FALSE)

  inter <- x$shap_interaction
  rows <- list()
  k <- 0L
  fn <- x$feature_names
  n <- length(fn)
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        feature_1 = fn[i],
        feature_2 = fn[j],
        interaction_value = (inter[i, j] + inter[j, i]) / 2,
        stringsAsFactors = FALSE
      )
    }
  }
  tab <- if (k > 0L) {
    do.call(rbind, rows)
  } else {
    data.frame(
      feature_1 = character(),
      feature_2 = character(),
      interaction_value = numeric(),
      stringsAsFactors = FALSE
    )
  }
  if (nrow(tab) > 0) {
    tab <- tab[order(tab$interaction_value, decreasing = TRUE), ]
  }
  p2 <- file.path(dir, paste0(prefix, "_pairwise_interactions.csv"))
  utils::write.csv(tab, p2, row.names = FALSE)
  invisible(c(p1, p2))
}
