#' XGBoost + SHAP explainability for species distribution inputs
#'
#' Trains a gradient-boosted tree on `presence` ~ environmental variables (the
#' same columns used elsewhere in castSDM). Computes SHAP contributions and
#' SHAP interaction values on a held-out subset, then [plot.cast_shap()] can
#' reproduce the **Impact intensity** network and **Impact direction** waterfall
#' (single hold-out row by default: bias + feature SHAPs add to the logit margin;
#' not the same as column means across rows).
#'
#' @param data A `data.frame` with `presence` (0/1) and numeric predictors.
#' @param response Name of the binary response column. Default `"presence"`.
#' @param env_vars Predictor column names. If `NULL`, columns are chosen in
#'   order: (1) if `dag` is set, **`intersect(dag$nodes, names(data))`** -- the
#'   same convention as [cast_fit()] for traditional models; (2) else if
#'   `screen` is set, `intersect(get_env_vars(), screen$selected)`; (3) else
#'   [get_env_vars()].
#' @param dag Optional [cast_dag]; when supplied with `env_vars = NULL`,
#'   fixes SHAP inputs to **DAG nodes** so the XGBoost surrogate uses the **same
#'   raw environmental columns** as `cast_fit(..., models = "rf", ...)` (not
#'   `cast_features()` interaction columns).
#' @param screen Optional [cast_screen]; used only when `dag` is `NULL` and
#'   `env_vars` is `NULL` (legacy screening of numeric env names).
#' @param nrounds Maximum boosting rounds. Default `400`.
#' @param max_depth,eta,subsample,colsample_bytree Passed to \pkg{xgboost}
#'   (see [xgboost::xgb.train()]).
#' @param test_fraction Fraction held out for SHAP visualization. Default `0.2`.
#' @param seed Optional integer seed.
#' @param verbose Passed to [xgboost::xgb.train()].
#'
#' @return A `cast_shap` object: `shap` (matrix, rows = test samples, columns =
#'   features only), `bias_shap` (numeric vector, length `nrow(shap)`: TreeSHAP
#'   bias / baseline margin contribution per row, sums with `shap` to the
#'   logit margin), `shap_interaction` (averaged absolute interaction matrix
#'   `p x p`), `feature_names`, `base_score` (mean of `bias_shap`), `xgb_model`,
#'   `expected_value` (mean predicted probability on test rows), `train_matrix`,
#'   `label`, `engine` (`"xgboost"`), `method` (`"xgboost"`), `fitted_model`
#'   (`NULL`), `shap_engine_note` (`NULL`). Objects from [cast_shap_fit()]
#'   reuse the same slots with `engine`/`method` describing the source model.
#'
#' @seealso [plot.cast_shap()], [cast_fit()], [cast_screen()]
#'
#' @export
cast_shap_xgb <- function(data,
                          response = "presence",
                          env_vars = NULL,
                          screen = NULL,
                          dag = NULL,
                          nrounds = 400L,
                          max_depth = 6L,
                          eta = 0.05,
                          subsample = 0.8,
                          colsample_bytree = 0.8,
                          test_fraction = 0.2,
                          seed = NULL,
                          verbose = FALSE) {
  check_suggested("xgboost", "for SHAP explainability")

  if (!is.null(env_vars)) {
    env_vars <- intersect(env_vars, names(data))
  } else if (!is.null(dag)) {
    env_vars <- intersect(dag$nodes, names(data))
  } else if (!is.null(screen)) {
    env_vars <- intersect(
      get_env_vars(data, response = response),
      screen$selected
    )
  } else {
    env_vars <- get_env_vars(data, response = response)
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

  bias_per_obs <- if (length(bias_col)) {
    as.numeric(shap_mat[, bias_col[1L]])
  } else {
    rep(NA_real_, nrow(shap_only))
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
    mean(bias_per_obs, na.rm = TRUE)
  } else {
    mean(predict(fit, dtest))
  }

  cap <- paste0(
    "Separate XGBoost classifier (not cast_fit RF/CAST). ",
    "TreeSHAP on p=", length(env_vars), " raw env columns",
    if (!is.null(dag)) " (= intersect(dag$nodes, names(data)))" else "",
    ". Waterfall: y = logit margin; P(presence) from plogis in subtitle."
  )

  structure(
    list(
      shap = shap_only,
      bias_shap = bias_per_obs,
      shap_interaction = mean_inter,
      feature_names = colnames(X_te),
      base_score = base_score,
      xgb_model = fit,
      expected_value = mean(predict(fit, dtest)),
      train_matrix = X_te,
      label = y_te,
      response = response,
      engine = "xgboost",
      method = "xgboost",
      fitted_model = NULL,
      shap_engine_note = NULL,
      feature_space = "raw_env",
      shap_scale = "logit_margin",
      n_features = length(env_vars),
      n_interactions = 0L,
      shap_plot_caption = cap
    ),
    class = "cast_shap"
  )
}


#' Plot SHAP diagnostics (interaction network or waterfall)
#'
#' @param x Object from \code{\link{cast_shap_xgb}} or \code{\link{cast_shap_fit}}.
#' @param type `"interaction_network"` (circular impact-intensity style) or
#'   `"waterfall"` (additive SHAP decomposition, impact-direction style).
#' @param top_n Number of features in the figure. Default `20`.
#' @param pos_color,neg_color Waterfall bar colours. Defaults match the
#'   reference (`#f1696d`, `#44c3df`).
#' @param waterfall_row Integer in `1:nrow(x$shap)` (a **hold-out test row**
#'   used only for SHAP). Used when `type = "waterfall"` and
#'   `waterfall_aggregate` is `FALSE`. The waterfall then matches the usual
#'   **single-observation** SHAP path (bias + feature SHAPs = logit margin for
#'   that row). Default `1`.
#' @param waterfall_aggregate If `TRUE`, bars use **column means** of `x$shap`
#'   and `base_score` as baseline (a summary across the hold-out set, **not** a
#'   single prediction path). Default `FALSE`.
#' @param ... Reserved.
#'
#' @return A `ggplot` object.
#' @export
plot.cast_shap <- function(x,
                           type = c("interaction_network", "waterfall"),
                           top_n = 20L,
                           pos_color = "#f1696d",
                           neg_color = "#44c3df",
                           waterfall_row = 1L,
                           waterfall_aggregate = FALSE,
                           ...) {
  type <- match.arg(type)
  check_suggested("ggplot2", "for plotting")

  if (type == "waterfall") {
    return(.plot_shap_waterfall(
      x, top_n, pos_color, neg_color,
      waterfall_row = waterfall_row,
      waterfall_aggregate = isTRUE(waterfall_aggregate)
    ))
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
    cli::cli_abort(
      "Interaction matrix missing or invalid in {.cls cast_shap} object."
    )
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

  # Keep labels close to nodes while giving larger nodes a bit more clearance.
  # This reduces node-label gap but avoids overlap on dominant nodes.
  max_imp <- max(imp, 1e-8)
  label_r <- 1.18 + 0.08 * (imp / max_imp)
  pos$label_x <- label_r * cos(ang - pi / 2)
  pos$label_y <- label_r * sin(ang - pi / 2)

  # hjust/vjust depend on angular position for neat radial placement
  cos_v <- cos(ang - pi / 2)
  sin_v <- sin(ang - pi / 2)
  pos$label_hjust <- ifelse(
    cos_v > 0.15, 0,
    ifelse(cos_v < -0.15, 1, 0.5)
  )
  pos$label_vjust <- ifelse(
    abs(cos_v) <= 0.15,
    ifelse(sin_v > 0, -0.2, 1.2),
    0.5
  )

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
        label = .data$feature,
        hjust = .data$label_hjust,
        vjust = .data$label_vjust
      ),
      inherit.aes = TRUE,
      size = 4.2,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_gradientn(
      colours = greens,
      limits = c(0, max_imp),
      name = "Vimp",
      oob = scales::squish
    ) +
    ggplot2::scale_size(range = c(2, 12), guide = "none") +
    ggplot2::coord_fixed(
      xlim = c(-2.55, 2.55),
      ylim = c(-2.05, 2.05),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::labs(
      title = sprintf(
        "Impact intensity (Top %d of %d features)",
        n,
        ncol(as.matrix(x$shap))
      ),
      subtitle = x$shap_plot_caption %||% NULL,
      x = NULL, y = NULL
    )

  if (!is.null(seg) && nrow(seg) > 0) {
    p <- p + ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title = "Importance\nLow \u2500\u2500\u25b6 High\nVimp",
        barwidth = ggplot2::unit(5.0, "cm"),
        barheight = ggplot2::unit(0.4, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        order = 1L
      ),
      colour = ggplot2::guide_colorbar(
        title = "Interaction Intensity\nWeak \u2500\u2500\u25b6 Strong\nVint",
        barwidth = ggplot2::unit(5.0, "cm"),
        barheight = ggplot2::unit(0.4, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        order = 2L
      )
    )
  } else {
    p <- p + ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        title = "Importance\nLow \u2500\u2500\u25b6 High\nVimp",
        barwidth = ggplot2::unit(5.0, "cm"),
        barheight = ggplot2::unit(0.4, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom"
      )
    )
  }

  p <- p +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0, size = 12
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0, size = 7.5, color = "gray35", lineheight = 1.15
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.spacing = ggplot2::unit(0.3, "cm"),
      legend.spacing.x = ggplot2::unit(1.8, "cm"),
      legend.justification = "center",
      legend.title = ggplot2::element_text(size = 10, face = "bold", lineheight = 1.1),
      legend.text = ggplot2::element_text(size = 8.5),
      plot.margin = ggplot2::margin(10, 26, 10, 26)
    )

  p
}


#' @keywords internal
#' @noRd
.shap_waterfall_arrow_poly <- function(x_c, y0, y1, xhw = 0.18, tip_frac = 0.14) {
  h <- y1 - y0
  if (abs(h) < 1e-10) {
    return(data.frame(
      x = c(x_c - xhw, x_c + xhw, x_c + xhw, x_c - xhw),
      y = rep(y0, 4)
    ))
  }
  tip_len <- min(abs(h) * tip_frac, abs(h) * 0.42)
  if (h > 0) {
    yt <- y1 - tip_len
    data.frame(
      x = c(x_c - xhw, x_c + xhw, x_c + xhw, x_c, x_c - xhw),
      y = c(y0, y0, yt, y1, yt)
    )
  } else {
    yt <- y1 + tip_len
    data.frame(
      x = c(x_c - xhw, x_c + xhw, x_c + xhw, x_c, x_c - xhw),
      y = c(y0, y0, yt, y1, yt)
    )
  }
}


#' @keywords internal
#' @noRd
.plot_shap_waterfall <- function(x, top_n, pos_color, neg_color,
                                 waterfall_row = 1L,
                                 waterfall_aggregate = FALSE) {
  sh <- x$shap
  feats <- colnames(sh)
  n_obs <- nrow(sh)

  if (isTRUE(waterfall_aggregate)) {
    ms <- colMeans(sh)
    base <- as.numeric(x$base_score)[1]
    row_lab <- "mean across hold-out rows"
  } else {
    wi <- as.integer(waterfall_row)[1L]
    if (is.na(wi) || wi < 1L || wi > n_obs) {
      cli::cli_abort(
        "{.arg waterfall_row} must be an integer in 1..{n_obs} (hold-out SHAP rows)."
      )
    }
    ms <- as.numeric(sh[wi, , drop = TRUE])
    names(ms) <- colnames(sh)
    bias <- x$bias_shap
    if (is.null(bias) || length(bias) != n_obs) {
      bias <- rep(as.numeric(x$base_score)[1], n_obs)
    }
    base <- as.numeric(bias[wi])
    row_lab <- sprintf("hold-out row %d", wi)
  }

  ord <- order(abs(ms), decreasing = TRUE)[seq_len(min(top_n, length(ms)))]
  ms <- ms[ord]
  feats <- feats[ord]

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
  final_margin <- utils::tail(ends, 1L)
  df$val_lab <- ifelse(
    df$val >= 0,
    sprintf("+%.4f", df$val),
    sprintf("%.4f", df$val)
  )

  poly_parts <- vector("list", nrow(df))
  for (k in seq_len(nrow(df))) {
    pr <- .shap_waterfall_arrow_poly(df$x[k], df$start[k], df$end[k])
    pr$piece <- k
    pr$direction <- df$direction[k]
    poly_parts[[k]] <- pr
  }
  poly_df <- do.call(rbind, poly_parts)

  y_rng <- range(c(df$ymin, df$ymax, base, final_margin))
  y_pad <- max(diff(y_rng) * 0.06, 1e-6)
  y_top <- max(df$ymax) + y_pad * 2.6
  y_bot <- min(df$ymin) - y_pad * 0.65

  df$lab_y <- ifelse(
    df$val >= 0,
    df$ymax + y_pad * 0.55,
    df$ymin - y_pad * 0.55
  )

  eng <- x$engine %||% "xgboost"
  if (identical(eng, "xgboost")) {
    sub_txt <- sprintf(
      "XGBoost margin (log-odds); %s | P(presence): %.4f \u2192 %.4f",
      row_lab,
      stats::plogis(base),
      stats::plogis(final_margin)
    )
    ylab_main <- "SHAP value"
  } else {
    note <- x$shap_engine_note %||% ""
    sub_txt <- sprintf(
      "%s | %s | endpoint pred.: %.4f (baseline + sum SHAP).",
      note,
      row_lab,
      final_margin
    )
    ylab_main <- "SHAP value"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x))

  # Baseline dashed horizontal
  p <- p +
    ggplot2::geom_hline(
      yintercept = base,
      linetype = "dotted",
      color = "gray40",
      linewidth = 0.35
    ) +
    ggplot2::geom_hline(
      yintercept = final_margin,
      linetype = "dotted",
      color = "gray40",
      linewidth = 0.35
    )

  # Vertical dotted connectors from bar tips
  for (k in seq_len(nrow(df))) {
    conn_y <- df$end[k]
    p <- p + ggplot2::annotate(
      "segment",
      x = df$x[k] + 0.22, xend = min(df$x[k] + 1, max(df$x) + 0.45),
      y = conn_y, yend = conn_y,
      linetype = "dotted", color = "grey60", linewidth = 0.3
    )
  }

  # Arrow-shaped bar polygons
  p <- p +
    ggplot2::geom_polygon(
      data = poly_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$piece,
        fill = .data$direction
      ),
      inherit.aes = FALSE,
      colour = "#2a2a2a",
      linewidth = 0.28
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = .data$lab_y,
        label = .data$val_lab
      ),
      size = 2.8,
      color = "#1a1a1a",
      inherit.aes = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = "Impact direction",
      values = c("Positive" = pos_color, "Negative" = neg_color),
      breaks = c("Positive", "Negative"),
      labels = c(
        "Positive" = "Positive",
        "Negative" = "Negative"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = df$x,
      labels = df$feature,
      expand = ggplot2::expansion(mult = 0.05)
    ) +
    ggplot2::scale_y_continuous(
      name = NULL,
      sec.axis = ggplot2::sec_axis(~ ., name = ylab_main)
    ) +
    ggplot2::expand_limits(
      y = c(y_bot, y_top),
      x = c(0.2, max(df$x) + 0.5)
    )

  # f(x) label on LEFT side, vertical (rotated 90)
  p <- p +
    ggplot2::annotate(
      "text",
      x = 0.15,
      y = final_margin,
      label = sprintf("italic(f(x))==%.4f", final_margin),
      parse = TRUE,
      hjust = 0.5,
      vjust = 1.3,
      size = 3.3,
      color = "gray20",
      angle = 90
    )

  # E[f(x)] label on RIGHT side, vertical (rotated 90)
  p <- p +
    ggplot2::annotate(
      "text",
      x = max(df$x) + 0.55,
      y = base,
      label = sprintf("italic(E)*'['*italic(f(x))*']'==%.4f", base),
      parse = TRUE,
      hjust = 0.5,
      vjust = -0.3,
      size = 3.3,
      color = "gray20",
      angle = 90
    )

  p <- p +
    ggplot2::labs(
      title = NULL,
      subtitle = sub_txt,
      caption = x$shap_plot_caption %||% NULL,
      x = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      # Outer box frame
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.6
      ),
      # Remove all grid lines
      panel.grid = ggplot2::element_blank(),
      # Subtitle (small, informational)
      plot.subtitle = ggplot2::element_text(
        hjust = 0, size = 7.5, color = "gray45",
        margin = ggplot2::margin(b = 4)
      ),
      plot.caption = ggplot2::element_text(
        size = 7, color = "gray40", hjust = 0, lineheight = 1.2
      ),
      # X-axis labels rotated
      axis.text.x = ggplot2::element_text(
        angle = 45, vjust = 1, hjust = 1,
        size = 8.5, color = "gray15"
      ),
      axis.ticks.x = ggplot2::element_line(
        color = "grey50", linewidth = 0.3
      ),
      axis.ticks.y = ggplot2::element_line(
        color = "grey50", linewidth = 0.3
      ),
      # Right-side Y axis title
      axis.title.y.right = ggplot2::element_text(
        color = "gray25", size = 10, face = "bold",
        angle = 270, vjust = 1.5
      ),
      axis.title.y.left = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 14, 10, 10),
      # Legend at bottom with arrow labels
      legend.position = "bottom",
      legend.title = ggplot2::element_text(
        face = "bold", size = 9.5
      ),
      legend.text = ggplot2::element_text(size = 9)
    )

  p
}


#' Write SHAP summary tables to CSV
#'
#' @param x Object from \code{\link{cast_shap_xgb}} or \code{\link{cast_shap_fit}}.
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
