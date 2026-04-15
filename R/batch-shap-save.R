# Shared SHAP figure export for cast_batch() and post-hoc RDS replays.
# Requires cfg$do_shap; cfg$response; cfg$shap_* fields matching cast_batch cfg list.

.cast_shap_one_pair <- function(sh_obj, stem, topn, fig_dir, fig_dpi) {
  if (is.null(sh_obj)) {
    return(invisible(NULL))
  }
  p_net <- tryCatch(
    plot(sh_obj, type = "interaction_network", top_n = topn),
    error = function(e) NULL
  )
  if (!is.null(p_net)) {
    tryCatch(
      ggplot2::ggsave(
        file.path(fig_dir, paste0(stem, "_interaction_network.png")),
        p_net, width = 9, height = 9, dpi = fig_dpi,
        bg = "white", limitsize = FALSE
      ),
      error = function(e) NULL
    )
  }
  bv <- sh_obj$bias_shap
  if (is.null(bv) || length(bv) != nrow(sh_obj$shap)) {
    bv <- rep(as.numeric(sh_obj$base_score)[1], nrow(sh_obj$shap))
  }
  marg <- bv + rowSums(as.matrix(sh_obj$shap))
  wi <- which.min(abs(marg - stats::median(marg)))[1L]
  p_wf <- tryCatch(
    plot(
      sh_obj,
      type = "waterfall",
      top_n = topn,
      waterfall_row = wi
    ),
    error = function(e) NULL
  )
  if (!is.null(p_wf)) {
    tryCatch(
      ggplot2::ggsave(
        file.path(fig_dir, paste0(stem, "_waterfall.png")),
        p_wf, width = 10, height = 6, dpi = fig_dpi,
        bg = "white", limitsize = FALSE
      ),
      error = function(e) NULL
    )
  }
  invisible(NULL)
}

.cast_shap_blank_panel <- function(subtitle) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.5, label = subtitle,
      size = 3.2, color = "grey45"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

.cast_shap_panel_plot <- function(sh_obj, type, topn, title) {
  if (is.null(sh_obj)) {
    return(.cast_shap_blank_panel(title))
  }
  tryCatch(
    {
      if (identical(type, "waterfall")) {
        bv <- sh_obj$bias_shap
        if (is.null(bv) || length(bv) != nrow(sh_obj$shap)) {
          bv <- rep(as.numeric(sh_obj$base_score)[1], nrow(sh_obj$shap))
        }
        marg <- bv + rowSums(as.matrix(sh_obj$shap))
        wi <- which.min(abs(marg - stats::median(marg)))[1L]
        plot(
          sh_obj,
          type = "waterfall",
          top_n = topn,
          waterfall_row = wi
        ) +
          ggplot2::labs(title = title) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"))
      } else {
        plot(sh_obj, type = "interaction_network", top_n = topn) +
          ggplot2::labs(title = title) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"))
      }
    },
    error = function(e) {
      .cast_shap_blank_panel(paste0(title, "\n(error)"))
    }
  )
}

.cast_shap_save_panel_2x3 <- function(sh_xgb, sh_rf, sh_cast, topn, fig_dir, fig_dpi) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(NULL))
  }
  p11 <- .cast_shap_panel_plot(sh_xgb, "interaction_network", topn, "XGBoost - network")
  p12 <- .cast_shap_panel_plot(sh_rf, "interaction_network", topn, "RF - network")
  p13 <- .cast_shap_panel_plot(sh_cast, "interaction_network", topn, "CAST - network")
  p21 <- .cast_shap_panel_plot(sh_xgb, "waterfall", topn, "XGBoost - waterfall")
  p22 <- .cast_shap_panel_plot(sh_rf, "waterfall", topn, "RF - waterfall")
  p23 <- .cast_shap_panel_plot(sh_cast, "waterfall", topn, "CAST - waterfall")
  comb <- patchwork::wrap_plots(p11, p12, p13, p21, p22, p23, nrow = 2L, ncol = 3L)
  tryCatch(
    ggplot2::ggsave(
      file.path(fig_dir, "shap_panel_2x3.png"),
      comb, width = 16, height = 11, dpi = fig_dpi,
      bg = "white", limitsize = FALSE
    ),
    error = function(e) NULL
  )
  invisible(NULL)
}

#' Save SHAP figures for one species (XGB / RF / CAST + 2x3 panel)
#'
#' Writes `shap_*_interaction_network.png`, `shap_*_waterfall.png`, and
#' `shap_panel_2x3.png` (patchwork: row1 = networks, row2 = waterfalls).
#' Used by [cast_batch()] when `do_shap = TRUE`; can be called manually after
#' reading `cast_result.rds` if you replay [cast_prepare()] with the same
#' `train_fraction` and `seed` as the original batch.
#'
#' @param train_df Training data (e.g. `cast_prepare(...)$train`).
#' @param fit,dag,screen Objects from a saved batch result.
#' @param cfg List with `do_shap`, `response`, and `shap_*` fields (same names
#'   as in [cast_batch()]).
#' @param fig_dir Output directory for PNGs.
#' @param fig_dpi Integer DPI.
#' @param seed_i Optional seed passed to SHAP routines.
#'
#' @return `NULL` invisibly.
#' @export
save_cast_batch_shap_outputs <- function(train_df,
                                         fit,
                                         dag,
                                         screen,
                                         cfg,
                                         fig_dir,
                                         fig_dpi,
                                         seed_i = NULL) {
  if (!isTRUE(cfg$do_shap)) {
    return(invisible(NULL))
  }
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  topn <- cfg$shap_plot_top_n

  sh_xgb <- NULL
  sh_rf <- NULL
  sh_cast <- NULL

  if (requireNamespace("xgboost", quietly = TRUE)) {
    sh_xgb <- tryCatch(
      cast_shap_xgb(
        train_df,
        response = cfg$response,
        env_vars = NULL,
        screen = screen,
        dag = dag,
        nrounds = cfg$shap_nrounds,
        max_depth = cfg$shap_max_depth,
        eta = cfg$shap_eta,
        subsample = cfg$shap_subsample,
        colsample_bytree = cfg$shap_colsample_bytree,
        test_fraction = cfg$shap_test_fraction,
        seed = seed_i,
        verbose = cfg$shap_verbose
      ),
      error = function(e) NULL
    )
    .cast_shap_one_pair(sh_xgb, "shap_xgb", topn, fig_dir, fig_dpi)
  }

  if (requireNamespace("fastshap", quietly = TRUE) &&
      requireNamespace("ranger", quietly = TRUE) &&
      "rf" %in% names(fit$models) && !is.null(fit$models$rf$model)) {
    sh_rf <- tryCatch(
      cast_shap_fit(
        fit = fit,
        which = "rf",
        data = train_df,
        response = cfg$response,
        test_fraction = cfg$shap_test_fraction,
        seed = seed_i,
        fastshap_nsim = cfg$shap_fastshap_nsim,
        max_explain_rows = cfg$shap_max_explain_rows,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    .cast_shap_one_pair(sh_rf, "shap_rf", topn, fig_dir, fig_dpi)
  }

  torch_ok <- requireNamespace("torch", quietly = TRUE) &&
    tryCatch(torch::torch_is_installed(), error = function(e) FALSE)
  if (requireNamespace("fastshap", quietly = TRUE) && isTRUE(torch_ok) &&
      "cast" %in% names(fit$models) && !is.null(fit$models$cast$model)) {
    sh_cast <- tryCatch(
      cast_shap_fit(
        fit = fit,
        which = "cast",
        data = train_df,
        response = cfg$response,
        test_fraction = cfg$shap_test_fraction,
        seed = seed_i,
        fastshap_nsim = cfg$shap_fastshap_nsim,
        max_explain_rows = cfg$shap_max_explain_rows,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    .cast_shap_one_pair(sh_cast, "shap_cast", topn, fig_dir, fig_dpi)
  }

  .cast_shap_save_panel_2x3(sh_xgb, sh_rf, sh_cast, topn, fig_dir, fig_dpi)
  invisible(NULL)
}
