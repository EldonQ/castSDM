# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
# Saved predictions only: replotting never fits a model or computes a PDP.
case_figures <- function(comparison, out_dir, scenario) {
  dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)
  draw <- function() plot(comparison, scenario = scenario,
    layers = c("current", "joint", "joint_change", "interaction"))
  grDevices::png(file.path(out_dir, "figures", "scenario_comparison.png"),
    width = 2400, height = 2000, res = 250)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::pdf(file.path(out_dir, "figures", "scenario_comparison.pdf"), width = 9.6, height = 8)
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(out_dir)
}

case_validation_figure <- function(run, out_dir) {
  info <- run$cv$fold_info
  metrics <- run$cv$fold_metrics
  if (is.null(info) || is.null(metrics) || !nrow(metrics)) return(invisible(out_dir))
  if ("replicate" %in% names(metrics)) metrics <- metrics[metrics$replicate == 1, ]
  colors <- grDevices::hcl.colors(run$cv$k, "Dark 3")
  draw <- function() {
    old <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    on.exit(graphics::par(old), add = TRUE)
    graphics::plot(info$lon, info$lat, col = colors[info$fold],
      pch = ifelse(info$obs == 1, 16, 46), cex = ifelse(info$obs == 1, .45, .7),
      ylim = range(info$lat, finite = TRUE) + c(0, .15*diff(range(info$lat, finite = TRUE))),
      xlab = "X (input coordinate units)", ylab = "Y (input coordinate units)",
      main = "Training spatial folds: repeat 1")
    graphics::legend("topleft", legend = paste("F", seq_len(run$cv$k)),
      col = colors, pch = 16, bty = "n", cex = .75, ncol = run$cv$k)
    models <- unique(metrics$model)
    values <- vapply(models, function(m) {
      z <- metrics[metrics$model == m, ]
      z$auc[match(seq_len(run$cv$k), z$fold)]
    }, numeric(run$cv$k))
    graphics::matplot(seq_len(run$cv$k), values, type = "b", pch = seq_along(models),
      col = seq_along(models), lty = 1, ylim = c(0, 1), xaxt = "n",
      xlab = "Spatial fold", ylab = "ROC AUC", main = "OOF discrimination: descriptive")
    graphics::axis(1, at = seq_len(run$cv$k))
    graphics::abline(h = .5, lty = 3, col = "grey60")
    graphics::legend("bottomright", legend = models, col = seq_along(models),
      pch = seq_along(models), lty = 1, bty = "n", cex = .8)
  }
  for (format in c("png", "pdf")) {
    path <- file.path(out_dir, "figures", paste0("spatial_validation.", format))
    if (format == "png") grDevices::png(path, width = 2400, height = 1200, res = 250) else
      grDevices::pdf(path, width = 9.6, height = 4.8)
    tryCatch(draw(), finally = grDevices::dev.off())
  }
  invisible(out_dir)
}
