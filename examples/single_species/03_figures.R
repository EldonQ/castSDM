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
