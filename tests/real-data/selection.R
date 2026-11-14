# Opt-in integration check using an existing actual castSDM case only.
# Usage: Rscript tests/real-data/selection.R <case_directory> [maximum_rows]
# No synthetic measurements, labels, pseudo-effects or national projections.
args <- commandArgs(TRUE)
case <- if (length(args)) args[1L] else Sys.getenv("CAST_REAL_CASE", "")
if (!nzchar(case)) stop("Supply an existing real case directory.")
suppressPackageStartupMessages(library(castSDM))
p <- readRDS(file.path(case, "training", "stages", "01_prepare.rds"))
r <- readRDS(file.path(case, "training", "result.rds"))
stopifnot(identical(p$train$lon, r$cv$fold_info$lon),
          identical(p$train$lat, r$cv$fold_info$lat))
rows <- seq_len(nrow(p$train))
if (length(args) > 1L) rows <- unique(as.integer(round(seq(1, nrow(p$train),
  length.out = min(nrow(p$train), as.integer(args[2L]))))))
x <- p$train[rows, , drop = FALSE]
timing <- system.time(s <- cast_select(x, seed = 42))
stopifnot(identical(s$diagnostics$algorithm, "vif"),
  all(s$selected %in% s$scores$variable[s$scores$eligible]),
  !any(s$scores$role[s$scores$selected] %in% c("identifier", "nominal", "cyclic")),
  all(s$scores$final_vif[s$scores$selected] <= 10))
# Re-express actual BIO04 from degC x100 to millidegC x100 in memory.
# Only the unit representation changes; membership and risk must be stable.
if ("bio04" %in% names(x)) {
  converted <- x
  converted$bio04 <- converted$bio04 * 1000
  u <- cast_select(converted, seed = 42, verbose = FALSE)
  stopifnot(identical(s$selected, u$selected),
    isTRUE(all.equal(s$scores$final_vif, u$scores$final_vif, tolerance = 1e-6)))
}
# Independent actual validation rows are used to check engine/prediction
# contracts only, not to select a policy or claim predictive superiority.
if (length(s$selected)) {
  fit <- cast_fit(x, screen = s, models = c("rf", "maxent", "gam"), seed = 42,
    num_threads = 1L, maxent_classes = "lq", verbose = FALSE)
  hold <- p$test[seq_len(min(500L, nrow(p$test))), , drop = FALSE]
  pred <- cast_predict(fit, hold, clamp = TRUE, extrapolation = FALSE)
  stopifnot(inherits(pred, "cast_predict"),
    nrow(pred$predictions) == nrow(hold),
    all(c("HSS_rf", "HSS_maxent", "HSS_gam") %in% names(pred$predictions)),
    all(is.finite(as.matrix(pred$predictions[, c("HSS_rf", "HSS_maxent", "HSS_gam")]))),
    all(as.matrix(pred$predictions[, c("HSS_rf", "HSS_maxent", "HSS_gam")]) >= 0),
    all(as.matrix(pred$predictions[, c("HSS_rf", "HSS_maxent", "HSS_gam")]) <= 1))
  print(fit)
  print(head(pred$predictions))
}
# Exercise both public plotting interfaces against actual selection results.
invisible(ggplot2::ggplot_build(plot(s)))
invisible(ggplot2::ggplot_build(plot(cast_importance(s))))
print(s); print(timing)
cat("Actual-data selection integration checks passed. No full workflow run.\n")
