# AI-assisted maintenance (2026-11-07); scope and verification: inst/AI_ASSISTANCE.md.
#' Diagnose a Fitted Workflow without Refitting
#'
#' Reports independent hold-out performance, the frozen ensemble, spatial-fold
#' composition, selection stability, and training OOF score-scale comparisons.
#' [cast_run()] saves these diagnostics automatically, including when no map
#' was requested. No policy is selected from hold-out labels.
#'
#' @param object A `cast_result` from [cast_run()].
#' @param comparison Optional `cast_comparison` from [cast_compare()].
#' @param test_data Optional original hold-out frame for older results that
#'   did not save evaluation predictions. Never supply training/refit data.
#' @param output_dir Optional directory for diagnostic CSVs and an RDS object.
#' @return A `cast_diagnostics` list of tables and interpretation notes.
#' @details
#' Score comparisons use first-repeat training OOF predictions with the
#' declared ensemble weights held fixed. Native means the same fitted engine
#' before its logit shift; it is not a fit with different class weights.
#' These comparisons are descriptive, as are CV metrics at a pooled OOF
#' threshold. A new policy needs new validation; this function never chooses it.
#'
#' Brier scores, log loss and reliability bins concern the sampled binary
#' response. With presence/background data they do not measure calibration of
#' ecological occurrence probabilities. Background points are not known absences.
#' Average precision (AP) uses the non-interpolated precision-recall step
#' integral; a constant score has AP equal to sampled prevalence. Historical
#' `pr_auc` uses trapezoidal interpolation and is retained with that definition,
#' except in the constant-baseline row where it is omitted.
#' Support is a marginal training-quantile diagnostic, distinct from multivariate
#' AOA. Model sign agreement is descriptive, not statistical confidence.
#' @examples
#' # After a run: diagnostics <- cast_diagnose(result)
#' @export
cast_diagnose <- function(object, comparison = NULL, test_data = NULL,
                          output_dir = NULL) {
  if (!inherits(object, "cast_result")) cli::cli_abort("object must be a cast_result.")
  if (!is.null(comparison) && !inherits(comparison, "cast_comparison"))
    cli::cli_abort("comparison must be a cast_comparison or NULL.")
  cv <- object$cv
  ev <- object$eval
  if (is.null(ev$predictions) && !is.null(test_data))
    ev <- cast_evaluate(object$fit, test_data, thresholds = cv$thresholds,
      threshold_method = object$protocol$threshold_method %||% "max_tss",
      clamp = object$eval$clamp %||% FALSE)
  policy <- .cast_diagnostic_policy(object)
  holdout <- list()
  p <- ev$predictions
  train_y <- object$fit$scaling$response
  prevalence <- if (length(train_y)) mean(train_y) else NA_real_
  if (is.data.frame(p) && "obs" %in% names(p)) {
    for (m in names(object$fit$models)) {
      col <- paste0("HSS_", m)
      if (col %in% names(p)) holdout[[m]] <- .cast_diagnostic_metrics(
        p[[col]], p$obs, unname(ev$thresholds[m]) %||% NA_real_, m, "declared",
        ev$threshold_source %||% "unavailable",
        reference = if (is.data.frame(cv$oof)) list(pred = cv$oof[[col]], obs = cv$oof$obs))
    }
    if (!is.null(policy)) {
      pred <- .cast_diagnostic_combine(p, policy)
      threshold <- .cast_diagnostic_threshold(cv$oof, cv, policy)
      holdout$ensemble <- .cast_diagnostic_metrics(pred, p$obs, threshold,
        "ensemble", "declared", if (is.finite(threshold)) "out_of_fold_cv" else "unavailable",
        reference = list(pred = .cast_diagnostic_combine(cv$oof, policy), obs = cv$oof$obs))
    }
    if (is.finite(prevalence)) holdout$constant_training_prevalence <-
      .cast_diagnostic_metrics(rep(prevalence, nrow(p)), p$obs, NA_real_,
        "constant_training_prevalence", "baseline", "not_applicable")
  }
  score_views <- list(); reliability <- list()
  for (view in c("declared", "native_same_fit")) {
    oof <- if (view == "declared") cv$oof else cv$native_oof
    if (!is.data.frame(oof) || !"obs" %in% names(oof)) next
    for (m in names(object$fit$models)) {
      col <- paste0("HSS_", m)
      if (!col %in% names(oof)) next
      threshold <- tryCatch(cast_threshold(oof[[col]], oof$obs,
        method = cv$threshold_method %||% "max_tss"), error = function(e) NA_real_)
      score_views[[length(score_views)+1L]] <- .cast_diagnostic_metrics(
        oof[[col]], oof$obs, threshold, m, view, "pooled_training_oof_descriptive")
      reliability[[length(reliability)+1L]] <- .cast_diagnostic_bins(oof[[col]], oof$obs, m, view)
    }
    if (!is.null(policy)) {
      view_policy <- policy
      if (view == "native_same_fit" && policy$method == "committee")
        view_policy$thresholds <- vapply(names(policy$weights), function(m)
          tryCatch(cast_threshold(oof[[paste0("HSS_", m)]], oof$obs,
            method = cv$threshold_method %||% "max_tss"), error = function(e) NA_real_), numeric(1))
      pred <- .cast_diagnostic_combine(oof, view_policy)
      threshold <- .cast_diagnostic_threshold(oof, cv, view_policy)
      score_views[[length(score_views)+1L]] <- .cast_diagnostic_metrics(pred, oof$obs,
        threshold, "ensemble", view, "pooled_training_oof_descriptive")
      reliability[[length(reliability)+1L]] <- .cast_diagnostic_bins(pred, oof$obs, "ensemble", view)
    }
  }
  bind <- function(x) if (length(x)) do.call(rbind, x) else NULL
  folds <- .cast_diagnostic_folds(cv)
  definitions <- data.frame(metric = c("average_precision", "pr_auc", "brier/logloss", "omission_at_threshold",
    "omission_5/10 (legacy)", "mpa (legacy)", "native_same_fit", "common_support", "model_sd"),
    definition = c("Non-interpolated precision-recall step integral (AP); constant-score baseline equals sampled positive fraction.",
    "Trapezoidal precision-recall area; not AP. A constant score gives (1 + prevalence)/2; omitted for the constant baseline row.",
    "Loss for sampled presence/background labels; not ecological probability calibration.",
    "Fraction of evaluation presences below the frozen training OOF threshold.",
    "Presence omission below the top 5/10 percent cutoff of all evaluation scores; not training-presence omission.",
    "Fraction of evaluation sample above its presence-score 10th percentile; not geographic area.",
    "Same fitted engine before score shift; fitting weights and variables are unchanged.",
    "Finite and within marginal training bounds in current and ALL requested futures; not multivariate AOA.",
    "Sample SD of paired model contrasts; not SE, confidence interval or significance."))
  notes <- data.frame(code = c("sampled_response", "holdout_policy", "oof_scope", "support_scope", "score_shift"),
    detail = c("Background labels are not verified absences; sampled prevalence changes probability losses and PR-AUC.",
      "Hold-out labels evaluate the declared policy only. Do not select a new policy on these labels.",
      "OOF score views, fold locations and AOA use repeat 1; fold_metrics and selection frequency cover all scheduled repeats.",
      "Scenario-supported and all-scenario common-supported summaries have different denominators; quantify both.",
      "Declared logit shifts standardize scores; they do not establish probability calibration. Native same-fit OOF views are descriptive, not a new fitted policy."))
  if (is.null(p)) notes <- rbind(notes, data.frame(code = "missing_holdout_predictions",
    detail = "Older evaluation lacks predictions. Supply the original test_data to add ensemble diagnostics without refitting."))
  if (is.data.frame(cv$oof) && is.null(cv$native_oof))
    notes <- rbind(notes, data.frame(code = "missing_native_oof",
      detail = "This older result did not save native OOF scores or fold offsets; they are not reconstructed from the final fit."))
  if (!identical(object$protocol$clamp %||% FALSE, ev$clamp %||% FALSE))
    notes <- rbind(notes, data.frame(code = "legacy_clamp_mismatch",
      detail = "This older result evaluated unclamped predictors while its map protocol clamps. Historical evaluation is preserved; a new declared run is needed for matching validation."))
  if (nrow(folds %||% data.frame()) && any(folds$n_presence == 0 | folds$n_background == 0))
    notes <- rbind(notes, data.frame(code = "single_class_fold", detail = "At least one spatial fold lacks a response class; inspect coverage before using thresholds."))
  d <- structure(list(holdout = bind(holdout), score_views_oof = bind(score_views),
    reliability_oof = bind(reliability), folds = folds, fold_metrics = cv$fold_metrics,
    selection_frequency = cv$selection_freq, score_offsets = cv$offsets,
    ensemble_policy = if (!is.null(policy)) data.frame(model = names(policy$weights),
      weight = as.numeric(policy$weights), method = policy$method),
    comparison_support = comparison$summary_common_support,
    comparison_binary_common_support = comparison$binary_common_support,
    model_agreement = comparison$agreement, threshold_sensitivity = comparison$threshold_sensitivity,
    metric_definitions = definitions, notes = notes), class = "cast_diagnostics")
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (nm in names(d)) {
      path <- file.path(output_dir, paste0(nm, ".csv"))
      if (is.data.frame(d[[nm]])) utils::write.csv(d[[nm]], path, row.names = FALSE) else
        if (file.exists(path)) unlink(path)
    }
    .cast_atomic_rds(d, file.path(output_dir, "diagnostics.rds"))
  }
  d
}

.cast_diagnostic_policy <- function(object) {
  cv <- object$cv
  if (is.null(cv$metrics) || !nrow(cv$metrics)) return(NULL)
  p <- object$protocol
  models <- intersect(names(object$fit$models), cv$metrics$model)
  if (!length(models)) return(NULL)
  method <- p$ensemble_method %||% "weighted"
  if (!is.null(object$ensemble$weights)) {
    weights <- object$ensemble$weights
    weights <- weights[is.finite(weights) & weights > 0]
    if (length(weights)) return(list(weights = weights/sum(weights), method = method,
                                     thresholds = cv$thresholds))
  }
  scores <- .cast_ensemble_scores(cv$metrics, models,
    min_metric_folds = p$ensemble_min_metric_folds %||% 2L)
  pw <- .cast_power_scores(scores, p$ensemble_decay %||% 1, method,
    p$ensemble_min_score %||% .5)
  weights <- .cast_ensemble_weights(pw$scores, method, pw$min_score,
    p$ensemble_fallback %||% "best")
  weights <- weights[is.finite(weights) & weights > 0]
  if (!length(weights)) return(NULL)
  list(weights = weights/sum(weights), method = method, thresholds = cv$thresholds)
}

.cast_diagnostic_combine <- function(predictions, policy) {
  cols <- paste0("HSS_", names(policy$weights))
  if (!is.data.frame(predictions) || !all(cols %in% names(predictions)))
    return(rep(NA_real_, nrow(predictions) %||% 0L))
  x <- as.matrix(predictions[, cols, drop = FALSE])
  ok <- .cast_finite_rows(x)
  out <- rep(NA_real_, nrow(x))
  if (policy$method == "committee") {
    thr <- policy$thresholds[names(policy$weights)]
    if (length(thr) == ncol(x) && all(is.finite(thr)) && ncol(x) >= 2L)
      out[ok] <- rowMeans(sweep(x[ok, , drop = FALSE], 2L, thr, `>=`))
  } else if (policy$method == "median") {
    out[ok] <- .ensemble_rowmedian(x[ok, , drop = FALSE])
  } else out[ok] <- as.numeric(x[ok, , drop = FALSE] %*% policy$weights)
  out
}

.cast_diagnostic_threshold <- function(oof, cv, policy) {
  if (!is.data.frame(oof) || !all(.cast_oof_fold_coverage(cv, names(policy$weights))$ok)) return(NA_real_)
  pred <- .cast_diagnostic_combine(oof, policy)
  ok <- is.finite(pred) & is.finite(oof$obs)
  if (sum(ok) < 10L || length(unique(oof$obs[ok])) < 2L) return(NA_real_)
  if (policy$method == "committee") return(.5)
  tryCatch(cast_threshold(pred[ok], oof$obs[ok], method = cv$threshold_method %||% "max_tss"),
    error = function(e) NA_real_)
}

.cast_diagnostic_metrics <- function(pred, obs, threshold, model, view, threshold_source,
                                      reference = NULL) {
  n_total <- length(obs)
  ok <- is.finite(pred) & is.finite(obs) & obs %in% c(0, 1)
  pred <- pred[ok]; obs <- obs[ok]
  # Passing a fixed finite threshold prevents evaluate_model_full from
  # optimizing on these labels. Binary fields are removed if unavailable.
  binary_ok <- length(threshold) == 1L && is.finite(threshold)
  ev <- evaluate_model_full(pred, obs, threshold = if (binary_ok) threshold else .5)
  if (!binary_ok) ev[intersect(names(ev), c("tss", "sedi", "kappa", "tss_threshold"))] <- NA_real_
  train_presence <- reference$pred[is.finite(reference$pred) & reference$obs == 1]
  q <- if (length(train_presence)) unname(stats::quantile(train_presence, c(.05, .1),
    type = 1, names = FALSE)) else rep(NA_real_, 2L)
  data.frame(model = model, view = view, n = length(obs), n_total = n_total,
    n_missing_score = n_total-length(obs), n_presence = sum(obs == 1),
    sampled_prevalence = if (length(obs)) mean(obs) else NA_real_,
    mean_score = if (length(pred)) mean(pred) else NA_real_,
    auc = unname(ev["auc"]), average_precision = compute_average_precision(obs, pred),
    pr_auc = if (identical(view, "baseline")) NA_real_ else unname(ev["pr_auc"]), tss = unname(ev["tss"]),
    brier = unname(ev["brier"]), logloss = unname(ev["logloss"]),
    threshold = if (binary_ok) threshold else NA_real_, threshold_source = threshold_source,
    omission_at_threshold = if (binary_ok && any(obs == 1)) mean(pred[obs == 1] < threshold) else NA_real_,
    training_oof_presence05_threshold = q[1L], training_oof_presence10_threshold = q[2L],
    omission_at_presence05 = if (is.finite(q[1L]) && any(obs == 1)) mean(pred[obs == 1] < q[1L]) else NA_real_,
    omission_at_presence10 = if (is.finite(q[2L]) && any(obs == 1)) mean(pred[obs == 1] < q[2L]) else NA_real_)
}

.cast_diagnostic_bins <- function(pred, obs, model, view) {
  ok <- is.finite(pred) & is.finite(obs)
  bin <- cut(pred[ok], seq(0, 1, .1), include.lowest = TRUE, labels = FALSE)
  do.call(rbind, lapply(seq_len(10L), function(i) {
    idx <- which(bin == i)
    data.frame(model = model, view = view, bin = i, n = length(idx),
      mean_score = if (length(idx)) mean(pred[ok][idx]) else NA_real_,
      sampled_positive_fraction = if (length(idx)) mean(obs[ok][idx]) else NA_real_)
  }))
}

.cast_diagnostic_folds <- function(cv) {
  info <- cv$fold_info
  if (is.null(info) && length(cv$folds) && nrow(cv$oof) == length(cv$folds))
    info <- data.frame(fold = cv$folds, obs = cv$oof$obs, lon = NA_real_, lat = NA_real_)
  if (!is.data.frame(info)) return(NULL)
  di <- cv$aoa$di_oof
  thr <- cv$aoa$threshold %||% NA_real_
  do.call(rbind, lapply(seq_len(cv$k), function(i) {
    idx <- which(info$fold == i)
    finite_range <- function(x, fun) if (any(is.finite(x))) fun(x[is.finite(x)]) else NA_real_
    vals <- if (length(di) == nrow(info)) di[idx] else numeric()
    vals <- vals[is.finite(vals)]
    data.frame(fold = i, n = length(idx), n_presence = sum(info$obs[idx] == 1),
      n_background = sum(info$obs[idx] == 0),
      xmin = finite_range(info$lon[idx], min), xmax = finite_range(info$lon[idx], max),
      ymin = finite_range(info$lat[idx], min), ymax = finite_range(info$lat[idx], max),
      n_aoa = length(vals), share_outside_aoa = if (length(vals) && is.finite(thr)) mean(vals > thr) else NA_real_)
  }))
}
