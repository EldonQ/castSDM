# Threshold selection and presence-background evaluation metrics ------------

#' Select a Binary Suitability Threshold
#'
#' Selects a threshold using training/calibration predictions. In validation
#' workflows, pass a threshold selected on the outer-training data and apply it
#' to the held-out fold; selecting the threshold on the test fold is optimistic.
#'
#' @param pred Numeric suitability scores.
#' @param obs Binary response (1 = presence, 0 = background/absence).
#' @param method Threshold rule: `"max_tss"`, `"equal_sens_spec"`,
#'   `"max_jaccard"`, `"max_kappa"`, `"no_omission"`, `"p10_presence"`,
#'   `"equal_prevalence"`, or `"sensitivity"`.
#' @param sensitivity Target sensitivity for `method = "sensitivity"`.
#' @param omission Allowed omission fraction for `method = "p10_presence"`.
#'
#' @return A scalar numeric threshold.
#' @export
cast_threshold <- function(pred, obs,
                           method = c("max_tss", "equal_sens_spec",
                                      "max_jaccard", "max_kappa",
                                      "no_omission", "p10_presence",
                                      "equal_prevalence", "sensitivity"),
                           sensitivity = 0.9, omission = 0.1) {
  method <- match.arg(method)
  if (!is.numeric(pred) || !is.numeric(obs) || length(pred) != length(obs)) {
    cli::cli_abort("{.arg pred} and {.arg obs} must be numeric vectors of equal length.")
  }
  ok <- is.finite(pred) & !is.na(obs)
  pred <- as.numeric(pred[ok])
  obs <- as.integer(obs[ok])
  if (!length(pred) || !all(obs %in% c(0L, 1L)) ||
      !all(c(0L, 1L) %in% unique(obs))) {
    cli::cli_abort("Threshold selection needs finite predictions and both response classes.")
  }
  if (!is.numeric(sensitivity) || length(sensitivity) != 1L ||
      !is.finite(sensitivity) || sensitivity <= 0 || sensitivity > 1) {
    cli::cli_abort("{.arg sensitivity} must be in (0, 1].")
  }
  if (!is.numeric(omission) || length(omission) != 1L ||
      !is.finite(omission) || omission < 0 || omission >= 1) {
    cli::cli_abort("{.arg omission} must be in [0, 1).")
  }

  n1 <- sum(obs == 1L)
  n0 <- sum(obs == 0L)
  ord <- order(pred, decreasing = TRUE)
  p <- pred[ord]
  y <- obs[ord]
  ends <- c(which(diff(p) != 0), length(p))
  tp <- cumsum(y)[ends]
  fp <- ends - tp
  fn <- n1 - tp
  tn <- n0 - fp
  thr <- p[ends]
  # Include both extreme classifiers (all-negative and all-positive).
  eps <- .Machine$double.eps * max(1, abs(range(pred)))
  thr <- c(max(pred) + eps, thr, min(pred) - eps)
  tp <- c(0, tp, n1); fp <- c(0, fp, n0)
  fn <- c(n1, fn, 0); tn <- c(n0, tn, 0)
  tpr <- tp / n1
  tnr <- tn / n0
  tss <- tpr + tnr - 1
  jaccard <- tp / pmax(1, tp + fp + fn)
  po <- (tp + tn) / (n1 + n0)
  pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / (n1 + n0)^2
  kappa <- (po - pe) / pmax(.Machine$double.eps, 1 - pe)

  if (method == "max_tss") return(unname(thr[which.max(tss)]))
  if (method == "equal_sens_spec") return(unname(thr[which.min(abs(tpr - tnr))]))
  if (method == "max_jaccard") return(unname(thr[which.max(jaccard)]))
  if (method == "max_kappa") return(unname(thr[which.max(kappa)]))
  if (method == "equal_prevalence") {
    predicted_prevalence <- (tp + fp) / (n1 + n0)
    return(unname(thr[which.min(abs(predicted_prevalence - n1 / (n1 + n0)))]))
  }
  if (method == "sensitivity") {
    eligible <- which(tpr >= sensitivity)
    return(unname(thr[eligible[which.max(thr[eligible])]]))
  }
  if (method == "no_omission") return(unname(min(pred[obs == 1L])))
  unname(stats::quantile(pred[obs == 1L], probs = omission, names = FALSE, type = 1))
}

#' Precision-Recall Area for Imbalanced Binary Data
#'
#' Trapezoidal area under the empirical precision-recall curve. Unlike ROC-AUC,
#' this metric makes the sampled prevalence explicit: always report the
#' presence:background ratio and the prevalence baseline alongside it.
#'
#' @param obs Binary response vector.
#' @param pred Numeric suitability scores.
#' @return Scalar PR-AUC in [0, 1], or `NA` if a class is absent.
#' @keywords internal
#' @noRd
compute_pr_auc <- function(obs, pred) {
  ok <- !is.na(obs) & is.finite(pred)
  obs <- as.integer(obs[ok]); pred <- as.numeric(pred[ok])
  n_pos <- sum(obs == 1L)
  n_neg <- sum(obs == 0L)
  if (!n_pos || !n_neg || !all(obs %in% c(0L, 1L))) return(NA_real_)
  ord <- order(pred, decreasing = TRUE)
  p <- pred[ord]; y <- obs[ord]
  ends <- c(which(diff(p) != 0), length(p))
  tp <- cumsum(y)[ends]
  fp <- ends - tp
  recall <- c(0, tp / n_pos)
  precision <- c(1, tp / pmax(1, tp + fp))
  sum(diff(recall) * (head(precision, -1L) + tail(precision, -1L)) / 2)
}

#' Symmetric Extremal Dependence Index
#'
#' @param obs Binary response vector.
#' @param pred Numeric suitability scores.
#' @param threshold Numeric threshold used to form the contingency table.
#' @return Scalar SEDI in [-1, 1], or `NA` for degenerate inputs.
#' @keywords internal
#' @noRd
compute_sedi <- function(obs, pred, threshold) {
  ok <- !is.na(obs) & is.finite(pred)
  obs <- as.integer(obs[ok]); pred <- as.numeric(pred[ok])
  if (!length(obs) || !all(obs %in% c(0L, 1L)) ||
      !all(c(0L, 1L) %in% unique(obs)) || !is.finite(threshold)) return(NA_real_)
  pos <- obs == 1L; hit <- pred >= threshold
  h <- mean(hit[pos])
  f <- mean(hit[!pos])
  eps <- 1 / (2 * length(obs))
  h <- min(1 - eps, max(eps, h))
  f <- min(1 - eps, max(eps, f))
  num <- log(f) - log(h) - log1p(-f) + log1p(-h)
  den <- log(f) + log(h) + log1p(-f) + log1p(-h)
  if (!is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
  num / den
}

#' Moving-Window Continuous Boyce Index
#'
#' @param pred Numeric predictions for presences and background/absence rows.
#' @param obs Binary response vector.
#' @param window_fraction Fraction of the prediction range used as moving
#'   window width. Default `0.1`.
#' @param n_windows Number of window positions. Default `100`.
#' @return Scalar Boyce index, or `NA` if too few usable windows.
#' @keywords internal
#' @noRd
compute_boyce <- function(pred, obs, window_fraction = 0.1, n_windows = 100L) {
  ok <- is.finite(pred) & !is.na(obs)
  pred <- as.numeric(pred[ok]); obs <- as.integer(obs[ok])
  if (!length(pred) || !all(obs %in% c(0L, 1L)) ||
      sum(obs == 1L) < 5L || sum(obs == 0L) < 2L) return(NA_real_)
  if (!is.numeric(window_fraction) || length(window_fraction) != 1L ||
      !is.finite(window_fraction) || window_fraction <= 0 || window_fraction >= 1) {
    cli::cli_abort("{.arg window_fraction} must be in (0, 1).")
  }
  n_windows <- as.integer(n_windows)
  if (is.na(n_windows) || n_windows < 2L) cli::cli_abort("{.arg n_windows} must be at least 2.")
  lo <- min(pred); hi <- max(pred); width <- (hi - lo) * window_fraction
  if (!is.finite(width) || width <= 0) return(NA_real_)
  starts <- seq(lo, hi - width, length.out = n_windows)
  pres <- pred[obs == 1L]; bg <- pred[obs == 0L]
  pe <- vapply(starts, function(a) {
    b <- a + width
    fp <- mean(pres >= a & pres <= b)
    fb <- mean(bg >= a & bg <= b)
    if (fp <= 0 || fb <= 0) NA_real_ else fp / fb
  }, numeric(1))
  centers <- starts + width / 2
  keep <- is.finite(pe) & pe > 0
  if (sum(keep) < 2L) return(NA_real_)
  pe <- pe[keep]; centers <- centers[keep]
  if (length(pe) > 1L) {
    distinct <- c(TRUE, diff(pe) != 0)
    pe <- pe[distinct]; centers <- centers[distinct]
  }
  if (length(pe) < 2L) return(NA_real_)
  suppressWarnings(as.numeric(stats::cor(centers, pe, method = "spearman")))
}
