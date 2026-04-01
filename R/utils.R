# Internal Utility Functions ---------------------------------------------------

#' Validate Species Data Input
#'
#' Checks that input data has required columns and correct types.
#'
#' @param data A `data.frame` to validate.
#' @param required_cols Character vector of required column names.
#' @param call Caller environment for error reporting.
#'
#' @return `data` (invisibly), or aborts with informative error.
#'
#' @keywords internal
#' @noRd
validate_species_data <- function(data,
                                  required_cols = c("lon", "lat", "presence"),
                                  call = rlang::caller_env()) {
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "{.arg data} must be a data.frame, not {.obj_type_friendly {data}}.",
      call = call
    )
  }
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    cli::cli_abort(
      "{.arg data} is missing required column{?s}: {.val {missing}}.",
      call = call
    )
  }
  invisible(data)
}


#' Extract Environmental Variable Names
#'
#' Returns column names excluding coordinates, response, and non-numeric
#' columns.
#'
#' @param data A `data.frame`.
#' @param response Character. Response column name.
#' @param coords Character vector. Coordinate column names.
#'
#' @return Character vector of environmental variable names.
#'
#' @keywords internal
#' @noRd
get_env_vars <- function(data, response = "presence",
                         coords = c("lon", "lat"),
                         meta = NULL) {
  # Default metadata columns commonly found in species data
  default_meta <- c(
    "HID", "species", "sid", "family", "category", "fraction",
    "id", "ID", "site", "cell_id", "grid_id",
    # disdat-specific metadata columns
    "group", "spid", "siteid", "occ", "fold"
  )
  exclude <- unique(c(response, coords, default_meta, meta))
  nms <- setdiff(names(data), exclude)
  nms <- nms[vapply(data[nms], is.numeric, logical(1))]

  # Exclude near-zero-variance columns (e.g., occ = all 1s, constant flags)
  # These carry no ecological signal and often indicate metadata slippage
  nms[vapply(nms, function(v) {
    vals <- data[[v]]
    vals <- vals[!is.na(vals)]
    length(vals) > 0 && stats::var(vals) > 1e-10
  }, logical(1))]
}


#' Check if torch is Available
#'
#' @return Logical.
#' @keywords internal
#' @noRd
has_torch <- function() {
  requireNamespace("torch", quietly = TRUE) &&
    torch::torch_is_installed()
}


#' Normalize a Numeric Vector to [0, 1]
#'
#' @param x Numeric vector.
#' @return Numeric vector scaled to [0, 1].
#' @keywords internal
#' @noRd
normalize01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}


#' Evaluate Model Predictions (AUC + TSS)
#'
#' @param pred Numeric vector of predicted probabilities.
#' @param obs Integer/numeric vector of observed binary outcomes.
#' @return Named numeric vector with `auc` and `tss`.
#' @keywords internal
#' @noRd
evaluate_model <- function(pred, obs) {
  m <- evaluate_model_full(pred, obs)
  c(auc = m["auc"], tss = m["tss"])
}


#' Full Model Evaluation: AUC, TSS, CBI, SEDI, Kappa, PRAUC
#'
#' @param pred Numeric vector of predicted probabilities [0,1].
#' @param obs  Integer/numeric binary observed outcomes (0/1).
#' @return Named numeric vector with all metrics.
#' @keywords internal
#' @noRd
evaluate_model_full <- function(pred, obs) {
  pred <- pmin(pmax(as.numeric(pred), 1e-7), 1 - 1e-7)
  obs  <- as.integer(obs)

  # ── AUC (ROC) ───────────────────────────────────────────────────────────────
  auc_val <- tryCatch({
    as.numeric(pROC::auc(pROC::roc(obs, pred, quiet = TRUE)))
  }, error = function(e) NA_real_)

  # ── TSS (at threshold maximising sensitivity+specificity) ────────────────────
  tss_val <- tryCatch({
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
    coords  <- pROC::coords(roc_obj, "best",
                            ret = c("sensitivity", "specificity"))
    as.numeric(coords$sensitivity[1] + coords$specificity[1] - 1)
  }, error = function(e) NA_real_)

  # ── PRAUC (Precision-Recall AUC) ────────────────────────────────────────────
  # More informative than ROC-AUC under strong class imbalance (many backgrounds)
  prauc_val <- tryCatch({
    ord   <- order(pred, decreasing = TRUE)
    o_ord <- obs[ord]
    p_ord <- pred[ord]
    tp    <- cumsum(o_ord)
    fp    <- cumsum(1L - o_ord)
    fn    <- sum(obs) - tp
    prec  <- tp / (tp + fp)
    rec   <- tp / (tp + fn)
    # Remove NaN at zero-recall boundary
    keep  <- !is.nan(prec) & !is.nan(rec)
    prec  <- prec[keep]; rec <- rec[keep]
    # Trapezoidal integration
    if (length(rec) < 2) return(NA_real_)
    abs(sum(diff(rec) * (prec[-length(prec)] + prec[-1]) / 2))
  }, error = function(e) NA_real_)

  # ── CBI (Continuous Boyce Index) ─────────────────────────────────────────────
  # Measures habitat preference along the suitability gradient;
  # positive = model predictions align with presence distribution;
  # range (-1, 1), values > 0.5 indicate good predictive power.
  cbi_val <- tryCatch({
    compute_cbi(pred, obs)
  }, error = function(e) NA_real_)

  # ── SEDI (Symmetric Extremal Dependence Index) ───────────────────────────────
  # Log-scale measure robust to prevalence; -1 (worst) to +1 (perfect).
  # Recommended for rare-species SDMs (Ferro & Stephenson 2011 MWR).
  sedi_val <- tryCatch({
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
    coords  <- pROC::coords(roc_obj, "best",
                            ret = c("sensitivity", "specificity"))
    sens <- as.numeric(coords$sensitivity[1])
    spec <- as.numeric(coords$specificity[1])
    fp_rate <- 1 - spec
    fn_rate <- 1 - sens
    # Guard against log(0)
    if (fp_rate <= 0 || fp_rate >= 1 || fn_rate <= 0 || fn_rate >= 1) {
      return(NA_real_)
    }
    (log(fp_rate) - log(sens) - log(fn_rate) + log(spec)) /
      (log(fp_rate) + log(sens) + log(fn_rate) + log(spec))
  }, error = function(e) NA_real_)

  # ── Cohen's Kappa ────────────────────────────────────────────────────────────
  kappa_val <- tryCatch({
    roc_obj  <- pROC::roc(obs, pred, quiet = TRUE)
    coords   <- pROC::coords(roc_obj, "best",
                             ret = c("sensitivity", "specificity", "threshold"))
    thr      <- as.numeric(coords$threshold[1])
    pred_bin <- as.integer(pred >= thr)
    n  <- length(obs)
    tp <- sum(pred_bin == 1L & obs == 1L)
    tn <- sum(pred_bin == 0L & obs == 0L)
    fp <- sum(pred_bin == 1L & obs == 0L)
    fn <- sum(pred_bin == 0L & obs == 1L)
    po <- (tp + tn) / n
    pe <- ((tp + fp) * (tp + fn) + (tn + fn) * (tn + fp)) / (n^2)
    if (abs(1 - pe) < 1e-10) return(NA_real_)
    (po - pe) / (1 - pe)
  }, error = function(e) NA_real_)

  c(auc   = auc_val,
    tss   = tss_val,
    cbi   = cbi_val,
    sedi  = sedi_val,
    kappa = kappa_val,
    prauc = prauc_val)
}


#' Continuous Boyce Index (CBI)
#'
#' Moving-window Spearman correlation between predicted suitability and
#' observed presence frequency. Implementation follows Hirzel et al. (2006).
#'
#' @param pred Numeric predicted probabilities.
#' @param obs  Binary 0/1 observed.
#' @param n_bins Integer. Number of moving window bins. Default 101.
#' @return Scalar CBI in [-1, 1].
#' @keywords internal
#' @noRd
compute_cbi <- function(pred, obs, n_bins = 101L) {
  pres_pred <- pred[obs == 1L]
  if (length(pres_pred) < 5) return(NA_real_)

  bins  <- seq(0, 1, length.out = n_bins + 1L)
  width <- bins[2] - bins[1]
  mids  <- (bins[-1] + bins[-(n_bins + 1L)]) / 2

  # Expected: fraction of all predictions in bin (random expectation)
  exp_f <- vapply(seq_len(n_bins), function(i) {
    mean(pred >= bins[i] & pred < bins[i + 1L])
  }, numeric(1))

  # Predicted: fraction of presence predictions in bin
  pred_f <- vapply(seq_len(n_bins), function(i) {
    mean(pres_pred >= bins[i] & pres_pred < bins[i + 1L])
  }, numeric(1))

  # Remove zero-expectation bins
  keep <- exp_f > 0
  if (sum(keep) < 5) return(NA_real_)

  ratio <- pred_f[keep] / exp_f[keep]
  as.numeric(stats::cor(mids[keep], ratio, method = "spearman"))
}


#' Check and Require a Suggested Package
#'
#' @param pkg Package name.
#' @param reason Why the package is needed.
#' @param call Caller environment.
#' @keywords internal
#' @noRd
check_suggested <- function(pkg, reason = NULL, call = rlang::caller_env()) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package {.pkg %s} is required", pkg)
    if (!is.null(reason)) msg <- paste0(msg, " ", reason)
    cli::cli_abort(
      c(paste0(msg, "."), i = "Install it with {.code install.packages(\"{pkg}\")}"),
      call = call
    )
  }
}
