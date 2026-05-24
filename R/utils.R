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
                                  call = parent.frame()) {
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
#' Returns column names excluding coordinates, response, metadata, and
#' non-numeric columns. Useful for identifying which columns in your
#' data are environmental predictors.
#'
#' @param data A `data.frame`.
#' @param response Character. Response column name. Default `"presence"`.
#' @param coords Character vector. Coordinate column names. Default
#'   `c("lon", "lat")`.
#' @param meta Optional character vector of additional metadata column
#'   names to exclude.
#'
#' @return Character vector of environmental variable names.
#'
#' @export
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


#' Compute AUC via Wilcoxon-Mann-Whitney Estimator
#'
#' A lightweight AUC computation that does not require pROC. Used
#' internally by [cast_esm()] for bivariate model weighting.
#'
#' @param y Binary integer vector (0/1).
#' @param pred Numeric vector of predicted probabilities.
#' @return Scalar AUC in \[0, 1\].
#' @keywords internal
#' @noRd
compute_auc <- function(y, pred) {
  ok <- !is.na(y) & !is.na(pred)
  y <- y[ok]; pred <- pred[ok]
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(pred)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
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


#' Compute Out-degree and In-degree from DAG Edges
#'
#' @param edges A `data.frame` with columns `from` and `to`.
#' @param variables Character vector of variable names.
#' @return A `data.frame` with columns `variable`, `out_degree`, `in_degree`.
#' @keywords internal
#' @noRd
compute_edge_degrees <- function(edges, variables) {
  if (nrow(edges) > 0) {
    out_agg <- stats::aggregate(to ~ from, data = edges, FUN = length)
    names(out_agg) <- c("variable", "out_degree")
    in_agg <- stats::aggregate(from ~ to, data = edges, FUN = length)
    names(in_agg) <- c("variable", "in_degree")
  } else {
    out_agg <- data.frame(variable = character(0), out_degree = integer(0),
                          stringsAsFactors = FALSE)
    in_agg <- data.frame(variable = character(0), in_degree = integer(0),
                         stringsAsFactors = FALSE)
  }
  df <- data.frame(variable = variables, stringsAsFactors = FALSE)
  df <- merge(df, out_agg, by = "variable", all.x = TRUE)
  df <- merge(df, in_agg, by = "variable", all.x = TRUE)
  df$out_degree[is.na(df$out_degree)] <- 0
  df$in_degree[is.na(df$in_degree)] <- 0
  df
}


#' Full Model Evaluation: AUC, TSS, CBI
#'
#' @param pred Numeric vector of predicted probabilities [0,1].
#' @param obs  Integer/numeric binary observed outcomes (0/1).
#' @return Named numeric vector with AUC, TSS, CBI.
#' @keywords internal
#' @noRd
evaluate_model_full <- function(pred, obs) {
  pred <- pmin(pmax(as.numeric(pred), 1e-7), 1 - 1e-7)
  obs  <- as.integer(obs)

  # -- AUC (ROC) -------------------------------------------------------------
  auc_val <- tryCatch({
    as.numeric(pROC::auc(pROC::roc(obs, pred, quiet = TRUE)))
  }, error = function(e) NA_real_)

  # -- TSS (at threshold maximising sensitivity+specificity) ------------------
  tss_val <- tryCatch({
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
    coords  <- pROC::coords(roc_obj, "best",
                            ret = c("sensitivity", "specificity"))
    as.numeric(coords$sensitivity[1] + coords$specificity[1] - 1)
  }, error = function(e) NA_real_)

  # -- CBI (Continuous Boyce Index) -------------------------------------------
  cbi_val <- tryCatch({
    compute_cbi(pred, obs)
  }, error = function(e) NA_real_)

  c(auc = auc_val,
    tss = tss_val,
    cbi = cbi_val)
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


#' Null-Coalescing Operator
#'
#' Returns `b` if `a` is `NULL`, else `a`. Defined locally to support R < 4.4.
#'
#' @param a,b Values to compare.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Check and Require a Suggested Package
#'
#' @param pkg Package name.
#' @param reason Why the package is needed.
#' @param call Caller environment.
#' @keywords internal
#' @noRd
check_suggested <- function(pkg, reason = NULL, call = parent.frame()) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package {.pkg %s} is required", pkg)
    if (!is.null(reason)) msg <- paste0(msg, " ", reason)
    cli::cli_abort(
      c(paste0(msg, "."), i = "Install it with {.code install.packages(\"{pkg}\")}"),
      call = call
    )
  }
}
