# Internal Utility Functions ---------------------------------------------------

#' Validate Species Data Input
#'
#' Checks that input data has required columns and correct types.
#'
#' @param data A `data.frame` to validate.
#' @param required_cols Character vector of required column names.
#' @param response Response column name; when present in `data` its values are
#'   checked to be binary (0/1 or logical). Rows with `NA` trigger a warning.
#'   Default `"presence"`.
#' @param call Caller environment for error reporting.
#'
#' @return `data` (invisibly), or aborts with informative error.
#'
#' @keywords internal
#' @noRd
validate_species_data <- function(data,
                                  required_cols = c("lon", "lat", "presence"),
                                  response = "presence",
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
  if (length(response) && response %in% names(data)) {
    .cast_check_response(data[[response]], response, call = call)
  }
  invisible(data)
}


#' Validate a Binary Response Vector
#'
#' Aborts unless every non-missing value is 0/1 (logical is accepted: it
#' coerces cleanly to 0/1). Missing rows only warn - dropping them is the
#' caller's decision.
#'
#' @param y Response vector.
#' @param response Column name used in messages.
#' @param call Caller environment for error reporting.
#' @keywords internal
#' @noRd
.cast_check_response <- function(y, response = "presence",
                                 call = parent.frame()) {
  n_na <- sum(is.na(y))
  if (n_na > 0L) {
    cli::cli_warn(
      "Response {.val {response}} has {n_na} missing row{?s}; model fitting will fail or drop them downstream.",
      call = call
    )
  }
  if (is.logical(y)) return(invisible(y))
  vals <- unique(y[!is.na(y)])
  if (!is.numeric(y) || !all(vals %in% c(0, 1))) {
    cli::cli_abort(c(
      "Response {.val {response}} must be binary 0/1; found value{?s}: {.val {head(vals, 5)}}.",
      i = "Recode presence/background as 1/0 before modelling."
    ), call = call)
  }
  invisible(y)
}


#' Reject Non-Numeric Predictor Columns
#'
#' Shared guard for the fit / evaluate / CV / predict pathways: silently
#' coercing a factor with `as.numeric()` would model its level codes, not its
#' values, so non-numeric predictors are rejected everywhere.
#'
#' @param x A `data.frame` of predictor columns.
#' @param arg Argument name used in messages.
#' @param call Caller environment for error reporting.
#' @keywords internal
#' @noRd
.cast_check_numeric_predictors <- function(x, arg = "data",
                                           call = parent.frame()) {
  bad <- names(x)[!vapply(x, is.numeric, logical(1))]
  if (length(bad)) {
    cli::cli_abort(c(
      "Non-numeric predictor{?s}: {.val {bad}}.",
      i = "{.arg {arg}} must hold numeric predictors only; convert factors/characters to numeric before modelling."
    ), call = call)
  }
  invisible(x)
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
#' @details
#' Columns named `HID`, `species`, `sid`, `family`, `category`, `fraction`,
#' `id`, `ID`, `site`, `cell_id`, `grid_id`, `group`, `spid`, `siteid`,
#' `occ`, or `fold` are treated as metadata and excluded from the returned
#' predictor set, together with any names in `meta`. If one of your real
#' predictors carries one of these names, rename the column or pass an
#' explicit predictor vector to the modelling functions instead of relying
#' on auto-detection.
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
  # These carry no ecological signal and often indicate metadata slippage.
  # `length(vals) <= 1L` guards var() of a single non-NA value (which returns
  # NA); a near-constant column is dropped instead of injecting an NA name.
  nms[vapply(nms, function(v) {
    vals <- data[[v]]
    vals <- vals[!is.na(vals)]
    if (length(vals) <= 1L) return(FALSE)
    vv <- stats::var(vals)
    is.finite(vv) && vv > 1e-10
  }, logical(1))]
}


#' Compute AUC via Wilcoxon-Mann-Whitney Estimator
#'
#' A lightweight AUC computation that does not require pROC. Used internally
#' where a Wilcoxon-Mann-Whitney AUC estimate is needed without the pROC
#' dependency.
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


#' Full Model Evaluation: AUC, TSS, CBI
#'
#' @param pred Numeric vector of predicted probabilities [0,1].
#' @param obs  Integer/numeric binary observed outcomes (0/1).
#' @return Named numeric vector with AUC, TSS, CBI.
#' @details The AUC fixes `pROC::roc(direction = "<")` so a predictor that
#'   ranks absences above presences correctly reports AUC < 0.5 instead of
#'   being silently flipped by `direction = "auto"`. The TSS threshold is the
#'   max-Youden point chosen *on the evaluation set itself* (the usual
#'   max-Youden convention for SDM evaluation), so TSS is mildly optimistic
#'   when the same threshold is reused elsewhere.
#' @keywords internal
#' @noRd
evaluate_model_full <- function(pred, obs) {
  pred <- pmin(pmax(as.numeric(pred), 1e-7), 1 - 1e-7)
  obs  <- as.integer(obs)

  # -- AUC (ROC), direction fixed: worse-than-random predictors must report
  #    AUC < 0.5, not be mirrored to 1 - AUC by direction = "auto".
  auc_val <- tryCatch({
    as.numeric(pROC::auc(pROC::roc(obs, pred, quiet = TRUE, direction = "<")))
  }, error = function(e) NA_real_)

  # -- TSS (at threshold maximising sensitivity+specificity) ------------------
  tss_val <- tryCatch({
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE, direction = "<")
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
#' Fixed equal-width binning (101 bins): Spearman correlation between the
#' bin centres of predicted suitability and observed presence frequency.
#' This is a fixed-bin implementation and deliberately does NOT reproduce
#' the moving-window variant of Hirzel et al. (2006); results are therefore
#' comparable across implementations only approximately.
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


#' Coerce Selected Columns to a Clean Numeric Matrix
#'
#' Converts the requested columns to numeric and imputes non-finite values with
#' the column median (falling back to zero). Shared by the variable-selection
#' and conditional-effect routines.
#'
#' @param data A `data.frame`.
#' @param vars Character vector of column names.
#' @return A numeric matrix with one column per `vars`.
#' @keywords internal
#' @noRd
.cast_numeric_matrix <- function(data, vars) {
  x <- as.data.frame(data[, vars, drop = FALSE])
  for (nm in names(x)) {
    x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
    med <- stats::median(x[[nm]], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    x[[nm]][!is.finite(x[[nm]]) | is.na(x[[nm]])] <- med
  }
  as.matrix(x)
}


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


#' Content Digest for Cache Keys
#'
#' Returns a stable hash of an arbitrary R object, used to parameterise
#' step- and tile-level checkpoints so that stale caches are never replayed
#' after a configuration or data change. Uses \pkg{digest} when available
#' and falls back to a fast FNV-1a hash over the serialised object.
#'
#' @param x Any R object.
#' @param len Number of hash characters to return (default all).
#' @return Character scalar.
#' @keywords internal
#' @noRd
.cast_digest <- function(x, len = NULL) {
  if (requireNamespace("digest", quietly = TRUE)) {
    out <- digest::digest(x)
  } else {
    # Polynomial rolling hash on the serialised bytes; double arithmetic
    # avoids the 32-bit overflow that breaks a bitwXor-based FNV-1a in R.
    raw <- serialize(x, NULL, version = 3)
    h <- 0
    for (b in as.integer(raw)) {
      h <- (h * 1000003 + (b + 2^31)) %% 2^52
    }
    out <- sprintf("%013.0f", h)
  }
  if (!is.null(len)) out <- substr(out, 1L, len)
  out
}


#' Configure Plot Fonts
#'
#' Sets the font family used by castSDM figures. The package default is the
#' portable \code{"sans"} family (available in every R graphics device); on
#' Windows a user can switch to \code{"Arial"}, which is registered with
#' [grDevices::windowsFont()] when available. The choice is stored in the
#' session option `castSDM.font_family` and read by every plot method; the
#' global ggplot2 theme is never modified. Use [cast_safe_ggsave()] for
#' export so PNG/PDF devices handle fonts consistently.
#'
#' @param family Character. Font family passed to ggplot2. Default
#'   \code{"sans"}.
#'
#' @return Invisibly returns \code{family}.
#' @export
cast_set_plot_defaults <- function(family = "sans") {
  if (.Platform$OS.type == "windows" && identical(tolower(family), "arial")) {
    tryCatch(
      grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")),
      error = function(e) NULL
    )
  }
  options(castSDM.font_family = family)
  invisible(family)
}


#' Safely Save a ggplot Object
#'
#' Wrapper around [ggplot2::ggsave()] that applies castSDM's configured
#' plotting font and chooses conservative graphics devices. PNG output uses
#' \pkg{ragg} when available and otherwise falls back to base PNG. PDF output
#' uses Cairo when available and otherwise base PDF with an explicit family.
#'
#' @param filename Output filename.
#' @param plot Plot object.
#' @param ... Additional arguments passed to [ggplot2::ggsave()]. When not
#'   supplied, \code{dpi} defaults to \code{600} and \code{bg} defaults to
#'   \code{"transparent"}.
#'
#' @return Invisibly returns \code{filename}.
#' @export
cast_safe_ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  check_suggested("ggplot2", "for plot export")
  family <- getOption("castSDM.font_family", "sans")
  if (.Platform$OS.type == "windows" && identical(tolower(family), "arial")) {
    tryCatch(
      grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")),
      error = function(e) NULL
    )
  }
  if (inherits(plot, c("gg", "ggplot"))) {
    plot <- plot + ggplot2::theme(text = ggplot2::element_text(family = family))
  } else if (inherits(plot, "patchwork")) {
    plot <- plot & ggplot2::theme(text = ggplot2::element_text(family = family))
  }

  ext <- tolower(tools::file_ext(filename))
  device <- NULL
  extra <- list(...)
  if (is.null(extra$dpi)) extra$dpi <- 600L
  if (is.null(extra$bg)) extra$bg <- "transparent"
  if (identical(ext, "png")) {
    device <- if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png
    } else {
      "png"
    }
  } else if (identical(ext, "pdf")) {
    device <- if (isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf
    } else {
      grDevices::pdf
    }
    if (is.null(extra$family)) {
      extra$family <- family
    }
  }

  args <- c(
    list(filename = filename, plot = plot),
    if (!is.null(device)) list(device = device) else list(),
    extra
  )
  do.call(ggplot2::ggsave, args)
  invisible(filename)
}
