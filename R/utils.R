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
    "id", "ID", "site", "cell_id", "grid_id"
  )
  exclude <- unique(c(response, coords, default_meta, meta))
  nms <- setdiff(names(data), exclude)
  nms[vapply(data[nms], is.numeric, logical(1))]
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
  pred <- pmin(pmax(pred, 1e-7), 1 - 1e-7)
  auc_val <- tryCatch({
    as.numeric(pROC::auc(pROC::roc(obs, pred, quiet = TRUE)))
  }, error = function(e) NA_real_)
  tss_val <- tryCatch({
    roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
    coords <- pROC::coords(roc_obj, "best", ret = c("sensitivity", "specificity"))
    as.numeric(coords$sensitivity + coords$specificity - 1)
  }, error = function(e) NA_real_)
  c(auc = auc_val, tss = tss_val)
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
