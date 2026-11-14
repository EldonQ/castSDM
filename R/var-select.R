# Conventional training-fold VIF screening; no predictor-count target.
# AI-assisted refactor (2026-11-13); see inst/AI_ASSISTANCE.md.

#' Select Environmental Variables by VIF
#'
#' Iteratively removes the highest-VIF predictor until all retained VIFs
#' meet the declared threshold. Repeated within each spatial training fold.
#' No response ranking, predictor-count target or category quota is used.
#' @param data Predictor frame with a binary presence/background response.
#' @param response Binary response name, default `"presence"`.
#' @param method `"vif"` (default) or `"full"` (all eligible predictors).
#' @param seed Optional seed for interface consistency; VIF is deterministic.
#' @param verbose Print the screening result.
#' @param env_vars Optional predictor names, excluding coordinates and response.
#' @param roles Optional named role metadata: numeric, binary, ordinal,
#'   identifier, nominal or cyclic. Identifiers and unencoded codes are excluded.
#' @param var_id,category Optional aligned character metadata, recorded only.
#' @param vif_threshold Maximum retained VIF, default `10`.
#' @return A `cast_select` with qualification reasons, initial/final VIF and
#'   the elimination log. No measurements are imputed.
#' @details VIF measures linear redundancy, not ecological importance or
#'   causality. Variables may substitute; inspect spatial-fold inclusion and
#'   fitted responses. The threshold is not a prescribed predictor count.
#' @seealso [cast_vif()], [cast_cv()]
#' @export
cast_select <- function(data, response = "presence", method = c("vif", "full"),
    seed = NULL, verbose = TRUE, env_vars = NULL, roles = NULL,
    var_id = NULL, category = NULL, vif_threshold = 10) {
  method <- match.arg(method)
  vif_threshold <- .cast_check_scalar(vif_threshold, "vif_threshold", lower = 1)
  if (!is.null(seed)) seed <- .cast_check_scalar(seed, "seed", integer = TRUE)
  if (!is.data.frame(data) || !response %in% names(data)) cli::cli_abort("Supply a predictor frame and binary response.")
  .cast_check_response(data[[response]], response)
  if (anyNA(data[[response]])) cli::cli_abort("The selection response must be complete.")
  vars <- env_vars %||% get_env_vars(data, response)
  if (!is.character(vars) || anyNA(vars) || anyDuplicated(vars) ||
      any(!vars %in% setdiff(names(data), c(response, "lon", "lat"))))
    cli::cli_abort("env_vars must name unique predictor columns, excluding coordinates and response.")
  qualification <- .cast_predictor_qualification(data, vars, roles)
  eligible <- qualification$variable[qualification$eligible]
  metadata <- function(x, name) {
    if (is.null(x)) return(rep(NA_character_, length(vars)))
    if (!is.character(x) || anyNA(x)) cli::cli_abort("{name} must be character predictor metadata.")
    if (!is.null(names(x))) x <- x[match(vars, names(x))]
    if (length(x) != length(vars) || anyNA(x)) cli::cli_abort("{name} must cover every predictor.")
    unname(x)
  }
  qualification$construct <- metadata(var_id, "var_id")
  qualification$category <- metadata(category, "category")
  qualification$selected <- FALSE
  qualification$initial_vif <- qualification$final_vif <- NA_real_
  diagnostics <- list(algorithm = method, status = "no_eligible_predictors",
    count_policy = "no prescribed predictor count or protected categories")
  if (!length(eligible)) return(new_cast_select(character(), qualification, method, diagnostics))
  x <- data[, eligible, drop = FALSE]
  .cast_check_numeric_predictors(x)
  if (any(!is.finite(as.matrix(x)))) cli::cli_abort("Selection requires complete finite predictors; declare missing-data QC first. No measurements are imputed.")
  if (method == "full") {
    selected <- eligible
    qualification$reason[qualification$eligible] <- "prespecified_full_set"
  } else {
    result <- cast_vif(x, threshold = vif_threshold, verbose = FALSE)
    selected <- result$selected
    ix <- match(eligible, qualification$variable)
    qualification$initial_vif[ix] <- unname(result$initial_vif[eligible])
    qualification$final_vif[ix] <- unname(result$final_vif[eligible])
    qualification$reason[qualification$eligible] <- "vif_redundancy_removed"
    qualification$reason[qualification$variable %in% selected] <- "vif_at_or_below_threshold"
    diagnostics$vif_threshold <- vif_threshold
    diagnostics$vif_log <- result$vif_log
    diagnostics$n_training_rows <- nrow(x)
  }
  qualification$selected <- qualification$variable %in% selected
  diagnostics$status <- "selected"
  if (verbose) cli::cli_inform(if (method == "vif")
    "VIF screening: {length(eligible)} eligible predictors -> {length(selected)} retained; VIF <= {vif_threshold}; no count target or inner CV."
    else "Full predictor set: {length(selected)} eligible predictors retained.")
  new_cast_select(selected, qualification, method, diagnostics)
}
