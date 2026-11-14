# Elastic-net screening with spatial validation; no predictor-count target.
# AI-assisted refactor (2026-11-13); see inst/AI_ASSISTANCE.md.

#' Select Environmental Variables by Spatially Validated Elastic Net
#'
#' Fits a binomial elastic net (equal L1/L2 mixing) to linear and quadratic
#' environmental terms. Spatial validation selects the strongest penalty
#' within one fold standard error of the minimum mean log loss. An original
#' variable is retained when either of its terms has a nonzero coefficient.
#' There is no predictor-count target, coefficient ranking or category quota.
#'
#' @param data Predictor frame with a binary presence/background response.
#' @param response Binary response name, default `"presence"`.
#' @param method `"screen"` (spatial elastic net) or `"full"` (all eligible
#'   varying numeric predictors without statistical screening).
#' @param seed Random seed for spatial folds.
#' @param verbose Print progress and the resulting variable count.
#' @param env_vars Optional explicit predictor names. Coordinates and response
#'   are never predictors.
#' @param roles Optional named predictor-role vector: `numeric`, `binary`,
#'   `ordinal`, `identifier`, `nominal`, or `cyclic`. Identifiers are excluded;
#'   nominal/cyclic codes require a declared numeric encoding before use.
#' @param var_id,category Optional aligned character metadata, recorded only.
#'   Multiple scales of a construct are not forcibly reduced to one variable.
#' @param inner_k Number of inner spatial folds, default `5`.
#' @param block_method Spatial fold builder, `"blockcv"` or `"grid"`.
#' @param block_size Optional spatial block size in metres.
#' @param fold_ids Optional externally supplied spatial fold assignments.
#' @param coords_crs Coordinate CRS, required for projected coordinates.
#' @return A `cast_select` with selected variables, qualification/reason table,
#'   spatial risk curve, fold losses and provenance. An empty selected set is
#'   a valid intercept-only result; no variable is forced into the model.
#' @details Linear and quadratic terms accommodate monotonic and simple
#'   unimodal responses without spline knots or a separate complexity search.
#'   Binary variables use one linear term. Centering/scaling and term creation
#'   are rebuilt using only the training rows of each inner fold. Elastic-net
#'   mixing is fixed at 0.5; the regularisation strength is spatially validated.
#'   Fold paths are aligned by their penalty/max-penalty ratio. The null model
#'   is evaluated explicitly. A minimum at the weak-penalty boundary extends
#'   the search. Failed fits/folds are reported, never silently dropped.
#'   The one-SE rule is a parsimony heuristic, not a confidence interval.
#'   Predictors must be complete and finite; no measurements or response labels
#'   are generated or imputed. Loss concerns sampled presence/background
#'   labels, not ecological occurrence probabilities. This simple additive
#'   screen can miss complex responses and pure interactions; correlated
#'   variables can substitute. Inspect outer-fold selection frequencies.
#'   Selection does not identify causal drivers.
#' @references Friedman J, Hastie T, Tibshirani R (2010). Regularization paths
#'   for generalized linear models via coordinate descent. Journal of
#'   Statistical Software 33(1):1-22. doi:10.18637/jss.v033.i01.
#' @seealso [cast_cv()], [cast_importance()]
#' @export
cast_select <- function(data, response = "presence", method = c("screen", "full"),
    seed = NULL, verbose = TRUE, env_vars = NULL, roles = NULL,
    var_id = NULL, category = NULL, inner_k = 5L,
    block_method = c("blockcv", "grid"), block_size = NULL, fold_ids = NULL,
    coords_crs = NULL) {
  if (missing(method)) method <- "screen"
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("screen", "full")) cli::cli_abort("method must be screen or full.")
  block_method <- match.arg(block_method)
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
  if (!length(eligible)) return(new_cast_select(character(), qualification, method,
    list(status = "no_eligible_predictors", count_policy = "no prescribed predictor count")))
  x <- data[, eligible, drop = FALSE]
  .cast_check_numeric_predictors(x)
  if (any(!is.finite(as.matrix(x)))) cli::cli_abort("Selection requires complete finite predictors; qualify or clean missing inputs explicitly.")
  if (method == "full") {
    qualification$selected <- qualification$eligible
    qualification$reason[qualification$selected] <- "prespecified_full_set"
    return(new_cast_select(eligible, qualification, method,
      list(status = "selected", count_policy = "all eligible predictors")))
  }
  check_suggested("glmnet", "for spatial elastic-net variable selection")
  inner_k <- .cast_check_scalar(inner_k, "inner_k", lower = 2, integer = TRUE)
  folds <- .cast_selection_folds(data, response, fold_ids, inner_k, block_method, block_size, seed, coords_crs)
  y <- as.integer(data[[response]])
  full_basis <- .cast_selection_basis(x)
  # glmnet's documented path floors depend on rows versus design columns.
  floor <- if (nrow(x) < ncol(full_basis$matrix)) 0.01 else 0.0001
  repeat {
    if (verbose) cli::cli_inform("Selection: elastic-net path ({nrow(x)} rows; {length(eligible)} eligible predictors).")
    fractions <- exp(seq(0, log(floor), length.out = 100L))
    full <- .cast_selection_path(full_basis$matrix, y, fractions)
    losses <- matrix(NA_real_, length(unique(folds)), length(fractions) + 1L,
      dimnames = list(paste0("fold_", seq_len(length(unique(folds)))),
        paste0("candidate_", 0:length(fractions))))
    fold_details <- vector("list", nrow(losses))
    for (i in seq_len(nrow(losses))) {
      train <- folds != i; test <- !train
      if (length(unique(y[train])) != 2L)
        cli::cli_abort("Inner spatial fold {i} has a single-class training set. Change the spatial design; no random-fold fallback is used.")
      basis <- .cast_selection_basis(x[train, , drop = FALSE])
      losses[i, 1L] <- .cast_selection_loss(rep(stats::qlogis(mean(y[train])), sum(test)), y[test])
      if (!length(basis$variables)) {
        losses[i, -1L] <- losses[i, 1L]
      } else {
        path <- .cast_selection_path(basis$matrix, y[train], fractions)
        design <- .cast_selection_design(basis, x[test, , drop = FALSE])
        prediction <- as.matrix(stats::predict(path, newx = design, type = "link"))
        losses[i, -1L] <- vapply(seq_along(fractions), function(j)
          .cast_selection_loss(prediction[, j], y[test]), numeric(1))
      }
      fold_details[[i]] <- data.frame(fold = i, n_train = sum(train), n_test = sum(test),
        train_presence = sum(y[train]), test_presence = sum(y[test]),
        varying_predictors = length(basis$variables))
      if (verbose) cli::cli_inform("Selection: inner spatial fold {i}/{nrow(losses)} evaluated.")
    }
    mean_loss <- colMeans(losses)
    se <- apply(losses, 2L, stats::sd)/sqrt(nrow(losses))
    best <- which.min(mean_loss)
    curve <- data.frame(candidate = seq_along(mean_loss)-1L,
      penalty_fraction = c(Inf, fractions), lambda = c(Inf, full$lambda),
      mean_logloss = mean_loss, se = se, selected = FALSE)
    if (best == nrow(curve)) {
      floor <- floor/10
      if (floor <= .Machine$double.eps) cli::cli_abort("The selection penalty boundary remains unresolved at numerical precision.")
      if (verbose) cli::cli_inform("Selection: extending the weak-penalty search; no variable count is imposed.")
      next
    }
    choice <- which(mean_loss <= mean_loss[best] + se[best])[1L]
    curve$selected[choice] <- TRUE
    beta <- if (choice == 1L) rep(0, ncol(full_basis$matrix)) else
      as.numeric(stats::coef(full, s = curve$lambda[choice]))[-1L]
    selected <- full_basis$variables[vapply(seq_along(full_basis$variables), function(j)
      any(beta[full_basis$group == j] != 0), logical(1))]
    break
  }
  qualification$selected <- qualification$variable %in% selected
  qualification$reason[qualification$eligible] <- "zero_terms_at_spatial_1se"
  qualification$reason[qualification$selected] <- "nonzero_term_at_spatial_1se"
  qualification$group_signal_sd <- NA_real_
  for (j in seq_along(full_basis$variables)) {
    at <- full_basis$group == j
    signal <- as.numeric(full_basis$matrix[, at, drop = FALSE] %*% beta[at])
    qualification$group_signal_sd[match(full_basis$variables[j], qualification$variable)] <- stats::sd(signal)
  }
  diagnostics <- list(status = if (length(selected)) "selected" else "null_model_selected",
    algorithm = "elastic_net", alpha = 0.5, terms = "linear + quadratic (binary: linear only)",
    path_floor = floor, inner_k = nrow(losses), folds = folds,
    fold_design = attr(folds, "design"), fold_details = do.call(rbind, fold_details),
    risk_curve = curve, fold_losses = losses, choice = choice-1L,
    lambda = curve$lambda[choice], penalty_fraction = curve$penalty_fraction[choice],
    rule = "strongest penalty within one fold-SE of minimum mean spatial log loss; null included",
    count_policy = "no prescribed predictor count; no protected category slots",
    limitations = "Simple additive screen; complex responses and pure interactions may be missed, correlated variables may substitute; selection is not causality.",
    software = c(glmnet = as.character(utils::packageVersion("glmnet"))))
  if (verbose) cli::cli_inform("Spatial elastic net: {length(eligible)} eligible predictors; {length(selected)} retained by spatial 1-SE (no count target).")
  new_cast_select(selected, qualification, method, diagnostics)
}

.cast_selection_loss <- function(eta, y) {
  if (any(!is.finite(eta))) cli::cli_abort("A selection model returned non-finite validation linear predictors.")
  mean(pmax(eta, 0)-y*eta+log1p(exp(-abs(eta))))
}

.cast_selection_path <- function(x, y, fractions) {
  args <- list(x = x, y = y, family = "binomial", alpha = 0.5, standardize = TRUE)
  pilot <- do.call(glmnet::glmnet, c(args, list(nlambda = 5L, lambda.min.ratio = 0.1)))
  args$lambda <- pilot$lambda[1L] * fractions
  fit <- do.call(glmnet::glmnet, args)
  if (fit$jerr != 0L || length(fit$lambda) != length(fractions) || any(!is.finite(fit$beta@x)))
    cli::cli_abort("Elastic-net solver did not complete the requested penalty path; no partially fitted path is used.")
  fit
}

.cast_selection_folds <- function(data, response, ids, k, method, size, seed, coords_crs = NULL) {
  if (is.null(ids)) {
    if (!all(c("lon", "lat") %in% names(data)))
      cli::cli_abort("Automatic selection requires spatial coordinates or externally supplied spatial fold_ids.")
    ids <- make_spatial_folds(data$lon, data$lat, k = k, method = method,
      seed = seed, y = data[[response]], block_size = size, coords_crs = coords_crs)
  }
  if (!is.numeric(ids) || length(ids) != nrow(data) || any(!is.finite(ids)) ||
      any(ids != as.integer(ids)) || length(unique(ids)) < 2L)
    cli::cli_abort("fold_ids must assign every row to one of at least two spatial folds.")
  design <- attr(ids, "design") %||% list(requested = "external", effective = "external", block_size_m = NA_real_)
  out <- as.integer(factor(ids, levels = sort(unique(ids))))
  attr(out, "design") <- design
  out
}

# Transform actual measurements, using training rows only; no imputation/jitter.
.cast_selection_basis <- function(x) {
  variables <- names(x)[vapply(x, .cast_is_varying, logical(1))]
  specs <- matrices <- vector("list", length(variables))
  for (j in seq_along(variables)) {
    v <- x[[variables[j]]]
    center <- mean(v); scale <- stats::sd(v)
    z <- (v-center)/scale
    quadratic <- length(unique(v)) > 2L
    matrices[[j]] <- if (quadratic) cbind(z, z^2) else matrix(z, ncol = 1L)
    specs[[j]] <- list(center = center, scale = scale, quadratic = quadratic)
  }
  names(specs) <- variables
  design <- if (length(matrices)) do.call(cbind, matrices) else matrix(numeric(), nrow(x), 0L)
  group <- rep(seq_along(variables), vapply(matrices, ncol, integer(1)))
  # glmnet requires two design columns. The constant is an intercept basis,
  # excluded internally for zero variance; it is never an environmental axis.
  if (ncol(design) == 1L) { design <- cbind(design, 1); group <- c(group, 0L) }
  colnames(design) <- paste0("term_", seq_len(ncol(design)))
  list(matrix = design, group = group, variables = variables, specs = specs)
}

.cast_selection_design <- function(basis, x) {
  parts <- lapply(basis$variables, function(v) {
    spec <- basis$specs[[v]]; z <- (x[[v]]-spec$center)/spec$scale
    if (spec$quadratic) cbind(z, z^2) else matrix(z, ncol = 1L)
  })
  out <- do.call(cbind, parts)
  if (ncol(out) == 1L) out <- cbind(out, 1)
  colnames(out) <- colnames(basis$matrix)
  out
}
