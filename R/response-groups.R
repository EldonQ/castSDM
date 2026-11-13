# AI-assisted development (2026-11-11); see inst/AI_ASSISTANCE.md.
# Summaries of a frozen model's responses, never new ecological observations.

.cast_response_nodes <- function(tree, data) {
  nodes <- rep(1L, nrow(data))
  frame <- tree$frame
  ids <- as.integer(row.names(frame))
  offset <- 1L
  for (i in seq_len(nrow(frame))) {
    if (frame$var[i] == "<leaf>") next
    split <- tree$splits[offset, ]
    at <- which(nodes == ids[i])
    if (length(at)) {
      left <- data[[frame$var[i]]][at] < split["index"]
      if (split["ncat"] > 0) left <- !left
      nodes[at] <- 2L * ids[i] + as.integer(!left)
    }
    offset <- offset + 1L + frame$ncompete[i] + frame$nsurrogate[i]
  }
  nodes
}

.cast_response_frame <- function(data, variables) {
  if (!is.data.frame(data) || !all(variables %in% names(data)))
    cli::cli_abort("{.arg data} must contain all response-group predictors.")
  x <- data[, variables, drop = FALSE]
  .cast_check_numeric_predictors(x, arg = "data")
  names(x) <- paste0("v", seq_along(variables))
  x
}

#' Find Environmental Subgroups with Similar Model Responses
#'
#' Summarise paired changes in a frozen SDM's predictions using a small
#' regression tree. Supplied spatial folds select pruning complexity by the
#' one-standard-error rule. A constant response or unsupported subdivision
#' can correctly produce a single group.
#'
#' @param data Data frame containing only numeric explanatory predictors.
#'   Coordinates, presence/background labels and row identifiers should not
#'   be included as predictors. Rows must match `response` and `folds` exactly.
#' @param response Numeric vector of paired predicted suitability changes,
#'   for example `attribution$data$total_change` from [cast_attribution()].
#'   It is a model output, not an observed treatment effect.
#' @param folds Vector of supplied spatial fold identifiers, one per row.
#'   At least two folds must contain supported, complete rows. No random folds
#'   or invented ecological responses are generated.
#' @param supported Optional logical vector identifying rows on which the
#'   response may be summarised. Unknown support is excluded. Default all rows.
#' @param maxdepth Maximum tree depth, default 3 (at most eight leaves).
#'   This controls explanatory complexity, not an ecological threshold.
#' @param minbucket Optional minimum leaf size; `NULL` uses the `rpart` default.
#' @return A `cast_response_groups` object containing the tree, rules, group
#'   summaries, row membership, spatial tuning losses and the settings used.
#' @details The method follows finite-change / conditional average marginal
#'   effect explanations. It does not estimate CATE, establish causality or
#'   validate the SDM. Cross-validation selects a descriptive tree of the
#'   supplied model responses; its error is not independent ecological
#'   validation. Environmental groups need not be geographically contiguous.
#'   Supplied folds do not remove spatial dependence within folds.
#'   Prediction assigns complete new rows to existing leaves; it does not
#'   assess support. Apply the response's support mask to mapped predictions.
#' @references
#' Scholbeck et al. (2024). Marginal effects for non-linear prediction
#' functions. \doi{10.1007/s10618-023-00993-x}.
#' Lowe et al. (2024). fmeffects: An R Package for Forward Marginal Effects.
#' \doi{10.32614/RJ-2024-024}.
#' @export
cast_response_groups <- function(data, response, folds, supported = NULL,
                                 maxdepth = 3L, minbucket = NULL) {
  check_suggested("rpart", "for response grouping")
  if (!is.data.frame(data) || !nrow(data) || !ncol(data) ||
      is.null(names(data)) || anyDuplicated(names(data)) || anyNA(names(data)) ||
      any(!nzchar(names(data))))
    cli::cli_abort("{.arg data} must be a non-empty, uniquely named data frame.")
  n <- nrow(data)
  if (!is.numeric(response) || length(response) != n)
    cli::cli_abort("{.arg response} must be numeric with one value per data row.")
  if (length(folds) != n || !is.atomic(folds))
    cli::cli_abort("{.arg folds} must have one spatial fold identifier per row.")
  if (is.null(supported)) supported <- rep(TRUE, n)
  if (!is.logical(supported) || length(supported) != n)
    cli::cli_abort("{.arg supported} must be logical with one value per row.")
  .cast_check_scalar(maxdepth, "maxdepth", lower = 1, upper = 30, integer = TRUE)
  if (!is.null(minbucket))
    .cast_check_scalar(minbucket, "minbucket", lower = 1, integer = TRUE)
  variables <- names(data)
  x <- .cast_response_frame(data, variables)
  keep <- .cast_finite_rows(x) & is.finite(response) &
    !is.na(supported) & supported & !is.na(folds)
  ids <- unique(as.character(folds[keep]))
  if (sum(keep) < 4L || length(ids) < 2L)
    cli::cli_abort("Response grouping needs complete supported rows in at least two spatial folds.")
  x <- x[keep, , drop = FALSE]
  y <- response[keep]
  f <- as.character(folds[keep])
  controls <- list(cp = 0, xval = 0L, maxdepth = as.integer(maxdepth),
    maxsurrogate = 0L, usesurrogate = 0L)
  if (!is.null(minbucket)) controls$minbucket <- as.integer(minbucket)
  control <- do.call(rpart::rpart.control, controls)
  grow <- function(rows) {
    d <- x[rows, , drop = FALSE]
    d$.response <- y[rows]
    rpart::rpart(stats::as.formula(".response ~ ."), data = d,
      method = "anova", control = control, model = FALSE, x = FALSE, y = FALSE)
  }
  full <- grow(seq_len(nrow(x)))
  candidates <- sort(unique(c(1, full$cptable[, "CP"])), decreasing = TRUE)
  loss <- matrix(NA_real_, length(ids), length(candidates))
  baseline <- numeric(length(ids))
  for (i in seq_along(ids)) {
    train <- which(f != ids[i]); test <- which(f == ids[i])
    if (length(train) < 2L)
      cli::cli_abort("A spatial training fold has fewer than two supported rows.")
    tree <- grow(train)
    baseline[i] <- mean((y[test] - mean(y[train]))^2)
    for (j in seq_along(candidates)) {
      pruned <- rpart::prune(tree, cp = candidates[j])
      p <- stats::predict(pruned, x[test, , drop = FALSE])
      loss[i, j] <- mean((y[test] - p)^2)
    }
  }
  errors <- colMeans(loss)
  se <- apply(loss, 2L, stats::sd) / sqrt(nrow(loss))
  best <- which.min(errors)
  # candidates descend from the root: first eligible is the simplest tree.
  chosen <- which(errors <= errors[best] + se[best])[1L]
  tree <- rpart::prune(full, cp = candidates[chosen])
  # Numeric predictor terms need no caller environment. Retaining grow()'s
  # frame captures the data, intermediate trees and self-references, making
  # saved values larger and their content hashes unstable after publication.
  if (inherits(tree$terms, "formula") || inherits(tree$terms, "terms"))
    environment(tree$terms) <- baseenv()
  membership <- .cast_response_nodes(tree, x)
  leaves <- as.integer(row.names(tree$frame)[tree$frame$var == "<leaf>"])
  paths <- rpart::path.rpart(tree, nodes = leaves, print.it = FALSE)
  rules <- vapply(paths, function(path) {
    rule <- paste(path[-1L], collapse = " & ")
    if (!nzchar(rule)) return("All supported rows")
    for (j in rev(seq_along(variables)))
      rule <- gsub(paste0("\\bv", j, "\\b"), variables[j], rule)
    rule
  }, character(1))
  summaries <- lapply(seq_along(leaves), function(i) {
    at <- membership == leaves[i]; z <- y[at]
    data.frame(group = leaves[i], n = sum(at), mean_response = mean(z),
      mean_abs_response = mean(abs(z)), sd_response = stats::sd(z),
      positive_fraction = mean(z > 0), folds = length(unique(f[at])),
      rule = unname(rules[i]), row.names = NULL)
  })
  group <- rep(NA_integer_, n); group[keep] <- membership
  fitted <- rep(NA_real_, n)
  fitted[keep] <- as.numeric(stats::predict(tree, x))
  structure(list(tree = tree, variables = variables,
    summary = do.call(rbind, summaries), group = group, fitted = fitted,
    cv = data.frame(cp = candidates, leaves = vapply(candidates, function(cp)
      sum(rpart::prune(full, cp = cp)$frame$var == "<leaf>"), integer(1)),
      mean_mse = errors, se = se, selected = seq_along(candidates) == chosen),
    fold_losses = data.frame(fold = ids, mse = loss[, chosen], baseline_mse = baseline),
    settings = list(control = control, cp = candidates[chosen],
      selection = "Spatial-fold one-standard-error pruning", n_input = n,
      n_supported = sum(keep), folds = length(ids)),
    interpretation = "Descriptive subgroups of frozen-model responses; not causal effects or independent ecological validation"),
    class = "cast_response_groups")
}

#' @rdname cast_response_groups
#' @param object A `cast_response_groups` object.
#' @param newdata Data frame containing the group's explanatory predictors.
#' @param type Return `"group"` membership or group mean `"response"`.
#' @param ... Reserved for S3 compatibility.
#' @export
predict.cast_response_groups <- function(object, newdata, type = c("group", "response"), ...) {
  type <- match.arg(type)
  x <- .cast_response_frame(newdata, object$variables)
  ok <- .cast_finite_rows(x)
  out <- if (type == "group") rep(NA_integer_, nrow(x)) else rep(NA_real_, nrow(x))
  if (any(ok)) out[ok] <- if (type == "group")
    .cast_response_nodes(object$tree, x[ok, , drop = FALSE]) else
      as.numeric(stats::predict(object$tree, x[ok, , drop = FALSE]))
  out
}

#' @rdname cast_response_groups
#' @param x A `cast_response_groups` object.
#' @export
print.cast_response_groups <- function(x, ...) {
  cli::cli_h2("Environmental response groups")
  print(x$summary, row.names = FALSE)
  cli::cli_text("{x$interpretation}")
  invisible(x)
}
