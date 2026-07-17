#' Select Variables for Species Distribution Models
#'
#' `method = "stable"` is the CAST selector. It combines a fast response-aware
#' shortlist with a cross-environment conditional-invariance test and a greedy
#' minimum-set search. `method = "rf"` is retained as a conventional benchmark.
#'
#' The stable selector is causal-inspired: it tests whether the response still
#' depends on the data-derived environment after conditioning on a predictor
#' set. It does not claim automatic identification of ecological causes.
#'
#' @param data Data frame with response, coordinates, and predictors.
#' @param response Binary response column.
#' @param method `"stable"` (default), `"stable_no_invariance"` (ablation),
#'   or conventional `"rf"`.
#' @param domains Optional precomputed domain factor. It must have one value per
#'   row and should be constructed from training data only.
#' @param domain_method `"spatial"` or `"environment"`.
#' @param n_domains Requested number of domains.
#' @param num_trees Number of RF trees used for shortlisting/ranking.
#' @param min_vars Minimum retained variables.
#' @param max_vars Maximum stable-selector shortlist size or RF output size.
#' @param alpha Significance threshold for invariance tests.
#' @param loss_tolerance Maximum allowed increase in leave-domain-out log-loss
#'   during greedy removal.
#' @param cor_threshold Absolute-correlation threshold for the RF benchmark.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_select` object.
#' @export
cast_select <- function(data,
                        response = "presence",
                        method = c("stable", "stable_no_invariance", "rf"),
                        domains = NULL,
                        domain_method = c("spatial", "environment"),
                        n_domains = 4L,
                        num_trees = 300L,
                        min_vars = 3L,
                        max_vars = 12L,
                        alpha = 0.05,
                        loss_tolerance = 0.02,
                        cor_threshold = 0.8,
                        seed = NULL,
                        verbose = TRUE) {
  method <- match.arg(method)
  domain_method <- match.arg(domain_method)
  env_vars <- get_env_vars(data, response)
  if (length(env_vars) < 3L) cli::cli_abort("Need at least three predictors.")

  if (method %in% c("stable", "stable_no_invariance")) {
    out <- .cast_select_stable(
      data = data, env_vars = env_vars, response = response,
      domains = domains, domain_method = domain_method,
      n_domains = n_domains, max_candidates = max_vars,
      min_vars = min_vars, alpha = alpha, loss_tolerance = loss_tolerance,
      num_trees = num_trees,
      use_invariance = identical(method, "stable"),
      seed = seed, verbose = verbose
    )
    return(new_cast_select(
      selected = out$selected, scores = out$scores,
      method = method, diagnostics = out$diagnostics
    ))
  }

  check_suggested("ranger", "for RF permutation selection")
  x <- as.data.frame(.cast_numeric_matrix(data, env_vars), check.names = FALSE)
  names(x) <- make.names(env_vars, unique = TRUE)
  x$.response <- factor(data[[response]])
  rf <- ranger::ranger(
    .response ~ ., data = x, importance = "permutation",
    num.trees = as.integer(num_trees), seed = seed %||% 42L,
    num.threads = 1L
  )
  imp <- rf$variable.importance
  original <- env_vars[match(names(imp), make.names(env_vars, unique = TRUE))]
  imp <- stats::setNames(as.numeric(imp), original)
  imp <- imp[env_vars]
  imp[!is.finite(imp) | is.na(imp)] <- 0
  max_vars <- min(max(as.integer(min_vars), as.integer(max_vars)), length(env_vars))
  ranked <- env_vars[order(imp, decreasing = TRUE)]

  selected <- character()
  for (v in ranked) {
    if (length(selected) >= max_vars) break
    if (!length(selected)) {
      selected <- v
    } else {
      cors <- abs(stats::cor(
        .cast_numeric_matrix(data, v),
        .cast_numeric_matrix(data, selected), use = "pairwise.complete.obs"
      ))
      if (!any(cors >= cor_threshold, na.rm = TRUE)) selected <- c(selected, v)
    }
  }
  if (length(selected) < min_vars) {
    selected <- unique(c(selected, ranked))[seq_len(min(min_vars, length(ranked)))]
  }
  scores <- data.frame(
    variable = env_vars,
    rf_importance = unname(imp),
    selected = env_vars %in% selected,
    exclusion_reason = ifelse(env_vars %in% selected, "selected", "lower_or_redundant"),
    stringsAsFactors = FALSE
  )
  scores <- scores[order(!scores$selected, -scores$rf_importance), , drop = FALSE]
  new_cast_select(
    selected = selected, scores = scores, method = "rf",
    diagnostics = list(cor_threshold = cor_threshold)
  )
}
