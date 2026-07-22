#' Select Variables for Species Distribution Models
#'
#' `method = "cpi"` (default) is the castSDM causal selector. It tests each
#' predictor's conditional independence from occurrence given every other
#' predictor via the conditional predictive impact (Watson & Wright 2021):
#' a random forest is cross-fitted once, each predictor is replaced by a
#' Gaussian knockoff, and the increase in log-loss measures its conditional
#' contribution. Because the learner is a random forest, symmetric unimodal
#' (bell-shaped) niche responses are detected. Predictors surviving
#' Benjamini-Hochberg FDR control are kept.
#'
#' `method = "dml"` is an alternative causal selector estimating each
#' predictor's Neyman-orthogonal partially linear effect with
#' [DoubleML][DoubleML::DoubleMLPLR]; it assumes an approximately linear
#' conditional effect and is best paired with [cast_effect()] for a signed,
#' interpretable coefficient. `method = "rf"` is retained as a conventional
#' (associational) permutation-importance benchmark for comparison.
#'
#' @param data Data frame with response, coordinates, and predictors.
#' @param response Binary response column.
#' @param method `"cpi"` (default), `"dml"`, or the conventional `"rf"`
#'   benchmark.
#' @param alpha FDR level for the CPI/DML selectors. Default `0.05`.
#' @param max_candidates Predictors tested with CPI/DML; larger sets are
#'   pre-screened by RF importance for feasibility. Also the output ceiling for
#'   the RF benchmark. Default `30`.
#' @param dml_folds Cross-fitting folds for the CPI/DML selectors. Default `5`.
#' @param num_trees Trees for the RF nuisance/benchmark forests. Default `300`.
#' @param min_vars Minimum retained variables. Default `3`.
#' @param cor_threshold Absolute-correlation threshold for the RF benchmark.
#' @param num_threads Threads for the ranger learners/benchmark. Default `1`.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_select` object.
#' @seealso [cast_effect()], [cast_counterfactual()]
#' @export
cast_select <- function(data,
                        response = "presence",
                        method = c("cpi", "dml", "rf"),
                        alpha = 0.05,
                        max_candidates = 30L,
                        dml_folds = 5L,
                        num_trees = 300L,
                        min_vars = 3L,
                        cor_threshold = 0.8,
                        num_threads = 1L,
                        seed = NULL,
                        verbose = TRUE) {
  method <- match.arg(method)
  env_vars <- get_env_vars(data, response)
  if (length(env_vars) < 3L) cli::cli_abort("Need at least three predictors.")

  if (identical(method, "cpi")) {
    out <- .cast_select_cpi(
      data = data, env_vars = env_vars, response = response,
      alpha = alpha, max_candidates = max_candidates, n_folds = dml_folds,
      num_trees = min(as.integer(num_trees), 200L), min_vars = min_vars,
      seed = seed, verbose = verbose, num_threads = num_threads
    )
    return(new_cast_select(
      selected = out$selected, scores = out$scores,
      method = method, diagnostics = out$diagnostics
    ))
  }

  if (identical(method, "dml")) {
    out <- .cast_select_dml(
      data = data, env_vars = env_vars, response = response,
      alpha = alpha, max_candidates = max_candidates, n_folds = dml_folds,
      num_trees = min(as.integer(num_trees), 200L), min_vars = min_vars,
      seed = seed, verbose = verbose, num_threads = num_threads
    )
    return(new_cast_select(
      selected = out$selected, scores = out$scores,
      method = method, diagnostics = out$diagnostics
    ))
  }

  # -- Conventional RF permutation-importance benchmark ---------------------
  check_suggested("ranger", "for RF permutation selection")
  x <- as.data.frame(.cast_numeric_matrix(data, env_vars), check.names = FALSE)
  names(x) <- make.names(env_vars, unique = TRUE)
  x$.response <- factor(data[[response]])
  rf <- ranger::ranger(
    .response ~ ., data = x, importance = "permutation",
    num.trees = as.integer(num_trees), seed = seed %||% 42L,
    num.threads = as.integer(num_threads), verbose = FALSE
  )
  imp <- rf$variable.importance
  original <- env_vars[match(names(imp), make.names(env_vars, unique = TRUE))]
  imp <- stats::setNames(as.numeric(imp), original)
  imp <- imp[env_vars]
  imp[!is.finite(imp) | is.na(imp)] <- 0
  max_vars <- min(max(as.integer(min_vars), as.integer(max_candidates)),
                  length(env_vars))
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
    diagnostics = list(engine = "RF permutation importance",
                       cor_threshold = cor_threshold)
  )
}
