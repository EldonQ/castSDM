#' Select Variables for Species Distribution Models
#'
#' `method = "dml"` is the castSDM causal selector. It estimates each
#' predictor's Neyman-orthogonal partially linear effect on occurrence with
#' [DoubleML][DoubleML::DoubleMLPLR] while flexibly controlling for every other
#' predictor, then keeps the predictors whose effect survives Benjamini-Hochberg
#' FDR control. The only ecological choice is the FDR level; cross-fitting folds
#' and the random-forest nuisance learner are method defaults. `method = "rf"`
#' is retained as a conventional (associational) permutation-importance
#' benchmark for comparison.
#'
#' @param data Data frame with response, coordinates, and predictors.
#' @param response Binary response column.
#' @param method `"dml"` (default) or the conventional `"rf"` benchmark.
#' @param alpha FDR level for the DML selector. Default `0.05`.
#' @param max_candidates Predictors tested with DML; larger sets are
#'   pre-screened by RF importance for feasibility. Also the output ceiling for
#'   the RF benchmark. Default `30`.
#' @param dml_folds Cross-fitting folds for the DML selector. Default `5`.
#' @param num_trees Trees for the RF nuisance/benchmark forests. Default `300`.
#' @param min_vars Minimum retained variables. Default `3`.
#' @param cor_threshold Absolute-correlation threshold for the RF benchmark.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `cast_select` object.
#' @seealso [cast_effect()], [cast_counterfactual()]
#' @export
cast_select <- function(data,
                        response = "presence",
                        method = c("dml", "rf"),
                        alpha = 0.05,
                        max_candidates = 30L,
                        dml_folds = 5L,
                        num_trees = 300L,
                        min_vars = 3L,
                        cor_threshold = 0.8,
                        seed = NULL,
                        verbose = TRUE) {
  method <- match.arg(method)
  env_vars <- get_env_vars(data, response)
  if (length(env_vars) < 3L) cli::cli_abort("Need at least three predictors.")

  if (identical(method, "dml")) {
    out <- .cast_select_dml(
      data = data, env_vars = env_vars, response = response,
      alpha = alpha, max_candidates = max_candidates, n_folds = dml_folds,
      num_trees = min(as.integer(num_trees), 200L), min_vars = min_vars,
      seed = seed, verbose = verbose
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
    num.threads = 1L
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
