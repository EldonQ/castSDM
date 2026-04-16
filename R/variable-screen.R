#' Adaptive Multi-Criteria Variable Screening
#'
#' Combines three complementary signals -- DAG centrality, ATE magnitude, and
#' RF permutation importance -- into a weighted composite score for variable
#' selection. Weights adapt based on DAG quality and ATE signal strength.
#'
#' @param dag A [cast_dag] object.
#' @param ate A [cast_ate] object.
#' @param data A `data.frame` used for RF importance calculation.
#' @param response Character. Response column name. Default `"presence"`.
#' @param min_vars Integer. Minimum variables to retain. Default `5`.
#' @param min_fraction Numeric. Minimum fraction of variables to retain.
#'   Default `0.5`.
#' @param num_trees Integer. Trees for RF importance. Default `300`.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A [cast_screen] object.
#'
#' @details
#' The composite score for each variable is:
#'
#' \deqn{S = w_{dag} \cdot s_{dag} + w_{ate} \cdot s_{ate} + w_{imp} \cdot s_{imp}}
#'
#' where weights adapt based on signal quality:
#' - \eqn{w_{dag} = 0.15 + 0.15 \times dag\_quality}
#' - \eqn{w_{ate} = 0.25 + 0.25 \times ate\_sig\_ratio}
#' - \eqn{w_{imp} = 1 - w_{dag} - w_{ate}}
#'
#' Variables are selected using k-means clustering on composite scores.
#'
#' @seealso [cast_dag()], [cast_ate()], [cast_roles()]
#'
#' @export
cast_screen <- function(dag,
                        ate,
                        data,
                        response = "presence",
                        min_vars = 5L,
                        min_fraction = 0.5,
                        num_trees = 300L,
                        seed = NULL,
                        verbose = TRUE) {
  check_suggested("ranger", "for RF importance")

  env_vars <- dag$nodes
  edges <- dag$edges

  # -- RF permutation importance --
  Y <- data[[response]]
  X <- as.data.frame(data[, env_vars, drop = FALSE])
  X[is.na(X)] <- 0
  if (!is.null(seed)) set.seed(seed)
  rf_imp <- ranger::ranger(
    y ~ .,
    data = cbind(y = as.factor(Y), X),
    num.trees = num_trees,
    importance = "permutation",
    verbose = FALSE
  )$variable.importance

  # -- DAG out-degree per variable --
  deg_df <- compute_edge_degrees(edges, env_vars)

  # -- DAG quality & ATE signal ratio --
  n_possible <- length(env_vars) * (length(env_vars) - 1) / 2
  dag_density <- nrow(edges) / max(n_possible, 1)
  dag_quality <- 1 - dag_density

  est <- ate$estimates
  n_tested <- nrow(est)
  ate_sig_ratio <- if (n_tested > 0) {
    sum(est$significant, na.rm = TRUE) / n_tested
  } else {
    0
  }

  # -- Adaptive weights --
  w_dag <- 0.15 + 0.15 * dag_quality
  w_ate <- 0.25 + 0.25 * ate_sig_ratio
  w_imp <- 1 - w_dag - w_ate

  # -- Build scoring data.frame --
  scr <- deg_df

  # Merge ATE estimates
  if (n_tested > 0) {
    ate_sub <- est[, c("variable", "coef", "p_value", "significant")]
    scr <- merge(scr, ate_sub, by = "variable", all.x = TRUE)
  } else {
    scr$coef <- NA_real_
    scr$p_value <- NA_real_
    scr$significant <- FALSE
  }
  scr$coef[is.na(scr$coef)] <- 0
  scr$p_value[is.na(scr$p_value)] <- 1
  scr$significant[is.na(scr$significant)] <- FALSE

  # RF importance
  scr$importance <- rf_imp[scr$variable]

  # -- Component scores --
  scr$score_dag <- normalize01(scr$out_degree)
  scr$score_ate_raw <- normalize01(abs(scr$coef))
  scr$ate_penalty <- pmin(
    1.0, -log10(pmax(scr$p_value, 1e-10)) / 3
  )
  scr$score_ate <- scr$score_ate_raw * scr$ate_penalty
  scr$score_imp <- normalize01(scr$importance)

  scr$score_total <- w_dag * scr$score_dag +
    w_ate * scr$score_ate +
    w_imp * scr$score_imp

  scr <- scr[order(-scr$score_total), ]
  rownames(scr) <- NULL

  # -- Score-gap selection (replaces k-means(k=2)) --
  # Rank variables by composite score and find the largest drop between
  # consecutive scores.  Variables above the largest gap are retained.
  # This is more stable than k-means because it does not depend on random
  # initialisation and does not assume a bimodal score distribution.
  min_keep <- max(min_vars, ceiling(length(env_vars) * min_fraction))
  max_keep <- nrow(scr)

  scores_sorted <- scr$score_total  # already sorted descending

  selected <- if (max_keep <= min_keep) {
    # Too few variables to search for a gap — keep all
    scr$variable
  } else {
    # Compute consecutive score differences (positive = big drop)
    gaps      <- scores_sorted[-max_keep] - scores_sorted[-1L]
    # Only consider cut-points that respect min/max bounds
    valid_idx <- seq_len(max_keep - 1L)
    valid_idx <- valid_idx[valid_idx >= min_keep & valid_idx < max_keep]

    if (length(valid_idx) == 0L) {
      scr$variable[seq_len(min_keep)]
    } else {
      cut_at <- valid_idx[which.max(gaps[valid_idx])]
      scr$variable[seq_len(cut_at)]
    }
  }

  weights <- c(w_dag = w_dag, w_ate = w_ate, w_imp = w_imp)

  if (verbose) {
    cli::cli_inform(
      "Screening: {length(selected)}/{length(env_vars)} variables selected."
    )
  }

  new_cast_screen(
    selected = selected,
    scores = scr,
    weights = weights
  )
}
