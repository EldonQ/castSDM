#' Estimate Average Treatment Effects via Double Machine Learning
#'
#' For each environmental variable, estimates the Average Treatment Effect
#' (ATE) on species presence using Double/Debiased Machine Learning
#' (Chernozhukov et al. 2018). Variables are binarized at their median,
#' and Random Forests are used as nuisance functions with K-fold
#' cross-fitting.
#'
#' @param data A `data.frame` with `presence` column and environmental
#'   variables.
#' @param response Character. Name of response column. Default `"presence"`.
#' @param variables Character vector of variable names to test. Default
#'   `NULL` (all numeric columns except `response`, `lon`, `lat`).
#' @param K Integer. Number of cross-fitting folds. Default `2`.
#' @param num_trees Integer. Number of trees for nuisance Random Forests.
#'   Default `300`.
#' @param alpha Numeric. Significance level. Default `0.05`.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A [cast_ate] object.
#'
#' @details
#' The DML procedure for each variable \eqn{X}:
#' 1. Binarize \eqn{X} at its median to create treatment
#'    \eqn{T \in \{0,1\}}.
#' 2. K-fold cross-fitting: train RF models for \eqn{E[Y|W]} and
#'    \eqn{E[T|W]} on training folds, compute residuals on held-out folds.
#' 3. Estimate ATE =
#'    \eqn{\sum(\tilde{T} \cdot \tilde{Y}) / \sum(\tilde{T}^2)}.
#' 4. Compute standard error via the influence function.
#'
#' @references
#' Chernozhukov, V., Chetverikov, D., Demirer, M., et al. (2018).
#' Double/debiased machine learning for treatment and structural parameters.
#' *The Econometrics Journal*, 21(1), C1-C68.
#'
#' @seealso [cast_dag()], [cast_screen()]
#'
#' @export
cast_ate <- function(data,
                     response = "presence",
                     variables = NULL,
                     K = 2L,
                     num_trees = 300L,
                     alpha = 0.05,
                     seed = NULL,
                     verbose = TRUE) {
  check_suggested("ranger", "for DML nuisance models")

  env_vars <- variables %||% get_env_vars(data, response = response)

  Y <- data[[response]]
  X_full <- as.data.frame(data[, env_vars, drop = FALSE])
  X_full[is.na(X_full)] <- 0

  if (verbose) {
    cli::cli_inform(
      "Estimating ATE for {length(env_vars)} variables (K={K})..."
    )
  }

  results <- vector("list", length(env_vars))
  for (i in seq_along(env_vars)) {
    v <- env_vars[i]
    results[[i]] <- tryCatch({
      T_bin <- as.integer(
        X_full[[v]] > stats::median(X_full[[v]], na.rm = TRUE)
      )
      W <- X_full[, setdiff(env_vars, v), drop = FALSE]
      if (!is.null(seed)) set.seed(seed)
      res <- dml_ate_internal(
        Y = Y, T_var = T_bin, W = W,
        K = K, num_trees = num_trees
      )
      data.frame(
        variable = v, coef = res$ate, se = res$se,
        p_value = res$p_value,
        significant = res$p_value < alpha,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        variable = v, coef = NA_real_, se = NA_real_,
        p_value = NA_real_, significant = FALSE,
        stringsAsFactors = FALSE
      )
    })
  }

  estimates <- do.call(rbind, results)
  n_sig <- sum(estimates$significant, na.rm = TRUE)
  if (verbose) {
    cli::cli_inform(
      "ATE: {n_sig}/{length(env_vars)} significant (p < {alpha})."
    )
  }

  new_cast_ate(estimates = estimates, K = K, alpha = alpha)
}


#' Double Machine Learning ATE (Internal)
#'
#' Core DML algorithm for a single variable.
#'
#' @param Y Numeric response vector.
#' @param T_var Binary treatment vector.
#' @param W data.frame of confounders.
#' @param K Number of cross-fitting folds.
#' @param num_trees Number of RF trees.
#' @return List with `ate`, `se`, `p_value`.
#' @keywords internal
#' @noRd
dml_ate_internal <- function(Y, T_var, W, K = 2, num_trees = 300) {
  n <- length(Y)
  folds <- sample(rep(seq_len(K), length.out = n))
  y_res <- numeric(n)
  t_res <- numeric(n)

  for (k in seq_len(K)) {
    train_idx <- which(folds != k)
    test_idx <- which(folds == k)
    W_train <- W[train_idx, , drop = FALSE]
    W_test <- W[test_idx, , drop = FALSE]

    # Nuisance model for outcome
    rf_y <- ranger::ranger(
      y ~ ., data = cbind(y = Y[train_idx], W_train),
      num.trees = num_trees, verbose = FALSE
    )
    y_res[test_idx] <- Y[test_idx] -
      stats::predict(rf_y, data = W_test)$predictions

    # Nuisance model for treatment
    rf_t <- ranger::ranger(
      y ~ ., data = cbind(y = T_var[train_idx], W_train),
      num.trees = num_trees, verbose = FALSE
    )
    t_res[test_idx] <- T_var[test_idx] -
      stats::predict(rf_t, data = W_test)$predictions
  }

  ate <- sum(t_res * y_res) / sum(t_res^2)
  residuals <- y_res - ate * t_res
  se <- sqrt(
    mean(residuals^2 * t_res^2) / (mean(t_res^2)^2) / n
  )
  p_value <- 2 * stats::pnorm(-abs(ate / se))
  list(ate = ate, se = se, p_value = p_value)
}
