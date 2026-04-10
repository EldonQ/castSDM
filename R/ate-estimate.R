#' Estimate Average Treatment Effects via Double Machine Learning
#'
#' For each environmental variable, estimates the Average Treatment Effect
#' (ATE) on species presence using Double/Debiased Machine Learning
#' (Chernozhukov et al. 2018). To handle non-linear ecological responses,
#' binarization is repeated at multiple quantile thresholds and the results
#' are averaged, greatly reducing the information loss of a single median cut.
#' Random Forests are used as nuisance functions with K-fold cross-fitting.
#'
#' @param data A `data.frame` with `presence` column and environmental
#'   variables.
#' @param response Character. Name of response column. Default `"presence"`.
#' @param variables Character vector of variable names to test. Default
#'   `NULL` (all numeric columns except `response`, `lon`, `lat`).
#' @param K Integer. Number of cross-fitting folds. Default `5`. The original
#'   default of 2 produced high-variance residual estimates; 5 is the standard
#'   recommendation from the DML literature.
#' @param num_trees Integer. Number of trees for nuisance Random Forests.
#'   Default `300`.
#' @param alpha Numeric. Significance level **after Bonferroni correction**.
#'   Default `0.05`.
#' @param quantile_cuts Numeric vector. Quantile thresholds used to binarize
#'   each treatment variable. ATEs from all cuts are averaged before
#'   significance testing. Default `c(0.25, 0.50, 0.75)`. Set to `0.5` to
#'   recover the original single-median behaviour.
#' @param bonferroni Logical. Apply Bonferroni correction for the number of
#'   variables tested. Default `TRUE`.
#' @param parallel Logical. If `TRUE`, estimate ATEs for each variable in
#'   parallel using \pkg{future.apply}. The user must set up a
#'   [future::plan()] beforehand (e.g.,
#'   `future::plan(future::multisession, workers = 4)`). Default `FALSE`.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A [cast_ate] object.
#'
#' @details
#' ## Multi-quantile DML procedure
#' For each variable \eqn{X} and each quantile threshold \eqn{q}:
#' 1. Binarize: \eqn{T^{(q)} = \mathbf{1}[X > Q_q(X)]}.
#' 2. K-fold cross-fitting: RF nuisances for \eqn{E[Y|W]} and
#'    \eqn{E[T^{(q)}|W]}; compute residuals on held-out folds.
#' 3. Estimate \eqn{\widehat{ATE}^{(q)}} via the DML moment condition.
#'
#' Final coefficient: \eqn{\bar{ATE} = \text{mean}_q(\widehat{ATE}^{(q)})}.
#' Final SE: \eqn{\bar{SE} = \text{mean}_q(\widehat{SE}^{(q)})}.
#' This captures both lower-tail and upper-tail effects, substantially
#' improving detection of non-linear (e.g. Gaussian, threshold) responses.
#'
#' ## Bonferroni correction
#' When `bonferroni = TRUE` the significance threshold is divided by the
#' number of variables tested, controlling the family-wise error rate.
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
                     response      = "presence",
                     variables     = NULL,
                     K             = 5L,
                     num_trees     = 300L,
                     alpha         = 0.05,
                     quantile_cuts = c(0.25, 0.50, 0.75),
                     bonferroni    = TRUE,
                     parallel      = FALSE,
                     seed          = NULL,
                     verbose       = TRUE) {
  check_suggested("ranger", "for DML nuisance models")

  env_vars <- variables %||% get_env_vars(data, response = response)

  Y      <- data[[response]]
  X_full <- as.data.frame(data[, env_vars, drop = FALSE])
  X_full[is.na(X_full)] <- 0

  # Bonferroni-corrected significance threshold
  alpha_adj <- if (bonferroni && length(env_vars) > 1) {
    alpha / length(env_vars)
  } else {
    alpha
  }

  if (verbose) {
    cli::cli_inform(c(
      "Estimating ATE for {length(env_vars)} variables",
      "i" = "K={K} folds, {length(quantile_cuts)} quantile cut(s), alpha_adj={round(alpha_adj,4)}"
    ))
  }

  # Worker function for a single variable
  ate_one_var <- function(v, Y, X_full, env_vars, quantile_cuts, K,
                          num_trees, alpha_adj, seed) {
    tryCatch({
      W <- X_full[, setdiff(env_vars, v), drop = FALSE]

      ate_vec <- numeric(length(quantile_cuts))
      se_vec  <- numeric(length(quantile_cuts))

      for (qi in seq_along(quantile_cuts)) {
        q_thr  <- stats::quantile(X_full[[v]], probs = quantile_cuts[qi],
                                  na.rm = TRUE)
        T_bin  <- as.integer(X_full[[v]] > q_thr)

        if (length(unique(T_bin)) < 2L) {
          ate_vec[qi] <- NA_real_
          se_vec[qi]  <- NA_real_
          next
        }

        if (!is.null(seed)) set.seed(seed + qi)
        res <- dml_ate_internal(
          Y = Y, T_var = T_bin, W = W,
          K = K, num_trees = num_trees
        )
        ate_vec[qi] <- res$ate
        se_vec[qi]  <- res$se
      }

      valid    <- !is.na(ate_vec)
      coef_avg <- if (any(valid)) mean(ate_vec[valid]) else NA_real_
      se_avg   <- if (any(valid)) mean(se_vec[valid])  else NA_real_

      p_val <- if (is.finite(coef_avg) && is.finite(se_avg) && se_avg > 0) {
        2 * stats::pnorm(-abs(coef_avg / se_avg))
      } else {
        NA_real_
      }

      data.frame(
        variable    = v,
        coef        = coef_avg,
        se          = se_avg,
        p_value     = p_val,
        significant = !is.na(p_val) && p_val < alpha_adj,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        variable = v, coef = NA_real_, se = NA_real_,
        p_value  = NA_real_, significant = FALSE,
        stringsAsFactors = FALSE
      )
    })
  }

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    if (verbose) cli::cli_inform("Running ATE in parallel ({length(env_vars)} variables)...")
    results <- future.apply::future_lapply(
      env_vars,
      ate_one_var,
      Y = Y, X_full = X_full, env_vars = env_vars,
      quantile_cuts = quantile_cuts, K = K,
      num_trees = num_trees, alpha_adj = alpha_adj, seed = seed,
      future.seed = TRUE
    )
  } else {
    results <- vector("list", length(env_vars))
    for (i in seq_along(env_vars)) {
      results[[i]] <- ate_one_var(
        env_vars[i], Y = Y, X_full = X_full, env_vars = env_vars,
        quantile_cuts = quantile_cuts, K = K,
        num_trees = num_trees, alpha_adj = alpha_adj, seed = seed
      )
    }
  }

  estimates <- do.call(rbind, results)
  n_sig <- sum(estimates$significant, na.rm = TRUE)
  if (verbose) {
    cli::cli_inform(
      "ATE: {n_sig}/{length(env_vars)} significant (Bonferroni p < {round(alpha_adj,4)})."
    )
  }

  new_cast_ate(estimates = estimates, K = K, alpha = alpha_adj)
}


#' Double Machine Learning ATE (Internal)
#'
#' Core DML algorithm for a single binary treatment variable.
#'
#' @param Y Numeric response vector.
#' @param T_var Binary treatment vector.
#' @param W data.frame of confounders.
#' @param K Number of cross-fitting folds.
#' @param num_trees Number of RF trees.
#' @return List with `ate`, `se`, `p_value`.
#' @keywords internal
#' @noRd
dml_ate_internal <- function(Y, T_var, W, K = 5, num_trees = 300) {
  n     <- length(Y)
  folds <- sample(rep(seq_len(K), length.out = n))
  y_res <- numeric(n)
  t_res <- numeric(n)

  for (k in seq_len(K)) {
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    W_train   <- W[train_idx, , drop = FALSE]
    W_test    <- W[test_idx,  , drop = FALSE]

    # Nuisance model for outcome E[Y | W]
    rf_y <- ranger::ranger(
      y ~ ., data = cbind(y = Y[train_idx], W_train),
      num.trees = num_trees, verbose = FALSE
    )
    y_res[test_idx] <- Y[test_idx] -
      stats::predict(rf_y, data = W_test)$predictions

    # Nuisance model for treatment E[T | W]
    rf_t <- ranger::ranger(
      y ~ ., data = cbind(y = T_var[train_idx], W_train),
      num.trees = num_trees, verbose = FALSE
    )
    t_res[test_idx] <- T_var[test_idx] -
      stats::predict(rf_t, data = W_test)$predictions
  }

  # DML moment estimator
  denom <- sum(t_res^2)
  if (denom < 1e-10) return(list(ate = NA_real_, se = NA_real_, p_value = NA_real_))

  ate       <- sum(t_res * y_res) / denom
  residuals <- y_res - ate * t_res
  se        <- sqrt(
    mean(residuals^2 * t_res^2) / (mean(t_res^2)^2) / n
  )
  p_value <- if (is.finite(se) && se > 0) {
    2 * stats::pnorm(-abs(ate / se))
  } else {
    NA_real_
  }
  list(ate = ate, se = se, p_value = p_value)
}
