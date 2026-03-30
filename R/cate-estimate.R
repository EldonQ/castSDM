#' Estimate Spatially Heterogeneous Treatment Effects (CATE)
#'
#' Uses causal forests ([grf::causal_forest()]) to estimate Conditional
#' Average Treatment Effects at each spatial location. Reveals how
#' environmental impacts on species presence vary across geographic space.
#'
#' @param data A `data.frame` with `presence`, environmental variables, and
#'   optionally `lon`, `lat` coordinates.
#' @param variables Character vector of treatment variables to estimate CATE
#'   for. Default `NULL` uses top significant variables from ATE.
#' @param ate A [cast_ate] object (used to select top variables if
#'   `variables` is `NULL`). Default `NULL`.
#' @param screen A [cast_screen] object (used to select top variables if
#'   `variables` is `NULL`). Default `NULL`.
#' @param response Character. Response column name. Default `"presence"`.
#' @param top_n Integer. Number of top variables if `variables` is `NULL`.
#'   Default `3`.
#' @param n_trees Integer. Number of causal forest trees. Default `1000`.
#' @param predict_data Optional `data.frame` for out-of-sample CATE
#'   prediction. If `NULL`, predicts on the training data. Default `NULL`.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A `cast_cate` object containing a `effects` data.frame with
#'   columns: `lon`, `lat`, `variable`, `cate`.
#'
#' @details
#' For each treatment variable, a causal forest partitions the feature space
#' to estimate location-specific effects:
#'
#' \deqn{\tau(x) = E[Y(1) - Y(0) | X = x]}
#'
#' This reveals spatial heterogeneity in how environmental drivers affect
#' species presence -- e.g., "temperature sensitivity is stronger at range
#' margins."
#'
#' @references
#' Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random
#' forests. *The Annals of Statistics*, 47(2), 1148-1178.
#'
#' @seealso [grf::causal_forest()], [cast_ate()]
#'
#' @export
cast_cate <- function(data,
                      variables = NULL,
                      ate = NULL,
                      screen = NULL,
                      response = "presence",
                      top_n = 3L,
                      n_trees = 1000L,
                      predict_data = NULL,
                      seed = NULL,
                      verbose = TRUE) {
  check_suggested("grf", "for causal forest estimation")

  env_vars <- get_env_vars(data, response = response)
  Y <- data[[response]]
  X_full <- as.data.frame(data[, env_vars, drop = FALSE])
  X_full[is.na(X_full)] <- 0

  # Determine which variables to estimate CATE for
  if (is.null(variables)) {
    if (!is.null(screen) && !is.null(ate)) {
      # Intersection: screened + significant
      sig_vars <- ate$estimates$variable[
        ate$estimates$significant == TRUE
      ]
      cate_vars <- intersect(screen$selected, sig_vars)
      if (length(cate_vars) == 0) {
        cate_vars <- screen$selected[seq_len(
          min(top_n, length(screen$selected))
        )]
      } else {
        cate_vars <- cate_vars[seq_len(min(top_n, length(cate_vars)))]
      }
    } else if (!is.null(ate)) {
      sig <- ate$estimates[ate$estimates$significant == TRUE, ]
      sig <- sig[order(abs(sig$coef), decreasing = TRUE), ]
      cate_vars <- sig$variable[seq_len(min(top_n, nrow(sig)))]
    } else {
      cate_vars <- env_vars[seq_len(min(top_n, length(env_vars)))]
    }
  } else {
    cate_vars <- variables
  }

  if (length(cate_vars) == 0) {
    cli::cli_abort("No variables available for CATE estimation.")
  }
  if (verbose) {
    cli::cli_inform("Estimating CATE for: {.val {cate_vars}}")
  }

  # Prediction data
  pred_X <- if (!is.null(predict_data)) {
    pd <- as.data.frame(predict_data[, env_vars, drop = FALSE])
    pd[is.na(pd)] <- 0
    pd
  } else {
    X_full
  }
  pred_source <- predict_data %||% data

  # Extract coords from prediction data
  has_coords <- all(c("lon", "lat") %in% names(pred_source))

  # Estimate CATE for each variable
  cate_list <- list()
  for (cv in cate_vars) {
    if (verbose) cli::cli_inform("  Causal forest: {.val {cv}}...")

    T_cont <- X_full[[cv]]
    W_covs <- as.matrix(
      X_full[, setdiff(env_vars, cv), drop = FALSE]
    )

    cf <- tryCatch({
      if (!is.null(seed)) set.seed(seed)
      grf::causal_forest(
        X = W_covs, Y = Y, W = T_cont,
        num.trees = n_trees, seed = seed %||% 42L
      )
    }, error = function(e) {
      if (verbose) cli::cli_warn("  Failed for {cv}: {e$message}")
      NULL
    })

    if (!is.null(cf)) {
      W_pred <- as.matrix(
        pred_X[, setdiff(env_vars, cv), drop = FALSE]
      )
      cate_preds <- stats::predict(
        cf, W_pred, estimate.variance = FALSE
      )$predictions

      df <- data.frame(
        variable = cv,
        cate = as.numeric(cate_preds),
        stringsAsFactors = FALSE
      )
      if (has_coords) {
        df$lon <- pred_source$lon
        df$lat <- pred_source$lat
      }
      cate_list[[cv]] <- df
    }
  }

  if (length(cate_list) == 0) {
    cli::cli_abort("All CATE estimations failed.")
  }

  effects <- do.call(rbind, cate_list)
  rownames(effects) <- NULL

  new_cast_cate(
    effects = effects,
    variables = names(cate_list),
    n_trees = n_trees
  )
}
