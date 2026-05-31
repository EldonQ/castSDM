#' Estimate Spatially Heterogeneous Treatment Effects (CATE)
#'
#' Uses causal forests ([grf::causal_forest()]) to estimate Conditional
#' Average Treatment Effects at each spatial location. Reveals how
#' environmental impacts on species presence vary across geographic space.
#'
#' When a DAG is provided, confounders for each treatment variable are
#' identified as the remaining Markov Blanket variables (DAG-guided).
#' Without a DAG, all other selected variables are used as confounders.
#'
#' @param data A `data.frame` with `presence`, environmental variables, and
#'   optionally `lon`, `lat` coordinates.
#' @param variables Character vector of treatment variables to estimate CATE
#'   for. Default `NULL` uses top variables from the screen (by role priority:
#'   parents first).
#' @param dag A [cast_dag] object (used for DAG-guided confounder
#'   selection). Default `NULL`.
#' @param screen A `cast_select` object from [cast_select()] (used to select top variables if
#'   `variables` is `NULL`). Default `NULL`.
#' @param response Character. Response column name. Default `"presence"`.
#' @param top_n Integer. Number of top variables if `variables` is `NULL`.
#'   Default `3`.
#' @param n_trees Integer. Number of causal forest trees. Default `1000`.
#' @param predict_data Optional `data.frame` for out-of-sample CATE
#'   prediction. Default `NULL`.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A `cast_cate` object containing an `effects` data.frame with
#'   columns: `lon`, `lat`, `variable`, `cate`.
#'
#' @references
#' Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random
#' forests. *The Annals of Statistics*, 47(2), 1148-1178.
#'
#' @seealso [grf::causal_forest()], [cast_select()]
#'
#' @export
cast_cate <- function(data,
                      variables = NULL,
                      dag = NULL,
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
  for (col in names(X_full)) X_full[[col]] <- as.numeric(X_full[[col]])
  X_full[is.na(X_full)] <- 0

  # Determine which variables to estimate CATE for
  if (is.null(variables)) {
    if (!is.null(screen) && !is.null(screen$roles)) {
      # Prioritise parents > children > co_parents > predictive
      role_priority <- c("parent", "child", "co_parent", "predictive")
      roles_df <- screen$roles
      roles_df$priority <- match(roles_df$role, role_priority)
      roles_df <- roles_df[order(roles_df$priority), ]
      cate_vars <- roles_df$variable[seq_len(
        min(top_n, nrow(roles_df))
      )]
    } else if (!is.null(screen)) {
      cate_vars <- screen$selected[seq_len(
        min(top_n, length(screen$selected))
      )]
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

  # Extract Markov Blanket for DAG-guided confounder selection
  mb <- NULL
  if (!is.null(dag) && !is.null(dag$response_node)) {
    mb <- extract_markov_blanket(dag, dag$response_node)
  }

  # Prediction data
  pred_X <- if (!is.null(predict_data)) {
    pd <- as.data.frame(predict_data[, env_vars, drop = FALSE])
    for (col in names(pd)) pd[[col]] <- as.numeric(pd[[col]])
    pd[is.na(pd)] <- 0
    pd
  } else {
    X_full
  }
  pred_source <- predict_data %||% data

  has_coords <- all(c("lon", "lat") %in% names(pred_source))

  # Estimate CATE for each variable
  cate_list <- list()
  for (cv in cate_vars) {
    if (verbose) cli::cli_inform("  Causal forest: {.val {cv}}...")

    T_cont <- X_full[[cv]]

    # DAG-guided confounder selection: use MB(presence) minus treatment var
    if (!is.null(mb) && length(mb$all) > 0) {
      confounders <- intersect(setdiff(mb$all, cv), env_vars)
      if (length(confounders) < 2) {
        # Fallback: use all env vars minus treatment
        confounders <- setdiff(env_vars, cv)
      }
    } else {
      confounders <- setdiff(env_vars, cv)
    }

    W_covs <- as.matrix(X_full[, confounders, drop = FALSE])

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
        pred_X[, confounders, drop = FALSE]
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
