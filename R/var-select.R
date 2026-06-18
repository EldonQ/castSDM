# Internal helpers for invariant causal screening.

#' @keywords internal
#' @noRd
fast_auc01 <- function(y, score) {
  y <- as.integer(y)
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' @keywords internal
#' @noRd
make_invariance_blocks <- function(data, block_var = NULL, n_blocks = 5L) {
  if (!is.null(block_var) && block_var %in% names(data)) {
    return(as.factor(data[[block_var]]))
  }
  if (all(c("lon", "lat") %in% names(data))) {
    n_blocks <- max(2L, as.integer(n_blocks))
    nx <- max(2L, ceiling(sqrt(n_blocks)))
    ny <- max(2L, ceiling(n_blocks / nx))
    xq <- unique(stats::quantile(data$lon, probs = seq(0, 1, length.out = nx + 1),
                                 na.rm = TRUE, type = 8))
    yq <- unique(stats::quantile(data$lat, probs = seq(0, 1, length.out = ny + 1),
                                 na.rm = TRUE, type = 8))
    if (length(xq) > 2L && length(yq) > 2L) {
      xb <- cut(data$lon, breaks = xq, include.lowest = TRUE, labels = FALSE)
      yb <- cut(data$lat, breaks = yq, include.lowest = TRUE, labels = FALSE)
      return(interaction(xb, yb, drop = TRUE))
    }
  }
  as.factor(rep(1L, nrow(data)))
}

#' @keywords internal
#' @noRd
compute_invariance_metrics <- function(data, env_vars, response, block_var = NULL,
                                       n_blocks = 5L) {
  y <- as.integer(data[[response]])
  blocks <- make_invariance_blocks(data, block_var = block_var, n_blocks = n_blocks)
  block_levels <- levels(blocks)
  out <- lapply(env_vars, function(v) {
    signs <- numeric(0)
    aucs <- numeric(0)
    effects <- numeric(0)
    for (b in block_levels) {
      idx <- which(blocks == b)
      if (length(idx) < 20L || length(unique(y[idx])) < 2L) next
      x <- suppressWarnings(as.numeric(data[[v]][idx]))
      yy <- y[idx]
      ok <- is.finite(x) & !is.na(yy)
      if (sum(ok) < 20L || stats::sd(x[ok]) == 0) next
      xz <- as.numeric(scale(x[ok]))
      fit <- tryCatch(
        suppressWarnings(stats::glm(yy[ok] ~ xz, family = stats::binomial())),
        error = function(e) NULL
      )
      if (is.null(fit) || length(stats::coef(fit)) < 2L) next
      beta <- unname(stats::coef(fit)[2])
      p <- stats::predict(fit, type = "response")
      auc <- fast_auc01(yy[ok], p)
      if (!is.finite(beta)) next
      signs <- c(signs, sign(beta))
      effects <- c(effects, abs(beta))
      aucs <- c(aucs, auc)
    }
    nz <- signs[signs != 0]
    sign_consistency <- if (length(nz) == 0L) {
      0
    } else {
      max(mean(nz > 0), mean(nz < 0))
    }
    effect_stability <- if (length(effects) <= 1L || mean(effects) <= 0) {
      if (length(effects) == 1L) 1 else 0
    } else {
      max(0, 1 - stats::sd(effects) / (mean(effects) + 1e-8))
    }
    auc_stability <- if (length(aucs) <= 1L) {
      if (length(aucs) == 1L) 1 else 0
    } else {
      max(0, 1 - stats::sd(aucs, na.rm = TRUE) / (abs(mean(aucs, na.rm = TRUE) - 0.5) + 1e-4))
    }
    data.frame(
      variable = v,
      n_blocks_used = length(effects),
      sign_consistency = sign_consistency,
      effect_stability = effect_stability,
      auc_stability = auc_stability,
      mean_block_auc = if (length(aucs)) mean(aucs, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' @keywords internal
#' @noRd
correlation_clusters <- function(X, threshold = 0.8) {
  vars <- names(X)
  if (length(vars) <= 1L) return(stats::setNames(seq_along(vars), vars))
  cors <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
  cors[!is.finite(cors)] <- 0
  adj <- abs(cors) >= threshold
  diag(adj) <- TRUE
  seen <- stats::setNames(rep(FALSE, length(vars)), vars)
  cluster <- stats::setNames(rep(NA_integer_, length(vars)), vars)
  cid <- 0L
  for (v in vars) {
    if (seen[[v]]) next
    cid <- cid + 1L
    queue <- v
    seen[[v]] <- TRUE
    while (length(queue)) {
      cur <- queue[[1]]
      queue <- queue[-1]
      cluster[[cur]] <- cid
      nb <- vars[adj[cur, vars]]
      nb <- nb[!seen[nb]]
      if (length(nb)) {
        seen[nb] <- TRUE
        queue <- c(queue, nb)
      }
    }
  }
  cluster
}

#' @keywords internal
#' @noRd
greedy_invariant_select <- function(scores, X, min_vars, max_vars = NULL,
                                    cor_threshold = 0.8,
                                    max_per_cluster = 1L,
                                    score_threshold = NULL) {
  scores <- scores[order(-scores$invariant_score, -scores$importance), , drop = FALSE]
  if (is.null(score_threshold) || !is.finite(score_threshold)) {
    sorted <- scores$invariant_score
    sorted <- sorted[is.finite(sorted)]
    if (length(sorted) <= min_vars) {
      score_threshold <- -Inf
    } else {
      gaps <- sorted[-length(sorted)] - sorted[-1L]
      valid_idx <- which(seq_along(gaps) >= min_vars)
      cut_at <- if (length(valid_idx)) valid_idx[which.max(gaps[valid_idx])] else min_vars
      score_threshold <- sorted[cut_at]
    }
  }
  max_vars <- if (is.null(max_vars) || !is.finite(max_vars)) {
    Inf
  } else {
    max(min_vars, as.integer(max_vars))
  }
  clusters <- correlation_clusters(X[, scores$variable, drop = FALSE], cor_threshold)
  scores$cor_cluster <- as.integer(clusters[scores$variable])
  selected <- character(0)
  selected_cluster_counts <- list()
  redundant <- character(0)
  for (v in scores$variable) {
    if (length(selected) >= max_vars) break
    score_v <- scores$invariant_score[scores$variable == v][1]
    if (length(selected) >= min_vars && is.finite(score_v) && score_v < score_threshold) next
    cl <- as.character(scores$cor_cluster[scores$variable == v][1])
    n_cl <- if (!is.null(selected_cluster_counts[[cl]])) selected_cluster_counts[[cl]] else 0L
    if (n_cl >= max_per_cluster && length(selected) >= min_vars) {
      redundant <- c(redundant, v)
      next
    }
    if (n_cl >= max_per_cluster) {
      redundant <- c(redundant, v)
      next
    }
    selected <- c(selected, v)
    selected_cluster_counts[[cl]] <- n_cl + 1L
  }
  if (length(selected) < min_vars) {
    add <- setdiff(scores$variable, selected)
    selected <- unique(c(selected, add[seq_len(min(min_vars - length(selected), length(add)))]))
  }
  list(
    selected = selected,
    clusters = clusters,
    redundant = setdiff(redundant, selected),
    score_threshold = score_threshold
  )
}

#' Causal-Aware Variable Selection
#'
#' Selects predictor variables with a response-focused causal screening
#' workflow. The default method, `"invariant_screen"`, combines RF permutation
#' importance, bootstrap selection frequency, spatial-block effect-direction
#' consistency, and correlation-cluster redundancy control. The legacy
#' `"mb_rf"` path preserves Markov Blanket + RF screening for audit and
#' comparison.
#'
#' @param dag A [cast_dag] object, ideally learned with
#'   `include_response = TRUE`.
#' @param data A `data.frame` with species occurrence and environmental
#'   variables.
#' @param response Character. Name of the response column. Default
#'   `"presence"`.
#' @param method Character. `"invariant_screen"` (default), `"mb_rf"` legacy
#'   Markov Blanket + RF fusion, or `"rf"` pure RF/stability screening.
#' @param num_trees Integer. Number of trees for the RF importance step.
#' @param min_vars Integer. Minimum variables to retain.
#' @param min_fraction Numeric in `[0, 1]`. Minimum fraction of candidate
#'   variables for the legacy `"mb_rf"` path.
#' @param seed Integer or `NULL`. Random seed.
#' @param verbose Logical. Print progress.
#' @param stability_reps Integer. Bootstrap repetitions for lightweight
#'   selection stability diagnostics.
#' @param stability_threshold Numeric in `[0, 1]`. Minimum bootstrap
#'   selection frequency.
#' @param max_vars Optional integer safety ceiling for variables retained by
#'   `"invariant_screen"` and `"rf"`. Default `NULL` means no fixed cap; the
#'   number of selected variables is determined by the adaptive score break and
#'   redundancy control.
#' @param score_threshold Optional numeric threshold on the invariant score.
#'   Default `NULL` estimates an adaptive threshold from the largest score gap
#'   after `min_vars`.
#' @param cor_threshold Numeric. Absolute correlation above which predictors
#'   are treated as redundant proxies.
#' @param max_per_cluster Integer. Maximum selected variables allowed from one
#'   correlation cluster.
#' @param block_var Optional column name defining spatial/environmental
#'   blocks. If `NULL`, blocks are built from `lon`/`lat` quantiles when
#'   available.
#' @param n_blocks Integer target number of blocks when `block_var = NULL`.
#'
#' @return A `cast_select` object with `selected`, `scores`, and `roles`.
#' @seealso [cast_dag()], [cast_fit()]
#' @export
cast_select <- function(dag,
                        data,
                        response = "presence",
                        method = c("invariant_screen", "mb_rf", "rf"),
                        num_trees = 300L,
                        min_vars = 5L,
                        min_fraction = 0.3,
                        seed = NULL,
                        verbose = TRUE,
                        stability_reps = 0L,
                        stability_threshold = 0.6,
                        max_vars = NULL,
                        cor_threshold = 0.8,
                        max_per_cluster = 1L,
                        score_threshold = NULL,
                        block_var = NULL,
                        n_blocks = 5L) {
  check_suggested("ranger", "for RF importance")
  method <- match.arg(method)

  # --- Determine candidate env vars (exclude response, lon, lat) -----------
  response_node <- dag$response_node
  env_vars <- setdiff(dag$nodes, c(response, "lon", "lat"))
  if (length(env_vars) < 3) {
    cli::cli_abort("Need at least 3 environmental variables for selection.")
  }

  # --- Extract Markov Blanket (if response was in DAG) ---------------------
  mb <- NULL
  if (!is.null(response_node)) {
    mb <- response_markov_blanket(dag, response_node, env_vars)
    if (verbose && length(mb$all) > 0) {
      cli::cli_inform(
        "Markov Blanket of {.val {response_node}}: {length(mb$all)} variables ({length(mb$parents)} direct MB, {length(mb$children) + length(mb$co_parents)} associated MB)."
      )
    }
  }

  # --- RF permutation importance -------------------------------------------
  Y <- data[[response]]
  X <- as.data.frame(data[, env_vars, drop = FALSE])
  for (col in names(X)) X[[col]] <- as.numeric(X[[col]])
  X[is.na(X)] <- 0
  if (!is.null(seed)) set.seed(seed)
  rf_imp <- ranger::ranger(
    y ~ .,
    data = cbind(y = as.factor(Y), X),
    num.trees = num_trees,
    importance = "permutation",
    verbose = FALSE
  )$variable.importance

  # --- Build scoring data.frame --------------------------------------------
  deg_df <- compute_edge_degrees(dag$edges, env_vars)
  scores_df <- deg_df
  scores_df$importance <- rf_imp[scores_df$variable]
  scores_df$importance[is.na(scores_df$importance)] <- 0

  # Mark MB membership
  scores_df$in_mb <- scores_df$variable %in% (mb$all %||% character())
  scores_df$mb_role <- vapply(scores_df$variable, function(v) {
    if (v %in% (mb$parents %||% character())) return("mb_direct")
    if (v %in% (mb$children %||% character())) return("mb_associated")
    if (v %in% (mb$co_parents %||% character())) return("mb_associated")
    "none"
  }, character(1))

  # Normalize importance to [0, 1]
  scores_df$imp_norm <- normalize01(scores_df$importance)

  # --- Selection logic -----------------------------------------------------
  min_keep <- max(min_vars, ceiling(length(env_vars) * min_fraction))
  select_by_importance <- function(scores, min_keep) {
    scores <- scores[order(-scores$importance), ]
    imp_sorted <- scores$importance

    if (length(imp_sorted) <= min_keep) {
      return(scores$variable)
    }

    gaps <- imp_sorted[-length(imp_sorted)] - imp_sorted[-1L]
    valid_idx <- which(seq_along(gaps) >= min_keep &
                         seq_along(gaps) < length(imp_sorted))
    if (length(valid_idx) == 0) {
      scores$variable[seq_len(min_keep)]
    } else {
      cut_at <- valid_idx[which.max(gaps[valid_idx])]
      scores$variable[seq_len(cut_at)]
    }
  }

  estimate_stability <- function(X, Y, scores, reps, min_keep) {
    reps <- as.integer(reps)
    if (reps <= 0L) return(rep(NA_real_, nrow(scores)))
    counts <- stats::setNames(rep(0L, nrow(scores)), scores$variable)
    trees <- max(50L, min(150L, as.integer(ceiling(num_trees / 3))))
    n <- nrow(X)
    for (ii in seq_len(reps)) {
      if (!is.null(seed)) set.seed(seed + 1000L + ii)
      idx <- sample.int(n, size = n, replace = TRUE)
      imp <- tryCatch(
        ranger::ranger(
          y ~ .,
          data = cbind(y = as.factor(Y[idx]), X[idx, , drop = FALSE]),
          num.trees = trees,
          importance = "permutation",
          verbose = FALSE
        )$variable.importance,
        error = function(e) NULL
      )
      if (is.null(imp)) next
      rep_scores <- scores
      rep_scores$importance <- imp[rep_scores$variable]
      rep_scores$importance[is.na(rep_scores$importance)] <- 0
      keep <- select_by_importance(rep_scores, min_keep)
      counts[keep] <- counts[keep] + 1L
    }
    as.numeric(counts[scores$variable]) / reps
  }

  scores_df$selection_frequency <- estimate_stability(
    X, Y, scores_df, stability_reps, min_keep
  )
  scores_df$stable <- ifelse(
    is.na(scores_df$selection_frequency),
    NA,
    scores_df$selection_frequency >= stability_threshold
  )

  mb_selective <- TRUE

  if (method %in% c("invariant_screen", "rf")) {
    min_keep_inv <- max(1L, as.integer(min_vars))
    X_numeric <- X
    for (col in names(X_numeric)) X_numeric[[col]] <- as.numeric(X_numeric[[col]])

    inv <- if (identical(method, "invariant_screen")) {
      compute_invariance_metrics(
        data = data,
        env_vars = env_vars,
        response = response,
        block_var = block_var,
        n_blocks = n_blocks
      )
    } else {
      data.frame(
        variable = env_vars,
        n_blocks_used = NA_integer_,
        sign_consistency = 0,
        effect_stability = 0,
        auc_stability = 0,
        mean_block_auc = NA_real_,
        stringsAsFactors = FALSE
      )
    }

    scores_df <- merge(scores_df, inv, by = "variable", all.x = TRUE, sort = FALSE)
    for (nm in c("sign_consistency", "effect_stability", "auc_stability")) {
      scores_df[[nm]][is.na(scores_df[[nm]])] <- 0
    }
    stability_component <- scores_df$selection_frequency
    if (all(is.na(stability_component))) stability_component <- scores_df$imp_norm
    stability_component[is.na(stability_component)] <- 0

    scores_df$invariant_score <- if (identical(method, "invariant_screen")) {
      0.35 * scores_df$imp_norm +
        0.20 * stability_component +
        0.25 * scores_df$sign_consistency +
        0.15 * scores_df$effect_stability +
        0.05 * scores_df$auc_stability
    } else {
      0.75 * scores_df$imp_norm + 0.25 * stability_component
    }

    pick <- greedy_invariant_select(
      scores = scores_df,
      X = X_numeric,
      min_vars = min_keep_inv,
      max_vars = max_vars,
      cor_threshold = cor_threshold,
      max_per_cluster = max_per_cluster,
      score_threshold = score_threshold
    )
    selected <- pick$selected
    scores_df$cor_cluster <- as.integer(pick$clusters[scores_df$variable])
    scores_df$adaptive_score_threshold <- pick$score_threshold
    scores_df$redundant_with_selected <- FALSE
    if (length(selected)) {
      cmat <- suppressWarnings(stats::cor(X_numeric[, env_vars, drop = FALSE],
                                          use = "pairwise.complete.obs"))
      cmat[!is.finite(cmat)] <- 0
      for (v in setdiff(env_vars, selected)) {
        scores_df$redundant_with_selected[scores_df$variable == v] <-
          any(abs(cmat[v, selected]) >= cor_threshold, na.rm = TRUE)
      }
    }

    roles_df <- data.frame(
      variable = selected,
      role = vapply(selected, function(v) {
        row <- scores_df[scores_df$variable == v, , drop = FALSE]
        if (identical(method, "rf")) return("stable_predictive")
        if (isTRUE(row$sign_consistency[1] >= 0.75) &&
            isTRUE(row$effect_stability[1] >= 0.4) &&
            isTRUE(row$n_blocks_used[1] >= 2)) {
          return("invariant_driver")
        }
        if (isTRUE((row$selection_frequency[1] %||% 0) >= stability_threshold) ||
            isTRUE(row$imp_norm[1] >= 0.5)) {
          return("stable_predictive")
        }
        "predictive_rescue"
      }, character(1)),
      stringsAsFactors = FALSE
    )
    roles_df$selection_frequency <- scores_df$selection_frequency[
      match(roles_df$variable, scores_df$variable)
    ]
    roles_df$causal_role <- roles_df$role
    scores_df$causal_role <- "unstable_rejected"
    scores_df$causal_role[scores_df$redundant_with_selected] <- "redundant_proxy"
    scores_df$causal_role[match(roles_df$variable, scores_df$variable)] <- roles_df$causal_role
    scores_df$selected <- scores_df$variable %in% selected
    scores_df$screening_method <- method
    roles_df$screening_method <- method

    scores_df <- scores_df[order(!scores_df$selected, -scores_df$invariant_score,
                                 -scores_df$importance), ]
    rownames(scores_df) <- NULL
    roles_df$score <- scores_df$invariant_score[match(roles_df$variable, scores_df$variable)]
    roles_df <- roles_df[order(-roles_df$score), ]
    roles_df$score <- NULL
    rownames(roles_df) <- NULL

    if (verbose) {
      n_inv <- sum(roles_df$role == "invariant_driver")
      n_stable <- sum(roles_df$role == "stable_predictive")
      n_rescue <- sum(roles_df$role == "predictive_rescue")
      cli::cli_inform(
        "Invariant screen: {length(selected)}/{length(env_vars)} variables ({n_inv} invariant drivers, {n_stable} stable predictive, {n_rescue} predictive rescue; cor threshold={cor_threshold})."
      )
    }

    return(new_cast_select(
      selected = selected,
      scores = scores_df,
      roles = roles_df
    ))
  }

  if (!is.null(mb) && length(mb$all) > 0) {
    # -- Fusion: MB vars auto-selected + importance-gated non-MB vars --
    mb_vars <- intersect(mb$all, env_vars)
    mb_dense <- length(mb_vars) >= ceiling(0.8 * length(env_vars))

    if (isTRUE(mb_dense)) {
      if (verbose) {
        cli::cli_warn(c(
          "Dense Markov Blanket: {length(mb_vars)}/{length(env_vars)} predictors are in MB.",
          "i" = "Treating MB membership as non-selective and falling back to RF importance screening."
        ))
      }
      mb_selective <- FALSE
      if (stability_reps > 0L && any(scores_df$stable %in% TRUE)) {
        selected <- scores_df$variable[scores_df$stable %in% TRUE]
        if (length(selected) < min_keep) {
          remaining <- setdiff(scores_df$variable[order(-scores_df$importance)], selected)
          selected <- unique(c(
            selected,
            remaining[seq_len(min(min_keep - length(selected), length(remaining)))]
          ))
        }
      } else {
        selected <- select_by_importance(scores_df, min_keep)
      }
    } else {
      mb_importance <- scores_df$importance[scores_df$variable %in% mb_vars]
      imp_threshold <- if (length(mb_importance) > 0) {
        stats::median(mb_importance)
      } else {
        0
      }

      # Non-MB vars that pass the threshold
      non_mb <- scores_df[!scores_df$in_mb, , drop = FALSE]
      extra_vars <- non_mb$variable[non_mb$importance > imp_threshold]

      selected <- unique(c(mb_vars, extra_vars))

      # Enforce minimum
      if (length(selected) < min_keep) {
        remaining <- setdiff(env_vars, selected)
        rem_imp <- scores_df$importance[match(remaining, scores_df$variable)]
        top_rem <- remaining[order(-rem_imp)][seq_len(
          min(min_keep - length(selected), length(remaining))
        )]
        selected <- unique(c(selected, top_rem))
      }
    }
  } else {
    # -- Fallback: pure RF importance with score-gap heuristic --
    if (stability_reps > 0L && any(scores_df$stable %in% TRUE)) {
      selected <- scores_df$variable[scores_df$stable %in% TRUE]
      if (length(selected) < min_keep) {
        remaining <- setdiff(scores_df$variable[order(-scores_df$importance)], selected)
        selected <- unique(c(
          selected,
          remaining[seq_len(min(min_keep - length(selected), length(remaining)))]
        ))
      }
    } else {
      selected <- select_by_importance(scores_df, min_keep)
    }
  }

  # --- Assign screening roles ----------------------------------------------
  roles_df <- data.frame(
    variable = selected,
    role = vapply(selected, function(v) {
      if (!isTRUE(mb_selective)) return("importance_screened")
      if (v %in% (mb$parents %||% character())) return("mb_direct")
      if (v %in% (mb$children %||% character())) return("mb_associated")
      if (v %in% (mb$co_parents %||% character())) return("mb_associated")
      "importance_added"
    }, character(1)),
    stringsAsFactors = FALSE
  )
  roles_df$selection_frequency <- scores_df$selection_frequency[
    match(roles_df$variable, scores_df$variable)
  ]
  roles_df$causal_role <- vapply(seq_len(nrow(roles_df)), function(i) {
    v <- roles_df$variable[i]
    legacy <- roles_df$role[i]
    stable <- roles_df$selection_frequency[i]
    stable_ok <- is.na(stable) || stable >= stability_threshold
    if (legacy == "mb_direct" && stable_ok) return("causal_core")
    if (legacy == "mb_associated" && stable_ok) return("causal_adjuster")
    if (!isTRUE(mb_selective) && v %in% (mb$all %||% character()) && stable_ok) {
      return("causal_core")
    }
    if (legacy %in% c("importance_added", "importance_screened")) {
      return("predictive_rescue")
    }
    "unstable_rejected"
  }, character(1))

  scores_df$causal_role <- roles_df$causal_role[
    match(scores_df$variable, roles_df$variable)
  ]
  scores_df$causal_role[is.na(scores_df$causal_role)] <- "unstable_rejected"

  # Sort: MB roles first, then by importance
  role_order <- c(
    mb_direct = 1,
    mb_associated = 2,
    importance_added = 3,
    importance_screened = 4,
    parent = 1,
    child = 2,
    co_parent = 2,
    predictive = 3
  )
  roles_df$role_rank <- role_order[roles_df$role]
  roles_df$imp_val <- scores_df$importance[match(roles_df$variable,
                                                  scores_df$variable)]
  roles_df <- roles_df[order(roles_df$role_rank, -roles_df$imp_val), ]
  roles_df$role_rank <- NULL
  roles_df$imp_val <- NULL
  rownames(roles_df) <- NULL

  # Order scores_df: selected first, then by importance
  scores_df <- scores_df[order(-scores_df$in_mb, -scores_df$importance), ]
  rownames(scores_df) <- NULL

  if (verbose) {
    n_mb <- sum(roles_df$role %in% c("mb_direct", "mb_associated"))
    n_pred <- sum(roles_df$role %in% c("importance_added", "importance_screened"))
    cli::cli_inform(
      "Selection: {length(selected)}/{length(env_vars)} variables ({n_mb} from MB, {n_pred} importance-screened)."
    )
    if (stability_reps > 0L) {
      cli::cli_inform(
        "Stability: {sum(scores_df$stable %in% TRUE)} variables above frequency >= {stability_threshold} across {stability_reps} bootstraps."
      )
    }
  }

  new_cast_select(
    selected = selected,
    scores = scores_df,
    roles = roles_df
  )
}
