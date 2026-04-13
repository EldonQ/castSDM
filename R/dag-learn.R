#' Learn Causal DAG (Multiple Structure Learning Methods)
#'
#' Discovers the causal structure among environmental variables. The default
#' (`structure_method = "bootstrap_hc"`) uses bootstrap-aggregated
#' Hill-Climbing as in the original CAST workflow. Alternative **constraint-based**
#' (`"pc"`, `"fci"`), **Bayesian MAP search** via \pkg{BiDAG} (`"bidag_bge"`,
#' recommended for many variables), and **continuous optimization** linear NOTEARS
#' (`"notears_linear"`, requires \pkg{torch}) are also available.
#'
#' The DAG is always learned on **environmental variables only** (presence
#' excluded); see Details in previous versions for rationale.
#'
#' @param data A `data.frame` containing the `presence` column and
#'   environmental variables. Must be numeric (no factors).
#' @param response Character. Name of the response column. Default
#'   `"presence"`.
#' @param env_vars Character vector or `NULL`. Environmental variable names
#'   to include in DAG learning. When `NULL` (default), all numeric columns
#'   in `data` excluding `response`, `lon`, and `lat` are used.
#' @param R Integer. Number of bootstrap replicates for `structure_method =
#'   "bootstrap_hc"`. Ignored otherwise (stored as `NA` for documentation).
#' @param algorithm Character. Score-based algorithm passed to
#'   \code{bnlearn::boot.strength()} when `structure_method = "bootstrap_hc"`.
#'   Default `"hc"`.
#' @param score Character. Scoring criterion for bootstrap HC. Default `"bic-g"`.
#' @param strength_threshold Numeric. Minimum edge strength (bootstrap HC) or
#'   auxiliary threshold metadata for other methods. Default `0.7`.
#' @param direction_threshold Numeric. Minimum direction consistency (bootstrap
#'   HC only). Default `0.6`.
#' @param max_rows Integer. Maximum rows; subsample if exceeded. Default `8000`.
#' @param seed Integer or `NULL`. Random seed. Default `NULL`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#' @param structure_method Character. One of `"bootstrap_hc"` (default),
#'   `"pc"`, `"fci"`, `"bidag_bge"`, `"notears_linear"`.
#' @param pc_alpha Significance level for \code{bnlearn::pc()} (Gaussian CI test).
#'   Default `0.05`.
#' @param fci_alpha Significance level for \code{bnlearn::fci()}. Default `0.05`.
#' @param bidag_algorithm Passed to \code{BiDAG::learnBN()}: `"order"` or
#'   `"orderIter"`. Default `"order"`.
#' @param bidag_iterations Optional integer MCMC iterations for BiDAG. `NULL`
#'   uses BiDAG defaults (data-driven).
#' @param notears_lambda L1 penalty weight for `notears_linear`. Default `0.03`.
#' @param notears_max_iter Maximum optimization steps. Default `2000`.
#' @param notears_lr Adam learning rate. Default `0.02`.
#' @param notears_tol Stop when acyclicity metric `|h(W)| < notears_tol`.
#'   Default `1e-3`.
#' @param notears_rho_init Initial augmented-Lagrangian rho. Default `0.1`.
#' @param notears_alpha_mult Multiplier on rho each 100 steps. Default `1.01`.
#'
#' @return A `cast_dag` object (`structure_method` is stored on the object).
#'
#' @references
#' Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package.
#' *Journal of Statistical Software*, 35(3), 1-22.
#'
#' Zheng, X., Aragam, B., Ravikumar, P., & Xing, E. P. (2018). DAGs with NO
#' TEARS: Continuous Optimization for Structure Learning. *NeurIPS*.
#'
#' Suter, P., Kuipers, J., Moffa, G., & Beerenwinkel, N. (2023). Bayesian
#' structure learning and sampling of Bayesian networks with the R package
#' BiDAG. *Journal of Statistical Software*, 105(9), 1-32.
#'
#' @seealso \pkg{bnlearn} (\code{boot.strength}, \code{pc}, \code{fci}),
#'   \pkg{BiDAG} (\code{learnBN}).
#'
#' @export
cast_dag <- function(data,
                     response = "presence",
                     env_vars = NULL,
                     R = 100L,
                     algorithm = "hc",
                     score = "bic-g",
                     strength_threshold = 0.7,
                     direction_threshold = 0.6,
                     max_rows = 8000L,
                     seed = NULL,
                     verbose = TRUE,
                     structure_method = c(
                       "bootstrap_hc", "pc", "fci",
                       "bidag_bge", "notears_linear"
                     ),
                     pc_alpha = 0.05,
                     fci_alpha = 0.05,
                     bidag_algorithm = c("order", "orderIter"),
                     bidag_iterations = NULL,
                     notears_lambda = 0.03,
                     notears_max_iter = 2000L,
                     notears_lr = 0.02,
                     notears_tol = 1e-3,
                     notears_rho_init = 0.1,
                     notears_alpha_mult = 1.01) {
  structure_method <- match.arg(structure_method)
  bidag_algorithm <- match.arg(bidag_algorithm)

  env_vars <- env_vars %||% get_env_vars(data, response = response)
  if (length(env_vars) < 3) {
    cli::cli_abort("Need at least 3 environmental variables for DAG learning.")
  }

  dag_df <- as.data.frame(data[, env_vars, drop = FALSE])
  for (col in names(dag_df)) {
    dag_df[[col]] <- as.numeric(dag_df[[col]])
  }
  dag_df <- stats::na.omit(dag_df)

  if (nrow(dag_df) < 10) {
    cli::cli_abort(
      "Fewer than 10 complete cases for DAG learning ({nrow(dag_df)} available)."
    )
  }

  if (nrow(dag_df) > max_rows) {
    if (!is.null(seed)) set.seed(seed)
    dag_df <- dag_df[sample.int(nrow(dag_df), max_rows), ]
  }

  if (verbose) {
    cli::cli_inform(
      "Learning DAG ({structure_method}): {length(env_vars)} env vars, {nrow(dag_df)} obs..."
    )
  }

  env_edges <- switch(structure_method,
    bootstrap_hc = .dag_bootstrap_hc(
      dag_df, R = R, algorithm = algorithm, score = score,
      strength_threshold = strength_threshold,
      direction_threshold = direction_threshold,
      seed = seed, verbose = verbose
    ),
    pc = .dag_pc_edges(dag_df, alpha = pc_alpha, seed = seed, verbose = verbose),
    fci = .dag_fci_edges(dag_df, alpha = fci_alpha, seed = seed, verbose = verbose),
    bidag_bge = .dag_bidag_edges(
      dag_df,
      algorithm = bidag_algorithm,
      iterations = bidag_iterations,
      seed = seed, verbose = verbose
    ),
    notears_linear = notears_learn_edges(
      dag_df,
      lambda = notears_lambda,
      max_iter = notears_max_iter,
      lr = notears_lr,
      tol = notears_tol,
      rho = notears_rho_init,
      alpha_mult = notears_alpha_mult,
      seed = seed,
      verbose = verbose
    ),
    cli::cli_abort("Unknown {.arg structure_method}.")
  )

  boot_R_out <- if (structure_method == "bootstrap_hc") R else NA_integer_
  score_out <- if (structure_method == "bootstrap_hc") {
    score
  } else if (structure_method == "bidag_bge") {
    "bge"
  } else if (structure_method == "notears_linear") {
    "notears_linear"
  } else {
    structure_method
  }

  if (verbose) {
    cli::cli_inform("DAG: {nrow(env_edges)} edges ({structure_method}).")
  }

  new_cast_dag(
    edges = env_edges,
    nodes = env_vars,
    boot_R = boot_R_out,
    strength_threshold = strength_threshold,
    direction_threshold = direction_threshold,
    score = score_out,
    structure_method = structure_method
  )
}


#' @keywords internal
#' @noRd
.dag_bootstrap_hc <- function(dag_df, R, algorithm, score,
                              strength_threshold, direction_threshold,
                              seed, verbose) {
  check_suggested("bnlearn", "for DAG structure learning")
  if (!is.null(seed)) set.seed(seed)
  boot_strength <- utils::getFromNamespace("boot.strength", "bnlearn")
  boot_str <- boot_strength(
    dag_df, R = R, algorithm = algorithm,
    algorithm.args = list(score = score)
  )
  strong <- boot_str[
    boot_str$strength >= strength_threshold &
      boot_str$direction >= direction_threshold, ,
    drop = FALSE
  ]
  as.data.frame(strong)
}


#' @keywords internal
#' @noRd
.dag_pc_edges <- function(dag_df, alpha, seed, verbose) {
  check_suggested("bnlearn", "for PC algorithm")
  if (!is.null(seed)) set.seed(seed)
  pc_fun <- utils::getFromNamespace("pc", "bnlearn")
  learned <- pc_fun(dag_df, test = "zf", alpha = alpha, debug = FALSE)
  .bn_arcs_to_cast_edges(learned)
}


#' @keywords internal
#' @noRd
.dag_fci_edges <- function(dag_df, alpha, seed, verbose) {
  check_suggested("bnlearn", "for FCI algorithm")
  if (!is.null(seed)) set.seed(seed)
  fci_fun <- utils::getFromNamespace("fci", "bnlearn")
  learned <- fci_fun(dag_df, test = "zf", alpha = alpha, debug = FALSE)
  .bn_arcs_to_cast_edges(learned)
}


#' @keywords internal
#' @noRd
.bn_arcs_to_cast_edges <- function(bn_or_fitted) {
  arcs_fun <- utils::getFromNamespace("arcs", "bnlearn")
  ar <- tryCatch(
    arcs_fun(bn_or_fitted),
    error = function(e) matrix(nrow = 0, ncol = 2)
  )
  if (is.null(ar) || nrow(ar) == 0) {
    return(data.frame(
      from = character(), to = character(),
      strength = numeric(), direction = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    from = ar[, 1],
    to = ar[, 2],
    strength = 1,
    direction = 1,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
#' @noRd
.dag_bidag_edges <- function(dag_df, algorithm, iterations, seed, verbose) {
  check_suggested("BiDAG", "for BiDAG / high-dimensional Bayesian DAG search")
  if (!is.null(seed)) set.seed(seed)

  scoreparameters <- utils::getFromNamespace("scoreparameters", "BiDAG")
  learn_bn <- utils::getFromNamespace("learnBN", "BiDAG")
  get_dag <- utils::getFromNamespace("getDAG", "BiDAG")

  scorepar <- scoreparameters("bge", dag_df)
  it <- iterations
  learn_args <- list(
    scorepar = scorepar,
    algorithm = algorithm,
    verbose = verbose,
    chainout = FALSE
  )
  if (!is.null(it)) learn_args$iterations <- it

  fit <- do.call(learn_bn, learn_args)
  Amat <- as.matrix(get_dag(fit))
  if (is.null(rownames(Amat))) {
    rownames(Amat) <- colnames(Amat) <- names(dag_df)
  }
  dn <- colnames(Amat)
  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(Amat))) {
    for (j in seq_len(ncol(Amat))) {
      if (Amat[i, j] != 0L) {
        k <- k + 1L
        rows[[k]] <- data.frame(
          from = dn[i],
          to = dn[j],
          strength = 1,
          direction = 1,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (k == 0L) {
    return(data.frame(
      from = character(), to = character(),
      strength = numeric(), direction = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# Linear NOTEARS (torch). See dag-notears.R history; kept here for load order.
#' @keywords internal
#' @noRd
notears_learn_edges <- function(dag_df,
                                  lambda = 0.03,
                                  max_iter = 2000L,
                                  lr = 0.02,
                                  tol = 1e-3,
                                  mu_init = 0.0,
                                  rho = 0.1,
                                  alpha_mult = 1.01,
                                  seed = NULL,
                                  verbose = FALSE) {
  check_suggested("torch", "for NOTEARS linear structure learning")

  X <- as.matrix(dag_df)
  storage.mode(X) <- "double"
  X <- scale(X)
  if (any(!is.finite(X))) {
    cli::cli_abort("NOTEARS: non-finite values after scaling (check inputs).")
  }

  p <- ncol(X)
  if (p > 60L) {
    cli::cli_abort(c(
      "NOTEARS linear is limited to moderate dimension (p <= 60).",
      i = "Try {.code structure_method = 'bidag_bge'} for high-dimensional data."
    ))
  }

  if (!is.null(seed)) {
    set.seed(seed)
    torch::torch_manual_seed(as.integer(seed))
  }

  Xt <- torch::torch_tensor(X, dtype = torch::torch_float32())

  W <- torch::torch_zeros(p, p, requires_grad = TRUE)
  opt <- torch::optim_adam(list(W), lr = lr)

  torch_ns <- asNamespace("torch")
  mxexp <- if (exists("linalg_matrix_exp", torch_ns, inherits = FALSE)) {
    get("linalg_matrix_exp", envir = torch_ns)
  } else {
    cli::cli_abort(c(
      "NOTEARS requires a recent {.pkg torch} with matrix exponential.",
      i = "Update torch or use {.code structure_method = 'bidag_bge'}."
    ))
  }

  mu <- mu_init
  h_np <- NA_real_

  for (it in seq_len(max_iter)) {
    opt$zero_grad()

    pred <- torch::torch_matmul(Xt, W)
    mse <- torch::torch_mean((Xt - pred)^2)
    l1 <- torch::torch_mean(torch::torch_abs(W))
    loss <- mse + lambda * l1

    M <- W * W
    E <- mxexp(M)
    h <- torch::torch_trace(E) - torch::torch_tensor(
      p, dtype = torch::torch_float32()
    )
    aug <- loss + mu * h + (rho / 2.0) * (h^2)

    aug$backward()
    opt$step()

    torch::with_no_grad({
      W$clamp_(-1.0, 1.0)
    })
    h_np <- as.numeric(h$detach())

    if (it %% 100L == 0L) {
      mu <- mu + rho * h_np
      rho <- rho * alpha_mult
      if (verbose) {
        cli::cli_inform("  NOTEARS iter {it}: h={round(h_np, 6)} rho={round(rho, 4)}")
      }
    }

    if (is.finite(h_np) && abs(h_np) < tol && it > 200L) break
  }

  Wm <- torch::as_array(W$detach())
  if (!is.matrix(Wm)) Wm <- matrix(Wm, nrow = p, ncol = p)
  diag(Wm) <- 0
  thr <- stats::quantile(abs(Wm), probs = 0.85, na.rm = TRUE)
  thr <- max(thr, 1e-3)

  edges <- list()
  ei <- 0L
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i == j) next
      if (abs(Wm[i, j]) >= thr) {
        ei <- ei + 1L
        edges[[ei]] <- data.frame(
          from = colnames(X)[j],
          to = colnames(X)[i],
          strength = min(1, abs(Wm[i, j])),
          direction = 1.0,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(edges) == 0) {
    return(data.frame(
      from = character(), to = character(),
      strength = numeric(), direction = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, edges)
  rownames(out) <- NULL
  out
}
