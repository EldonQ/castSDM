#' Learn Causal DAG (Multiple Structure Learning Methods)
#'
#' Discovers the causal structure among environmental variables and
#' (optionally) the species response.
#'
#' When `include_response = TRUE` (the default), the `response` column is
#' included as a factor node in the DAG so that its **Markov Blanket** can be
#' extracted later by [cast_select()].
#' Because the network is then *mixed* (one discrete node among continuous
#' variables), the conditional-independence test / scoring criterion is
#' automatically adjusted:
#'
#' * **PC** (`structure_method = "pc"`): uses `test = "mi-cg"` (mutual
#'   information for conditional Gaussian data).
#' * **Bootstrap HC** (`structure_method = "bootstrap_hc"`): uses
#'   `score = "bic-cg"`.
#' * **BiDAG** (`structure_method = "bidag_bge"`): the BGe score requires
#'   continuous data, so the response is kept as numeric 0/1 (an acceptable
#'   approximation).
#'
#' The default method is `"mb_first"` (two-stage: fast MB discovery + local PC).
#'
#' @param data A `data.frame` containing the `presence` column and
#'   environmental variables.
#' @param response Character. Name of the response column. Default
#'   `"presence"`.
#' @param env_vars Character vector or `NULL`. Environmental variable names
#'   to include in DAG learning. When `NULL` (default), all numeric columns
#'   in `data` excluding `response`, `lon`, and `lat` are used.
#' @param include_response Logical. Include the `response` column as a node
#'   in DAG learning so that its Markov Blanket can be extracted. Default
#'   `TRUE`.
#' @param response_as_sink Logical. If `TRUE` (default), forbids directed
#'   edges from the response node to environmental predictors. This encodes
#'   the SDM assumption that environmental variables may explain occurrence,
#'   but occurrence does not cause static environmental predictors.
#' @param R Integer. Number of bootstrap replicates for `structure_method =
#'   "bootstrap_hc"`. Ignored otherwise (stored as `NA` for documentation).
#' @param algorithm Character. **Only used** when `structure_method =
#'   "bootstrap_hc"`: score-based learner passed to
#'   \code{bnlearn::boot.strength()} (e.g. `"hc"`, `"tabu"`, `"mmhc"`,
#'   `"pc.stable"`). Ignored for `"pc"` and `"bidag_bge"`.
#' @param score Character. **Only used** when `structure_method =
#'   "bootstrap_hc"`: scoring criterion in `algorithm.args` for
#'   \code{bnlearn::boot.strength()}. Default `"bic-cg"` (conditional
#'   Gaussian) when `include_response = TRUE`, `"bic-g"` otherwise.
#' @param strength_threshold Numeric. Minimum edge strength (bootstrap HC) or
#'   auxiliary threshold metadata for other methods. Default `0.7`.
#' @param direction_threshold Numeric. Minimum direction consistency (bootstrap
#'   HC only). Default `0.6`.
#' @param max_rows Integer. Maximum rows; subsample if exceeded. Default `8000`.
#' @param seed Integer or `NULL`. Random seed. Default `NULL`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#' @param structure_method Character. One of `"mb_first"` (default),
#'   `"pc"`, `"bootstrap_hc"`, `"bidag_bge"`.
#'
#'   The `"mb_first"` method uses a **two-stage** approach that is
#'   dramatically faster than `"pc"` for high-dimensional data (p > 20):
#'
#'   * **Stage 1** - Direct Markov Blanket discovery via IAMB-family
#'     algorithms (\code{bnlearn::learn.mb()}).  Complexity: O(|MB|^2 x p).
#'   * **Stage 2** - Full DAG learning via PC, restricted to the MB
#'     subset plus the response node.  Because |MB| is typically 5-15,
#'     this is near-instantaneous.
#'
#'   Under standard Markov Blanket and faithfulness assumptions, the two-stage
#'   result targets the response-relevant local structure while avoiding a
#'   full high-dimensional graph search (Pellet & Elisseeff 2008).
#' @param pc_alpha Significance level for constraint-based PC. Default `0.05`.
#'   Used when `structure_method` is `"pc"` or `"mb_first"` (Stage 2).
#' @param pc_test Conditional-independence test passed as `test` to
#'   \code{bnlearn::pc.stable()}. Default `"mi-cg"` (mutual information for
#'   conditional Gaussian data when `include_response = TRUE`; `"zf"` for
#'   pure-numeric networks). Used when `structure_method` is `"pc"` or
#'   `"mb_first"`.
#' @param mb_method Character.  Algorithm for Stage 1 Markov Blanket
#'   discovery (only used when `structure_method = "mb_first"`).  One of
#'   `"fast.iamb"` (default), `"iamb"`, `"inter.iamb"`, `"gs"`,
#'   `"iamb.fdr"`.  See \code{bnlearn::learn.mb}.
#' @param mb_alpha Numeric.  Significance level for Stage 1 MB discovery
#'   (only used when `structure_method = "mb_first"`).  Default `0.05`.
#'   This may differ from `pc_alpha` if you want a more liberal (larger)
#'   threshold in the MB screening phase.
#' @param mb_pc_max_nodes Integer. Maximum number of local nodes for Stage 2
#'   PC when `structure_method = "mb_first"`. If the discovered MB plus
#'   response exceeds this value and `mb_dense_action = "star"`, local PC is
#'   skipped. Default `15`.
#' @param mb_dense_action Character. Action when MB discovery is dense.
#'   Default `"star"` returns a response-focused MB graph instead of running
#'   local PC; `"pc"` forces local PC.
#' @param bidag_algorithm Passed to \code{BiDAG::learnBN()}: `"order"` or
#'   `"orderIter"`. Default `"order"`.
#' @param bidag_iterations Optional integer MCMC iterations for BiDAG.
#'   `NULL` uses a **bounded** default derived from the number of nodes.
#' @param blacklist A `data.frame` with columns `from` and `to` specifying
#'   **forbidden** directed edges. Default `NULL`.
#' @param whitelist A `data.frame` with columns `from` and `to` specifying
#'   **required** directed edges. Default `NULL`.
#'
#' @return A `cast_dag` object with a `response_node` field when
#'   `include_response = TRUE`.
#'
#' @references
#' Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package.
#' *Journal of Statistical Software*, 35(3), 1-22.
#'
#' Suter, P., Kuipers, J., Moffa, G., & Beerenwinkel, N. (2023). Bayesian
#' structure learning and sampling of Bayesian networks with the R package
#' BiDAG. *Journal of Statistical Software*, 105(9), 1-32.
#'
#' Tsamardinos, I., Aliferis, C. F., & Statnikov, A. (2003). Algorithms for
#' Large Scale Markov Blanket Discovery. *Proceedings of the Sixteenth
#' International Florida Artificial Intelligence Research Society Conference*,
#' 376-381.
#'
#' Pellet, J.-P. & Elisseeff, A. (2008). Using Markov Blankets for Causal
#' Structure Learning. *Journal of Machine Learning Research*, 9(43),
#' 1295-1342.
#'
#' @seealso \pkg{bnlearn} (\code{boot.strength}, \code{pc.stable}),
#'   \pkg{BiDAG} (\code{learnBN}).
#'
#' @export
cast_dag <- function(data,
                     response = "presence",
                     env_vars = NULL,
                     structure_method = c(
                       "mb_first", "pc", "bootstrap_hc", "bidag_bge"
                     ),
                     include_response = TRUE,
                     response_as_sink = TRUE,
                     R = 100L,
                     algorithm = "hc",
                     score = NULL,
                     strength_threshold = 0.7,
                     direction_threshold = 0.6,
                     max_rows = 8000L,
                     seed = NULL,
                     verbose = TRUE,
                     pc_alpha = 0.05,
                     pc_test = NULL,
                     mb_method = c(
                       "fast.iamb", "iamb", "inter.iamb", "gs", "iamb.fdr"
                     ),
                     mb_alpha = 0.05,
                     mb_pc_max_nodes = 15L,
                     mb_dense_action = c("star", "pc"),
                     bidag_algorithm = c("order", "orderIter"),
                     bidag_iterations = NULL,
                     blacklist = NULL,
                     whitelist = NULL) {
  structure_method <- match.arg(structure_method)
  bidag_algorithm <- match.arg(bidag_algorithm)
  mb_method <- match.arg(mb_method)
  mb_dense_action <- match.arg(mb_dense_action)

  # mb_first requires response in the DAG

  if (structure_method == "mb_first" && !isTRUE(include_response)) {
    cli::cli_abort(c(
      "{.code structure_method = \"mb_first\"} requires {.code include_response = TRUE}.",
      "i" = "The two-stage MB approach discovers the Markov Blanket of the response node."
    ))
  }

  # --- Resolve smart defaults for mixed-data settings ----------------------
  has_mixed <- isTRUE(include_response) && response %in% names(data)

  if (is.null(score)) {
    score <- if (has_mixed) "bic-cg" else "bic-g"
  }
  if (is.null(pc_test)) {
    pc_test <- if (has_mixed) "mi-cg" else "zf"
  }

  # --- bootstrap_hc only: valid bnlearn learners for boot.strength() -------
  boot_learners <- c(
    "gs", "iamb", "fast.iamb", "inter.iamb", "iamb.fdr", "pc.stable", "mmpc",
    "si.hiton.pc", "hpc", "hc", "tabu", "rsmax2", "mmhc", "h2pc",
    "chow.liu", "aracne", "naive.bayes", "tree.bayes", "structural.em"
  )
  if (structure_method == "bootstrap_hc") {
    algo_lc <- tolower(trimws(as.character(algorithm)))
    if (algo_lc == "pc") {
      cli::cli_abort(c(
        "{.arg algorithm} = {.val {algorithm}} is not valid for {.code structure_method = \"bootstrap_hc\"}.",
        "i" = "Use {.code structure_method = \"pc\"} for the constraint-based PC algorithm.",
        "i" = "For PC *inside* each bootstrap replicate, use {.code algorithm = \"pc.stable\"} instead."
      ))
    }
    if (!algorithm %in% boot_learners) {
      cli::cli_abort(c(
        "{.arg algorithm} = {.val {algorithm}} is not a valid learner for {.fn bnlearn::boot.strength}.",
        "i" = "Examples: {.val hc}, {.val tabu}, {.val mmhc}, {.val pc.stable}."
      ))
    }
  } else if (isTRUE(verbose) && (algorithm != "hc")) {
    cli::cli_inform(c(
      "i" = "Ignoring {.arg algorithm} (only applies to {.code structure_method = \"bootstrap_hc\"})."
    ))
  }

  env_vars <- env_vars %||% get_env_vars(data, response = response)
  if (length(env_vars) < 3) {
    cli::cli_abort("Need at least 3 environmental variables for DAG learning.")
  }

  if (isTRUE(include_response) && isTRUE(response_as_sink) &&
      response %in% names(data)) {
    sink_blacklist <- data.frame(
      from = response,
      to = env_vars,
      stringsAsFactors = FALSE
    )
    blacklist <- if (is.null(blacklist)) {
      sink_blacklist
    } else {
      unique(rbind(blacklist[, c("from", "to"), drop = FALSE],
                   sink_blacklist))
    }
  }

  # --- Build DAG data frame ------------------------------------------------
  # Environmental columns: always numeric
  dag_df <- as.data.frame(data[, env_vars, drop = FALSE])
  for (col in names(dag_df)) {
    dag_df[[col]] <- as.numeric(dag_df[[col]])
  }

  # Include response as a node (for Markov Blanket extraction)
  response_node <- NULL
  if (has_mixed) {
    resp_col <- data[[response]]
    if (structure_method == "bidag_bge") {
      # BiDAG BGe requires all-continuous; keep as numeric 0/1
      dag_df[[response]] <- as.numeric(resp_col)
    } else {
      # bnlearn CG networks require discrete nodes as factor
      dag_df[[response]] <- factor(resp_col, levels = c(0, 1))
    }
    response_node <- response
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

  # All node names in DAG
  all_nodes <- names(dag_df)

  if (verbose) {
    n_nodes <- length(all_nodes)
    resp_msg <- if (!is.null(response_node)) {
      " (incl. response)"
    } else {
      ""
    }
    cli::cli_inform(
      "Learning DAG ({structure_method}): {n_nodes} nodes{resp_msg}, {nrow(dag_df)} obs..."
    )
  }

  # Validate blacklist / whitelist
  if (!is.null(blacklist)) {
    if (!is.data.frame(blacklist) || !all(c("from", "to") %in% names(blacklist))) {
      cli::cli_abort("{.arg blacklist} must be a data.frame with columns {.val from} and {.val to}.")
    }
    blacklist <- blacklist[, c("from", "to"), drop = FALSE]
  }
  if (!is.null(whitelist)) {
    if (!is.data.frame(whitelist) || !all(c("from", "to") %in% names(whitelist))) {
      cli::cli_abort("{.arg whitelist} must be a data.frame with columns {.val from} and {.val to}.")
    }
    whitelist <- whitelist[, c("from", "to"), drop = FALSE]
  }

  dag_metadata <- list()
  env_edges <- switch(structure_method,
    pc = .dag_pc_edges(
      dag_df, alpha = pc_alpha, test = pc_test,
      seed = seed, verbose = verbose,
      blacklist = blacklist, whitelist = whitelist
    ),
    bootstrap_hc = .dag_bootstrap_hc(
      dag_df, R = R, algorithm = algorithm, score = score,
      strength_threshold = strength_threshold,
      direction_threshold = direction_threshold,
      seed = seed, verbose = verbose,
      blacklist = blacklist, whitelist = whitelist
    ),
    mb_first = .dag_mb_first(
      dag_df, response = response_node, mb_method = mb_method,
      mb_alpha = mb_alpha, pc_alpha = pc_alpha, test = pc_test,
      mb_pc_max_nodes = mb_pc_max_nodes,
      dense_action = mb_dense_action,
      seed = seed, verbose = verbose,
      blacklist = blacklist, whitelist = whitelist
    ),
    bidag_bge = .dag_bidag_edges(
      dag_df,
      algorithm = bidag_algorithm,
      iterations = bidag_iterations,
      seed = seed, verbose = verbose,
      blacklist = blacklist, whitelist = whitelist
    ),
    cli::cli_abort("Unknown {.arg structure_method}.")
  )

  if (structure_method == "mb_first" && isTRUE(attr(env_edges, "dense_mb"))) {
    dag_metadata$dense_mb <- TRUE
    dag_metadata$graph_role <- "response_focused_screening"
    dag_metadata$mb_size <- attr(env_edges, "mb_size")
    dag_metadata$mb_total <- attr(env_edges, "mb_total")
  }
  if (structure_method == "mb_first") {
    dag_metadata$mb_vars_stage1 <- attr(env_edges, "mb_vars_stage1")
    dag_metadata$response_as_sink <- isTRUE(response_as_sink)
  }

  boot_R_out <- if (structure_method == "bootstrap_hc") R else NA_integer_
  score_out <- if (structure_method == "bootstrap_hc") {
    score
  } else if (structure_method == "bidag_bge") {
    "bge"
  } else if (structure_method == "mb_first") {
    paste0("mb_first:", mb_method)
  } else {
    structure_method
  }

  if (verbose) {
    cli::cli_inform("DAG: {nrow(env_edges)} edges ({structure_method}).")
  }

  new_cast_dag(
    edges = env_edges,
    nodes = all_nodes,
    boot_R = boot_R_out,
    strength_threshold = strength_threshold,
    direction_threshold = direction_threshold,
    score = score_out,
    structure_method = structure_method,
    response_node = response_node,
    metadata = dag_metadata
  )
}


# ---- Internal: Extract Markov Blanket from DAG ---------------------------

#' Extract Markov Blanket of a Node from a cast_dag
#'
#' `MB(Y) = parents(Y) U children(Y) U parents_of_children(Y)`.
#' Under standard Markov Blanket assumptions, this is the local set targeted
#' for response-focused variable screening.
#'
#' @param dag A `cast_dag` object.
#' @param node Character. Target node name (typically the response).
#'
#' @return A list with components `parents`, `children`, `co_parents`, and
#'   `all` (the union).
#'
#' @keywords internal
#' @noRd
extract_markov_blanket <- function(dag, node) {
  edges <- dag$edges
  if (is.null(edges) || nrow(edges) == 0) {
    return(list(
      parents = character(), children = character(),
      co_parents = character(), all = character()
    ))
  }

  # Parents: nodes with an edge TO the target node
  parents <- unique(edges$from[edges$to == node])

  # Children: nodes that the target node has an edge TO
  children <- unique(edges$to[edges$from == node])

  # Co-parents (spouses): other parents of children, excluding node itself
  co_parents <- character()
  for (ch in children) {
    ch_parents <- edges$from[edges$to == ch]
    co_parents <- c(co_parents, ch_parents)
  }
  co_parents <- unique(setdiff(co_parents, c(node, parents, children)))

  list(
    parents = parents,
    children = children,
    co_parents = co_parents,
    all = unique(c(parents, children, co_parents))
  )
}


#' Response-Focused Markov Blanket Variables
#'
#' Internal helper used by variable screening and CATE so both stages use the
#' same MB definition. For `mb_first` DAGs, Stage 1 MB membership is preferred
#' because it is the screening target even when local PC edge orientation is
#' sparse or unstable.
#'
#' @keywords internal
#' @noRd
response_markov_blanket <- function(dag, response_node, env_vars = NULL) {
  mb <- extract_markov_blanket(dag, response_node)
  env_vars <- env_vars %||% setdiff(dag$nodes, c(response_node, "lon", "lat"))
  stage1_mb <- dag$metadata$mb_vars_stage1 %||% character()

  if (length(stage1_mb) > 0) {
    mb$all <- unique(intersect(stage1_mb, env_vars))
    if (isTRUE(dag$metadata$response_as_sink)) {
      mb$parents <- mb$all
      mb$children <- character()
      mb$co_parents <- character()
    }
  }

  mb
}


# ---- Structure learning backends -----------------------------------------

#' @keywords internal
#' @noRd
.dag_bootstrap_hc <- function(dag_df, R, algorithm, score,
                              strength_threshold, direction_threshold,
                              seed, verbose,
                              blacklist = NULL, whitelist = NULL) {
  check_suggested("bnlearn", "for DAG structure learning")
  if (!is.null(seed)) set.seed(seed)
  boot_strength <- utils::getFromNamespace("boot.strength", "bnlearn")
  boot_args <- list(
    data = dag_df, R = R, algorithm = algorithm,
    algorithm.args = list(score = score)
  )
  if (!is.null(blacklist)) boot_args$algorithm.args$blacklist <- blacklist
  if (!is.null(whitelist)) boot_args$algorithm.args$whitelist <- whitelist
  boot_str <- suppressWarnings(do.call(boot_strength, boot_args))
  strong <- boot_str[
    boot_str$strength >= strength_threshold &
      boot_str$direction >= direction_threshold, ,
    drop = FALSE
  ]
  # Re-inject whitelisted edges that may not pass threshold
  if (!is.null(whitelist)) {
    for (r in seq_len(nrow(whitelist))) {
      wl_from <- whitelist$from[r]
      wl_to   <- whitelist$to[r]
      found <- any(strong$from == wl_from & strong$to == wl_to)
      if (!found) {
        strong <- rbind(strong, data.frame(
          from = wl_from, to = wl_to,
          strength = 1.0, direction = 1.0,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  as.data.frame(strong)
}


#' @keywords internal
#' @noRd
.dag_pc_edges <- function(dag_df, alpha, test, seed, verbose,
                          blacklist = NULL, whitelist = NULL) {
  check_suggested("bnlearn", "for PC algorithm")
  if (!is.null(seed)) set.seed(seed)
  pc_fun <- tryCatch(
    utils::getFromNamespace("pc.stable", "bnlearn"),
    error = function(e) {
      cli::cli_abort(c(
        "Could not load {.fn bnlearn::pc.stable} from {.pkg bnlearn}.",
        "i" = "Update {.pkg bnlearn} to a recent CRAN version."
      ))
    }
  )
  pc_args <- list(x = dag_df, test = test, alpha = alpha, debug = FALSE)
  if (!is.null(blacklist)) pc_args$blacklist <- blacklist
  if (!is.null(whitelist)) pc_args$whitelist <- whitelist
  learned <- suppressWarnings(do.call(pc_fun, pc_args))
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


# ---- MB-First: Two-stage MB discovery + local DAG --------------------------

#' Two-Stage MB-Guided Structure Learning
#'
#' Stage 1: Uses [bnlearn::learn.mb()] (IAMB-family) to discover the
#' Markov Blanket of the response - complexity O(|MB|^2 x p).
#' Stage 2: Runs [bnlearn::pc.stable()] on the MB subset + response to
#' recover edge directions - near-instantaneous for typical |MB| ~ 5-15.
#'
#' Under standard Markov Blanket and faithfulness assumptions this targets the
#' same response-relevant local structure as full-DAG PC while avoiding a full
#' high-dimensional graph search (Pellet & Elisseeff 2008, JMLR).
#'
#' @keywords internal
#' @noRd
.dag_mb_first <- function(dag_df, response, mb_method, mb_alpha,
                          pc_alpha, test, mb_pc_max_nodes, dense_action,
                          seed, verbose,
                          blacklist = NULL, whitelist = NULL) {
  check_suggested("bnlearn", "for Markov Blanket discovery")
  if (!is.null(seed)) set.seed(seed)

  # --- Stage 1: Direct MB discovery -----------------------------------------
  learn_mb_fn <- utils::getFromNamespace("learn.mb", "bnlearn")

  if (is.null(response) || !response %in% names(dag_df)) {
    cli::cli_abort(c(
      "{.code structure_method = \"mb_first\"} requires a response node in the DAG.",
      "i" = "Set {.code include_response = TRUE} (the default)."
    ))
  }

  if (verbose) {
    p <- ncol(dag_df) - 1L
    cli::cli_inform(
      "Stage 1: MB discovery ({mb_method}) for {.val {response}} among {p} variables..."
    )
  }

  mb_args <- list(x = dag_df, node = response, method = mb_method,
                  test = test, alpha = mb_alpha)
  mb_vars <- tryCatch(
    suppressWarnings(do.call(learn_mb_fn, mb_args)),
    error = function(e) {
      cli::cli_warn("MB discovery ({mb_method}) failed: {e$message}")
      character()
    }
  )

  if (verbose) {
    cli::cli_inform(
      "  MB({response}): {length(mb_vars)} variable{?s}: {.val {mb_vars}}"
    )
  }

  # --- Fallback: empty MB -> full PC -----------------------------------------
  if (length(mb_vars) == 0L) {
    if (verbose) {
      cli::cli_warn(
        "MB discovery returned empty set; falling back to full PC on all {ncol(dag_df)} nodes."
      )
    }
    out <- .dag_pc_edges(
      dag_df, alpha = pc_alpha, test = test,
      seed = seed, verbose = verbose,
      blacklist = blacklist, whitelist = whitelist
    )
    attr(out, "mb_vars_stage1") <- character()
    return(out)
  }

  # --- Stage 2: Local DAG on MB subset --------------------------------------
  local_vars <- unique(c(mb_vars, response))
  local_df <- dag_df[, local_vars, drop = FALSE]

  if (identical(dense_action, "star") &&
      length(local_vars) > as.integer(mb_pc_max_nodes)) {
    if (verbose) {
      cli::cli_warn(c(
        "Dense MB: {length(mb_vars)}/{ncol(dag_df) - 1L} predictors selected; skipping local PC on {length(local_vars)} nodes.",
        "i" = "Returning a response-focused MB graph instead. Increase {.arg mb_pc_max_nodes} or set {.arg mb_dense_action = \"pc\"} to force local PC."
      ))
    }
    out <- .dag_mb_star_edges(
      mb_vars = mb_vars,
      response = response,
      total_vars = ncol(dag_df) - 1L,
      blacklist = blacklist,
      whitelist = whitelist
    )
    attr(out, "mb_vars_stage1") <- mb_vars
    return(out)
  }

  if (verbose) {
    cli::cli_inform(
      "Stage 2: Local DAG (PC) on {length(local_vars)} nodes..."
    )
  }

  # Filter blacklist/whitelist to local vars only
  local_bl <- .filter_edge_list(blacklist, local_vars)
  local_wl <- .filter_edge_list(whitelist, local_vars)

  out <- .dag_pc_edges(
    local_df, alpha = pc_alpha, test = test,
    seed = seed, verbose = FALSE,
    blacklist = local_bl, whitelist = local_wl
  )

  if (!any(out$from == response | out$to == response)) {
    if (verbose) {
      cli::cli_warn(c(
        "Local PC did not retain any edge adjacent to the response node.",
        "i" = "Returning a response-focused screening graph from Stage-1 MB membership instead."
      ))
    }
    out <- .dag_mb_star_edges(
      mb_vars = mb_vars,
      response = response,
      total_vars = ncol(dag_df) - 1L,
      blacklist = blacklist,
      whitelist = whitelist
    )
  }

  attr(out, "mb_vars_stage1") <- mb_vars
  out
}


#' @keywords internal
#' @noRd
.dag_mb_star_edges <- function(mb_vars, response,
                               total_vars = length(mb_vars),
                               blacklist = NULL, whitelist = NULL) {
  mb_vars <- unique(setdiff(as.character(mb_vars), response))
  if (length(mb_vars) == 0L) {
    out <- data.frame(
      from = character(), to = character(),
      strength = numeric(), direction = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      from = mb_vars,
      to = response,
      strength = 1,
      direction = 1,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(blacklist) && nrow(blacklist) > 0L && nrow(out) > 0L) {
    bad <- paste(blacklist$from, blacklist$to, sep = "->")
    out <- out[!paste(out$from, out$to, sep = "->") %in% bad, , drop = FALSE]
  }

  if (!is.null(whitelist) && nrow(whitelist) > 0L) {
    wl <- whitelist[
      whitelist$from %in% c(mb_vars, response) &
        whitelist$to %in% c(mb_vars, response),
      c("from", "to"),
      drop = FALSE
    ]
    if (nrow(wl) > 0L) {
      wl$strength <- 1
      wl$direction <- 1
      out <- unique(rbind(out, wl))
    }
  }

  rownames(out) <- NULL
  attr(out, "dense_mb") <- TRUE
  attr(out, "mb_size") <- length(mb_vars)
  attr(out, "mb_total") <- total_vars
  out
}


#' Filter an edge list (blacklist/whitelist) to only contain edges
#' between nodes in `node_set`.
#' @keywords internal
#' @noRd
.filter_edge_list <- function(edge_df, node_set) {
  if (is.null(edge_df) || nrow(edge_df) == 0L) return(NULL)
  keep <- edge_df$from %in% node_set & edge_df$to %in% node_set
  if (!any(keep)) return(NULL)
  edge_df[keep, , drop = FALSE]
}


# ---- BiDAG backend -------------------------------------------------------

#' @keywords internal
#' @noRd
.bidag_stabilize_df <- function(dag_df, max_rows = 1000L, seed = NULL) {
  dag_df <- as.data.frame(dag_df, stringsAsFactors = FALSE)
  for (nm in names(dag_df)) {
    dag_df[[nm]] <- suppressWarnings(as.numeric(dag_df[[nm]]))
  }
  dag_df <- stats::na.omit(dag_df)
  if (nrow(dag_df) < 10L) {
    cli::cli_abort("BiDAG: need at least 10 complete rows after NA removal.")
  }
  for (nm in names(dag_df)) {
    x <- dag_df[[nm]]
    rep_val <- stats::median(x[is.finite(x)], na.rm = TRUE)
    if (!is.finite(rep_val)) rep_val <- 0
    x[!is.finite(x)] <- rep_val
    dag_df[[nm]] <- x
  }
  sdv <- vapply(dag_df, stats::sd, numeric(1), na.rm = TRUE)
  keep <- names(dag_df)[is.finite(sdv) & sdv >= 1e-12]
  if (length(keep) < 3L) {
    cli::cli_abort("BiDAG: need at least 3 non-constant numeric columns.")
  }
  dag_df <- dag_df[, keep, drop = FALSE]
  nr <- nrow(dag_df)
  if (nr > max_rows) {
    if (!is.null(seed)) set.seed(seed)
    dag_df <- dag_df[sample.int(nr, max_rows), , drop = FALSE]
  }
  sc <- as.data.frame(
    lapply(dag_df, function(z) {
      z <- as.numeric(z)
      s <- stats::sd(z, na.rm = TRUE)
      if (!is.finite(s) || s < 1e-12) {
        return(z)
      }
      as.numeric(scale.default(z, center = TRUE, scale = TRUE))
    }),
    stringsAsFactors = FALSE
  )
  colnames(sc) <- names(dag_df)
  stats::na.omit(sc)
}


#' @keywords internal
#' @noRd
.bidag_safe_iterations <- function(p, iterations) {
  logp <- if (p <= 1L) 0 else log(as.numeric(p))
  n_mc <- 6 * (as.numeric(p)^2) * logp
  if (!is.finite(n_mc) || n_mc < 1) {
    n_mc <- 200
  }
  it_default <- max(200L, min(4000L, as.integer(round(n_mc))))
  it_use <- if (is.null(iterations)) {
    it_default
  } else {
    max(100L, min(5000L, suppressWarnings(as.integer(iterations))))
  }
  if (length(it_use) != 1L || is.na(it_use) || it_use < 50L) {
    it_use <- 500L
  }
  stepsave <- max(1L, min(it_use, as.integer(it_use / 200L)))
  list(iterations = it_use, stepsave = stepsave)
}


#' @keywords internal
#' @noRd
.dag_bidag_edges <- function(dag_df, algorithm, iterations, seed, verbose,
                             blacklist = NULL, whitelist = NULL) {
  check_suggested("BiDAG", "for BiDAG / high-dimensional Bayesian DAG search")
  if (!is.null(seed)) set.seed(seed)

  dag_df <- .bidag_stabilize_df(dag_df, max_rows = 1000L, seed = seed)

  scoreparameters <- utils::getFromNamespace("scoreparameters", "BiDAG")
  learn_bn <- utils::getFromNamespace("learnBN", "BiDAG")
  get_dag <- utils::getFromNamespace("getDAG", "BiDAG")

  scorepar <- scoreparameters("bge", dag_df)
  p <- ncol(dag_df)
  it_ctl <- .bidag_safe_iterations(p, iterations)

  learn_args <- list(
    scorepar = scorepar,
    algorithm = algorithm,
    verbose = verbose,
    chainout = FALSE,
    iterations = it_ctl$iterations,
    stepsave = it_ctl$stepsave
  )

  # BiDAG uses a binary adjacency matrix as startspace for constraints
  if (!is.null(blacklist) || !is.null(whitelist)) {
    vnames <- names(dag_df)
    ss <- matrix(1L, nrow = p, ncol = p,
                 dimnames = list(vnames, vnames))
    diag(ss) <- 0L
    if (!is.null(blacklist)) {
      for (r in seq_len(nrow(blacklist))) {
        fi <- match(blacklist$from[r], vnames)
        ti <- match(blacklist$to[r], vnames)
        if (!is.na(fi) && !is.na(ti)) ss[fi, ti] <- 0L
      }
    }
    learn_args$startspace <- ss
  }

  fit <- tryCatch(
    do.call(learn_bn, learn_args),
    error = function(e) {
      cli::cli_abort(c(
        "{.pkg BiDAG}::{.fn learnBN} failed: {e$message}",
        "i" = "Try fewer variables, smaller {.arg max_rows} before BiDAG, or {.code structure_method = \"pc\"}."
      ))
    }
  )
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
