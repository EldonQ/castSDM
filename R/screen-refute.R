#' Refute a castSDM Variable Screen
#'
#' Runs lightweight DoWhy-style refutation diagnostics for a response-focused
#' variable screen. These diagnostics are not a proof of causality; they check
#' whether selected variables are stable under common negative controls and
#' resampling perturbations.
#'
#' @param dag A [cast_dag] object.
#' @param screen A [cast_select] object.
#' @param data A `data.frame` with response and environmental variables.
#' @param response Character. Response column. Default `"presence"`.
#' @param reps Integer. Number of subset/bootstrap perturbations.
#' @param subset_fraction Numeric. Fraction of rows retained in subset tests.
#' @param num_trees Integer. Number of trees for refutation RF screens.
#' @param seed Integer or `NULL`.
#' @param verbose Logical.
#'
#' @return A `cast_refute` object with `tests` and per-variable `summary`
#'   tables.
#' @export
cast_refute_screen <- function(dag,
                               screen,
                               data,
                               response = "presence",
                               reps = 20L,
                               subset_fraction = 0.8,
                               num_trees = 100L,
                               seed = NULL,
                               verbose = TRUE) {
  check_suggested("ranger", "for screen refutation")

  env_vars <- setdiff(dag$nodes, c(response, "lon", "lat"))
  env_vars <- intersect(env_vars, names(data))
  selected <- screen$selected %||% character()
  min_keep <- max(1L, length(selected))
  reps <- max(1L, as.integer(reps))

  rf_screen <- function(df, y_col = response, extra_vars = character()) {
    vars <- unique(c(env_vars, extra_vars))
    vars <- intersect(vars, names(df))
    X <- as.data.frame(df[, vars, drop = FALSE])
    for (nm in names(X)) X[[nm]] <- as.numeric(X[[nm]])
    X[is.na(X)] <- 0
    Y <- df[[y_col]]
    imp <- ranger::ranger(
      y ~ .,
      data = cbind(y = as.factor(Y), X),
      num.trees = num_trees,
      importance = "permutation",
      verbose = FALSE
    )$variable.importance
    imp[is.na(imp)] <- 0
    names(sort(imp, decreasing = TRUE))[seq_len(min(min_keep, length(imp)))]
  }

  subset_counts <- stats::setNames(rep(0L, length(selected)), selected)
  block_counts <- stats::setNames(rep(0L, length(selected)), selected)
  test_rows <- list()

  for (ii in seq_len(reps)) {
    if (!is.null(seed)) set.seed(seed + ii)
    idx <- sample.int(nrow(data), max(10L, floor(nrow(data) * subset_fraction)))
    sub_df <- data[idx, , drop = FALSE]
    keep <- tryCatch(
      cast_select(
        dag, sub_df, response = response,
        min_vars = min_keep, min_fraction = 0,
        num_trees = num_trees, seed = if (!is.null(seed)) seed + ii else NULL,
        verbose = FALSE, stability_reps = 0L
      )$selected,
      error = function(e) character()
    )
    subset_counts[intersect(selected, keep)] <-
      subset_counts[intersect(selected, keep)] + 1L
    test_rows[[paste0("subset_", ii)]] <- data.frame(
      test = "data_subset",
      replicate = ii,
      overlap_fraction = length(intersect(selected, keep)) / max(1L, length(selected)),
      placebo_selected = NA,
      stringsAsFactors = FALSE
    )
  }

  placebo_hits <- 0L
  random_hits <- 0L
  for (ii in seq_len(reps)) {
    if (!is.null(seed)) set.seed(seed + 100L + ii)
    perm_df <- data
    perm_df[[response]] <- sample(perm_df[[response]])
    keep_perm <- tryCatch(rf_screen(perm_df), error = function(e) character())
    placebo_hits <- placebo_hits + length(intersect(selected, keep_perm))
    test_rows[[paste0("placebo_", ii)]] <- data.frame(
      test = "permuted_response",
      replicate = ii,
      overlap_fraction = length(intersect(selected, keep_perm)) / max(1L, length(selected)),
      placebo_selected = NA,
      stringsAsFactors = FALSE
    )

    rand_df <- data
    rand_df$cast_random_common_cause <- stats::rnorm(nrow(rand_df))
    keep_rand <- tryCatch(
      rf_screen(rand_df, extra_vars = "cast_random_common_cause"),
      error = function(e) character()
    )
    random_hits <- random_hits + as.integer("cast_random_common_cause" %in% keep_rand)
    test_rows[[paste0("random_", ii)]] <- data.frame(
      test = "random_common_cause",
      replicate = ii,
      overlap_fraction = length(intersect(selected, keep_rand)) / max(1L, length(selected)),
      placebo_selected = "cast_random_common_cause" %in% keep_rand,
      stringsAsFactors = FALSE
    )
  }

  if (all(c("lon", "lat") %in% names(data))) {
    block_id <- interaction(
      data$lon > stats::median(data$lon, na.rm = TRUE),
      data$lat > stats::median(data$lat, na.rm = TRUE),
      drop = TRUE
    )
    blocks <- levels(block_id)
    jj <- 0L
    for (b in blocks) {
      block_df <- data[block_id == b, , drop = FALSE]
      if (nrow(block_df) < 30L || length(unique(block_df[[response]])) < 2L) next
      jj <- jj + 1L
      keep_block <- tryCatch(
        cast_select(
          dag, block_df, response = response,
          min_vars = min_keep, min_fraction = 0,
          num_trees = num_trees, seed = if (!is.null(seed)) seed + 500L + jj else NULL,
          verbose = FALSE, stability_reps = 0L
        )$selected,
        error = function(e) character()
      )
      block_counts[intersect(selected, keep_block)] <-
        block_counts[intersect(selected, keep_block)] + 1L
      test_rows[[paste0("spatial_block_", jj)]] <- data.frame(
        test = "spatial_block",
        replicate = jj,
        overlap_fraction = length(intersect(selected, keep_block)) / max(1L, length(selected)),
        placebo_selected = NA,
        stringsAsFactors = FALSE
      )
    }
  }

  tests <- if (length(test_rows)) do.call(rbind, test_rows) else data.frame()
  rownames(tests) <- NULL
  block_den <- max(1L, sum(tests$test == "spatial_block", na.rm = TRUE))
  summary <- data.frame(
    variable = selected,
    subset_frequency = as.numeric(subset_counts[selected]) / reps,
    spatial_block_frequency = as.numeric(block_counts[selected]) / block_den,
    stringsAsFactors = FALSE
  )
  if (!is.null(screen$roles) && nrow(screen$roles) > 0) {
    summary <- merge(
      screen$roles,
      summary,
      by = "variable",
      all.y = TRUE,
      sort = FALSE
    )
  }

  settings <- list(
    reps = reps,
    subset_fraction = subset_fraction,
    num_trees = num_trees,
    permuted_response_mean_overlap =
      mean(tests$overlap_fraction[tests$test == "permuted_response"], na.rm = TRUE),
    random_common_cause_hit_rate =
      random_hits / reps
  )

  if (verbose) {
    cli::cli_inform(
      "Screen refutation: permuted-response overlap={round(settings$permuted_response_mean_overlap, 3)}, random-covariate hit rate={round(settings$random_common_cause_hit_rate, 3)}."
    )
  }

  new_cast_refute(tests = tests, summary = summary, settings = settings)
}
