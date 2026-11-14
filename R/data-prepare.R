#' VIF-Based Collinearity Screening
#'
#' Iteratively removes environmental variables with Variance Inflation Factor
#' (VIF) above the specified threshold. Optionally applies expert-guided
#' pre-filtering before VIF elimination.
#'
#' @param data A `data.frame` of environmental variables (numeric columns
#'   only). Rows are sites, columns are variables.
#' @param threshold Numeric. VIF threshold above which variables are removed.
#'   Default `10` (Zuur et al. 2010).
#' @param exclude Character vector of column names to exclude from screening
#'   (e.g., coordinates, presence). Default `NULL`.
#' @param expert_filter Character vector of variable names to remove before
#'   VIF screening (domain-knowledge redundancies). Default `NULL`.
#' @param verbose Logical. Print iteration details. Default `TRUE`.
#'
#' @details
#' Elimination continues until the declared VIF criterion is met or only one
#' predictor remains; there is no environmental-variable count floor.
#'
#' @return A list with components:
#' \describe{
#'   \item{`selected`}{Character vector of retained variable names.}
#'   \item{`removed`}{Character vector of removed variable names (in order).}
#'   \item{`vif_log`}{A `data.frame` with iteration-by-iteration VIF values;
#'   `removed` is `NA` on the converged final iteration (nothing removed).}
#'   \item{`data`}{Filtered `data.frame` with only retained variables.}
#' }
#'
#' @references
#' Zuur, A.F., Ieno, E.N., & Elphick, C.S. (2010). A protocol for data
#' exploration to avoid common statistical problems. *Methods in Ecology and
#' Evolution*, 1(1), 3-14.
#'
#' @export
cast_vif <- function(data,
                     threshold = 10,
                     exclude = NULL,
                     expert_filter = NULL,
                     verbose = TRUE) {
  threshold <- .cast_check_scalar(threshold, "threshold", lower = 1)
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data.frame.")
  }

  # Identify numeric columns only
  all_cols <- names(data)
  exclude <- exclude %||% character(0)
  num_cols <- setdiff(all_cols, exclude)
  num_cols <- num_cols[vapply(data[num_cols], is.numeric, logical(1))]

  # Stage 1: Expert pre-filter
  removed <- character(0)
  if (!is.null(expert_filter)) {
    hit <- intersect(expert_filter, num_cols)
    if (length(hit) > 0 && verbose) {
      cli::cli_inform(
        "Expert pre-filter removes {length(hit)} variable{?s}."
      )
    }
    removed <- hit
    num_cols <- setdiff(num_cols, hit)
  }

  # Complete cases for lm
  dat_clean <- stats::na.omit(data[, num_cols, drop = FALSE])
  if (nrow(dat_clean) < 10) {
    cli::cli_abort("Fewer than 10 complete cases after NA removal.")
  }
  nonfinite <- names(dat_clean)[vapply(
    dat_clean, function(v) any(!is.finite(v)), logical(1)
  )]
  if (length(nonfinite)) {
    cli::cli_abort(c(
      "Non-finite values (NA/NaN/Inf) in predictor{?s}: {.val {nonfinite}}.",
      i = "Clean or impute these columns before VIF screening."
    ))
  }

  # Stage 2: Iterative VIF elimination
  current_vars <- num_cols
  vif_log <- list()
  iteration <- 0L

  repeat {
    iteration <- iteration + 1L
    if (length(current_vars) < 2) break

    # Standard VIF: regress each observed predictor on the other predictors.
    # No manufactured response, row-index regression or absolute-unit cutoff.
    tmp <- scale(as.matrix(dat_clean[, current_vars, drop = FALSE]))
    varying <- vapply(dat_clean[, current_vars, drop = FALSE], .cast_is_varying, logical(1))
    if (any(!varying)) {
      removed <- c(removed, current_vars[!varying])
      current_vars <- current_vars[varying]
      next
    }
    vif_vals <- vapply(seq_along(current_vars), function(j) {
      residual <- stats::lm.fit(cbind(1, tmp[, -j, drop = FALSE]), tmp[, j])$residuals
      rss <- sum(residual^2); tss <- sum(tmp[, j]^2)
      if (rss <= .Machine$double.eps*tss) Inf else tss/rss
    }, numeric(1))
    names(vif_vals) <- current_vars

    max_vif <- max(vif_vals)
    max_var <- names(which.max(vif_vals))

    converged <- max_vif <= threshold
    vif_log[[iteration]] <- data.frame(
      iteration = iteration,
      n_vars = length(current_vars),
      max_vif = round(max_vif, 2),
      # NA on the converged final iteration: nothing was removed there.
      removed = if (converged) NA_character_ else max_var,
      stringsAsFactors = FALSE
    )

    if (converged) {
      if (verbose) {
        n_remain <- length(current_vars)
        cli::cli_inform(
          "VIF converged: all VIF <= {threshold} ({n_remain} vars)."
        )
      }
      break
    }

    if (verbose) {
      cli::cli_inform(
        "Iter {iteration}: VIF = {round(max_vif, 1)} -> remove {.val {max_var}}"
      )
    }
    removed <- c(removed, max_var)
    current_vars <- setdiff(current_vars, max_var)

  }

  log_df <- if (length(vif_log) > 0) {
    do.call(rbind, vif_log)
  } else {
    data.frame(
      iteration = integer(0), n_vars = integer(0),
      max_vif = numeric(0), removed = character(0)
    )
  }

  list(
    selected = current_vars,
    removed = removed,
    vif_log = log_df,
    data = data[, current_vars, drop = FALSE]
  )
}


#' Prepare Species Data for Modeling
#'
#' Validates and prepares species occurrence data for the CAST pipeline.
#' Checks required columns, handles missing values, and splits into
#' training/test sets.
#'
#' Auto-detects environmental variable columns by excluding coordinates,
#' response, known metadata columns, and constant columns. If the
#' detection is incorrect, pass `env_vars` explicitly to override.
#'
#' The hold-out split defaults to a **spatial block** split (whole spatial
#' blocks are assigned to the test set), which gives a more honest estimate of
#' transferable performance than a random split by reducing spatial
#' autocorrelation leakage between train and test (Roberts et al. 2017). When
#' coordinates are unavailable or a spatial split cannot preserve both
#' classes, it stops with a design diagnostic. A non-spatial split must be
#' requested explicitly; it is never substituted silently.
#'
#' @param data A `data.frame` with columns: `lon`, `lat`, `presence` (0/1),
#'   and environmental variables.
#' @param train_fraction Numeric in (0, 1), strictly. Fraction for the
#'   training set. Values outside the open interval are rejected.
#'   Default `0.7`.
#' @param seed Integer or `NULL`. Random seed for reproducible splitting.
#' @param env_vars Character vector or `NULL`. Environmental variable names
#'   to use. If `NULL` (default), auto-detected from `data`. Useful when
#'   auto-detection picks up metadata columns (e.g., occurrence indicators).
#' @param split Hold-out strategy: `"spatial"` (default), `"stratified"`, or
#'   `"random"`. Degenerate spatial/stratified splits fall back automatically.
#' @param block_method Spatial blocking for `split = "spatial"`: `"blockcv"`
#'   (default; \pkg{blockCV} square blocks with autocorrelation-range sizing),
#'   `"grid"` (k-means on grid-cell centroids; dependency-free fallback), or
#'   `"env"` (k-means on the scaled predictor matrix, so blocks are contiguous
#'   in environmental space; Roberts et al. 2017-style environmental
#'   blocking). See [cast_cv()].
#' @param n_blocks Number of groups for explicit grid/environmental blocking,
#'   default `20`. blockCV hold-out uses its calibrated physical blocks.
#' @param coords_crs Coordinate CRS, required for projected coordinates.
#' @param block_size Optional physical block size in metres.
#' @param verbose Logical. Print detected variables and excluded columns.
#'   Default `TRUE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`train`}{Training `data.frame`.}
#'   \item{`test`}{Test `data.frame`.}
#'   \item{`env_vars`}{Character vector of environmental variable names.}
#'   \item{`split`}{Character describing the hold-out strategy actually used.}
#' }
#'
#' @references
#' Roberts, D. R. et al. (2017). Cross-validation strategies for data with
#' temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*,
#' 40(8), 913-929.
#'
#' @export
cast_prepare <- function(data, train_fraction = 0.7, seed = NULL,
                         env_vars = NULL,
                         split = c("spatial", "stratified", "random"),
                         block_method = c("blockcv", "grid", "env"),
                         n_blocks = 20L, verbose = TRUE, coords_crs = NULL, block_size = NULL) {
  split <- match.arg(split)
  block_method <- match.arg(block_method)
  # Documented as a fraction; unenforced, 0 built an empty training set and 1
  # produced all-NA evaluation metrics instead of an error.
  train_fraction <- .cast_check_scalar(
    train_fraction, "train_fraction", lower = 0, upper = 1,
    min_inclusive = FALSE, max_inclusive = FALSE
  )
  n_blocks <- .cast_check_scalar(n_blocks, "n_blocks", lower = 1, integer = TRUE)
  validate_species_data(data)

  if (!is.null(env_vars)) {
    # User-supplied: validate they all exist
    missing_cols <- setdiff(env_vars, names(data))
    if (length(missing_cols) > 0) {
      cli::cli_abort(
        "Columns specified in {.arg env_vars} not found in data: {.val {missing_cols}}."
      )
    }
    non_numeric <- env_vars[!vapply(
      data[env_vars], is.numeric, logical(1)
    )]
    if (length(non_numeric) > 0) {
      cli::cli_warn(
        "Non-numeric columns in {.arg env_vars} will be ignored: {.val {non_numeric}}."
      )
      env_vars <- setdiff(env_vars, non_numeric)
    }
    declared_vars <- env_vars
    final_vars <- env_vars
    if (verbose) {
      cli::cli_inform(c(
        "v" = "Using {length(final_vars)} user-specified environmental variable{?s}:",
        " " = "{.val {final_vars}}"
      ))
    }
  } else {
    # Auto-detect
    all_numeric <- names(data)[vapply(data, is.numeric, logical(1))]
    final_vars  <- get_env_vars(data)

    # Compute what was excluded and why
    non_numeric_cols <- setdiff(names(data), all_numeric)
    meta_excluded    <- setdiff(all_numeric, c(final_vars, "presence"))
    const_excluded   <- meta_excluded[vapply(meta_excluded, function(v) {
      vals <- data[[v]][!is.na(data[[v]])]
      !.cast_is_varying(vals)
    }, logical(1))]
    name_excluded <- setdiff(meta_excluded, const_excluded)

    if (length(final_vars) == 0) {
      # Show diagnostics to help the user fix it
      cli::cli_abort(c(
        "No environmental variables detected in {.arg data}.",
        "i" = "All numeric columns: {.val {all_numeric}}",
        "i" = "Excluded (non-numeric): {.val {non_numeric_cols}}",
        "i" = "Excluded (known metadata names): {.val {name_excluded}}",
        "i" = "Excluded (constant or all missing): {.val {const_excluded}}",
        ">" = "Pass {.arg env_vars} to override: {.code cast_prepare(data, env_vars = c(...))}"
      ))
    }

    if (verbose) {
      cli::cli_inform(c(
        "v" = "Detected {length(final_vars)} environmental variable{?s}:",
        " " = "{.val {final_vars}}"
      ))
      if (length(name_excluded) > 0) {
        cli::cli_inform(c(
          "i" = "Excluded (metadata name): {.val {name_excluded}}"
        ))
      }
      if (length(const_excluded) > 0) {
        cli::cli_inform(c(
          "i" = "Excluded (constant in the supplied sample): {.val {const_excluded}}"
        ))
      }
      cli::cli_inform(c(
        "i" = "If incorrect, override with: {.code cast_prepare(data, env_vars = c(...))}"
      ))
    }
  }

  n <- nrow(data)
  sp <- .cast_with_seed(seed,
    .cast_holdout_split(data, train_fraction, split, block_method,
                        as.integer(n_blocks), "presence", verbose, coords_crs, block_size))
  train_idx <- sp$train
  if (verbose) {
    cli::cli_inform(
      "Hold-out split: {sp$method} ({length(train_idx)}/{n} in training)."
    )
  }

  # Variance screening runs on the TRAINING rows, after the split: a layer can
  # vary across the full study area yet be constant inside the hold-out
  # training sample (e.g. land-cover patch statistics over a uniform plateau
  # window). Screening before the split would keep such a column here and then
  # lose it again in get_env_vars() downstream, so var_id/category vectors
  # built for the declared set would no longer match the predictors seen by
  # cast_select(). One rule for both branches: declared = upper bound,
  # final_vars = what the training rows actually support.
  const_in_train <- final_vars[vapply(final_vars, function(v) {
    vals <- data[[v]][train_idx]
    vals <- vals[!is.na(vals)]
    !.cast_is_varying(vals)
  }, logical(1))]
  if (length(const_in_train) > 0L) {
    if (verbose) {
      cli::cli_inform(c(
        "i" = "{length(const_in_train)} declared variable{?s} constant in the training sample and cannot be fitted; dropping:",
        " " = "{.val {const_in_train}}"
      ))
    }
    final_vars <- setdiff(final_vars, const_in_train)
  }

  # env_vars is a contract, not a label: when the caller names the predictors,
  # the returned splits carry exactly those that survived (plus coordinates
  # and the response). Before this, `train`/`test` were row-slices of the FULL
  # data frame, so a caller that pre-filtered the pool (e.g. temporally
  # matched layers) passed var_id/category vectors that no longer matched
  # get_env_vars(split$train).
  keep_cols <- unique(c(
    final_vars,
    intersect(c("lon", "lat", "presence"), names(data))
  ))

  list(
    train    = data[train_idx, keep_cols, drop = FALSE],
    test     = data[-train_idx, keep_cols, drop = FALSE],
    train_index = train_idx,
    test_index = setdiff(seq_len(nrow(data)), train_idx),
    env_vars = final_vars,
    split    = sp$method,
    split_type = sp$type,
    fold_design = sp$design
  )
}

#' Hold-out index construction with spatial / stratified fallbacks
#' @keywords internal
#' @noRd
.cast_holdout_split <- function(data, train_fraction, split, block_method,
                                n_blocks, response, verbose, coords_crs = NULL, block_size = NULL) {
  n <- nrow(data)
  has_coords <- all(c("lon", "lat") %in% names(data)) &&
    is.numeric(data$lon) && is.numeric(data$lat)
  y <- if (response %in% names(data)) data[[response]] else NULL
  two_class <- !is.null(y) && length(unique(y[!is.na(y)])) >= 2L

  make_random <- function() sample.int(n, size = round(train_fraction * n))
  make_stratified <- function() {
    if (!two_class) return(make_random())
    idx <- integer(0)
    for (cl in unique(y[!is.na(y)])) {
      pos <- which(y == cl)
      take <- max(1L, round(train_fraction * length(pos)))
      idx <- c(idx, if (length(pos) == 1L) pos else sample(pos, size = take))
    }
    sort(idx)
  }

  if (identical(split, "random")) {
    return(list(train = make_random(), method = "random", type = "random"))
  }
  if (identical(split, "stratified")) {
    return(list(train = make_stratified(), method = "stratified", type = "stratified"))
  }

  # split == "spatial" (block_method "env" partitions in environmental space)
  if (!has_coords) {
    cli::cli_abort("Spatial hold-out requires coordinates. Supply lon/lat or explicitly choose a non-spatial split.")
  }
  env_block <- if (identical(block_method, "env")) {
    cand <- setdiff(names(data), c("lon", "lat", response))
    cand <- cand[vapply(data[cand], is.numeric, logical(1))]
    as.matrix(data[, cand, drop = FALSE])
  } else NULL
  blocks <- make_spatial_folds(data$lon, data$lat, k = if (block_method == "blockcv") 2L else n_blocks,
                       method = block_method, env = env_block, y = y,
                       block_size = block_size, coords_crs = coords_crs)
  if (is.null(blocks) || length(unique(blocks)) < 2L) {
    cli::cli_abort("Spatial blocking failed; no random/stratified hold-out is substituted. Review coordinates and the spatial design.")
  }
  design <- attr(blocks, "design")
  # Hold out physical blocks, not arbitrary groups formed only to reach k.
  # blockCV's two-fold assignment is used only to construct its block geometry.
  if (block_method == "blockcv") blocks <- design$block_ids
  ub <- sample(unique(blocks))
  target_test <- round((1 - train_fraction) * n)
  # Greedy block pick minimising |test size - target|: a block is added only
  # when it brings the running total closer to the target. Whole blocks can
  # make the achieved fraction differ from the requested one.
  test_blocks <- integer(0)
  cur_test <- 0L
  for (b in ub) {
    b_n <- sum(blocks == b)
    if (abs((cur_test + b_n) - target_test) < abs(cur_test - target_test)) {
      test_blocks <- c(test_blocks, b)
      cur_test <- cur_test + b_n
    }
  }
  test_idx <- which(blocks %in% test_blocks)
  train_idx <- setdiff(seq_len(n), test_idx)
  ok_classes <- if (two_class) {
    length(unique(y[train_idx])) >= 2L && length(unique(y[test_idx])) >= 2L
  } else TRUE
  if (!length(train_idx) || !length(test_idx) || !ok_classes) {
    cli::cli_abort("The spatial hold-out cannot preserve both response classes. Review the spatial design; no random split is substituted.")
  }
  if (inherits(design$blocks, "sf")) {
    design$blocks$folds <- NULL
    design$blocks$split <- ifelse(seq_len(nrow(design$blocks)) %in% test_blocks, "holdout", "training")
  }
  list(train = sort(train_idx),
       method = sprintf("%s block (%s)",
                        if (identical(block_method, "env")) "environmental" else "spatial",
                        (design$effective %||% block_method)),
       design = c(design, list(test_blocks = test_blocks)), type = "spatial")
}
