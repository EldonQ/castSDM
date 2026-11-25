# Spatial designs are inspected before fitting and reused without changing
# their scale, membership or requested fold count to improve model scores.
.cast_cv_plan <- function(data, k, block_method, block_size, coords_crs,
                          buffer, buffer_unit, n_repeat, seed,
                          response = "presence", fold_ids = NULL) {
  env <- NULL
  if (is.null(fold_ids) && identical(block_method, "env")) {
    vars <- setdiff(names(data), c("lon", "lat", response))
    vars <- vars[vapply(data[vars], is.numeric, logical(1))]
    if (length(vars) < 2L) cli::cli_abort("Environmental blocking requires two numeric predictors.")
    env <- as.matrix(data[, vars, drop = FALSE])
  }
  if (!is.null(fold_ids)) {
    if (!is.list(fold_ids)) {
      if (n_repeat != 1L) cli::cli_abort("A fixed fold vector requires n_repeat = 1; supply one vector per repeat otherwise.")
      fold_ids <- list(fold_ids)
    }
    if (length(fold_ids) != n_repeat)
      cli::cli_abort("Supply one fold_ids vector per repeat.")
  }
  folds <- vector("list", n_repeat)
  rows <- list()
  train_indices <- vector("list", n_repeat)
  for (ri in seq_len(n_repeat)) {
    rs <- if (is.null(seed)) NULL else seed + (ri - 1L) * 7919L
    ids <- if (is.null(fold_ids)) make_spatial_folds(data$lon, data$lat, k,
      block_method, rs, env = env, y = data[[response]], block_size = block_size,
      coords_crs = coords_crs) else fold_ids[[ri]]
    if (!is.numeric(ids) || length(ids) != nrow(data) || any(!is.finite(ids)) ||
        any(ids != as.integer(ids)) || !setequal(unique(ids), seq_len(k)))
      cli::cli_abort("Each fold_ids vector must assign every row to 1:k, with all folds present.")
    design <- attr(ids, "design")
    if (is.null(design)) design <- list(method = "external", size_source = "external")
    ids <- as.integer(ids); attr(ids, "design") <- design
    folds[[ri]] <- ids
    train_indices[[ri]] <- lapply(seq_len(k), function(i)
      .cast_buffer_train_idx(data$lon, data$lat, which(ids == i), buffer, buffer_unit))
    for (i in seq_len(k)) {
      tr <- data[[response]][train_indices[[ri]][[i]]]
      te <- data[[response]][ids == i]
      train_both <- all(c(0L, 1L) %in% unique(tr))
      test_both <- all(c(0L, 1L) %in% unique(te))
      rows[[length(rows) + 1L]] <- data.frame(replicate = ri, fold = i,
        train_presence = sum(tr == 1L), train_background = sum(tr == 0L),
        test_presence = sum(te == 1L), test_background = sum(te == 0L),
        can_fit = train_both, discrimination_defined = test_both,
        status = if (!train_both) "single_class_training" else
          if (!test_both) "single_class_test" else "two_class_test",
        block_size_m = design$block_size_m %||% NA_real_,
        size_source = design$size_source %||% design$effective %||% "external",
        occupied_blocks = design$occupied_blocks %||% NA_integer_,
        range_warnings = paste(design$range_warnings, collapse = "; "))
    }
  }
  list(folds = folds, train_indices = train_indices, balance = do.call(rbind, rows))
}

.cast_cv_check_plan <- function(plan, verbose = TRUE) {
  if (verbose) for (ri in seq_along(plan$folds)) {
    z <- plan$balance[plan$balance$replicate == ri, , drop = FALSE]
    scale <- if (is.finite(z$block_size_m[1L]))
      paste0(round(z$block_size_m[1L]), " m (", z$size_source[1L], ")") else z$size_source[1L]
    occupied <- if (is.na(z$occupied_blocks[1L])) "" else paste0("; ", z$occupied_blocks[1L], " occupied blocks")
    cli::cli_inform("CV design repeat {ri}: {scale}{occupied}; test presences [{paste(z$test_presence, collapse = ', ')}]; backgrounds [{paste(z$test_background, collapse = ', ')}].")
  }
  bad <- plan$balance[!plan$balance$can_fit, , drop = FALSE]
  if (nrow(bad)) cli::cli_abort(
    "Spatial CV cannot fit {nrow(bad)} fold(s): the training side lacks a response class after buffering. Inspect the saved design; no fold count or scale is changed automatically.",
    class = "cast_cv_invalid_design", fold_balance = plan$balance)
  single <- plan$balance[!plan$balance$discrimination_defined, , drop = FALSE]
  if (nrow(single) && verbose) cli::cli_inform(
    "CV design: {nrow(single)} single-class test fold(s) will be predicted. AUC/TSS are undefined there; Brier/log loss remain available. Metric summaries report their own valid fold counts.")
  invisible(plan)
}

.cast_cv_fold_complete <- function(value, models) {
  identical(value$status, "evaluated") && all(vapply(models, function(m) {
    x <- value$updates[[m]]
    !is.null(x) && length(x$idx) > 0L && length(x$pred) == length(x$idx) &&
      all(is.finite(x$pred)) && length(x$native) == length(x$idx) &&
      all(is.finite(x$native))
  }, logical(1)))
}
