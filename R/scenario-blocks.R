# One content-verified, bounded-memory copy of the observed current inputs.
.cast_current_blocks_key <- function(inputs, rows, variables) {
  .cast_digest(list(inputs, rows, variables, "current_blocks/2",
    .cast_code_identity(.cast_raster_cells), .cast_code_identity(.cast_raster_block)))
}

.cast_write_current_block <- function(value, path, key) {
  # Keep the regular receipt format readable by earlier releases. Wide
  # numeric blocks use faster lossless compression without changing keys.
  payload <- serialize(value, NULL, version = 3L)
  receipt <- list(schema = "cast_cache/3", key = key, payload = payload,
    digest = digest::digest(payload, algo = "sha256", serialize = FALSE), files = character())
  temp <- tempfile(".cast-block-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temp), add = TRUE)
  con <- gzfile(temp, "wb", compression = 1L)
  tryCatch(saveRDS(receipt, con, version = 3L), finally = close(con))
  .cast_publish_raster(temp, path, overwrite = TRUE)
}

.cast_compare_current_blocks <- function(plan, block_size, directory, key, verbose) {
  rows <- max(1L, as.integer(block_size/terra::ncol(plan$current)))
  starts <- seq.int(1L, terra::nrow(plan$current), by = rows)
  key <- .cast_current_blocks_key(key, rows, plan$vars)
  receipt <- file.path(directory, "manifest.rds")
  old <- .cast_cache_read(receipt, key)
  if (!is.null(old)) {
    old$directory <- directory
    if (verbose) cli::cli_inform("Current environmental blocks: verified cache; no source-layer reread.")
    attr(old, "cache_hit") <- TRUE
    return(old)
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(directory, paste0("rows_", starts, ".rds"))
  old_context <- options(castSDM.read_context = "prepare current environmental blocks")
  on.exit(options(old_context), add = TRUE)
  for (i in seq_along(starts)) {
    first <- starts[i]; count <- min(rows, terra::nrow(plan$current)-first+1L)
    block_key <- .cast_digest(list(key, first, count))
    if (is.null(.cast_cache_read(paths[i], block_key))) {
      active <- .cast_mask_rows(plan$mask, first, count, count*terra::ncol(plan$current))
      indices <- which(active)
      values <- as.data.frame(stats::setNames(rep(list(numeric()), length(plan$vars)), plan$vars))
      if (length(indices)) {
        values <- .cast_raster_cells(plan$current, first, count, indices)
        finite <- .cast_finite_rows(values)
        indices <- indices[finite]; values <- values[finite, , drop = FALSE]
      }
      .cast_write_current_block(list(indices = indices, values = values), paths[i], block_key)
    }
    if (verbose && (i == 1L || i %% 10L == 0L || i == length(starts)))
      cli::cli_inform("Current environmental blocks: {i}/{length(starts)} ready; reused across all scenarios.")
  }
  value <- list(directory = directory, key = key, rows = rows, starts = starts,
    paths = basename(paths), variables = plan$vars)
  .cast_cache_write(value, receipt, key, paths)
  stale <- setdiff(list.files(directory, pattern = "^rows_[0-9]+[.]rds$", full.names = TRUE), paths)
  if (length(stale)) unlink(stale)
  attr(value, "cache_hit") <- FALSE
  value
}

.cast_compare_current_block <- function(blocks, first, count, current = NULL, mask = NULL) {
  if (is.null(blocks)) {
    indices <- which(.cast_mask_rows(mask, first, count, count*terra::ncol(current)))
    values <- as.data.frame(stats::setNames(rep(list(numeric()), terra::nlyr(current)), names(current)))
    if (length(indices)) {
      values <- .cast_raster_cells(current, first, count, indices)
      finite <- .cast_finite_rows(values)
      indices <- indices[finite]; values <- values[finite, , drop = FALSE]
    }
    return(list(indices = indices, values = values))
  }
  i <- match(first, blocks$starts)
  if (is.na(i)) cli::cli_abort("Current environmental block layout does not match the scenario layout.")
  path <- file.path(blocks$directory, blocks$paths[i])
  value <- .cast_cache_read(path, .cast_digest(list(blocks$key, first, count)))
  if (is.null(value) || !identical(names(value$values), blocks$variables) ||
      length(value$indices) != nrow(value$values))
    cli::cli_abort("Current environmental block failed integrity checks: {.path {path}}. Rerun the unchanged comparison to rebuild it.")
  value
}

.cast_compare_progress <- function(path, scenario, phase, first = NA_integer_, total = NA_integer_) {
  if (is.null(path)) return(invisible(NULL))
  .cast_atomic_rds(list(scenario = scenario, phase = phase, row = first,
    total_rows = total, updated = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")), path)
  invisible(NULL)
}

.cast_compare_input_stage <- function(expr, stage, directory) {
  path <- if (!is.null(directory)) file.path(directory, "failures", paste0("stage_", stage, ".rds"))
  value <- tryCatch(force(expr), error = function(e) {
    if (!is.null(path)) {
      record <- list(stage = stage, message = conditionMessage(e),
        row = e$row, nrows = e$nrows, sources = e$sources, warnings = e$warnings,
        finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
      # A logging failure must not replace the original I/O error.
      try({
        .cast_atomic_rds(record, path)
        if (is.data.frame(record$sources)) utils::write.csv(record$sources,
          sub("[.]rds$", ".csv", path), row.names = FALSE)
      }, silent = TRUE)
    }
    stop(e)
  })
  if (!is.null(path)) {
    stale <- c(path, sub("[.]rds$", ".csv", path))
    unlink(stale[file.exists(stale)])
  }
  value
}
