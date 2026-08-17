#' Run a Pipeline Step with Checkpoint, Timing, and Peak RAM Logging
#'
#' Internal helper used inside `.cast_batch_run_one_species()` to make the
#' per-species pipeline resumable. Results are cached to
#' `<output_dir>/<species>/.steps/<step>.rds`, so a re-invocation of
#' [cast_batch()] over the same `output_dir` skips already-finished steps.
#' Elapsed seconds and peak RAM (in MiB) are appended to the per-species
#' resource log.
#'
#' Cache entries store the step's `params` next to its value. A cache hit
#' requires the stored `params` to be identical to the current call, so
#' changing any step input (data, configuration, seed) silently invalidates
#' the cache instead of replaying a stale result. `NULL` results are cached
#' too (as a sentinel entry with an empty value), so a step that legitimately
#' produced nothing is not recomputed on every resume. Bare-value cache files
#' written by castSDM < 0.7.0 are treated as misses and overwritten.
#'
#' @param step_name Character. Step identifier (filename-safe).
#' @param output_dir Top-level batch directory.
#' @param species Species name (used as subdirectory under `output_dir`).
#' @param expr Lazy expression producing the step's result.
#' @param params Named list describing the step inputs; the cache key.
#' @param verbose Logical. Print a `cli_inform` line on cache hit / miss.
#'
#' @return The step's result (either freshly computed or loaded from cache).
#' @keywords internal
#' @noRd
cast_run_step <- function(step_name, output_dir, species, expr,
                          params = list(), verbose = FALSE) {
  ckpt_dir <- file.path(output_dir, species, ".steps")
  dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)
  ckpt <- file.path(ckpt_dir, paste0(step_name, ".rds"))

  if (file.exists(ckpt)) {
    val <- tryCatch(readRDS(ckpt), error = function(e) NULL)
    if (!is.null(val) && is.list(val) &&
        all(c("params", "value") %in% names(val)) &&
        identical(val$params, params)) {
      if (verbose) {
        cli::cli_inform("  [{species}] {step_name}: cache hit; skipped.")
      }
      return(val$value)
    }
  }

  has_pram <- requireNamespace("peakRAM", quietly = TRUE)
  t0 <- Sys.time()
  val <- NULL
  expr_err <- NULL
  if (has_pram) {
    pram <- tryCatch(
      peakRAM::peakRAM({ val <- expr }),
      error = function(e) { expr_err <<- e; NULL }
    )
    if (!is.null(expr_err)) {
      cli::cli_abort("Step {.val {step_name}} failed: {conditionMessage(expr_err)}")
    }
    peak_mb <- if (is.null(pram)) NA_real_ else
      suppressWarnings(as.numeric(pram$Peak_RAM_Used_MiB[1]))
  } else {
    invisible(gc(reset = TRUE))
    val <- expr
    g <- tryCatch(gc(), error = function(e) NULL)
    peak_mb <- if (is.null(g)) NA_real_ else tryCatch({
      idx <- which(colnames(g) == "(Mb)" & seq_len(ncol(g)) > 4L)
      if (length(idx)) sum(g[, idx[1]], na.rm = TRUE) else NA_real_
    }, error = function(e) NA_real_)
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Cache even NULL results (sentinel entry): a step that produced nothing
  # should not be recomputed on every resume.
  saveRDS(list(params = params, value = val), ckpt)
  cast_log_resource(output_dir, species, step_name, elapsed, peak_mb)
  if (verbose) {
    cli::cli_inform(
      "  [{species}] {step_name}: {round(elapsed,2)}s, peak {round(peak_mb,1)} MiB."
    )
  }
  val
}


#' Append a Row to the Per-Species Batch Resource Log
#'
#' Writes a single CSV row to `<output_dir>/resource_log_<species>.csv`.
#' One file per species avoids the header/append race that parallel
#' workers would create when sharing a single log file. The per-species
#' files are merged into a single `<output_dir>/resource_log.csv` by
#' `.cast_merge_resource_logs()` after all workers have finished (called
#' by [cast_batch()]).
#'
#' @param output_dir Top-level batch directory.
#' @param species Species name.
#' @param step Step identifier.
#' @param elapsed_sec Numeric. Elapsed wall-clock time in seconds.
#' @param peak_mb Numeric. Peak RAM (MiB), or `NA_real_` if unmeasured.
#'
#' @return Invisibly `TRUE` on success, `FALSE` on failure.
#' @keywords internal
#' @noRd
cast_log_resource <- function(output_dir, species, step,
                              elapsed_sec, peak_mb = NA_real_) {
  safe_sp <- gsub("[^A-Za-z0-9_-]", "_", species)
  log_path <- file.path(output_dir, paste0("resource_log_", safe_sp, ".csv"))
  row <- data.frame(
    timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    species     = species,
    step        = step,
    elapsed_sec = round(as.numeric(elapsed_sec), 3),
    peak_mb     = if (is.na(peak_mb)) NA_real_ else round(peak_mb, 1),
    stringsAsFactors = FALSE
  )
  ok <- tryCatch({
    has_header <- file.exists(log_path)
    utils::write.table(
      row, log_path, append = has_header,
      sep = ",", row.names = FALSE,
      col.names = !has_header, quote = FALSE
    )
    TRUE
  }, error = function(e) FALSE)
  invisible(ok)
}


#' Merge Per-Species Resource Logs into One CSV
#'
#' Combines all `<output_dir>/resource_log_<species>.csv` files written by
#' [cast_log_resource()] into `<output_dir>/resource_log.csv`, sorted by
#' timestamp. Safe to call repeatedly; existing per-species files are kept
#' so interrupted runs can re-merge without losing rows.
#'
#' @param output_dir Top-level batch directory.
#' @return Invisibly `TRUE` if a merged file was written, `FALSE` otherwise.
#' @keywords internal
#' @noRd
.cast_merge_resource_logs <- function(output_dir) {
  files <- list.files(output_dir, pattern = "^resource_log_.*\\.csv$",
                      full.names = TRUE)
  files <- files[!grepl("/resource_log\\.csv$", files)]
  if (!length(files)) return(invisible(FALSE))
  rows <- lapply(files, function(f) {
    tryCatch(
      utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(invisible(FALSE))
  merged <- do.call(rbind, rows)
  merged <- merged[order(merged$timestamp, merged$species, merged$step), ]
  rownames(merged) <- NULL
  utils::write.csv(merged, file.path(output_dir, "resource_log.csv"),
                   row.names = FALSE)
  invisible(TRUE)
}
