#' Run a Pipeline Step with Checkpoint, Timing, and Peak RAM Logging
#'
#' Internal helper used inside [.cast_batch_run_one_species()] to make the
#' per-species pipeline resumable. Results are cached to
#' `<output_dir>/<species>/.steps/<step>.rds`, so a re-invocation of
#' [cast_batch()] over the same `output_dir` skips already-finished steps.
#' Elapsed seconds and peak RAM (in MiB) are appended to
#' `<output_dir>/resource_log.csv`.
#'
#' Lazy semantics: `expr` is only evaluated when no cache hit exists. This
#' is what makes the wrapper non-invasive — wrapping `cast_run_step("dag",
#' ..., cast_dag(...))` does not call `cast_dag()` if its output is on disk.
#'
#' @param step_name Character. Step identifier (filename-safe).
#' @param output_dir Top-level batch directory.
#' @param species Species name (used as subdirectory under `output_dir`).
#' @param expr Lazy expression producing the step's result.
#' @param verbose Logical. Print a `cli_inform` line on cache hit / miss.
#'
#' @return The step's result (either freshly computed or loaded from cache).
#' @keywords internal
#' @noRd
cast_run_step <- function(step_name, output_dir, species, expr,
                          verbose = FALSE) {
  ckpt_dir <- file.path(output_dir, species, ".steps")
  dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)
  ckpt <- file.path(ckpt_dir, paste0(step_name, ".rds"))

  if (file.exists(ckpt)) {
    val <- tryCatch(readRDS(ckpt), error = function(e) NULL)
    if (!is.null(val)) {
      if (verbose)
        cli::cli_inform("  [{species}] {step_name}: cache hit; skipped.")
      return(val)
    }
  }

  has_pram <- requireNamespace("peakRAM", quietly = TRUE)
  t0 <- Sys.time()
  if (has_pram) {
    pram <- tryCatch(
      peakRAM::peakRAM({ val <- expr }),
      error = function(e) NULL
    )
    peak_mb <- if (is.null(pram)) NA_real_ else
      suppressWarnings(as.numeric(pram$Peak_RAM_Used_MiB[1]))
  } else {
    invisible(gc(reset = TRUE))
    val <- expr
    g <- tryCatch(gc(), error = function(e) NULL)
    peak_mb <- if (is.null(g)) NA_real_ else
      suppressWarnings(sum(g[, "max used (Mb)"], na.rm = TRUE))
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (!is.null(val)) {
    saveRDS(val, ckpt)
  }
  cast_log_resource(output_dir, species, step_name, elapsed, peak_mb)
  if (verbose)
    cli::cli_inform(
      "  [{species}] {step_name}: {round(elapsed,2)}s, peak {round(peak_mb,1)} MiB."
    )
  val
}


#' Append a Row to the Batch Resource Log
#'
#' Writes a single CSV row to `<output_dir>/resource_log.csv`. The file is
#' created with a header on first use; subsequent calls append. Rows
#' contain timestamp, species, step, elapsed seconds, and peak RAM (MiB).
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
  log_path <- file.path(output_dir, "resource_log.csv")
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
