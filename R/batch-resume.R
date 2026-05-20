#' Resume an Interrupted `cast_batch()` Run
#'
#' Scans `output_dir` for species that already have a saved
#' `cast_result.rds` and re-runs [cast_batch()] only on the remainder.
#' Combined with the per-step `.steps/<name>.rds` checkpoints written by
#' [cast_batch()], this lets you recover from crashes, OOM kills, or
#' hard reboots without recomputing finished work.
#'
#' Usage:
#' \itemize{
#'   \item Same `species_list` and `output_dir` you originally passed to
#'         [cast_batch()].
#'   \item Same `...` config arguments (changing them silently invalidates
#'         step caches; in that case clear `output_dir` first).
#' }
#'
#' @param output_dir Top-level directory used by the original
#'   [cast_batch()] call.
#' @param species_list Original named list of species data.frames.
#' @param env_data Optional shared environmental data.frame.
#' @param force Logical. If `TRUE`, ignore existing `cast_result.rds`
#'   files and re-run all species (still benefits from per-step caches).
#'   Default `FALSE`.
#' @param verbose Logical. Default `TRUE`.
#' @param ... Forwarded to [cast_batch()] (must match the original call).
#'
#' @return A `cast_batch` object covering the species that were re-run.
#'   Already-finished species are not re-loaded into the result; load
#'   them manually from `<output_dir>/<species>/cast_result.rds`.
#'
#' @seealso [cast_batch()]
#' @export
cast_batch_resume <- function(output_dir,
                              species_list,
                              env_data = NULL,
                              force    = FALSE,
                              verbose  = TRUE,
                              ...) {
  if (!dir.exists(output_dir))
    cli::cli_abort("Output directory does not exist: {.path {output_dir}}.")
  if (!is.list(species_list) || is.null(names(species_list)))
    cli::cli_abort("{.arg species_list} must be a named list of data.frames.")

  done <- character(0)
  for (sp in names(species_list)) {
    rds <- file.path(output_dir, sp, "cast_result.rds")
    if (file.exists(rds)) done <- c(done, sp)
  }
  pending <- if (force) names(species_list) else
    setdiff(names(species_list), done)

  if (verbose) {
    cli::cli_h1("cast_batch_resume")
    cli::cli_inform("Output dir : {.path {output_dir}}")
    cli::cli_inform("Complete   : {length(done)} / {length(species_list)}")
    cli::cli_inform("Pending    : {length(pending)}")
    if (length(done))
      cli::cli_inform("  done: {.val {done}}")
    if (length(pending))
      cli::cli_inform("  pending: {.val {pending}}")
  }

  if (length(pending) == 0L) {
    if (verbose)
      cli::cli_inform("All species already complete; nothing to resume.")
    return(invisible(NULL))
  }

  cast_batch(
    species_list = species_list[pending],
    env_data     = env_data,
    output_dir   = output_dir,
    verbose      = verbose,
    ...
  )
}
