#' Run a Full castSDM Pipeline from a YAML Config File
#'
#' One-stop entry point that reads a YAML configuration file, loads
#' per-species CSVs (and optionally a shared environmental grid),
#' constructs a [cast_worker_budget()], and invokes [cast_batch()] (or
#' [cast_batch_resume()] if `run.resume: true` and partial output exists).
#'
#' This is the recommended way to drive reproducible multi-species runs:
#' the YAML acts as a self-documenting record of every parameter, so a
#' run can be re-executed by anyone who has the same data files.
#'
#' Generate a starter YAML with [cast_config_template()].
#'
#' @section YAML schema:
#' Top-level keys recognised by `cast_run_from_config()`:
#' \describe{
#'   \item{`run.output_dir`}{Where outputs are written.}
#'   \item{`run.seed`}{Base seed (each species gets `seed + i`).}
#'   \item{`run.parallel`}{Whether to run species in parallel.}
#'   \item{`run.worker_budget.total_workers`}{Optional integer; passed to
#'     [cast_worker_budget()]. Omit to use `detectCores() - 1`.}
#'   \item{`run.resume`}{If `TRUE`, calls [cast_batch_resume()] when any
#'     species in `output_dir` already has a saved `cast_result.rds`.}
#'   \item{`species`}{Named map. Each entry has a `file:` path to a CSV.}
#'   \item{`env_data.file`}{Optional CSV with shared environmental data
#'     for spatial HSS prediction.}
#'   \item{`models`}{Character vector forwarded to [cast_batch()].}
#'   \item{`cast_batch_args`}{Map of any other named arg accepted by
#'     [cast_batch()].}
#' }
#'
#' All file paths are resolved relative to the directory containing the
#' YAML file.
#'
#' @param config_path Path to a YAML config file.
#' @param ... Additional named arguments that override values from the
#'   YAML (forwarded to [cast_batch()]).
#'
#' @return The `cast_batch` object returned by [cast_batch()] (or
#'   [cast_batch_resume()]).
#'
#' @seealso [cast_batch()], [cast_batch_resume()],
#'   [cast_worker_budget()], [cast_config_template()]
#' @export
cast_run_from_config <- function(config_path, ...) {
  check_suggested("yaml", "for parsing castSDM YAML configs")
  if (!file.exists(config_path))
    cli::cli_abort("Config file not found: {.path {config_path}}.")

  cfg <- yaml::read_yaml(config_path)
  base_dir <- dirname(normalizePath(config_path, winslash = "/", mustWork = TRUE))

  resolve <- function(p) {
    if (is.null(p) || !nzchar(p)) return(NULL)
    if (file.exists(p)) return(p)
    fp <- file.path(base_dir, p)
    if (!file.exists(fp))
      cli::cli_abort("Cannot find file referenced in config: {.path {p}}.")
    fp
  }

  # Species
  if (is.null(cfg$species) || length(cfg$species) == 0L)
    cli::cli_abort("Config must define at least one species under {.field species}.")

  species_list <- lapply(cfg$species, function(spec) {
    if (is.null(spec$file))
      cli::cli_abort("Each species entry needs a {.field file:} key.")
    utils::read.csv(resolve(spec$file), stringsAsFactors = FALSE)
  })
  names(species_list) <- names(cfg$species)

  env_data <- NULL
  if (!is.null(cfg$env_data) && !is.null(cfg$env_data$file)) {
    env_data <- utils::read.csv(resolve(cfg$env_data$file),
                                stringsAsFactors = FALSE)
  }

  # Run options
  run_opts <- cfg$run %||% list()
  output_dir <- run_opts$output_dir %||% "castSDM_output"
  output_dir <- if (file.exists(output_dir) ||
                    grepl("^([A-Za-z]:|/)", output_dir))
    output_dir else file.path(base_dir, output_dir)
  seed <- run_opts$seed
  parallel <- isTRUE(run_opts$parallel %||% TRUE)

  # Worker budget
  budget <- NULL
  if (!is.null(run_opts$worker_budget)) {
    budget <- cast_worker_budget(
      total_workers = run_opts$worker_budget$total_workers,
      n_species     = length(species_list)
    )
    if (parallel && requireNamespace("future", quietly = TRUE))
      cast_setup_species_plan(budget)
  }

  # cast_batch args
  ba <- cfg$cast_batch_args %||% list()
  if (!is.null(cfg$models))   ba$models   <- cfg$models
  ba$species_list <- species_list
  ba$env_data     <- env_data
  ba$output_dir   <- output_dir
  ba$seed         <- seed
  ba$parallel     <- parallel

  # Caller overrides
  ba <- utils::modifyList(ba, list(...))

  resume <- isTRUE(run_opts$resume %||% FALSE)
  any_done <- FALSE
  if (resume && dir.exists(output_dir)) {
    for (sp in names(species_list)) {
      if (file.exists(file.path(output_dir, sp, "cast_result.rds"))) {
        any_done <- TRUE; break
      }
    }
  }

  cli::cli_h1("cast_run_from_config")
  cli::cli_inform("Config     : {.path {config_path}}")
  cli::cli_inform("Species    : {length(species_list)}")
  cli::cli_inform("Output dir : {.path {output_dir}}")
  if (!is.null(budget))
    cli::cli_inform("Workers    : total={budget$total} species={budget$species} intra={budget$intra}")
  cli::cli_inform("Resume     : {resume && any_done}")

  if (resume && any_done) {
    resume_args <- ba
    resume_args$output_dir <- NULL
    do.call(cast_batch_resume, c(list(output_dir = output_dir), resume_args))
  } else {
    do.call(cast_batch, ba)
  }
}


#' Copy the Bundled YAML Config Template to Disk
#'
#' Writes the starter YAML schema (shipped under
#' `inst/templates/castSDM_config.yaml.tmpl`) to a destination of your
#' choosing. Edit the result, then pass it to
#' [cast_run_from_config()].
#'
#' @param destination Path to write the template. Default
#'   `"castSDM_config.yaml"` in the current directory.
#' @param overwrite Logical. Overwrite existing files. Default `FALSE`.
#'
#' @return The path to the written file (invisibly).
#' @export
cast_config_template <- function(destination = "castSDM_config.yaml",
                                 overwrite   = FALSE) {
  src <- system.file("templates", "castSDM_config.yaml.tmpl",
                     package = "castSDM")
  if (!nzchar(src))
    cli::cli_abort("Could not find bundled config template (re-install package).")
  if (file.exists(destination) && !overwrite)
    cli::cli_abort(c(
      "{.path {destination}} already exists.",
      i = "Pass {.code overwrite = TRUE} to replace it."
    ))
  ok <- file.copy(src, destination, overwrite = overwrite)
  if (!ok)
    cli::cli_abort("Failed to copy template to {.path {destination}}.")
  cli::cli_inform("Wrote template to {.path {destination}}.")
  invisible(destination)
}
