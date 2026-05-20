#' Construct a Multi-layer Worker Budget for `cast_batch()`
#'
#' Distributes a fixed total of parallel workers across three nested layers
#' (species -> intra-species step -> within-step folds/bootstraps) so that
#' the available cores on a Windows workstation are filled without
#' oversubscription. Used internally by [cast_batch()] and exposed for users
#' who want to inspect or override the allocation.
#'
#' Why two layers and not just one? When the batch contains many species,
#' the species-level fan-out alone saturates the CPU. When the batch is
#' small (1–3 species), nesting the spare cores into per-species DAG
#' bootstraps / ATE folds / CV folds avoids idle cores.
#'
#' @param total_workers Integer. Total parallel worker budget. Default
#'   `parallel::detectCores() - 1L`, capped at the available cores.
#' @param n_species Integer. Number of species in the batch.
#'
#' @return A list of class `cast_worker_budget` with components:
#' \describe{
#'   \item{total}{Effective total workers (after capping).}
#'   \item{species}{Number of species processed in parallel.}
#'   \item{intra}{Number of workers allocated to per-species steps.}
#' }
#'
#' @examples
#' cast_worker_budget(total_workers = 8L, n_species = 4L)
#'
#' @export
cast_worker_budget <- function(total_workers = NULL,
                               n_species     = 1L) {
  avail <- max(1L, parallel::detectCores() - 1L)
  total <- if (is.null(total_workers)) avail else
    max(1L, min(as.integer(total_workers), avail))
  n_species <- max(1L, as.integer(n_species))

  species_w <- min(total, n_species)
  intra_w   <- max(1L, total %/% max(1L, species_w))

  out <- list(total = total, species = species_w, intra = intra_w)
  class(out) <- "cast_worker_budget"
  out
}


#' @export
print.cast_worker_budget <- function(x, ...) {
  cat("<cast_worker_budget>\n")
  cat("  total workers     :", x$total, "\n")
  cat("  species in parallel:", x$species, "\n")
  cat("  intra-species cores:", x$intra,
      sprintf(" (= %d / %d)", x$total, x$species), "\n")
  invisible(x)
}


#' Set Up a `future` Plan for the Species Layer of a `cast_batch()` Run
#'
#' Convenience wrapper around `future::plan(multisession)` that is safe on
#' Windows (uses PSOCK clusters) and idempotent. Falls back silently to
#' `sequential` when `future` is not installed.
#'
#' @param budget A `cast_worker_budget` returned by [cast_worker_budget()].
#' @return Previous `future` plan (returned by `future::plan`), so callers
#'   can restore it via `on.exit()`.
#' @keywords internal
#' @noRd
cast_setup_species_plan <- function(budget) {
  if (!requireNamespace("future", quietly = TRUE)) return(NULL)
  if (budget$species <= 1L) {
    return(future::plan(future::sequential))
  }
  future::plan(future::multisession, workers = budget$species)
}
