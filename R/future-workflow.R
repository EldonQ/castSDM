#' Load Future Environmental Scenarios
#'
#' Reads future environmental data from an RDS file, a CSV file, or a directory
#' of CSV files. The return value is always a named list of data.frames suitable
#' for [cast_project()].
#'
#' @param path Path to an RDS file, CSV file, or directory of CSV files.
#' @return Named list of data.frames.
#' @export
cast_load_future_envs <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    cli::cli_abort("{.arg path} must point to an RDS file, CSV file, or CSV directory.")
  }
  if (!file.exists(path)) {
    cli::cli_abort("Future environment path not found: {.path {path}}.")
  }

  if (dir.exists(path)) {
    files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
    if (!length(files)) {
      cli::cli_abort("No CSV files found in future environment directory: {.path {path}}.")
    }
    out <- lapply(files, utils::read.csv, check.names = FALSE)
    names(out) <- tools::file_path_sans_ext(basename(files))
  } else if (tolower(tools::file_ext(path)) == "rds") {
    obj <- readRDS(path)
    if (is.data.frame(obj)) {
      out <- list(future = obj)
    } else if (is.list(obj) && all(vapply(obj, is.data.frame, logical(1)))) {
      out <- obj
    } else {
      cli::cli_abort("RDS future environment must be a data.frame or named list of data.frames.")
    }
  } else if (tolower(tools::file_ext(path)) == "rda") {
    e <- new.env(parent = emptyenv())
    nm <- load(path, envir = e)
    obj <- if (length(nm) == 1L) get(nm, envir = e) else as.list(e)
    if (is.data.frame(obj)) {
      out <- list(future = obj)
    } else if (is.list(obj) && all(vapply(obj, is.data.frame, logical(1)))) {
      out <- obj
    } else {
      cli::cli_abort("RDA future environment must contain a data.frame or named list of data.frames.")
    }
  } else if (tolower(tools::file_ext(path)) == "csv") {
    out <- list(future = utils::read.csv(path, check.names = FALSE))
    names(out) <- tools::file_path_sans_ext(basename(path))
  } else {
    cli::cli_abort("Unsupported future environment file extension: {.val {tools::file_ext(path)}}.")
  }

  if (is.null(names(out)) || any(!nzchar(names(out)))) {
    names(out) <- paste0("scenario_", seq_along(out))
  }
  out
}


#' Run Future Projection and Save Standard Outputs
#'
#' Runs [cast_project()], saves CSV/GeoTIFF outputs via `cast_project()`, and
#' saves current, future, and change-map figures for every scenario.
#'
#' @param fit A `cast_fit` object.
#' @param cv A `cast_cv` object.
#' @param current_env Current prediction grid.
#' @param future_envs Named list of future prediction grids.
#' @param save_dir Output directory.
#' @param basemap Basemap passed to plotting methods.
#' @param fig_dpi Figure DPI.
#' @param method Ensemble method.
#' @param threshold_method Threshold method.
#' @param models Optional models used by [cast_project()].
#' @param prefix Filename prefix.
#' @return A `cast_project` object.
#' @export
cast_save_future_projection <- function(fit, cv, current_env, future_envs,
                                        save_dir,
                                        basemap = "world",
                                        fig_dpi = 300L,
                                        method = "weighted",
                                        threshold_method = "maxTSS",
                                        models = NULL,
                                        prefix = "future") {
  if (is.null(cv)) {
    cli::cli_abort("Future projection requires a non-null {.cls cast_cv} object for ensemble weights.")
  }
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  proj <- cast_project(
    fit = fit,
    cv = cv,
    current_env = current_env,
    future_envs = future_envs,
    method = method,
    threshold_method = threshold_method,
    models = models,
    save_dir = save_dir
  )

  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("sf", quietly = TRUE)) {
    save_plot <- if (exists("cast_safe_ggsave", mode = "function")) {
      cast_safe_ggsave
    } else {
      ggplot2::ggsave
    }

    p_cur <- tryCatch(plot(proj$current, basemap = basemap), error = function(e) NULL)
    if (!is.null(p_cur)) {
      save_plot(
        file.path(save_dir, paste0(prefix, "_current_ensemble.png")),
        p_cur, width = 10, height = 7, dpi = fig_dpi,
        bg = "white", limitsize = FALSE
      )
    }

    for (scen in names(proj$future)) {
      p_fut <- tryCatch(plot(proj$future[[scen]], basemap = basemap), error = function(e) NULL)
      if (!is.null(p_fut)) {
        save_plot(
          file.path(save_dir, paste0(prefix, "_", scen, "_ensemble.png")),
          p_fut, width = 10, height = 7, dpi = fig_dpi,
          bg = "white", limitsize = FALSE
        )
      }

      p_chg <- tryCatch(plot(proj, scenario = scen, basemap = basemap), error = function(e) NULL)
      if (!is.null(p_chg)) {
        save_plot(
          file.path(save_dir, paste0(prefix, "_", scen, "_change.png")),
          p_chg, width = 10, height = 7, dpi = fig_dpi,
          bg = "white", limitsize = FALSE
        )
      }
    }
  }

  saveRDS(proj, file.path(save_dir, paste0(prefix, "_projection.rds")))
  proj
}
