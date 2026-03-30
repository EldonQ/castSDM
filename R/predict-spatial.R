#' Generate Spatial Habitat Suitability Predictions
#'
#' Predicts habitat suitability scores (HSS) for new environmental data using
#' fitted models. Each model produces a probability representing the
#' predicted suitability of each site.
#'
#' @param fit A [cast_fit] object.
#' @param new_data A `data.frame` with `lon`, `lat`, and the same
#'   environmental variables used in fitting.
#' @param models Character vector. Which models to predict with. Default
#'   `NULL` (all fitted models).
#'
#' @return A `cast_predict` object containing a `predictions` data.frame
#'   with `lon`, `lat`, and one `HSS_*` column per model.
#'
#' @seealso [cast_fit()], [cast_evaluate()], [cast_cate()]
#'
#' @export
cast_predict <- function(fit, new_data, models = NULL) {
  env_vars <- fit$env_vars
  scaling <- fit$scaling
  mdl_names <- models %||% names(fit$models)
  mdl_names <- intersect(mdl_names, names(fit$models))

  if (length(mdl_names) == 0) {
    cli::cli_abort("No matching fitted models found.")
  }

  # Prepare new data
  X_raw <- as.data.frame(new_data[, env_vars, drop = FALSE])
  X_raw[is.na(X_raw)] <- 0
  X_sc <- as.data.frame(
    scale(X_raw, center = scaling$means, scale = scaling$sds)
  )
  X_sc[is.na(X_sc)] <- 0

  # Extract coordinates if available
  has_coords <- all(c("lon", "lat") %in% names(new_data))
  pred_df <- if (has_coords) {
    data.frame(lon = new_data$lon, lat = new_data$lat)
  } else {
    data.frame(site = seq_len(nrow(new_data)))
  }

  for (mdl_name in mdl_names) {
    mdl_info <- fit$models[[mdl_name]]
    col_name <- paste0("HSS_", mdl_name)

    pred_df[[col_name]] <- tryCatch(
      predict_single_model(
        mdl_info, X_raw, X_sc,
        fit$screen, fit$dag, fit$ate
      ),
      error = function(e) {
        cli::cli_warn("Prediction failed for {.val {mdl_name}}: {e$message}")
        rep(NA_real_, nrow(new_data))
      }
    )
  }

  new_cast_predict(
    predictions = pred_df,
    models = mdl_names
  )
}
