#' Generate Spatial Habitat Suitability Predictions
#'
#' Predicts habitat suitability scores (HSS) for new environmental data using
#' fitted models. Each model produces a probability representing the
#' predicted suitability of each site.
#'
#' HSS is a **relative habitat suitability** score whose ranking is meaningful,
#' not a calibrated probability of occurrence: its absolute level depends on the
#' arbitrary presence:background prevalence set in [cast_background()]. Compare
#' and threshold HSS within a projection, and avoid reading it as an absolute
#' occurrence probability.
#'
#' When the fit carries a training reference (it does by default), each
#' prediction row is scored for **extrapolation** via the multivariate
#' environmental similarity surface (MESS; Elith et al. 2010). Negative MESS
#' marks sites outside the training envelope of at least one predictor, where
#' model output is an extrapolation and should be treated cautiously. Optional
#' `clamp`ing caps each predictor to its training range before prediction,
#' which curbs runaway extrapolation without hiding it (the MESS flag is still
#' computed on the unclamped input).
#'
#' @param fit A [cast_fit] object.
#' @param new_data A `data.frame` with `lon`, `lat`, and the same
#'   environmental variables used in fitting.
#' @param models Character vector. Which models to predict with. Default
#'   `NULL` (all fitted models).
#' @param clamp Logical. Clamp predictors to the training range before
#'   prediction. Default `FALSE`.
#' @param extrapolation Logical. Append `mess` and `extrapolating` columns
#'   flagging out-of-envelope sites. Default `TRUE` (skipped if the fit lacks
#'   a stored reference).
#'
#' @return A `cast_predict` object containing a `predictions` data.frame
#'   with `lon`, `lat`, one `HSS_*` column per model, and (when enabled)
#'   `mess` / `extrapolating` columns.
#'
#' @references
#' Elith, J., Kearney, M. & Phillips, S. (2010). The art of modelling
#' range-shifting species. *Methods in Ecology and Evolution*, 1(4), 330-342.
#'
#' @seealso [cast_fit()], [cast_evaluate()], [cast_ensemble()]
#'
#' @export
cast_predict <- function(fit, new_data, models = NULL,
                         clamp = FALSE, extrapolation = TRUE) {
  env_vars <- fit$env_vars
  missing_vars <- setdiff(env_vars, names(new_data))
  if (length(missing_vars)) {
    cli::cli_abort(
      "{.arg new_data} is missing fitted predictor{?s}: {.val {missing_vars}}."
    )
  }
  mdl_names <- models %||% names(fit$models)
  mdl_names <- intersect(mdl_names, names(fit$models))

  if (length(mdl_names) == 0) {
    cli::cli_abort("No matching fitted models found.")
  }

  # Prepare new data
  X_raw <- as.data.frame(new_data[, env_vars, drop = FALSE], check.names = FALSE)
  for (col in names(X_raw)) X_raw[[col]] <- as.numeric(X_raw[[col]])
  X_raw <- .cast_impute(X_raw, fit$scaling$impute)

  reference <- fit$scaling$reference
  # Extrapolation diagnostics on the (imputed, unclamped) input.
  mess_vals <- NULL
  if (isTRUE(extrapolation) && !is.null(reference)) {
    mess_vals <- tryCatch(.cast_mess(reference, X_raw),
                          error = function(e) NULL)
  }
  if (isTRUE(clamp) && !is.null(reference)) {
    X_raw <- .cast_clamp(X_raw, reference)
  }

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
      predict_single_model(mdl_info, X_raw),
      error = function(e) {
        cli::cli_warn("Prediction failed for {.val {mdl_name}}: {e$message}")
        rep(NA_real_, nrow(new_data))
      }
    )
  }

  if (!is.null(mess_vals) && length(mess_vals) == nrow(pred_df)) {
    pred_df$mess <- mess_vals
    pred_df$extrapolating <- mess_vals < 0
  }

  new_cast_predict(
    predictions = pred_df,
    models = mdl_names
  )
}

#' Clamp predictors to the training range
#' @keywords internal
#' @noRd
.cast_clamp <- function(X, reference) {
  for (col in intersect(names(X), names(reference))) {
    rng <- range(reference[[col]], na.rm = TRUE)
    if (all(is.finite(rng))) {
      X[[col]] <- pmin(pmax(X[[col]], rng[1]), rng[2])
    }
  }
  X
}

#' Multivariate Environmental Similarity Surface (MESS)
#'
#' Row-wise MESS of `newdata` relative to `reference` (Elith et al. 2010).
#' Each predictor's similarity is the interpolation percentile within the
#' reference distribution (negative outside the reference min/max); the point
#' MESS is the minimum across predictors. Negative values flag extrapolation.
#'
#' @keywords internal
#' @noRd
.cast_mess <- function(reference, newdata) {
  vars <- intersect(names(reference), names(newdata))
  if (!length(vars)) return(rep(NA_real_, nrow(newdata)))
  sim <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(vars))
  for (j in seq_along(vars)) {
    v <- vars[j]
    ref <- reference[[v]][is.finite(reference[[v]])]
    p <- as.numeric(newdata[[v]])
    if (!length(ref)) next
    mn <- min(ref); mx <- max(ref); rng <- mx - mn
    n <- length(ref)
    f <- vapply(p, function(pi) {
      if (!is.finite(pi)) return(NA_real_)
      100 * sum(ref < pi) / n
    }, numeric(1))
    s <- numeric(length(p))
    for (i in seq_along(p)) {
      pi <- p[i]; fi <- f[i]
      if (!is.finite(pi)) { s[i] <- NA_real_; next }
      if (rng <= 0) {
        s[i] <- if (pi == mn) 100 else -100
      } else if (fi == 0) {
        s[i] <- (pi - mn) / rng * 100
      } else if (fi <= 50) {
        s[i] <- 2 * fi
      } else if (fi < 100) {
        s[i] <- 2 * (100 - fi)
      } else {
        s[i] <- (mx - pi) / rng * 100
      }
    }
    sim[, j] <- s
  }
  apply(sim, 1L, function(r) if (all(is.na(r))) NA_real_ else min(r, na.rm = TRUE))
}
