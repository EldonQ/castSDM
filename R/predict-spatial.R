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
#' Predictors must be numeric: non-numeric columns (e.g. factor or
#' character) in `new_data` abort with an error instead of being silently
#' coerced via `as.numeric()`.
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
#' On top of the univariate MESS envelope check, `aoa = TRUE` scores each
#' prediction row against the **area of applicability** (AOA; Meyer & Pebesma
#' 2021): the dissimilarity index (DI) is the distance from the site to its
#' nearest training neighbour in the standardised predictor space, compared
#' against a threshold calibrated by [cast_cv()] with `aoa = TRUE` as the 95th
#' percentile of the held-out-fold DI values. Sites with `aoa = FALSE` (DI
#' above the threshold) are more dissimilar from training conditions than any
#' held-out training fold was from its own training data; model quality there
#' is not covered by the cross-validation estimate and predictions should not
#' be interpolated into those sites.
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
#' @param cv Optional `cast_cv` object. Required when `aoa = TRUE`; it must
#'   carry an AOA calibration (i.e. come from [cast_cv()] with `aoa = TRUE`).
#' @param aoa Logical. Append `aoa_di` (dissimilarity index) and `aoa`
#'   (logical, `DI <= threshold`) columns. Default `FALSE`.
#' @param vi Optional named numeric vector of variable importance (any non-negative scale, normalised internally)
#'   (any positive scale; normalised internally), used to weight the
#'   standardized predictor dimensions in the DI (Meyer & Pebesma 2021,
#'   weighted variant). Default `NULL` (uniform weights). A warning is issued
#'   when `vi` is supplied but the `cv` threshold was calibrated with uniform
#'   weights, because threshold and DI then use different metrics.
#'
#' @return A `cast_predict` object containing a `predictions` data.frame
#'   with `lon`, `lat`, one `HSS_*` column per model, and (when enabled)
#'   `mess` / `extrapolating` and `aoa_di` / `aoa` columns.
#'
#' @references
#' Elith, J., Kearney, M. & Phillips, S. (2010). The art of modelling
#' range-shifting species. *Methods in Ecology and Evolution*, 1(4), 330-342.
#'
#' Meyer, H. & Pebesma, E. (2021). Predicting into unknown space? Estimating
#' the area of applicability of spatial prediction models. *Methods in
#' Ecology and Evolution*, 12(9), 1620-1633.
#'
#' @seealso [cast_fit()], [cast_evaluate()], [cast_ensemble()], [cast_cv()]
#'
#' @export
cast_predict <- function(fit, new_data, models = NULL,
                         clamp = FALSE, extrapolation = TRUE,
                         cv = NULL, aoa = FALSE, vi = NULL) {
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
  .cast_check_numeric_predictors(X_raw)
  for (col in names(X_raw)) X_raw[[col]] <- as.numeric(X_raw[[col]])
  X_raw <- .cast_impute(X_raw, fit$scaling$impute)

  reference <- fit$scaling$reference
  # Extrapolation diagnostics on the (imputed, unclamped) input.
  mess_vals <- NULL
  if (isTRUE(extrapolation) && !is.null(reference)) {
    mess_vals <- tryCatch(.cast_mess(reference, X_raw),
                          error = function(e) NULL)
  }

  # Area of applicability on the same (imputed, unclamped) input: clamping
  # must never soften the DI, exactly as it never hides MESS extrapolation.
  aoa_di_vals <- NULL
  if (isTRUE(aoa)) {
    if (is.null(cv) || !isTRUE(cv$aoa$enabled) ||
        !is.finite(cv$aoa$threshold %||% NA_real_)) {
      cli::cli_abort(c(
        "{.code aoa = TRUE} requires {.arg cv} from {.code cast_cv(aoa = TRUE)}.",
        "i" = "The supplied {.arg cv} carries no AOA calibration."
      ))
    }
    if (is.null(reference)) {
      cli::cli_abort("The fit carries no stored training reference; AOA cannot be computed.")
    }
    if (!is.null(vi) && !isTRUE(cv$aoa$weighted)) {
      cli::cli_warn(c(
        "{.arg vi} weights the DI but the {.arg cv} AOA threshold was calibrated with uniform weights.",
        "i" = "Re-run {.code cast_cv(aoa = TRUE)} with the same importance source, or pass {.code vi = NULL}."
      ))
    }
    ref_m <- as.data.frame(reference, check.names = FALSE)
    ref_m <- ref_m[, intersect(fit$env_vars, names(ref_m)), drop = FALSE]
    prep <- tryCatch(
      .cast_aoa_prep(ref_m, X_raw[, fit$env_vars, drop = FALSE], vi = vi),
      error = function(e) e
    )
    if (inherits(prep, "error")) {
      cli::cli_warn("AOA computation failed ({prep$message}); aoa columns are skipped.")
    } else {
      aoa_di_vals <- .cast_aoa_di(prep$train, prep$new)
    }
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
  if (!is.null(aoa_di_vals) && length(aoa_di_vals) == nrow(pred_df)) {
    pred_df$aoa_di <- aoa_di_vals
    pred_df$aoa <- aoa_di_vals <= cv$aoa$threshold
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
#' Vectorised implementation: per-predictor reference quantiles come from a
#' single `findInterval()` against the sorted reference instead of a
#' per-cell `sum(ref < pi)` loop. Branch-for-branch identical to Elith et
#' al. (2010), including the degenerate `range <= 0` case and the
#' `NA`-in-`newdata` -> `NA` semantics.
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
    # f = 100 * (number of reference values strictly below p) / n.
    # left.open = TRUE makes findInterval() count ref < p (not ref <= p).
    f <- 100 * findInterval(p, sort(ref), left.open = TRUE) / n
    ok <- is.finite(p)
    s <- rep(NA_real_, length(p))
    if (any(ok)) {
      po <- p[ok]; fo <- f[ok]
      sj <- numeric(length(po))
      if (rng <= 0) {
        sj[] <- ifelse(po == mn, 100, -100)
      } else {
        below  <- fo == 0
        low    <- fo > 0 & fo <= 50
        high   <- fo > 50 & fo < 100
        above  <- fo >= 100
        sj[below] <- (po[below] - mn) / rng * 100
        sj[low]   <- 2 * fo[low]
        sj[high]  <- 2 * (100 - fo[high])
        sj[above] <- (mx - po[above]) / rng * 100
      }
      s[ok] <- sj
    }
    sim[, j] <- s
  }
  # Row-wise minimum over predictors; rows with no finite similarity -> NA.
  sim[!is.finite(sim)] <- Inf
  mess <- do.call(pmin, lapply(seq_len(ncol(sim)), function(j) sim[, j]))
  mess[!is.finite(mess)] <- NA_real_
  mess
}

#' Standardised, Optionally Importance-Weighted AOA Space
#'
#' Builds the predictor space in which the AOA dissimilarity index (DI) is
#' measured (Meyer & Pebesma 2021): columns are standardized on the training
#' means/SDs (training gaps imputed with column medians first), then each
#' column is multiplied by `sqrt(vi / sum(vi))` so Euclidean distances become
#' importance-weighted. Uniform weights reduce to plain standardization.
#'
#' `.cast_aoa_params()` estimates and applies the transform to the training
#' rows once; `.cast_aoa_apply()` projects new rows with the stored
#' parameters (used block by block on rasters); `.cast_aoa_prep()` is the
#' one-shot convenience wrapper.
#'
#' @param X_train,X_new Data frames or matrices with identical column names.
#' @param vi Optional named numeric vector of variable importance; matched
#'   to the columns of `X_train`.
#' @return `.cast_aoa_params()` returns the transform parameters with the
#'   standardized training matrix; `.cast_aoa_apply()` a standardized
#'   matrix; `.cast_aoa_prep()` a list with `train`, `new` and `weighted`.
#' @keywords internal
#' @noRd
.cast_aoa_params <- function(X_train, vi = NULL) {
  train_m <- as.matrix(X_train)
  storage.mode(train_m) <- "numeric"
  cn <- colnames(train_m)
  if (is.null(cn) || any(!nzchar(cn))) {
    cli::cli_abort("AOA training predictors must be a matrix with named columns.")
  }
  if (nrow(train_m) < 2L) {
    cli::cli_abort("AOA needs at least two training rows.")
  }
  w <- NULL
  if (!is.null(vi)) {
    if (!is.numeric(vi) || is.null(names(vi)) ||
        anyNA(names(vi)) || anyDuplicated(names(vi))) {
      cli::cli_abort("{.arg vi} must be a uniquely named numeric vector.")
    }
    if (any(vi < 0, na.rm = TRUE) || sum(vi, na.rm = TRUE) <= 0) {
      cli::cli_abort("{.arg vi} must be non-negative with a positive sum.")
    }
    w <- vi[match(cn, names(vi))]
    if (anyNA(w)) {
      cli::cli_abort("{.arg vi} is missing importance for: {.val {cn[is.na(w)]}}.")
    }
  }

  # Median-impute the training columns; the same medians are stored so new
  # rows (and raster blocks) are projected with identical parameters.
  med <- vapply(seq_len(ncol(train_m)), function(j) {
    stats::median(train_m[, j], na.rm = TRUE)
  }, numeric(1))
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(train_m))) {
    train_m[!is.finite(train_m[, j]), j] <- med[j]
  }

  cm <- colMeans(train_m)
  cs <- vapply(seq_len(ncol(train_m)), function(j) stats::sd(train_m[, j]),
               numeric(1))
  cs[!is.finite(cs) | cs <= 0] <- 1
  train_s <- sweep(sweep(train_m, 2L, cm), 2L, cs, "/")

  weighted <- FALSE
  if (!is.null(w)) {
    w <- w / sum(w)
    train_s <- sweep(train_s, 2L, sqrt(w), "*")
    weighted <- TRUE
  }
  list(cn = cn, med = med, cm = cm, cs = cs, w = w, weighted = weighted,
       train = train_s)
}

#' @keywords internal
#' @noRd
.cast_aoa_apply <- function(params, X_new) {
  new_m <- as.matrix(X_new)
  storage.mode(new_m) <- "numeric"
  if (is.null(colnames(new_m)) || !identical(colnames(new_m), params$cn)) {
    cli::cli_abort("AOA predictor matrices must share identically named columns.")
  }
  for (j in seq_len(ncol(new_m))) {
    new_m[!is.finite(new_m[, j]), j] <- params$med[j]
  }
  new_s <- sweep(sweep(new_m, 2L, params$cm), 2L, params$cs, "/")
  if (!is.null(params$w)) {
    new_s <- sweep(new_s, 2L, sqrt(params$w), "*")
  }
  new_s
}

#' @keywords internal
#' @noRd
.cast_aoa_prep <- function(X_train, X_new, vi = NULL) {
  p <- .cast_aoa_params(X_train, vi = vi)
  list(train = p$train, new = .cast_aoa_apply(p, X_new),
       weighted = p$weighted)
}

#' AOA Dissimilarity Index (DI)
#'
#' Distance from every query row to its nearest neighbour in the reference
#' (training) rows of the standardized AOA space, via \pkg{FNN}.
#'
#' @keywords internal
#' @noRd
.cast_aoa_di <- function(train_s, new_s) {
  if (!nrow(new_s)) return(numeric(0))
  check_suggested("FNN", "for the area of applicability (aoa = TRUE)")
  FNN::get.knnx(data = train_s, query = new_s, k = 1L)$nn.dist[, 1L]
}

#' Within-Training Nearest-Other-Neighbour Distances
#'
#' For every training row, the distance to its nearest *other* training row
#' (k = 2 skips the self-match at distance 0). The pooled upper quantile of
#' these distances is the AOA threshold in Meyer & Pebesma (2021); castSDM
#' calibrates it from held-out CV folds instead (see [cast_cv()] with `aoa = TRUE`),
#' so this helper is exported for diagnostics only.
#'
#' @keywords internal
#' @noRd
.cast_aoa_di_train <- function(train_s) {
  if (nrow(train_s) < 2L) return(rep(NA_real_, nrow(train_s)))
  check_suggested("FNN", "for the area of applicability (aoa = TRUE)")
  FNN::get.knn(train_s, k = 2L)$nn.dist[, 2L]
}
