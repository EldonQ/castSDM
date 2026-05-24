#' Ensemble of Small Models (ESM) for Rare Species
#'
#' Implements the Ensemble of Small Models approach of Breiner et al. (2015)
#' for species with low presence counts. Fits one bivariate sub-model
#' (GLM or GAM) for each pair of selected predictors and combines them via
#' AUC-weighted averaging on a held-out validation split.
#'
#' Designed as the recommended fallback when BRT would over-fit
#' on rare species (typically `n_pres < presence_threshold`). Complementary
#' to [cast_fit()] — see [cast_batch()]'s `auto_esm` argument for automatic
#' selection.
#'
#' @param data A `data.frame` with `presence` column and predictor columns.
#' @param vars Character vector of candidate predictor names. If `NULL`,
#'   uses [get_env_vars()].
#' @param top_k Integer. Maximum number of predictors retained before
#'   enumerating pairs (defends against combinatorial explosion). Default
#'   `8L`.
#' @param base_algo Character. `"glm"` (default) or `"gam"`. GAMs use
#'   thin-plate smooths; GLMs use linear + quadratic terms.
#' @param val_fraction Numeric. Stratified hold-out fraction for AUC weight
#'   estimation. Default `0.25`.
#' @param response Character. Response column. Default `"presence"`.
#' @param seed Integer. Random seed.
#' @param verbose Logical. Print progress.
#'
#' @return A list with class `"cast_esm"` containing:
#' \describe{
#'   \item{models}{List of fitted bivariate sub-models.}
#'   \item{pairs}{Matrix of column-name pairs used.}
#'   \item{weights}{AUC-derived weights on the held-out split.}
#'   \item{vars}{The retained predictors.}
#'   \item{base_algo}{`"glm"` or `"gam"`.}
#'   \item{val_auc}{Ensemble AUC on the held-out split.}
#' }
#'
#' @references
#' Breiner, F. T., Guisan, A., Bergamini, A., & Nobis, M. P. (2015).
#'   Overcoming limitations of modelling rare species by using ensembles of
#'   small models. *Methods in Ecology and Evolution*, 6, 1210-1218.
#'
#' @seealso [cast_fit()], [cast_batch()]
#' @export
cast_esm <- function(data,
                     vars         = NULL,
                     top_k        = 8L,
                     base_algo    = c("glm", "gam"),
                     val_fraction = 0.25,
                     response     = "presence",
                     seed         = NULL,
                     verbose      = TRUE) {
  base_algo <- match.arg(base_algo)
  if (base_algo == "gam") check_suggested("mgcv", "for ESM with GAM")

  Y <- data[[response]]
  if (is.null(vars)) vars <- get_env_vars(data, response)
  vars <- intersect(vars, names(data))
  if (length(vars) < 2L)
    cli::cli_abort("ESM requires at least 2 predictors; got {length(vars)}.")

  # Defend against combinatorial explosion — keep top_k by univariate AUC.
  if (length(vars) > top_k) {
    if (verbose) cli::cli_inform(
      "ESM: ranking {length(vars)} predictors and keeping top {top_k} by univariate AUC."
    )
    uauc <- vapply(vars, function(v) {
      x <- data[[v]]
      tryCatch(compute_auc(Y, x), error = function(e) NA_real_)
    }, numeric(1))
    # AUC and 1 - AUC are equivalent for ranking; take pmax to handle direction
    uauc <- pmax(uauc, 1 - uauc, na.rm = TRUE)
    keep <- order(uauc, decreasing = TRUE)[seq_len(top_k)]
    vars <- vars[keep]
  }

  pairs <- utils::combn(vars, 2L)  # 2 x P matrix

  # Stratified val split
  if (!is.null(seed)) set.seed(seed)
  pos_idx <- which(Y == 1L); neg_idx <- which(Y == 0L)
  val_pos <- sample(pos_idx, max(1L, round(val_fraction * length(pos_idx))))
  val_neg <- sample(neg_idx, max(1L, round(val_fraction * length(neg_idx))))
  val_idx <- c(val_pos, val_neg)
  tr_idx  <- setdiff(seq_len(nrow(data)), val_idx)
  tr_df   <- data[tr_idx, ]
  va_df   <- data[val_idx, ]
  y_va    <- Y[val_idx]

  fit_one <- function(v1, v2) {
    df_tr <- data.frame(
      presence = tr_df[[response]],
      x1 = tr_df[[v1]], x2 = tr_df[[v2]]
    )
    if (base_algo == "glm") {
      stats::glm(presence ~ x1 + I(x1^2) + x2 + I(x2^2) + x1:x2,
                 data = df_tr, family = stats::binomial())
    } else {
      mgcv::gam(presence ~ s(x1, k = 4) + s(x2, k = 4) + ti(x1, x2, k = c(3, 3)),
                data = df_tr, family = stats::binomial(), method = "REML")
    }
  }

  predict_one <- function(mdl, v1, v2, df) {
    nd <- data.frame(x1 = df[[v1]], x2 = df[[v2]])
    as.numeric(stats::predict(mdl, newdata = nd, type = "response"))
  }

  P <- ncol(pairs)
  models <- vector("list", P)
  weights <- numeric(P)

  for (j in seq_len(P)) {
    v1 <- pairs[1L, j]; v2 <- pairs[2L, j]
    mdl <- tryCatch(fit_one(v1, v2),
                    error = function(e) NULL)
    if (is.null(mdl)) {
      models[[j]] <- NULL; weights[j] <- 0; next
    }
    pv <- tryCatch(predict_one(mdl, v1, v2, va_df),
                   error = function(e) rep(NA_real_, nrow(va_df)))
    auc <- tryCatch(compute_auc(y_va, pv), error = function(e) 0.5)
    if (!is.finite(auc)) auc <- 0.5
    # Convert to weight: scale (AUC - 0.5) and clip below 0.
    w <- pmax(0, auc - 0.5)
    models[[j]] <- mdl
    weights[j] <- w
  }
  if (sum(weights) <= 0) weights[] <- 1 / length(weights)
  weights <- weights / sum(weights)

  # Validate ensemble
  ens_va <- rep(0, nrow(va_df))
  for (j in seq_len(P)) {
    if (is.null(models[[j]])) next
    pv <- tryCatch(predict_one(models[[j]], pairs[1L, j], pairs[2L, j], va_df),
                   error = function(e) rep(NA_real_, nrow(va_df)))
    pv[is.na(pv)] <- mean(pv, na.rm = TRUE)
    ens_va <- ens_va + weights[j] * pv
  }
  val_auc <- compute_auc(y_va, ens_va)

  if (verbose) cli::cli_inform(
    "ESM: fitted {sum(!vapply(models, is.null, logical(1)))} pairs from {length(vars)} predictors; val_AUC={signif(val_auc, 4)}."
  )

  out <- list(
    models    = models,
    pairs     = pairs,
    weights   = weights,
    vars      = vars,
    base_algo = base_algo,
    val_auc   = val_auc
  )
  class(out) <- "cast_esm"
  out
}


#' Predict from a Fitted Cast ESM
#'
#' @param object A `cast_esm` object.
#' @param newdata A `data.frame` containing all variables in `object$vars`.
#' @param ... Unused.
#' @return Numeric vector of probabilities in `[0, 1]`.
#' @keywords internal
#' @noRd
predict_cast_esm <- function(object, newdata, ...) {
  if (!inherits(object, "cast_esm"))
    cli::cli_abort("{.arg object} must be a cast_esm.")
  pairs <- object$pairs
  out <- rep(0, nrow(newdata))
  for (j in seq_len(ncol(pairs))) {
    if (is.null(object$models[[j]])) next
    v1 <- pairs[1L, j]; v2 <- pairs[2L, j]
    if (!(v1 %in% names(newdata)) || !(v2 %in% names(newdata))) next
    nd <- data.frame(x1 = newdata[[v1]], x2 = newdata[[v2]])
    pv <- tryCatch(
      as.numeric(stats::predict(object$models[[j]], newdata = nd,
                                type = "response")),
      error = function(e) rep(NA_real_, nrow(newdata))
    )
    pv[is.na(pv)] <- mean(pv, na.rm = TRUE)
    if (any(is.nan(pv))) pv[is.nan(pv)] <- 0.5
    out <- out + object$weights[j] * pv
  }
  pmin(pmax(out, 0), 1)
}


#' Print Method for Cast ESM
#' @param x A `cast_esm` object.
#' @param ... Ignored.
#' @export
print.cast_esm <- function(x, ...) {
  cat("<cast_esm>\n")
  cat("  base algorithm:    ", x$base_algo, "\n", sep = "")
  cat("  predictors kept:    ", length(x$vars), "\n", sep = "")
  cat("  bivariate sub-models:", ncol(x$pairs), "\n", sep = "")
  cat("  hold-out val AUC:   ",
      if (is.null(x$val_auc) || !is.finite(x$val_auc)) "NA"
      else format(round(x$val_auc, 4)), "\n", sep = "")
  invisible(x)
}
