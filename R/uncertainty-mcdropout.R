#' Prediction Uncertainty via MC Dropout
#'
#' Estimates prediction uncertainty for the CI-MLP model using Monte Carlo
#' Dropout (Gal & Ghahramani, 2016). Runs multiple stochastic forward passes
#' with dropout active, producing a distribution of predictions from which
#' mean, standard deviation, and coefficient of variation are computed per site.
#'
#' @param fit A [cast_fit] object that contains at least one fitted neural
#'   network model (`"cast"`, `"mlp_ate"`, or `"mlp"`). If multiple
#'   neural models are present, `"cast"` is used preferentially.
#' @param new_data A `data.frame` with `lon`, `lat`, and environmental
#'   variables matching those in `fit`.
#' @param n_forward Integer. Number of stochastic forward passes.
#'   Default `50`.
#' @param verbose Logical. Print progress. Default `TRUE`.
#'
#' @return A `data.frame` with columns:
#'   \describe{
#'     \item{lon, lat}{Coordinates from `new_data`.}
#'     \item{mean}{Mean predicted probability across MC samples.}
#'     \item{sd}{Standard deviation of predictions (epistemic uncertainty).}
#'     \item{cv}{Coefficient of variation (`sd / max(mean, 1e-6)`).}
#'     \item{q05, q95}{5th and 95th percentile of MC predictions.}
#'   }
#'
#' @details
#' Standard neural network predictions use a single deterministic forward pass
#' (dropout off). MC Dropout instead keeps dropout **on** during prediction,
#' producing different outputs each pass. The variance across passes is an
#' approximation to Bayesian model uncertainty, capturing regions where the
#' model is less confident.
#'
#' High uncertainty areas typically correspond to:
#' - Extrapolation beyond training data range
#' - Environmental conditions underrepresented in training
#' - Transition zones between suitable and unsuitable habitat
#'
#' @references
#' Gal, Y. & Ghahramani, Z. (2016). Dropout as a Bayesian Approximation:
#' Representing Model Uncertainty in Deep Learning. *ICML*.
#'
#' @seealso [cast_predict()], [cast_fit()]
#'
#' @export
cast_uncertainty <- function(fit, new_data, n_forward = 50L,
                             verbose = TRUE) {
  if (!inherits(fit, "cast_fit")) {
    cli::cli_abort("{.arg fit} must be a {.cls cast_fit} object.")
  }
  nn_priority <- c("cast", "mlp_ate", "mlp", "ci_mlp")
  nn_available <- nn_priority[nn_priority %in% names(fit$models)]
  if (length(nn_available) == 0L) {
    cli::cli_abort(
      "MC Dropout requires a fitted neural model ({.val cast}, {.val mlp_ate}, or {.val mlp}) in {.arg fit}."
    )
  }
  check_suggested("torch", "for MC Dropout uncertainty estimation")

  # Prepare data
  env_vars <- fit$env_vars
  scaling <- fit$scaling
  X_raw <- as.data.frame(new_data[, env_vars, drop = FALSE])
  X_raw[is.na(X_raw)] <- 0
  X_sc <- as.data.frame(
    scale(X_raw, center = scaling$means, scale = scaling$sds)
  )
  X_sc[is.na(X_sc)] <- 0

  model_name <- nn_available[1L]
  mdl_info <- fit$models[[model_name]]
  if (!identical(mdl_info$type, "nn") || is.null(mdl_info$model)) {
    cli::cli_abort(
      "Selected model {.val {model_name}} is not a fitted neural-network model."
    )
  }
  model <- mdl_info$model

  feature_type <- mdl_info$feature_type %||% "mlp"
  if (feature_type != "mlp") {
    if (is.null(fit$screen) || is.null(fit$dag) || is.null(fit$ate)) {
      cli::cli_abort(
        "Neural models with causal features require {.arg fit} to carry {.val screen}, {.val dag}, and {.val ate}."
      )
    }
    feat <- cast_features(
      X_sc,
      screen = fit$screen,
      dag = fit$dag,
      ate = fit$ate,
      model_type = feature_type
    )
    X_pred <- feat$features
  } else {
    X_pred <- X_sc
  }

  Xt <- torch::torch_tensor(
    as.matrix(X_pred),
    dtype = torch::torch_float32()
  )

  if (verbose) {
    cli::cli_inform(
      "MC Dropout ({model_name}): {n_forward} forward passes on {nrow(new_data)} sites..."
    )
  }

  # MC Dropout: keep model in TRAIN mode so dropout stays active
  model$train()

  mc_preds <- matrix(NA_real_, nrow = nrow(new_data), ncol = n_forward)

  torch::with_no_grad({
    for (k in seq_len(n_forward)) {
      logits <- model(Xt)
      probs <- as.numeric(
        torch::torch_sigmoid(logits)$squeeze()$cpu()
      )
      mc_preds[, k] <- probs
    }
  })

  # Switch back to eval mode
  model$eval()

  # Summarize
  mc_mean <- rowMeans(mc_preds, na.rm = TRUE)
  mc_sd <- apply(mc_preds, 1, stats::sd, na.rm = TRUE)
  mc_cv <- mc_sd / pmax(mc_mean, 1e-6)
  mc_q05 <- apply(mc_preds, 1, stats::quantile, probs = 0.05, na.rm = TRUE)
  mc_q95 <- apply(mc_preds, 1, stats::quantile, probs = 0.95, na.rm = TRUE)

  has_coords <- all(c("lon", "lat") %in% names(new_data))
  out <- if (has_coords) {
    data.frame(
      lon = new_data$lon, lat = new_data$lat,
      mean = mc_mean, sd = mc_sd, cv = mc_cv,
      q05 = mc_q05, q95 = mc_q95
    )
  } else {
    data.frame(
      site = seq_len(nrow(new_data)),
      mean = mc_mean, sd = mc_sd, cv = mc_cv,
      q05 = mc_q05, q95 = mc_q95
    )
  }

  if (verbose) {
    cli::cli_inform(c(
      "Uncertainty estimated for {nrow(out)} sites",
      "i" = "Mean SD: {round(mean(mc_sd), 4)}, Max SD: {round(max(mc_sd), 4)}"
    ))
  }

  out
}
