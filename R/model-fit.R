#' Fit Species Distribution Models
#'
#' Trains one or more SDM models on prepared data. Supports the causal CI-MLP
#' architecture (requires \pkg{torch}) and traditional methods (RF, MaxEnt,
#' BRT).
#'
#' @param data A `data.frame` with `presence` column and predictor variables.
#' @param screen A [cast_screen] object, or `NULL`.
#' @param dag A [cast_dag] object, or `NULL`.
#' @param ate A [cast_ate] object, or `NULL`.
#' @param models Character vector. Models to fit: `"cast"`, `"mlp_ate"`,
#'   `"mlp"`, `"rf"`, `"maxent"`, `"brt"`. Default `c("cast", "rf")`.
#' @param response Character. Response column name. Default `"presence"`.
#' @param n_epochs Integer. Training epochs for neural networks. Default `200`.
#' @param n_runs Integer. Number of random seeds for NN ensembling. Default `3`.
#' @param patience Integer. Early stopping patience (epochs). Default `20`.
#' @param val_fraction Numeric. Validation split for early stopping. Default
#'   `0.2`.
#' @param focal_gamma Numeric. Focal loss gamma parameter. Default `2.0`.
#' @param rf_ntree Integer. Number of RF trees. Default `300`.
#' @param brt_n_trees Integer. Number of BRT iterations. Default `500`.
#' @param brt_depth Integer. BRT tree depth. Default `5`.
#' @param hidden_size Integer or `NULL`. NN hidden layer size. `NULL` for auto.
#' @param dropout Numeric. NN dropout rate. Default `0.2`.
#' @param lr Numeric. Learning rate. Default `1e-3`.
#' @param batch_size Integer or `NULL`. `NULL` for auto-tuning.
#' @param seed Integer or `NULL`. Base random seed.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return A `cast_fit` object containing fitted models and metadata.
#'
#' @details
#' ## CI-MLP Architecture (CAST, MLP_ATE, MLP)
#' A 5-layer feedforward network with Layer Normalization, SiLU activation,
#' and dropout. Trained with AdamW optimizer, warmup + cosine annealing
#' schedule, and focal loss for class imbalance handling.
#'
#' Requires the \pkg{torch} package. If \pkg{torch} is not installed, neural
#' network models will be skipped with a warning.
#'
#' ## Traditional Models
#' - **RF**: [ranger::ranger()] with probability output.
#' - **MaxEnt**: [maxnet::maxnet()] with logistic output.
#' - **BRT**: [gbm::gbm()] with Bernoulli loss and 5-fold CV.
#'
#' @seealso [cast_features()], [cast_evaluate()], [cast_predict()]
#'
#' @export
cast_fit <- function(data,
                     screen = NULL,
                     dag = NULL,
                     ate = NULL,
                     models = c("cast", "rf"),
                     response = "presence",
                     n_epochs = 200L,
                     n_runs = 3L,
                     patience = 20L,
                     val_fraction = 0.2,
                     focal_gamma = 2.0,
                     rf_ntree = 300L,
                     brt_n_trees = 500L,
                     brt_depth = 5L,
                     hidden_size = NULL,
                     dropout = 0.2,
                     lr = 1e-3,
                     batch_size = NULL,
                     seed = NULL,
                     verbose = TRUE) {
  models <- tolower(models)
  nn_models <- intersect(models, c("cast", "mlp_ate", "mlp"))
  trad_models <- intersect(models, c("rf", "maxent", "brt"))

  # Validate NN prerequisites
  if (length(nn_models) > 0) {
    if (!has_torch()) {
      cli::cli_warn(
        "torch not available. Skipping NN models: {.val {nn_models}}."
      )
      nn_models <- character(0)
    }
    if ("cast" %in% nn_models || "mlp_ate" %in% nn_models) {
      if (is.null(screen) || is.null(dag) || is.null(ate)) {
        cli::cli_abort(
          "CAST/MLP_ATE models require {.arg screen}, {.arg dag}, and {.arg ate}."
        )
      }
    }
  }

  env_vars <- if (!is.null(dag)) dag$nodes else get_env_vars(data, response)
  Y <- data[[response]]
  X_raw <- as.data.frame(data[, env_vars, drop = FALSE])
  X_raw[is.na(X_raw)] <- 0

  # -- Standardize --
  X_means <- colMeans(X_raw, na.rm = TRUE)
  X_sds <- apply(X_raw, 2, stats::sd, na.rm = TRUE)
  X_sds[X_sds < 1e-10] <- 1
  X_sc <- as.data.frame(scale(X_raw, center = X_means, scale = X_sds))
  X_sc[is.na(X_sc)] <- 0

  # -- Validation split (stratified) --
  if (!is.null(seed)) set.seed(seed + 1000L)
  pos_idx <- which(Y == 1)
  neg_idx <- which(Y == 0)
  val_pos <- sample(pos_idx, round(val_fraction * length(pos_idx)))
  val_neg <- sample(neg_idx, round(val_fraction * length(neg_idx)))
  val_idx <- c(val_pos, val_neg)
  y_val <- Y[val_idx]
  y_tr <- Y[-val_idx]

  focal_alpha <- 1 - mean(y_tr)
  bs <- batch_size %||%
    min(128L, max(32L, as.integer(length(y_tr) / 100)))

  fitted_models <- list()
  cast_vars <- if (!is.null(screen)) screen$selected else env_vars

  # ======================================================================
  # Neural Network Models
  # ======================================================================
  seeds <- if (!is.null(seed)) seed + seq_len(n_runs) else
    sample.int(1e6, n_runs)

  for (mdl in nn_models) {
    if (verbose) cli::cli_inform("Training {.val {mdl}}...")

    # Build features
    feat_type <- switch(mdl, cast = "cast", mlp_ate = "mlp_ate", "mlp")
    if (feat_type != "mlp") {
      feat_all <- cast_features(X_sc, screen, dag, ate, model_type = feat_type)
      X_feat <- feat_all$features
    } else {
      X_feat <- X_sc
    }

    X_tr_nn <- X_feat[-val_idx, , drop = FALSE]
    X_val_nn <- X_feat[val_idx, , drop = FALSE]
    n_input <- ncol(X_tr_nn)
    hs <- hidden_size %||%
      max(32L, min(128L, as.integer(n_input * 4)))

    run_aucs <- numeric(n_runs)
    run_trained <- vector("list", n_runs)

    for (ri in seq_len(n_runs)) {
      tryCatch({
        torch::torch_manual_seed(seeds[ri])
        set.seed(seeds[ri])

        ds <- build_flat_dataset(X_tr_nn, y_tr)
        dl <- torch::dataloader(ds, batch_size = bs,
                                shuffle = TRUE, drop_last = TRUE)
        vt <- torch::torch_tensor(
          as.matrix(X_val_nn), dtype = torch::torch_float()
        )

        m <- build_ci_mlp(n_input, hs, dropout)
        res <- train_ci_mlp(
          m, dl, vt, y_val,
          epochs = n_epochs, lr = lr, patience = patience,
          focal_alpha = focal_alpha, focal_gamma = focal_gamma
        )

        run_aucs[ri] <- res$best_val_auc
        run_trained[[ri]] <- res$model
      }, error = function(e) {
        run_aucs[ri] <<- NA_real_
        if (verbose) cli::cli_warn("  Run {ri} failed: {e$message}")
      })
    }

    # Select best run
    best_ri <- which.max(run_aucs)
    best_model <- if (length(best_ri) > 0) run_trained[[best_ri]] else NULL

    fitted_models[[mdl]] <- list(
      type = "nn",
      model = best_model,
      auc_runs = run_aucs,
      feature_type = feat_type,
      n_input = n_input,
      hidden_size = hs
    )

    if (verbose && !is.null(best_model)) {
      cli::cli_inform(
        "  {mdl}: best val AUC = {round(max(run_aucs, na.rm = TRUE), 4)}"
      )
    }
  }

  # ======================================================================
  # Traditional Models
  # ======================================================================
  for (mdl in trad_models) {
    if (verbose) cli::cli_inform("Training {.val {mdl}}...")
    fitted_models[[mdl]] <- tryCatch(
      fit_traditional(mdl, X_raw, Y, rf_ntree, brt_n_trees, brt_depth, seed),
      error = function(e) {
        cli::cli_warn("{mdl} failed: {e$message}")
        list(type = "traditional", model = NULL, name = mdl)
      }
    )
  }

  new_cast_fit(
    models = fitted_models,
    cast_vars = cast_vars,
    env_vars = env_vars,
    scaling = list(means = X_means, sds = X_sds),
    dag = dag, ate = ate, screen = screen
  )
}


# ========================================================================
# Internal: CI-MLP torch module
# ========================================================================

#' Build CI-MLP torch Module
#' @keywords internal
#' @noRd
build_ci_mlp <- function(n_input, hidden = 64L, dropout = 0.2) {
  if (!has_torch()) cli::cli_abort("torch is required.")
  half_h <- as.integer(hidden %/% 2)
  torch::nn_module(
    "CI_MLP",
    initialize = function(n_in, h, dr) {
      self$net <- torch::nn_sequential(
        torch::nn_linear(n_in, h),
        torch::nn_layer_norm(h),
        torch::nn_silu(),
        torch::nn_dropout(dr),
        torch::nn_linear(h, h),
        torch::nn_layer_norm(h),
        torch::nn_silu(),
        torch::nn_dropout(dr),
        torch::nn_linear(h, h),
        torch::nn_layer_norm(h),
        torch::nn_silu(),
        torch::nn_dropout(dr),
        torch::nn_linear(h, as.integer(h %/% 2)),
        torch::nn_layer_norm(as.integer(h %/% 2)),
        torch::nn_silu(),
        torch::nn_dropout(dr * 0.5),
        torch::nn_linear(as.integer(h %/% 2), 1L)
      )
    },
    forward = function(x) {
      self$net(x)
    }
  )(n_input, hidden, dropout)
}


#' Focal Loss
#' @keywords internal
#' @noRd
focal_loss <- function(logits, targets, alpha = 0.25, gamma = 2.0) {
  bce <- torch::nn_bce_with_logits_loss(reduction = "none")(logits, targets)
  probs <- torch::torch_sigmoid(logits)
  pt <- targets * probs + (1 - targets) * (1 - probs)
  focal_w <- alpha * (1 - pt)^gamma
  (focal_w * bce)$mean()
}


#' Train CI-MLP
#' @keywords internal
#' @noRd
train_ci_mlp <- function(model, train_dl, val_tensor, y_val_vec,
                         epochs = 200, lr = 1e-3, wd = 1e-4,
                         patience = 30, warmup_epochs = 10,
                         focal_alpha = 0.25, focal_gamma = 2.0) {
  optimizer <- torch::optim_adamw(
    model$parameters, lr = lr, weight_decay = wd
  )
  best_auc <- 0
  best_state <- NULL
  no_imp <- 0L

  for (epoch in seq_len(epochs)) {
    # Learning rate schedule: warmup + cosine annealing
    current_lr <- if (epoch <= warmup_epochs) {
      lr * epoch / warmup_epochs
    } else {
      1e-5 + 0.5 * (lr - 1e-5) *
        (1 + cos(pi * (epoch - warmup_epochs) /
                   (epochs - warmup_epochs)))
    }
    for (pg in optimizer$param_groups) pg$lr <- current_lr

    # Training
    model$train()
    coro::loop(for (batch in train_dl) {
      optimizer$zero_grad()
      logits <- model(batch$x)
      loss <- focal_loss(
        logits, batch$y,
        alpha = focal_alpha, gamma = focal_gamma
      )
      if (is.nan(loss$item())) next
      loss$backward()
      torch::nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
      optimizer$step()
    })

    # Validation
    model$eval()
    torch::with_no_grad({
      vp <- as.numeric(
        torch::torch_sigmoid(model(val_tensor))$squeeze()$cpu()
      )
    })

    va <- if (any(is.nan(vp))) 0 else {
      tryCatch(
        as.numeric(pROC::auc(pROC::roc(y_val_vec, vp, quiet = TRUE))),
        error = function(e) 0
      )
    }

    if (va > best_auc + 1e-4) {
      best_auc <- va
      best_state <- lapply(
        model$state_dict(), function(p) p$clone()
      )
      no_imp <- 0L
    } else {
      no_imp <- no_imp + 1L
    }
    if (no_imp >= patience) break
  }

  if (!is.null(best_state)) model$load_state_dict(best_state)
  list(model = model, best_val_auc = best_auc)
}


#' Build Flat Dataset for torch
#' @keywords internal
#' @noRd
build_flat_dataset <- function(X, y) {
  torch::dataset("FlatDS",
    initialize = function(X, y) {
      self$x <- torch::torch_tensor(
        as.matrix(X), dtype = torch::torch_float()
      )
      self$y <- torch::torch_tensor(
        y, dtype = torch::torch_float()
      )$unsqueeze(2)
    },
    .getitem = function(i) {
      list(x = self$x[i, ], y = self$y[i, ])
    },
    .length = function() {
      self$y$size(1)
    }
  )(X, y)
}


#' Fit Traditional SDM
#' @keywords internal
#' @noRd
fit_traditional <- function(name, X, Y, rf_ntree, brt_n_trees,
                            brt_depth, seed) {
  switch(name,
    "rf" = {
      check_suggested("ranger", "for Random Forest")
      if (!is.null(seed)) set.seed(seed)
      m <- ranger::ranger(
        presence ~ .,
        data = cbind(presence = as.factor(Y), X),
        num.trees = rf_ntree, probability = TRUE, seed = seed %||% 42L
      )
      list(type = "traditional", model = m, name = "rf")
    },
    "maxent" = {
      check_suggested("maxnet", "for MaxEnt")
      m <- maxnet::maxnet(
        p = Y, data = X,
        maxnet::maxnet.formula(p = Y, data = X)
      )
      list(type = "traditional", model = m, name = "maxent")
    },
    "brt" = {
      check_suggested("gbm", "for BRT")
      if (!is.null(seed)) set.seed(seed)
      m <- gbm::gbm(
        presence ~ .,
        data = cbind(presence = Y, X),
        distribution = "bernoulli",
        n.trees = brt_n_trees,
        interaction.depth = brt_depth,
        shrinkage = 0.01,
        cv.folds = 5L,
        verbose = FALSE
      )
      bt <- gbm::gbm.perf(m, method = "cv", plot.it = FALSE)
      list(type = "traditional", model = m, name = "brt", best_trees = bt)
    }
  )
}
