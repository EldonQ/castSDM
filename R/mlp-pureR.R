# Pure-R Residual MLP backend for CI-MLP -------------------------------------
#
# A torch-free implementation of the CI-MLP architecture used in cast_fit().
# Mirrors the topology of build_ci_mlp() in model-fit.R:
#
#   input -> proj+LN+SiLU
#         -> 3x residual block: pre-LN -> linear -> SiLU -> LN -> linear -> SiLU -> dropout
#         -> bottleneck (linear half_h + SiLU)
#         -> head (linear -> 1)
#
# Trained with AdamW + warmup + cosine LR schedule + focal loss + early stop.
# All matrix ops use base R (BLAS-backed) for portability; no compiled code.
#
# This backend is the default for cast_fit(backend = "r"). The legacy torch
# backend is retained for users who want to keep it (`backend = "torch"`).
# -----------------------------------------------------------------------------


# ---- Activation & loss helpers ---------------------------------------------

silu <- function(x) x / (1 + exp(-x))
silu_grad <- function(x) {
  s <- 1 / (1 + exp(-x))
  s + x * s * (1 - s)
}
sigmoid <- function(x) 1 / (1 + exp(-x))


# ---- LayerNorm forward / backward ------------------------------------------

ln_forward <- function(X, gamma, beta, eps = 1e-5) {
  mu <- rowMeans(X)
  Xc <- X - mu
  var_ <- rowMeans(Xc^2)
  inv_std <- 1 / sqrt(var_ + eps)
  Xn <- Xc * inv_std
  Y  <- sweep(Xn, 2, gamma, "*")
  Y  <- sweep(Y,  2, beta,  "+")
  list(Y = Y, Xn = Xn, inv_std = inv_std)
}

ln_backward <- function(dY, cache, gamma) {
  Xn <- cache$Xn
  inv_std <- cache$inv_std
  N <- ncol(Xn)
  dgamma <- colSums(dY * Xn)
  dbeta  <- colSums(dY)
  dXn <- sweep(dY, 2, gamma, "*")
  mean_dXn <- rowMeans(dXn)
  mean_dXn_Xn <- rowMeans(dXn * Xn)
  dX <- (dXn - mean_dXn - Xn * mean_dXn_Xn) * inv_std
  list(dX = dX, dgamma = dgamma, dbeta = dbeta)
}


# ---- Initialise model parameters -------------------------------------------

#' Initialise pure-R residual MLP weights.
#' @keywords internal
#' @noRd
init_mlp_params <- function(n_in, hidden, dropout = 0.2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  half_h <- as.integer(hidden %/% 2L)

  # Kaiming-style init for SiLU (gain ~ sqrt(2)).
  he <- function(n_in, n_out) {
    matrix(stats::rnorm(n_in * n_out, sd = sqrt(2 / n_in)),
           nrow = n_in, ncol = n_out)
  }
  bzero <- function(n) numeric(n)
  ln_pair <- function(d) list(gamma = rep(1, d), beta = numeric(d))

  list(
    n_in = n_in, hidden = hidden, half_h = half_h, dropout = dropout,
    proj_W = he(n_in, hidden), proj_b = bzero(hidden),
    proj_ln = ln_pair(hidden),
    blocks = lapply(seq_len(3L), function(i) list(
      n1   = ln_pair(hidden),
      a_W  = he(hidden, hidden),
      a_b  = bzero(hidden),
      n2   = ln_pair(hidden),
      b_W  = he(hidden, hidden),
      b_b  = bzero(hidden)
    )),
    neck_ln = ln_pair(hidden),
    neck_W  = he(hidden, half_h),
    neck_b  = bzero(half_h),
    head_W  = matrix(stats::rnorm(half_h, sd = sqrt(2 / half_h)),
                     nrow = half_h, ncol = 1L),
    head_b  = 0
  )
}


# ---- Forward pass with optional caches for backprop ------------------------

forward_mlp <- function(X, p, training = TRUE, drop_mask_seed = NULL) {
  cache <- list(X = X)

  # Input projection
  z0_pre <- X %*% p$proj_W + matrix(p$proj_b, nrow = nrow(X), ncol = p$hidden,
                                    byrow = TRUE)
  ln <- ln_forward(z0_pre, p$proj_ln$gamma, p$proj_ln$beta)
  z0_act_in <- ln$Y
  z0 <- silu(z0_act_in)
  cache$proj <- list(z0_pre = z0_pre, ln = ln, z0_act_in = z0_act_in, z0 = z0)

  # 3 residual blocks
  blocks_cache <- list()
  z <- z0
  for (i in seq_len(3L)) {
    blk <- p$blocks[[i]]
    # pre-norm 1
    ln1 <- ln_forward(z, blk$n1$gamma, blk$n1$beta)
    a_in_act <- ln1$Y
    a_act <- silu(a_in_act)
    # linear a
    a_lin <- a_act %*% blk$a_W +
      matrix(blk$a_b, nrow = nrow(z), ncol = p$hidden, byrow = TRUE)
    # pre-norm 2
    ln2 <- ln_forward(a_lin, blk$n2$gamma, blk$n2$beta)
    b_in_act <- ln2$Y
    b_act <- silu(b_in_act)
    # linear b
    b_lin <- b_act %*% blk$b_W +
      matrix(blk$b_b, nrow = nrow(z), ncol = p$hidden, byrow = TRUE)
    # dropout
    if (training && p$dropout > 0) {
      if (!is.null(drop_mask_seed)) set.seed(drop_mask_seed + i)
      mask <- matrix(stats::runif(length(b_lin)) >= p$dropout,
                     nrow = nrow(b_lin), ncol = ncol(b_lin))
      r <- (b_lin * mask) / (1 - p$dropout)
    } else {
      mask <- NULL
      r <- b_lin
    }
    z_new <- z + r
    blocks_cache[[i]] <- list(
      z_in = z, ln1 = ln1, a_in_act = a_in_act, a_act = a_act,
      a_lin = a_lin, ln2 = ln2, b_in_act = b_in_act, b_act = b_act,
      b_lin = b_lin, mask = mask, z_out = z_new
    )
    z <- z_new
  }
  cache$blocks <- blocks_cache

  # Bottleneck
  ln_neck <- ln_forward(z, p$neck_ln$gamma, p$neck_ln$beta)
  neck_in_act <- ln_neck$Y
  neck_act_in <- silu(neck_in_act)
  neck_lin <- neck_act_in %*% p$neck_W +
    matrix(p$neck_b, nrow = nrow(z), ncol = p$half_h, byrow = TRUE)
  neck_out <- silu(neck_lin)
  cache$neck <- list(ln = ln_neck, neck_in_act = neck_in_act,
                     neck_act_in = neck_act_in,
                     neck_lin = neck_lin, neck_out = neck_out)

  # Head
  logits <- as.numeric(neck_out %*% p$head_W) + p$head_b
  cache$logits <- logits
  list(logits = logits, cache = cache)
}


# ---- Backward pass ---------------------------------------------------------

backward_mlp <- function(d_logits, p, cache) {
  grads <- list()
  N <- length(d_logits)

  # Head
  d_head_W <- t(cache$neck$neck_out) %*% matrix(d_logits, ncol = 1L)
  d_head_b <- sum(d_logits)
  d_neck_out <- matrix(d_logits, ncol = 1L) %*% t(p$head_W)

  # Bottleneck
  d_neck_lin <- d_neck_out * silu_grad(cache$neck$neck_lin)
  d_neck_W <- t(cache$neck$neck_act_in) %*% d_neck_lin
  d_neck_b <- colSums(d_neck_lin)
  d_neck_act_in <- d_neck_lin %*% t(p$neck_W)
  d_neck_in_act <- d_neck_act_in * silu_grad(cache$neck$neck_in_act)
  ln_back <- ln_backward(d_neck_in_act, cache$neck$ln, p$neck_ln$gamma)
  d_z_after_blocks <- ln_back$dX
  grads$neck_ln <- list(gamma = ln_back$dgamma, beta = ln_back$dbeta)
  grads$neck_W <- d_neck_W
  grads$neck_b <- d_neck_b
  grads$head_W <- d_head_W
  grads$head_b <- d_head_b

  # Residual blocks (reverse)
  d_z <- d_z_after_blocks
  blocks_grads <- vector("list", 3L)
  for (i in seq.int(3L, 1L)) {
    blk <- p$blocks[[i]]
    bc <- cache$blocks[[i]]
    # z_out = z_in + r   so dz_in_residual_path = d_z and dr = d_z
    d_r <- d_z
    if (!is.null(bc$mask)) {
      d_b_lin <- (d_r * bc$mask) / (1 - p$dropout)
    } else {
      d_b_lin <- d_r
    }
    d_b_W <- t(bc$b_act) %*% d_b_lin
    d_b_b <- colSums(d_b_lin)
    d_b_act <- d_b_lin %*% t(blk$b_W)
    d_b_in_act <- d_b_act * silu_grad(bc$b_in_act)
    ln2_back <- ln_backward(d_b_in_act, bc$ln2, blk$n2$gamma)
    d_a_lin <- ln2_back$dX
    d_a_W <- t(bc$a_act) %*% d_a_lin
    d_a_b <- colSums(d_a_lin)
    d_a_act <- d_a_lin %*% t(blk$a_W)
    d_a_in_act <- d_a_act * silu_grad(bc$a_in_act)
    ln1_back <- ln_backward(d_a_in_act, bc$ln1, blk$n1$gamma)
    d_z_in_through_branch <- ln1_back$dX
    d_z <- d_z + d_z_in_through_branch
    blocks_grads[[i]] <- list(
      n1 = list(gamma = ln1_back$dgamma, beta = ln1_back$dbeta),
      a_W = d_a_W, a_b = d_a_b,
      n2 = list(gamma = ln2_back$dgamma, beta = ln2_back$dbeta),
      b_W = d_b_W, b_b = d_b_b
    )
  }
  grads$blocks <- blocks_grads

  # Input projection
  d_z0 <- d_z
  d_z0_act_in <- d_z0 * silu_grad(cache$proj$z0_act_in)
  ln_proj_back <- ln_backward(d_z0_act_in, cache$proj$ln, p$proj_ln$gamma)
  d_z0_pre <- ln_proj_back$dX
  d_proj_W <- t(cache$X) %*% d_z0_pre
  d_proj_b <- colSums(d_z0_pre)
  grads$proj_W <- d_proj_W
  grads$proj_b <- d_proj_b
  grads$proj_ln <- list(gamma = ln_proj_back$dgamma, beta = ln_proj_back$dbeta)

  grads
}


# ---- Focal loss ------------------------------------------------------------

focal_loss_grad <- function(logits, y, alpha = 0.25, gamma = 2.0) {
  # Returns list(loss, d_logits) for binary focal loss.
  p_ <- sigmoid(logits)
  pt <- y * p_ + (1 - y) * (1 - p_)
  pt <- pmin(pmax(pt, 1e-7), 1 - 1e-7)
  loss <- -alpha * (1 - pt)^gamma * log(pt)

  # d_loss / d_logits
  # focal = -alpha * (1-pt)^gamma * log(pt)
  # d/dlogit [pt] = (2y - 1) * p_ * (1 - p_) ... but easier via chain on (p_)
  # derive numerically-stable form via standard trick:
  # ce = -[y*log(p)+(1-y)*log(1-p)],  d_ce/d_logit = p - y
  ce_grad <- p_ - y
  # d focal / d logit:
  # f = -alpha (1-pt)^g log(pt)
  # df/dpt = alpha[(1-pt)^g / pt + g (1-pt)^(g-1) log(pt)]    note sign care
  df_dpt <- alpha * ((1 - pt)^gamma / pt - gamma * (1 - pt)^(gamma - 1) * log(pt))
  # d pt / d logit = (2y-1) * p * (1-p)
  dpt_dlogit <- (2 * y - 1) * p_ * (1 - p_)
  d_logits <- df_dpt * dpt_dlogit
  list(loss = mean(loss), d_logits = d_logits / length(y))
}


# ---- AdamW optimiser state -------------------------------------------------

adamw_init <- function(p) {
  zeros_like <- function(x) {
    if (is.matrix(x)) matrix(0, nrow = nrow(x), ncol = ncol(x))
    else numeric(length(x))
  }
  list(
    proj_W = zeros_like(p$proj_W), proj_b = zeros_like(p$proj_b),
    proj_ln = list(gamma = zeros_like(p$proj_ln$gamma),
                   beta  = zeros_like(p$proj_ln$beta)),
    blocks = lapply(p$blocks, function(b) list(
      n1 = list(gamma = zeros_like(b$n1$gamma), beta = zeros_like(b$n1$beta)),
      a_W = zeros_like(b$a_W), a_b = zeros_like(b$a_b),
      n2 = list(gamma = zeros_like(b$n2$gamma), beta = zeros_like(b$n2$beta)),
      b_W = zeros_like(b$b_W), b_b = zeros_like(b$b_b)
    )),
    neck_ln = list(gamma = zeros_like(p$neck_ln$gamma),
                   beta  = zeros_like(p$neck_ln$beta)),
    neck_W  = zeros_like(p$neck_W),  neck_b = zeros_like(p$neck_b),
    head_W  = zeros_like(p$head_W),  head_b = 0
  )
}

adamw_step <- function(p, g, m, v, t, lr = 1e-3, wd = 1e-4,
                       b1 = 0.9, b2 = 0.999, eps = 1e-8) {
  bc1 <- 1 - b1^t
  bc2 <- 1 - b2^t
  upd <- function(theta, grad, mt, vt) {
    mt2 <- b1 * mt + (1 - b1) * grad
    vt2 <- b2 * vt + (1 - b2) * grad^2
    m_hat <- mt2 / bc1
    v_hat <- vt2 / bc2
    theta2 <- theta - lr * (m_hat / (sqrt(v_hat) + eps) + wd * theta)
    list(theta = theta2, m = mt2, v = vt2)
  }
  # proj
  s <- upd(p$proj_W, g$proj_W, m$proj_W, v$proj_W); p$proj_W <- s$theta
  m$proj_W <- s$m; v$proj_W <- s$v
  s <- upd(p$proj_b, g$proj_b, m$proj_b, v$proj_b); p$proj_b <- s$theta
  m$proj_b <- s$m; v$proj_b <- s$v
  s <- upd(p$proj_ln$gamma, g$proj_ln$gamma, m$proj_ln$gamma, v$proj_ln$gamma)
  p$proj_ln$gamma <- s$theta; m$proj_ln$gamma <- s$m; v$proj_ln$gamma <- s$v
  s <- upd(p$proj_ln$beta, g$proj_ln$beta, m$proj_ln$beta, v$proj_ln$beta)
  p$proj_ln$beta <- s$theta; m$proj_ln$beta <- s$m; v$proj_ln$beta <- s$v
  # blocks
  for (i in seq_len(3L)) {
    pb <- p$blocks[[i]]; gb <- g$blocks[[i]]; mb <- m$blocks[[i]]; vb <- v$blocks[[i]]
    s <- upd(pb$n1$gamma, gb$n1$gamma, mb$n1$gamma, vb$n1$gamma)
    pb$n1$gamma <- s$theta; mb$n1$gamma <- s$m; vb$n1$gamma <- s$v
    s <- upd(pb$n1$beta, gb$n1$beta, mb$n1$beta, vb$n1$beta)
    pb$n1$beta <- s$theta; mb$n1$beta <- s$m; vb$n1$beta <- s$v
    s <- upd(pb$a_W, gb$a_W, mb$a_W, vb$a_W); pb$a_W <- s$theta
    mb$a_W <- s$m; vb$a_W <- s$v
    s <- upd(pb$a_b, gb$a_b, mb$a_b, vb$a_b); pb$a_b <- s$theta
    mb$a_b <- s$m; vb$a_b <- s$v
    s <- upd(pb$n2$gamma, gb$n2$gamma, mb$n2$gamma, vb$n2$gamma)
    pb$n2$gamma <- s$theta; mb$n2$gamma <- s$m; vb$n2$gamma <- s$v
    s <- upd(pb$n2$beta, gb$n2$beta, mb$n2$beta, vb$n2$beta)
    pb$n2$beta <- s$theta; mb$n2$beta <- s$m; vb$n2$beta <- s$v
    s <- upd(pb$b_W, gb$b_W, mb$b_W, vb$b_W); pb$b_W <- s$theta
    mb$b_W <- s$m; vb$b_W <- s$v
    s <- upd(pb$b_b, gb$b_b, mb$b_b, vb$b_b); pb$b_b <- s$theta
    mb$b_b <- s$m; vb$b_b <- s$v
    p$blocks[[i]] <- pb; m$blocks[[i]] <- mb; v$blocks[[i]] <- vb
  }
  # neck + head
  s <- upd(p$neck_ln$gamma, g$neck_ln$gamma, m$neck_ln$gamma, v$neck_ln$gamma)
  p$neck_ln$gamma <- s$theta; m$neck_ln$gamma <- s$m; v$neck_ln$gamma <- s$v
  s <- upd(p$neck_ln$beta, g$neck_ln$beta, m$neck_ln$beta, v$neck_ln$beta)
  p$neck_ln$beta <- s$theta; m$neck_ln$beta <- s$m; v$neck_ln$beta <- s$v
  s <- upd(p$neck_W, g$neck_W, m$neck_W, v$neck_W); p$neck_W <- s$theta
  m$neck_W <- s$m; v$neck_W <- s$v
  s <- upd(p$neck_b, g$neck_b, m$neck_b, v$neck_b); p$neck_b <- s$theta
  m$neck_b <- s$m; v$neck_b <- s$v
  s <- upd(p$head_W, g$head_W, m$head_W, v$head_W); p$head_W <- s$theta
  m$head_W <- s$m; v$head_W <- s$v
  m$head_b <- b1 * m$head_b + (1 - b1) * g$head_b
  v$head_b <- b2 * v$head_b + (1 - b2) * g$head_b^2
  p$head_b <- p$head_b - lr * (m$head_b / bc1 / (sqrt(v$head_b / bc2) + eps) +
                                 wd * p$head_b)
  list(p = p, m = m, v = v)
}


# ============================================================================
# Public API
# ============================================================================

#' Fit a Pure-R Residual MLP (CPU-Native CI-MLP Backend)
#'
#' A torch-free, dependency-free implementation of the residual MLP used by
#' [cast_fit()] when `backend = "r"` (the default). Mirrors the architecture
#' of the legacy torch backend: input projection + LayerNorm + SiLU,
#' three residual blocks (pre-norm, two linear-SiLU layers, dropout),
#' bottleneck + head. Trained with AdamW + warmup + cosine LR schedule and
#' focal loss for class imbalance.
#'
#' @param X Numeric matrix or data.frame. Predictor matrix (already scaled).
#' @param y Numeric vector of 0/1 binary labels.
#' @param hidden Integer. Hidden width. Default `64`.
#' @param dropout Numeric in `[0, 1)`. Default `0.2`.
#' @param epochs Integer. Maximum epochs. Default `300`.
#' @param batch_size Integer. Mini-batch size. Default `64`.
#' @param lr Numeric. Peak learning rate. Default `1e-3`.
#' @param weight_decay Numeric. AdamW weight decay. Default `1e-4`.
#' @param warmup_epochs Integer. Linear warmup length. Default `10`.
#' @param patience Integer. Early-stopping patience (epochs without
#'   validation AUC improvement). Default `40`.
#' @param val_idx Integer vector. Indices of validation rows. If `NULL`,
#'   a stratified 20% split is built internally.
#' @param val_fraction Numeric. Used when `val_idx = NULL`. Default `0.2`.
#' @param focal_alpha,focal_gamma Numeric. Focal loss parameters.
#' @param seed Integer. Random seed for reproducibility.
#' @param verbose Logical. Print per-epoch progress every 25 epochs.
#'
#' @return A list with class `"cast_mlp"` containing fitted weights, the
#'   training history, and the configuration. Compatible with
#'   [cast_mlp_predict()].
#'
#' @seealso [cast_mlp_predict()], [cast_fit()]
#' @export
cast_mlp_fit <- function(X, y,
                         hidden        = 64L,
                         dropout       = 0.2,
                         epochs        = 300L,
                         batch_size    = 64L,
                         lr            = 1e-3,
                         weight_decay  = 1e-4,
                         warmup_epochs = 10L,
                         patience      = 40L,
                         val_idx       = NULL,
                         val_fraction  = 0.2,
                         focal_alpha   = 0.25,
                         focal_gamma   = 2.0,
                         seed          = NULL,
                         verbose       = FALSE) {
  X <- as.matrix(X)
  if (!is.numeric(y)) y <- as.numeric(y)
  storage.mode(X) <- "double"
  N <- nrow(X)
  if (length(y) != N) cli::cli_abort("length(y) must equal nrow(X).")

  # Validation split
  if (is.null(val_idx)) {
    if (!is.null(seed)) set.seed(seed + 1000L)
    pos <- which(y == 1); neg <- which(y == 0)
    val_pos <- sample(pos, max(1L, round(val_fraction * length(pos))))
    val_neg <- sample(neg, max(1L, round(val_fraction * length(neg))))
    val_idx <- c(val_pos, val_neg)
  }
  X_tr <- X[-val_idx, , drop = FALSE]
  y_tr <- y[-val_idx]
  X_va <- X[ val_idx, , drop = FALSE]
  y_va <- y[ val_idx]

  p <- init_mlp_params(ncol(X_tr), hidden, dropout, seed)
  m_state <- adamw_init(p); v_state <- adamw_init(p); v_state$head_b <- 0
  best_auc <- -Inf; best_p <- p; no_imp <- 0L
  hist <- list(epoch = integer(), train_loss = numeric(), val_auc = numeric())
  step_counter <- 0L

  for (epoch in seq_len(epochs)) {
    # Cosine warmup schedule
    cur_lr <- if (epoch <= warmup_epochs) {
      lr * epoch / warmup_epochs
    } else {
      1e-5 + 0.5 * (lr - 1e-5) *
        (1 + cos(pi * (epoch - warmup_epochs) / max(1, epochs - warmup_epochs)))
    }

    if (!is.null(seed)) set.seed(seed + 13L * epoch)
    perm <- sample(seq_len(nrow(X_tr)))
    train_loss_sum <- 0
    n_batches <- 0L

    for (start in seq(1, length(perm), by = batch_size)) {
      end <- min(start + batch_size - 1L, length(perm))
      bi <- perm[start:end]
      Xb <- X_tr[bi, , drop = FALSE]
      yb <- y_tr[bi]
      step_counter <- step_counter + 1L
      fwd <- forward_mlp(Xb, p, training = TRUE,
                         drop_mask_seed = if (!is.null(seed))
                           seed + 1000L * epoch + step_counter else NULL)
      fl <- focal_loss_grad(fwd$logits, yb,
                            alpha = focal_alpha, gamma = focal_gamma)
      train_loss_sum <- train_loss_sum + fl$loss
      n_batches <- n_batches + 1L
      grads <- backward_mlp(fl$d_logits, p, fwd$cache)
      step <- adamw_step(p, grads, m_state, v_state, step_counter,
                         lr = cur_lr, wd = weight_decay)
      p <- step$p; m_state <- step$m; v_state <- step$v
    }

    fwd_val <- forward_mlp(X_va, p, training = FALSE)
    val_p <- sigmoid(fwd_val$logits)
    auc <- compute_auc(y_va, val_p)
    hist$epoch <- c(hist$epoch, epoch)
    hist$train_loss <- c(hist$train_loss, train_loss_sum / max(1L, n_batches))
    hist$val_auc <- c(hist$val_auc, auc)

    if (verbose && epoch %% 25 == 0L) {
      cli::cli_inform(
        "  epoch {epoch}: lr={signif(cur_lr,3)}  loss={signif(train_loss_sum/max(1L,n_batches),4)}  val_AUC={signif(auc,4)}"
      )
    }

    if (is.finite(auc) && auc > best_auc + 1e-4) {
      best_auc <- auc; best_p <- p; no_imp <- 0L
    } else {
      no_imp <- no_imp + 1L
    }
    if (no_imp >= patience) break
  }

  out <- list(
    weights = best_p,
    history = hist,
    config = list(hidden = hidden, dropout = dropout,
                  epochs = epochs, batch_size = batch_size, lr = lr,
                  weight_decay = weight_decay, warmup_epochs = warmup_epochs,
                  patience = patience, focal_alpha = focal_alpha,
                  focal_gamma = focal_gamma, seed = seed),
    best_val_auc = best_auc,
    n_input = ncol(X)
  )
  class(out) <- "cast_mlp"
  out
}


#' Predict from a Fitted Pure-R Residual MLP
#'
#' Runs the forward pass without dropout and applies the sigmoid to obtain
#' calibrated probabilities.
#'
#' @param object A `cast_mlp` object from [cast_mlp_fit()].
#' @param newx Numeric matrix or data.frame of predictors (already scaled
#'   with the same parameters used in fitting).
#' @param ... Unused.
#' @return Numeric vector of probabilities in `[0, 1]`.
#'
#' @seealso [cast_mlp_fit()]
#' @export
cast_mlp_predict <- function(object, newx, ...) {
  if (!inherits(object, "cast_mlp"))
    cli::cli_abort("{.arg object} must be a cast_mlp.")
  X <- as.matrix(newx)
  storage.mode(X) <- "double"
  fwd <- forward_mlp(X, object$weights, training = FALSE)
  as.numeric(sigmoid(fwd$logits))
}


# ---- Internal: tiny AUC for early stopping (no pROC dep) -------------------

#' @keywords internal
#' @noRd
compute_auc <- function(y, p) {
  ord <- order(p)
  y_o <- y[ord]
  n_pos <- sum(y == 1); n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  ranks <- rank(p, ties.method = "average")
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
