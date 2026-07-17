# Internal implementation of CAST cross-environment stable selection.

#' @keywords internal
#' @noRd
.cast_shortlist <- function(data, env_vars, response, max_candidates,
                            num_trees, seed) {
  check_suggested("glmnet", "for CAST stable selection")
  check_suggested("ranger", "for CAST stable selection")
  x <- .cast_numeric_matrix(data, env_vars)
  y <- as.integer(data[[response]])
  if (!is.null(seed)) set.seed(seed)

  foldid <- sample(rep(seq_len(min(5L, max(3L, sum(y == 1L) %/% 5L))),
                       length.out = length(y)))
  en <- glmnet::cv.glmnet(
    x, y, family = "binomial", alpha = 0.5, foldid = foldid,
    type.measure = "deviance", standardize = TRUE
  )
  cf <- as.matrix(stats::coef(en, s = "lambda.1se"))
  en_vars <- setdiff(rownames(cf)[abs(cf[, 1L]) > 0], "(Intercept)")

  rf_dat <- as.data.frame(x, check.names = FALSE)
  names(rf_dat) <- make.names(env_vars, unique = TRUE)
  rf_dat$.response <- factor(y)
  rf <- ranger::ranger(
    .response ~ ., data = rf_dat, probability = TRUE,
    num.trees = as.integer(num_trees), importance = "permutation",
    seed = seed %||% 42L, num.threads = 1L
  )
  imp <- rf$variable.importance
  names(imp) <- env_vars[match(names(imp), make.names(env_vars, unique = TRUE))]
  imp <- imp[env_vars]
  imp[!is.finite(imp) | is.na(imp)] <- 0

  en_abs <- stats::setNames(rep(0, length(env_vars)), env_vars)
  coef_names <- intersect(env_vars, rownames(cf))
  en_abs[coef_names] <- abs(cf[coef_names, 1L])
  score <- normalize01(imp) + normalize01(en_abs)
  ranked <- env_vars[order(score, decreasing = TRUE)]
  pool <- unique(c(en_vars, ranked))
  pool <- pool[seq_len(min(length(pool), as.integer(max_candidates)))]

  list(
    candidates = pool,
    rf_importance = unname(imp),
    enet_coefficient = unname(en_abs),
    combined_score = unname(score)
  )
}

#' @keywords internal
#' @noRd
.cast_invariance_test <- function(data, vars, domains, response) {
  y <- as.integer(data[[response]])
  x <- as.data.frame(.cast_numeric_matrix(data, vars))
  names(x) <- paste0("x", seq_along(vars))
  d <- factor(domains)
  df <- data.frame(y = y, domain = d, x, check.names = FALSE)
  rhs <- paste(names(x), collapse = " + ")
  f0 <- stats::as.formula(paste("y ~", rhs))
  f1 <- stats::as.formula(paste("y ~", rhs, "+ domain"))
  f2 <- stats::as.formula(paste("y ~ (", rhs, ") * domain"))

  fits <- tryCatch(
    suppressWarnings(list(
      base = stats::glm(f0, data = df, family = stats::binomial()),
      domain = stats::glm(f1, data = df, family = stats::binomial()),
      interaction = stats::glm(f2, data = df, family = stats::binomial())
    )),
    error = function(e) NULL
  )
  if (is.null(fits)) {
    return(list(p_domain = 0, p_interaction = 0, min_p = 0))
  }
  p1 <- tryCatch(
    stats::anova(fits$base, fits$domain, test = "Chisq")$`Pr(>Chi)`[2L],
    error = function(e) NA_real_
  )
  p2 <- tryCatch(
    stats::anova(fits$domain, fits$interaction, test = "Chisq")$`Pr(>Chi)`[2L],
    error = function(e) NA_real_
  )
  p1 <- if (is.finite(p1)) p1 else 0
  p2 <- if (is.finite(p2)) p2 else 0
  list(p_domain = p1, p_interaction = p2, min_p = min(p1, p2))
}

#' @keywords internal
#' @noRd
.cast_domain_logloss <- function(data, vars, domains, response) {
  losses <- vapply(levels(factor(domains)), function(level) {
    test <- domains == level
    train <- !test
    y_train <- as.integer(data[[response]][train])
    y_test <- as.integer(data[[response]][test])
    if (length(unique(y_train)) < 2L || length(unique(y_test)) < 2L) return(NA_real_)
    x <- as.data.frame(.cast_numeric_matrix(data, vars))
    names(x) <- paste0("x", seq_along(vars))
    df <- data.frame(y = as.integer(data[[response]]), x)
    fit <- tryCatch(
      suppressWarnings(stats::glm(y ~ ., data = df[train, , drop = FALSE],
                                  family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    pred <- tryCatch(stats::predict(fit, df[test, , drop = FALSE], type = "response"),
                     error = function(e) rep(NA_real_, sum(test)))
    pred <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
    if (anyNA(pred)) return(NA_real_)
    -mean(y_test * log(pred) + (1 - y_test) * log(1 - pred))
  }, numeric(1))
  if (all(!is.finite(losses))) Inf else mean(losses, na.rm = TRUE)
}

#' @keywords internal
#' @noRd
.cast_select_stable <- function(data, env_vars, response, domains = NULL,
                                domain_method = "spatial", n_domains = 4L,
                                max_candidates = 12L, min_vars = 3L,
                                alpha = 0.05, loss_tolerance = 0.02,
                                num_trees = 300L, use_invariance = TRUE,
                                seed = NULL, verbose = TRUE) {
  domains <- domains %||% cast_domains(
    data, method = domain_method, k = n_domains, response = response,
    env_vars = env_vars, min_class = 2L, seed = seed
  )
  short <- .cast_shortlist(
    data, env_vars, response, max_candidates, num_trees, seed
  )
  current <- short$candidates
  min_vars <- min(max(2L, as.integer(min_vars)), length(current))

  evaluate <- function(vars) {
    inv <- if (use_invariance) {
      .cast_invariance_test(data, vars, domains, response)
    } else {
      list(p_domain = NA_real_, p_interaction = NA_real_, min_p = 1)
    }
    list(
      vars = vars,
      loss = .cast_domain_logloss(data, vars, domains, response),
      p_domain = inv$p_domain,
      p_interaction = inv$p_interaction,
      min_p = inv$min_p,
      pass = !use_invariance || inv$min_p >= alpha
    )
  }

  current_eval <- evaluate(current)
  path <- list(current_eval)
  best_pass <- if (current_eval$pass) current_eval else NULL
  while (length(current) > min_vars) {
    trials <- lapply(current, function(drop) evaluate(setdiff(current, drop)))
    losses <- vapply(trials, `[[`, numeric(1), "loss")
    min_ps <- vapply(trials, `[[`, numeric(1), "min_p")
    passes <- vapply(trials, `[[`, logical(1), "pass")
    allowed <- is.finite(losses) & losses <= current_eval$loss + loss_tolerance
    if (!any(allowed)) break
    ord <- order(!passes, -min_ps, losses)
    choice <- ord[which(allowed[ord])[1L]]
    candidate <- trials[[choice]]
    improves <- candidate$pass || candidate$min_p > current_eval$min_p + 1e-8
    if (!improves) break
    current <- candidate$vars
    current_eval <- candidate
    path[[length(path) + 1L]] <- current_eval
    if (current_eval$pass) best_pass <- current_eval
  }
  final <- best_pass %||% current_eval

  scores <- data.frame(
    variable = env_vars,
    rf_importance = short$rf_importance,
    enet_coefficient = short$enet_coefficient,
    combined_score = short$combined_score,
    candidate = env_vars %in% short$candidates,
    selected = env_vars %in% final$vars,
    exclusion_reason = ifelse(
      env_vars %in% final$vars, "selected",
      ifelse(env_vars %in% short$candidates, "removed_by_invariance_search", "not_shortlisted")
    ),
    stringsAsFactors = FALSE
  )
  scores <- scores[order(!scores$selected, -scores$combined_score), , drop = FALSE]
  rownames(scores) <- NULL
  diagnostics <- list(
    domains = factor(domains),
    n_domains = nlevels(factor(domains)),
    candidate_variables = short$candidates,
    p_domain = final$p_domain,
    p_interaction = final$p_interaction,
    invariance_pass = final$pass,
    domain_logloss = final$loss,
    search_path = path,
    use_invariance = use_invariance
  )
  if (verbose) {
    cli::cli_inform(
      "CAST stable selection: {length(final$vars)}/{length(env_vars)} variables; domains={diagnostics$n_domains}; invariance={diagnostics$invariance_pass}."
    )
  }
  list(selected = final$vars, scores = scores, diagnostics = diagnostics)
}
