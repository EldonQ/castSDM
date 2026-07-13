# Role-constrained causal-core selection.

.cast_prior_roles <- c(
  "direct_candidate", "proxy_only", "mediator", "adjustment_only",
  "sampling_bias", "unknown"
)

.cast_prior_confidence <- c("high", "medium", "low")

#' Build an Auditable Causal-Role Specification
#'
#' Creates a conservative role table for environmental predictors. Roles are
#' scientific assumptions used to constrain selection; they are not causal
#' directions learned from occurrence records. The default habitat template
#' treats climate, soil, topography, hydrology, vegetation and land-cover
#' variables as candidate ecological drivers, excludes anthropogenic sampling
#' proxies, and fails unrecognised variables closed to `"unknown"`.
#'
#' @param variables Character vector of predictor names.
#' @param template Character. Currently `"habitat"` (the validated default) or
#'   `"upstream"` (vegetation and land-cover are treated as mediators).
#' @param overrides Optional data frame containing `variable` and any of
#'   `role`, `confidence`, `rationale`, or `parent`.
#'
#' @return A data frame with one auditable role declaration per variable.
#' @export
cast_causal_spec <- function(variables,
                             template = c("habitat", "upstream"),
                             overrides = NULL) {
  template <- match.arg(template)
  variables <- unique(trimws(as.character(variables)))
  variables <- variables[nzchar(variables)]
  if (!length(variables)) cli::cli_abort("{.arg variables} must not be empty.")

  category <- rep("unknown", length(variables))
  category[grepl("^(bio|prec|sunp|tave|tmax|tmin)|temperature|precip", variables,
                 ignore.case = TRUE)] <- "bioclim"
  category[grepl("soil|soc|clay|silt|sand|ph(_|$)", variables,
                 ignore.case = TRUE)] <- "soil"
  category[grepl("elev|slope|aspect|rugged|topo|terrain|tri($|_)", variables,
                 ignore.case = TRUE)] <- "topography"
  category[grepl("hydro|river|stream|runoff|discharge|wetness|water", variables,
                 ignore.case = TRUE)] <- "hydrology"
  category[grepl("vege|veget|ndvi|evi|forest|pasture|crop|land.?cover|^lu_", variables,
                 ignore.case = TRUE)] <- "habitat"
  category[grepl("anthro|population|people|road|night.?light|gdp|urban", variables,
                 ignore.case = TRUE)] <- "anthropogenic"

  role <- ifelse(category == "anthropogenic", "sampling_bias",
                 ifelse(category == "unknown", "unknown", "direct_candidate"))
  if (template == "upstream") role[category == "habitat"] <- "mediator"
  confidence <- ifelse(category %in% c("bioclim", "soil", "topography", "habitat"),
                       "medium", "low")
  confidence[category == "anthropogenic"] <- "medium"

  rationale <- vapply(seq_along(variables), function(i) {
    switch(category[i],
      bioclim = "Climate is a plausible ecological constraint; this role does not establish a direct mechanism.",
      soil = "Soil conditions are plausible habitat constraints; this role remains species- and scale-dependent.",
      topography = "Topography is a plausible upstream habitat constraint; this declaration is an ecological prior.",
      hydrology = "Hydrology can constrain habitat availability, but its directness requires species-specific review.",
      habitat = if (template == "habitat") {
        "Vegetation or land cover is treated as a habitat mechanism under the validated habitat template."
      } else {
        "Vegetation or land cover is treated as a mediator under the conservative upstream template."
      },
      anthropogenic = "Human-pressure layers can encode observation effort or accessibility and are excluded from the default causal core.",
      "Unrecognised metadata category; excluded from the causal core until explicitly reviewed."
    )
  }, character(1))

  out <- data.frame(
    variable = variables,
    role = role,
    rationale = rationale,
    confidence = confidence,
    parent = "",
    metadata_category = category,
    specification_version = paste0("castSDM_", template, "_v1"),
    stringsAsFactors = FALSE
  )

  if (!is.null(overrides)) {
    if (!is.data.frame(overrides) || !"variable" %in% names(overrides)) {
      cli::cli_abort("{.arg overrides} must be a data frame with a {.field variable} column.")
    }
    if (anyDuplicated(overrides$variable)) {
      cli::cli_abort("{.arg overrides} must contain at most one row per variable.")
    }
    idx <- match(as.character(overrides$variable), out$variable)
    if (anyNA(idx)) {
      cli::cli_abort("Override variables not found in {.arg variables}: {.val {overrides$variable[is.na(idx)]}}.")
    }
    for (nm in intersect(c("role", "confidence", "rationale", "parent"), names(overrides))) {
      out[[nm]][idx] <- as.character(overrides[[nm]])
    }
    out$specification_version <- paste0(out$specification_version, "+reviewed")
  }
  .validate_cast_causal_spec(out, variables)
}

.validate_cast_causal_spec <- function(causal_spec, env_vars) {
  source <- "generated_default"
  spec <- if (is.null(causal_spec)) {
    cast_causal_spec(env_vars)
  } else if (is.character(causal_spec) && length(causal_spec) == 1L) {
    source <- normalizePath(causal_spec, winslash = "/", mustWork = TRUE)
    utils::read.csv(causal_spec, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(causal_spec)) {
    source <- "in_memory_specification"
    causal_spec
  } else {
    cli::cli_abort("{.arg causal_spec} must be NULL, a CSV path, or a data frame.")
  }

  required <- c("variable", "role", "rationale", "confidence")
  missing <- setdiff(required, names(spec))
  if (length(missing)) {
    cli::cli_abort("Causal specification is missing: {.field {missing}}.")
  }
  for (nm in required) spec[[nm]] <- trimws(as.character(spec[[nm]]))
  if (any(!nzchar(spec$variable)) || anyDuplicated(spec$variable)) {
    cli::cli_abort("The causal specification requires one unique non-empty row per variable.")
  }
  bad_roles <- setdiff(unique(spec$role), .cast_prior_roles)
  if (length(bad_roles)) cli::cli_abort("Unsupported causal roles: {.val {bad_roles}}.")
  bad_conf <- setdiff(unique(spec$confidence), .cast_prior_confidence)
  if (length(bad_conf)) cli::cli_abort("Unsupported confidence levels: {.val {bad_conf}}.")
  if (any(!nzchar(spec$rationale))) {
    cli::cli_abort("Every causal-role declaration requires a rationale.")
  }

  spec <- spec[match(intersect(env_vars, spec$variable), spec$variable), , drop = FALSE]
  absent <- setdiff(env_vars, spec$variable)
  if (length(absent)) {
    extra <- data.frame(
      variable = absent,
      role = "unknown",
      rationale = "Absent from the supplied specification; excluded from the causal core.",
      confidence = "low",
      stringsAsFactors = FALSE
    )
    for (nm in setdiff(names(spec), names(extra))) extra[[nm]] <- NA
    extra <- extra[, names(spec), drop = FALSE]
    spec <- rbind(spec, extra)
  }
  if (!"parent" %in% names(spec)) spec$parent <- ""
  if (!"specification_version" %in% names(spec)) {
    spec$specification_version <- "user_specification"
  }
  spec <- spec[match(env_vars, spec$variable), , drop = FALSE]
  spec$spec_source <- source
  rownames(spec) <- NULL
  spec
}

.prior_numeric_frame <- function(data, vars) {
  x <- as.data.frame(data[, vars, drop = FALSE])
  for (nm in names(x)) {
    x[[nm]] <- as.numeric(x[[nm]])
    if (anyNA(x[[nm]])) {
      med <- stats::median(x[[nm]], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      x[[nm]][is.na(x[[nm]])] <- med
    }
  }
  x
}

.prior_normalize <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x) | is.na(x)] <- 0
  if (!length(x)) return(x)
  r <- range(x)
  if (r[2] <= r[1]) return(if (r[2] > 0) rep(1, length(x)) else rep(0, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

.prior_univariate_scores <- function(data, vars, response) {
  y <- as.numeric(data[[response]])
  score <- vapply(vars, function(v) {
    x <- as.numeric(data[[v]])
    if (!is.finite(stats::sd(x, na.rm = TRUE)) || stats::sd(x, na.rm = TRUE) == 0) return(0)
    abs(stats::cor(x, y, use = "pairwise.complete.obs"))
  }, numeric(1))
  score[!is.finite(score)] <- 0
  data.frame(variable = vars, univ_score = score, stringsAsFactors = FALSE)
}

.prior_category <- function(vars) {
  out <- rep("other", length(vars))
  out[grepl("^(bio|prec|sunp|tave|tmax|tmin)", vars, ignore.case = TRUE)] <- "bioclim"
  out[grepl("soil|soc|clay|silt|sand", vars, ignore.case = TRUE)] <- "soil"
  out[grepl("vege|veget|ndvi|evi|forest|pasture|crop|land.?cover|^lu_", vars, ignore.case = TRUE)] <- "habitat"
  out[grepl("hydro|river|stream|runoff|discharge|water", vars, ignore.case = TRUE)] <- "hydrology"
  out[grepl("elev|slope|aspect|rugged|topo|terrain|tri($|_)", vars, ignore.case = TRUE)] <- "topography"
  out[grepl("anthro|population|people|road|night.?light|gdp|urban", vars, ignore.case = TRUE)] <- "anthropogenic"
  out
}

.prior_family <- function(vars) {
  fam <- gsub("_[0-9]+_[0-9]+cm", "", vars)
  fam <- gsub("_r[0-9]+cell", "", fam)
  fam <- gsub("_20[0-9]{2}", "", fam)
  gsub("_(mean|median|min|max|sd|valid_fraction)$", "", fam)
}

.prior_select_diverse <- function(ranked, data, target_n, cor_threshold = 0.7,
                                  family_cap = 2L) {
  ranked <- intersect(unique(ranked[!is.na(ranked)]), names(data))
  selected <- character()
  can_add <- function(v) {
    if (!length(selected)) return(TRUE)
    fam <- .prior_family(c(selected, v))
    if (sum(fam == tail(fam, 1)) > family_cap) return(FALSE)
    cors <- abs(stats::cor(as.numeric(data[[v]]),
                           .prior_numeric_frame(data, selected),
                           use = "pairwise.complete.obs"))
    !any(cors >= cor_threshold, na.rm = TRUE)
  }
  for (cat in unique(.prior_category(ranked))) {
    for (v in ranked[.prior_category(ranked) == cat]) {
      if (length(selected) >= target_n) break
      if (!v %in% selected && can_add(v)) {
        selected <- c(selected, v)
        break
      }
    }
  }
  for (v in ranked) {
    if (length(selected) >= target_n) break
    if (!v %in% selected && can_add(v)) selected <- c(selected, v)
  }
  if (length(selected) < target_n) {
    add <- setdiff(ranked, selected)
    selected <- c(selected, head(add, target_n - length(selected)))
  }
  selected
}

.prior_signal_pool <- function(data, vars, response, pool_max = 35L) {
  scores <- .prior_univariate_scores(data, vars, response)
  ranked <- scores$variable[order(-scores$univ_score)]
  .prior_select_diverse(ranked, data, min(pool_max, length(ranked)),
                        cor_threshold = 0.85, family_cap = 2L)
}

.prior_blocks <- function(data, vars, n_blocks = 5L, seed = NULL) {
  if (".causal_env" %in% names(data)) {
    out <- as.factor(data$.causal_env)
    if (nlevels(out) >= 2L) return(out)
  }
  if (!is.null(seed)) set.seed(seed)
  climate <- grep("^bio|^prec|^sunp|^tave|^tmax|^tmin", vars,
                  value = TRUE, ignore.case = TRUE)
  if (length(climate) < 3L) {
    sds <- vapply(vars, function(v) stats::sd(data[[v]], na.rm = TRUE), numeric(1))
    climate <- vars[order(sds, decreasing = TRUE)][seq_len(min(5L, length(vars)))]
  }
  x <- as.matrix(.prior_numeric_frame(data, climate))
  pca <- stats::prcomp(scale(x), rank. = min(3L, ncol(x)))
  pcs <- pca$x[, seq_len(min(3L, ncol(pca$x))), drop = FALSE]
  stats::kmeans(pcs, centers = min(n_blocks, nrow(pcs)), nstart = 5,
                iter.max = 50)$cluster |> as.factor()
}

.prior_invariance_scores <- function(data, pool, response, seed = NULL) {
  y <- data[[response]]
  blocks <- .prior_blocks(data, pool, seed = seed)
  rows <- lapply(pool, function(v) {
    beta <- se <- pval <- numeric()
    for (e in levels(blocks)) {
      idx <- which(blocks == e)
      x <- as.numeric(data[[v]][idx]); yy <- as.integer(y[idx])
      ok <- is.finite(x) & !is.na(yy)
      if (sum(ok) < 20L || length(unique(yy[ok])) < 2L || stats::sd(x[ok]) == 0) next
      fit <- tryCatch(suppressWarnings(stats::glm(yy[ok] ~ scale(x[ok]),
                                                   family = stats::binomial())),
                      error = function(e) NULL)
      if (is.null(fit) || length(stats::coef(fit)) < 2L) next
      b <- unname(stats::coef(fit)[2])
      s <- tryCatch(sqrt(stats::vcov(fit)[2, 2]), error = function(e) NA_real_)
      if (!is.finite(b) || !is.finite(s) || s <= 0) next
      beta <- c(beta, b); se <- c(se, s)
      pval <- c(pval, 2 * stats::pnorm(abs(b / s), lower.tail = FALSE))
    }
    if (length(beta) < 3L) {
      return(data.frame(variable = v, p_invar = 0, min_p_sig = 1,
                        sign_consistency = 0, mean_abs_coef = 0))
    }
    bbar <- stats::weighted.mean(beta, 1 / se^2)
    pinv <- stats::pchisq(sum(((beta - bbar) / se)^2), length(beta) - 1L,
                             lower.tail = FALSE)
    signs <- sign(beta); signs <- signs[signs != 0]
    data.frame(variable = v, p_invar = pinv, min_p_sig = min(pval),
               sign_consistency = if (length(signs)) max(mean(signs > 0), mean(signs < 0)) else 0,
               mean_abs_coef = mean(abs(beta)))
  })
  out <- do.call(rbind, rows)
  out$ics_score <- .prior_normalize(out$p_invar) * 0.40 +
    .prior_normalize(-out$min_p_sig) * 0.30 +
    .prior_normalize(out$sign_consistency) * 0.20 +
    .prior_normalize(out$mean_abs_coef) * 0.10
  out
}

.prior_partial_p <- function(y, x, z = NULL) {
  ok <- is.finite(y) & is.finite(x)
  if (!is.null(z) && ncol(as.data.frame(z))) {
    z <- as.data.frame(z)
    for (nm in names(z)) ok <- ok & is.finite(z[[nm]])
  }
  y <- y[ok]; x <- x[ok]
  if (length(y) < 10L || stats::sd(y) == 0 || stats::sd(x) == 0) return(1)
  if (is.null(z) || !ncol(as.data.frame(z))) {
    r <- stats::cor(y, x); df <- length(y) - 2L
  } else {
    z <- as.data.frame(z)[ok, , drop = FALSE]
    r <- stats::cor(stats::residuals(stats::lm(y ~ ., data = z)),
                    stats::residuals(stats::lm(x ~ ., data = z)))
    df <- length(y) - ncol(z) - 2L
  }
  if (!is.finite(r) || df <= 0) return(1)
  r <- max(min(r, 0.999999), -0.999999)
  2 * stats::pt(abs(r) * sqrt(df / max(1 - r^2, 1e-12)), df, lower.tail = FALSE)
}

.prior_conditional_scores <- function(data, pool, response, target_n, min_vars) {
  y <- as.numeric(data[[response]])
  selected <- character(); remaining <- pool; rows <- list()
  for (step in seq_len(min(target_n, length(pool)))) {
    scores <- do.call(rbind, lapply(remaining, function(v) {
      z <- if (length(selected)) .prior_numeric_frame(data, selected) else NULL
      p <- tryCatch(.prior_partial_p(y, as.numeric(data[[v]]), z),
                    error = function(e) 1)
      data.frame(variable = v, ci_pvalue = p,
                 ci_score = -log10(max(p, 1e-300)))
    }))
    best <- scores[which.max(scores$ci_score), , drop = FALSE]
    rows[[length(rows) + 1L]] <- best
    if (length(selected) >= min_vars && best$ci_pvalue > 0.05) break
    selected <- c(selected, best$variable)
    remaining <- setdiff(remaining, best$variable)
    if (!length(remaining)) break
  }
  out <- .prior_univariate_scores(data, pool, response)
  out$ci_pvalue <- NA_real_; out$ci_score <- 0
  if (length(rows)) {
    picked <- do.call(rbind, rows); idx <- match(picked$variable, out$variable)
    out$ci_pvalue[idx] <- picked$ci_pvalue; out$ci_score[idx] <- picked$ci_score
  }
  out
}

.select_causal_prior_rf <- function(data, env_vars, response, causal_spec,
                                    min_vars, max_vars, cor_threshold,
                                    num_trees, seed) {
  started <- proc.time()
  spec <- .validate_cast_causal_spec(causal_spec, env_vars)
  eligible <- spec$variable[spec$role == "direct_candidate"]
  if (!length(eligible)) cli::cli_abort("No direct_candidate variables are available.")
  target <- max(min_vars, min(max_vars, ceiling(log2(max(sum(data[[response]] == 1), 2L)))))
  target <- min(target, length(eligible))
  pool <- .prior_signal_pool(data, eligible, response, pool_max = 35L)

  if (!is.null(seed)) set.seed(seed)
  x <- .prior_numeric_frame(data, pool)
  rf <- ranger::ranger(y ~ ., data = cbind(y = as.factor(data[[response]]), x),
                       num.trees = num_trees, importance = "impurity",
                       verbose = FALSE)
  rf_score <- data.frame(variable = pool,
                         global_rf_importance = pmax(rf$variable.importance[pool], 0))
  inv <- .prior_invariance_scores(data, pool, response, seed)
  cond <- .prior_conditional_scores(data, pool, response, target,
                                    min(min_vars, target))
  scores <- Reduce(function(a, b) merge(a, b, by = "variable", all = TRUE),
                   list(rf_score, inv, cond))
  for (nm in c("global_rf_importance", "ics_score", "ci_score", "univ_score",
               "sign_consistency")) {
    if (!nm %in% names(scores)) scores[[nm]] <- 0
    scores[[nm]][!is.finite(scores[[nm]]) | is.na(scores[[nm]])] <- 0
  }
  scores$prior_rf_score <- .prior_normalize(scores$global_rf_importance) * 0.45 +
    .prior_normalize(scores$ics_score) * 0.20 +
    .prior_normalize(scores$ci_score) * 0.20 +
    .prior_normalize(scores$univ_score) * 0.10 +
    .prior_normalize(scores$sign_consistency) * 0.05
  ranked <- scores$variable[order(-scores$prior_rf_score)]
  selected <- .prior_select_diverse(ranked, data, target, cor_threshold, 2L)

  audit <- spec
  audit$specified_role <- audit$role
  audit$specified_confidence <- audit$confidence
  audit$prior_eligible <- audit$specified_role == "direct_candidate"
  audit$in_fast_pool <- audit$variable %in% pool
  for (nm in setdiff(names(scores), "variable")) {
    audit[[nm]] <- scores[[nm]][match(audit$variable, scores$variable)]
  }
  audit$selected <- audit$variable %in% selected
  audit$selection_tier <- ifelse(audit$selected, "causal_core",
                                  ifelse(audit$prior_eligible, "eligible_not_selected", "excluded"))
  audit$selection_reason <- ifelse(
    audit$selected, "selected_from_prior_eligible_pool",
    ifelse(!audit$prior_eligible, paste0("excluded_by_prior_role:", audit$specified_role),
           ifelse(audit$in_fast_pool, "ranked_below_core_threshold",
                  "eligible_but_not_in_fast_pool"))
  )
  audit$role <- ifelse(audit$selected, "prior_causal_core",
                       ifelse(audit$prior_eligible, "prior_eligible_unselected",
                              paste0("prior_", audit$specified_role)))
  audit$screening_method <- "causal_prior_rf"
  audit$selection_time <- unname((proc.time() - started)["elapsed"])
  roles <- audit[audit$selected,
                 c("variable", "role", "specified_role", "specified_confidence",
                   "selection_tier", "screening_method"), drop = FALSE]
  list(selected = selected, scores = audit, roles = roles,
       specification = spec,
       time = unname((proc.time() - started)["elapsed"]))
}
