# Vendored covsel core -------------------------------------------------------
#
# Source: covsel v1.1.1 (https://github.com/antadde/covsel, previously
# N-SDM/covsel), the reference implementation of the two-step embedded
# covariate selection procedure of
#
#   Adde, A., Rey, P.-L., Fopp, F., Petitpierre, B., Schweiger, A. K.,
#   Broennimann, O., Lehmann, A., Zimmermann, N. E., Altermatt, F.,
#   Pellissier, L., Guisan, A. (2023). Too many candidates: Embedded
#   covariate selection procedure for species distribution modelling with
#   the covsel R package. Ecological Informatics, 75, 102080.
#   \doi{10.1016/j.ecoinf.2023.102080}
#
# Copyright (c) 2023 Antoine Adde. MIT License (see inst/COPYRIGHTS).
#
# Vendored so the published procedure is available in-package: covsel is
# GitHub-only, its NAMESPACE declares no imports (Depends-only), and its
# glm branch is broken under data.table >= 1.18. Local modifications, all
# behaviour-preserving except where noted:
#
#   1. Bare calls to Depends packages are namespaced (glmnet::, RRF::,
#      mgcv::). covsel 1.1.1 lists those packages in Depends with an empty
#      NAMESPACE, so `covsel::` calls fail unless the caller's search path
#      happens to carry all of them.
#   2. The `data.table::setDT(...)[, .SD[which.max(coef)], by = covariate]`
#      de-duplication is replaced with the equivalent base-R
#      `!duplicated(covariate)` after ordering. Under data.table >= 1.18
#      the `coef` symbol in that expression resolves to `base::coef` and
#      the call errors.
#   3. RNG side effects removed (`RNGkind()`/`set.seed()` at entry). The
#      caller seeds locally with `.cast_with_seed()`; the algorithm itself
#      is otherwise unchanged.
#   4. `print()` progress replaced with `cli::cli_inform()` (silenced by
#      `verbose = FALSE`).
#   5. Group-level `force` handling intersects the forced names with each
#      split group's columns first; the official code indexes
#      `covdata[, force]` inside every group, which errors whenever a
#      forced covariate lives in a different group.
#   6. Intercept removal in the elastic-net branch is name-based;
#      the official `[-1, ]` dropped the first covariate whenever the
#      intercept coefficient was exactly zero.
#   7. The union-rank sequence in the aggregation starts at 1 when no
#      covariate was selected by every algorithm (the official
#      `max(integer(0) + 1)` produced a corrupt sequence on that path).
#   8. Elastic-net label extraction uses the exact inverse of the
#      `poly(<name>,2)<n>` construction. The official
#      `stringi::stri_sub(x, 6, -6)` drops the last character of every
#      covariate name, silently corrupting the cross-algorithm label
#      matching (and pulling in a stringi dependency).

#' Collinearity Filtering of a Candidate Covariate Set (vendored covsel)
#'
#' Step A of the covsel procedure: order covariates by univariate
#' quadratic-logistic p-value, then greedily keep the head of that order and
#' drop every covariate correlated with it beyond `corcut`, repeating on the
#' remainder. Covariates with fewer than 10 unique values are discarded
#' (forced covariates excepted). See the file header for the modifications
#' to the vendored source.
#'
#' @param covdata Data frame of candidate covariates.
#' @param pa Numeric 0/1 response vector.
#' @param weights Optional row weights.
#' @param force Optional character vector of covariates forced into the set
#'   (collinear partners are removed).
#' @param corcut Absolute Pearson correlation threshold.
#' @param verbose Print progress. Default `FALSE`.
#' @return The filtered covariate data frame.
#' @keywords internal
#' @noRd
.covsel_filteralgo <- function(covdata, pa, weights = NULL, force = NULL,
                               corcut = 0.7, verbose = FALSE) {
  candidates <- names(covdata)
  if (ncol(covdata) == 1L) {
    return(covdata)
  }
  pointless10 <- which(apply(covdata, 2, function(x) length(unique(x))) < 10)
  if (length(pointless10) > 0) {
    pointless <- pointless10[!names(pointless10) %in% force]
    if (length(pointless) > 0) {
      if (verbose) {
        cli::cli_inform("covsel Step A: covariate{?s} {.val {names(covdata)[pointless]}} have fewer than 10 unique values and are discarded.")
      }
      covdata <- covdata[, -pointless, drop = FALSE]
    }
    pointless_f <- pointless10[names(pointless10) %in% force]
    if (length(pointless_f) > 0 && verbose) {
      cli::cli_inform("covsel Step A: forced covariate{?s} {.val {names(covdata)[pointless_f]}} have fewer than 10 unique values.")
    }
  }
  force_dat <- data.frame()
  if (is.character(force) && any(force %in% names(covdata))) {
    force_in <- intersect(force, names(covdata))
    force_dat <- data.frame(covdata[, force_in, drop = FALSE])
    names(force_dat) <- force_in
    if (length(force_in) > 1) {
      ix <- unique(rownames(which(abs(cor(covdata))[, force_in] > corcut,
                                  arr.ind = TRUE)))
    } else {
      ix <- names(which(abs(cor(covdata))[, force_in] > corcut))
    }
    covdata <- covdata[, !(colnames(covdata) %in% ix), drop = FALSE]
  }
  covdata <- cbind(covdata, pa = as.factor(pa))
  res <- sort(apply(covdata[, setdiff(names(covdata), "pa"), drop = FALSE],
                    2, function(x) {
    min(summary(suppressWarnings(stats::glm(
      covdata$pa ~ stats::poly(x, degree = 2),
      family = "binomial", weights = weights
    )))$coefficients[2:3, 4])
  }))
  covranked <- data.frame(pval = res, row.names = names(res))
  covdata.ranked <- data.frame(covdata[, row.names(covranked), drop = FALSE],
                               check.names = FALSE)
  cor.mat <- abs(stats::cor(covdata.ranked, use = "pairwise.complete.obs"))
  sanctuaire <- NULL
  if (all(cor.mat[cor.mat != 1] < corcut)) {
    sanctuaire <- c(sanctuaire, row.names(cor.mat))
  } else {
    while (length(cor.mat) > 1) {
      tmp_cm <- data.frame(cor.mat[, 1])
      ix <- which(tmp_cm[, 1] <= corcut)
      var.name <- row.names(cor.mat)[1]
      sanctuaire <- c(sanctuaire, var.name)
      cor.mat <- cor.mat[ix, ix, drop = FALSE]
      if (all(cor.mat[cor.mat != 1] < corcut)) {
        sanctuaire <- c(sanctuaire, row.names(cor.mat))
        break
      }
    }
  }
  covdata.filter <- as.data.frame(covdata[, sanctuaire, drop = FALSE])
  if (ncol(force_dat) > 0) {
    covdata.filter <- cbind(covdata.filter, force_dat)
  }
  if (length(sanctuaire) == 1) colnames(covdata.filter) <- sanctuaire
  covdata.filter
}


#' Hierarchical Collinearity Filtering (vendored covsel)
#'
#' Step A applied sequentially at up to three target levels: one covariate
#' per variable-level identity (`corcut = 0` within each variable), then per
#' category, then on the remaining set. See the file header for the
#' modifications to the vendored source.
#'
#' @param covdata Data frame of candidate covariates.
#' @param pa Numeric 0/1 response vector.
#' @param weights Optional row weights.
#' @param force Optional character vector of forced covariates.
#' @param corcut Absolute Pearson correlation threshold for the category and
#'   global levels.
#' @param variables Optional character vector aligned with the columns of
#'   `covdata` giving each covariate's variable-level identity.
#' @param categories Optional character vector aligned with the columns of
#'   `covdata` giving each covariate's category.
#' @param verbose Print progress. Default `FALSE`.
#' @return The filtered covariate data frame.
#' @keywords internal
#' @noRd
.covsel_filter <- function(covdata, pa, weights = NULL, force = NULL,
                           corcut = 0.7, variables = NULL, categories = NULL,
                           verbose = FALSE) {
  strip_prefix <- function(x) sub("^.*\\.", "", x)
  force_in_group <- function(grp, force) {
    names(grp)[strip_prefix(names(grp)) %in% force]
  }
  if (length(variables) > 0) {
    covdata.candidates <- split.default(covdata, variables)
    covdata.variable.filter <- lapply(covdata.candidates, function(grp) {
      .covsel_filteralgo(grp, pa, weights = weights,
                         force = force_in_group(grp, force),
                         corcut = 0, verbose = verbose)
    })
    covdata <- do.call("cbind", covdata.variable.filter)
    names(covdata) <- strip_prefix(names(covdata))
    if (length(categories) > 0) {
      categories <- categories[match(names(covdata.variable.filter), variables)]
    }
  }
  if (length(categories) > 0) {
    covdata.candidates <- split.default(covdata, categories)
    covdata.category.filter <- lapply(covdata.candidates, function(grp) {
      .covsel_filteralgo(grp, pa, weights = weights,
                         force = force_in_group(grp, force),
                         corcut = corcut, verbose = verbose)
    })
    covdata <- do.call("cbind", covdata.category.filter)
    names(covdata) <- strip_prefix(names(covdata))
  }
  covdata.filter <- .covsel_filteralgo(covdata, pa, weights = weights,
                                       force = force, corcut = corcut,
                                       verbose = verbose)
  names(covdata.filter) <- strip_prefix(names(covdata.filter))
  covdata.filter
}


#' Model-Specific Embedded Selection (vendored covsel)
#'
#' Step B of the covsel procedure: elastic-net GLM, null-space-penalized GAM
#' and guided regularized random forest each rank the filtered covariates;
#' the aggregate ranking puts covariates selected by every algorithm first,
#' then the remaining union, and keeps the top `ncov`. See the file header
#' for the modifications to the vendored source.
#'
#' @param covdata Data frame of filtered covariates.
#' @param pa Numeric 0/1 response vector.
#' @param weights Optional row weights.
#' @param force Optional character vector of forced covariates.
#' @param algorithms Embedding learners: any of `"glm"`, `"gam"`, `"rf"`.
#' @param ncov Target covariate count.
#' @param maxncov Maximum covariate count.
#' @param nthreads Threads for the GAM embedding.
#' @param verbose Print progress. Default `FALSE`.
#' @return `NULL` when no algorithm selects anything, else a list with
#'   `covdata` (selected covariates), `ranks_1` (per-algorithm ranks),
#'   `ranks_2` (aggregate ranking) and `selected` (selected column names).
#' @keywords internal
#' @noRd
.covsel_embed <- function(covdata, pa, weights = NULL, force = NULL,
                          algorithms = c("glm", "gam", "rf"),
                          ncov = ceiling(log2(length(which(pa == 1)))),
                          maxncov = 12, nthreads = 1L, verbose = FALSE) {
  ranks_1 <- data.frame()
  if (!is.numeric(weights)) weights <- rep(1, length(pa))
  if ("glm" %in% algorithms) {
    form <- stats::as.formula(paste0(
      "as.factor(pa) ~ ",
      paste(paste0("poly(", names(covdata), ",2)"), collapse = " + "), "-1"
    ))
    x <- stats::model.matrix(form, covdata)
    mdl.glm <- suppressWarnings(glmnet::cv.glmnet(
      x, as.factor(pa), alpha = 0.5, weights = weights,
      family = "binomial", type.measure = "deviance"
    ))
    beta_at <- function(s) {
      b <- as.data.frame(as.matrix(stats::coef(mdl.glm, s = s)))
      b <- data.frame(covariate = row.names(b),
                      coef = as.numeric(abs(b[, 1])),
                      stringsAsFactors = FALSE)
      b[b$coef > 0 & b$covariate != "(Intercept)", , drop = FALSE]
    }
    glm.beta <- beta_at(mdl.glm$lambda.1se)
    if (nrow(glm.beta) < 1) glm.beta <- beta_at(mdl.glm$lambda.min)
    if (nrow(glm.beta) < 1) {
      if (verbose) {
        cli::cli_inform("covsel embed: no covariate selected after elastic-net regularization.")
      }
    } else {
      glm.beta <- data.frame(
        glm.beta[order(glm.beta$coef, decreasing = TRUE), , drop = FALSE],
        model = "glm"
      )
      # Robust label extraction. The official `stringi::stri_sub(x, 6, -6)`
      # truncates covariate names (the exact cut depends on the spacing
      # `model.matrix()` inserts), which corrupts the cross-algorithm label
      # matching; the regex reverses exactly the "poly(<name>,2)<n>"
      # construction, with or without the inserted space.
      glm.beta$covariate <- sub("^poly\\((.*),\\s*2\\)[0-9]+$", "\\1",
                                glm.beta$covariate)
      glm.beta <- glm.beta[!duplicated(glm.beta$covariate), , drop = FALSE]
      glm.beta$rank <- seq_len(nrow(glm.beta))
      ranks_1 <- rbind(ranks_1, glm.beta[, c("covariate", "rank", "model")])
    }
  }
  if ("gam" %in% algorithms) {
    force_gam <- intersect(force %||% character(0), names(covdata))
    pointless10 <- integer(0)
    if (length(force_gam)) {
      df_force <- data.frame(covdata[, force_gam, drop = FALSE])
      names(df_force) <- force_gam
      pointless10 <- which(apply(df_force, 2,
                                 function(x) length(unique(x))) < 10)
    }
    smooth_names <- names(covdata)[!names(covdata) %in% names(pointless10)]
    build_form <- function(bs) {
      stats::as.formula(paste0(
        "pa ~ ", paste(paste0("s(", smooth_names, ",bs='", bs, "')"),
                       collapse = " + ")
      ))
    }
    mdl.gam <- suppressWarnings(mgcv::bam(
      build_form("cr"), data = cbind(covdata, as.factor(pa)),
      weights = weights, family = "binomial", method = "fREML",
      select = TRUE, discrete = TRUE,
      control = list(nthreads = nthreads)
    ))
    t <- try(summary(mdl.gam), TRUE)
    if (inherits(t, "try-error")) {
      mdl.gam <- suppressWarnings(mgcv::bam(
        build_form("ts"), data = cbind(covdata, as.factor(pa)),
        weights = weights, family = "binomial", method = "fREML",
        select = TRUE, discrete = TRUE,
        control = list(nthreads = nthreads)
      ))
    }
    gam.beta <- data.frame(
      covariate = names(mdl.gam$model)[!names(mdl.gam$model) %in%
        c("(weights)", "pa")],
      summary(mdl.gam)$s.table,
      row.names = NULL
    )
    gam.beta <- gam.beta[gam.beta$p.value < 0.9, , drop = FALSE]
    if (nrow(gam.beta) < 1) {
      if (verbose) {
        cli::cli_inform("covsel embed: no covariate selected after GAM null-space penalization.")
      }
    } else {
      gam.beta <- data.frame(
        gam.beta[order(abs(gam.beta$Chi.sq), decreasing = TRUE), , drop = FALSE],
        rank = seq_len(nrow(gam.beta)),
        model = "gam"
      )
      ranks_1 <- rbind(ranks_1, gam.beta[, c("covariate", "rank", "model")])
    }
  }
  if ("rf" %in% algorithms) {
    rf <- RRF::RRF(covdata, as.factor(pa), flagReg = 0)
    impRF <- rf$importance[, "MeanDecreaseGini"]
    imp <- impRF / max(impRF)
    gamma <- 0.5
    coefReg <- (1 - gamma) + gamma * imp
    mdl.rf <- RRF::RRF(
      covdata, as.factor(pa),
      classwt = c(`0` = min(weights), `1` = max(weights)),
      coefReg = coefReg, flagReg = 1
    )
    rf.beta <- data.frame(covariate = row.names(mdl.rf$importance),
                          mdl.rf$importance, row.names = NULL)
    rf.beta <- rf.beta[which(rf.beta$MeanDecreaseGini > 0), , drop = FALSE]
    if (nrow(rf.beta) < 1) {
      if (verbose) {
        cli::cli_inform("covsel embed: no covariate selected after guided regularized random forest.")
      }
    } else {
      rf.beta <- data.frame(
        rf.beta[order(rf.beta$MeanDecreaseGini, decreasing = TRUE), ,
                drop = FALSE],
        rank = seq_len(nrow(rf.beta)),
        model = "rf"
      )
      ranks_1 <- rbind(ranks_1, rf.beta[, c("covariate", "rank", "model")])
    }
  }
  if (nrow(ranks_1) < 1) {
    if (verbose) {
      cli::cli_inform("covsel embed: no covariate selected after the embedding procedure.")
    }
    return(NULL)
  }
  n_models <- length(unique(ranks_1$model))
  in_all <- names(which(table(ranks_1$covariate) == n_models))
  ranks_2 <- NULL
  if (length(in_all)) {
    inter <- ranks_1[ranks_1$covariate %in% in_all, , drop = FALSE]
    inter <- stats::aggregate(inter[, "rank", drop = FALSE],
                              list(inter$covariate), sum)
    colnames(inter) <- c("covariate", "rank")
    inter.sel <- data.frame(
      inter[order(inter$rank, decreasing = FALSE), , drop = FALSE],
      rank.f = seq_len(nrow(inter))
    )
    ranks_2 <- inter.sel
  }
  rest <- ranks_1[!ranks_1$covariate %in% in_all, , drop = FALSE]
  if (nrow(rest) > 0) {
    rest <- stats::aggregate(rest[, "rank", drop = FALSE],
                             list(rest$covariate), sum)
    colnames(rest) <- c("covariate", "rank")
    start <- if (!is.null(ranks_2)) max(ranks_2$rank.f) + 1L else 1L
    union.sel <- data.frame(
      rest[order(rest$rank, decreasing = FALSE), , drop = FALSE],
      rank.f = seq.int(start, length.out = nrow(rest))
    )
    ranks_2 <- if (is.null(ranks_2)) union.sel else rbind(ranks_2, union.sel)
  }
  ncov <- min(ncov, maxncov)
  ncov <- min(ncov, nrow(ranks_2))
  ranks_2 <- ranks_2[seq_len(ncov), , drop = FALSE]
  if (is.character(force)) {
    tf <- force[which(!(force %in% ranks_2$covariate))]
    # Official condition `length(tf > 1)`, i.e. length(tf) > 1: forced
    # covariates are spliced into the tail slots of the ncov set. A single
    # missing forced covariate is not spliced here; the caller's
    # union(keep, selected) restores it.
    if (length(tf) > 1L) {
      toforce <- data.frame(covariate = tf, rank = "forced", rank.f = "forced")
      ranks_2[c(nrow(ranks_2) - nrow(toforce) + 1):c(nrow(ranks_2)), ] <- toforce
    }
  }
  ranks_2 <- ranks_2[, c("covariate", "rank.f")]
  covdata <- covdata[, unlist(ranks_2["covariate"]), drop = FALSE]
  list(covdata = covdata, ranks_1 = ranks_1, ranks_2 = ranks_2,
       selected = names(covdata))
}
