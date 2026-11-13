# AI-assisted development (2026-11-13); see inst/AI_ASSISTANCE.md.
# Compare actual scenario responses on the same supported locations.

.cast_response_direction <- function(only, retained) {
  a <- sign(only); b <- sign(retained)
  out <- rep(NA_integer_, length(a))
  ok <- is.finite(a) & is.finite(b)
  out[ok & a == 0 & b == 0] <- 0L
  out[ok & a == b & a > 0] <- 1L
  out[ok & a == b & a < 0] <- -1L
  out[ok & a * b < 0] <- 2L
  out[ok & xor(a == 0, b == 0)] <- 3L
  out
}

.cast_response_evidence <- function(data, drivers, current, future) {
  only <- as.matrix(data[paste0("effect_", drivers)])
  retained <- as.matrix(data[paste0("retained_", drivers)])
  support <- as.matrix(data[c("support_total", paste0("support_effect_", drivers),
    paste0("support_retained_", drivers))])
  finite <- is.finite(data$total_change) & rowSums(!is.finite(only) | !is.finite(retained)) == 0L
  common <- finite & rowSums(is.na(support) | support != 1) == 0L
  changed <- as.matrix(current[drivers]) != as.matrix(future[drivers])
  common <- common & rowSums(is.na(changed)) == 0L
  n <- sum(common)
  mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_
  summaries <- lapply(seq_along(drivers), function(j) {
    v <- drivers[j]; a <- only[common, j]; b <- retained[common, j]
    varied <- changed[, j]
    direction <- .cast_response_direction(a, b)
    agreement <- function(prefix, changed_only = FALSE) {
      field <- paste0(prefix, v)
      if (!field %in% names(data)) return(NA_real_)
      z <- data[[field]][common]
      if (changed_only) z <- z[varied[common]]
      if (any(!is.finite(z))) return(NA_real_)
      mean_or_na(z == 1)
    }
    data.frame(driver = v, n_common = n,
      n_changed = sum(varied, na.rm = TRUE),
      n_changed_common = sum(varied[common]),
      status = if (all(is.na(varied))) "missing_predictor_pairs" else
        if (!any(varied, na.rm = TRUE) && anyNA(varied)) "unchanged_known_pairs" else
        if (!any(varied, na.rm = TRUE)) "unchanged_scenario" else
        if (!n) "no_common_support" else if (!any(varied[common])) "unchanged_on_common_domain" else
          if (all(a == 0 & b == 0)) "zero_model_response" else "changed",
      mean_only = mean_or_na(a), mean_retained = mean_or_na(b),
      mean_abs_only = mean_or_na(abs(a)), mean_abs_retained = mean_or_na(abs(b)),
      mean_abs_context_gap = mean_or_na(abs(b - a)),
      context_reversal_fraction = mean_or_na(direction == 2L),
      one_background_zero_fraction = mean_or_na(direction == 3L),
      all_models_agree_only = agreement("effect_agreement_"),
      all_models_agree_retained = agreement("retained_agreement_"),
      all_models_agree_changed_only = agreement("effect_agreement_", TRUE),
      all_models_agree_changed_retained = agreement("retained_agreement_", TRUE))
  })
  summary <- do.call(rbind, summaries)
  summary$rank <- NA_integer_
  eligible <- summary$status == "changed" & is.finite(summary$mean_abs_only)
  summary$rank[eligible] <- as.integer(rank(-summary$mean_abs_only[eligible], ties.method = "min"))
  summary <- summary[order(summary$rank, na.last = TRUE), , drop = FALSE]
  row.names(summary) <- NULL
  winner <- function(x) {
    out <- rep(NA_integer_, nrow(x))
    if (n) {
      z <- abs(x[common, , drop = FALSE]); biggest <- apply(z, 1L, max)
      ties <- rowSums(z == biggest) > 1L
      id <- max.col(z, ties.method = "first")
      id[ties] <- -1L; id[biggest == 0] <- 0L
      out[common] <- id
    }
    out
  }
  first <- winner(only); last <- winner(retained)
  same <- rep(NA_integer_, nrow(data))
  resolved <- !is.na(first) & !is.na(last) & first > 0L & last > 0L
  same[resolved] <- as.integer(first[resolved] == last[resolved])
  reversals <- rowSums(sign(only) * sign(retained) < 0)
  reversals[!common] <- NA_integer_
  list(summary = summary,
    locations = data.frame(common_support = as.integer(common),
      largest_only = first, largest_retained = last, same_largest = same, context_reversals = reversals),
    drivers = data.frame(id = seq_along(drivers), driver = drivers),
    population = data.frame(n_input = nrow(data), n_finite = sum(finite), n_common = n,
      common_fraction = n / nrow(data)),
    interpretation = paste("Shared range-supported saved locations; rankings concern the supplied scenario changes.",
      "Changed-location agreement excludes predictors unchanged at those saved locations; unqualified agreement includes all common rows.",
      "Unchanged-scenario status refers to supplied locations, not the entire raster. Background dependence and model sign agreement are descriptive, not causal identification or confidence.",
      "Largest-response IDs: -1 exact tie, 0 no response, NA outside common support."))
}
