#' Sensitivity Analysis via E-values
#'
#' Computes E-values for each ATE estimate, quantifying the minimum strength
#' of unmeasured confounding (on the risk ratio scale) needed to explain away
#' the observed causal effect. Larger E-values indicate more robust findings.
#'
#' @param ate A [cast_ate] object.
#' @param transform Character. How to convert ATE coefficients (risk
#'   differences on a 0-1 outcome) to approximate risk ratios.
#'   \describe{
#'     \item{`"RR"` (default)}{Use \eqn{RR \approx 1 + |ATE| / p_0} where
#'       \eqn{p_0 = 0.5} is assumed baseline prevalence. Appropriate when ATE
#'       represents a difference in proportions.}
#'     \item{`"OR"`}{Convert via \eqn{OR = (p_0 + ATE) / p_0 \cdot
#'       (1-p_0) / (1-p_0-ATE)}, then approximate RR from OR using the
#'       square-root bound: \eqn{RR \approx \sqrt{OR}}.}
#'   }
#' @param p0 Numeric. Assumed baseline outcome prevalence for RR/OR
#'   approximation. Default `0.5`.
#' @param verbose Logical. Print summary. Default `TRUE`.
#'
#' @return A `cast_ate` object with the `estimates` data.frame augmented by
#'   columns `rr` (approximate risk ratio), `evalue_point` (E-value for the
#'   point estimate), and `evalue_ci` (E-value for the confidence interval
#'   bound closest to the null).
#'
#' @details
#' The E-value for a risk ratio \eqn{RR \geq 1} is:
#' \deqn{E = RR + \sqrt{RR \times (RR - 1)}}
#' For the confidence interval bound, the same formula is applied to the
#' lower CI bound's RR; if the CI includes the null (RR = 1), the E-value
#' is 1.
#'
#' @references
#' VanderWeele, T.J. and Ding, P. (2017). Sensitivity analysis in
#' observational research: Introducing the E-value.
#' *Annals of Internal Medicine*, 167(4), 268-274.
#'
#' @seealso [cast_ate()]
#'
#' @export
cast_evalue <- function(ate, transform = c("RR", "OR"), p0 = 0.5,
                        verbose = TRUE) {
  if (!inherits(ate, "cast_ate")) {
    cli::cli_abort("{.arg ate} must be a {.cls cast_ate} object.")
  }
  transform <- match.arg(transform)

  est <- ate$estimates
  n <- nrow(est)

  rr_vec <- numeric(n)
  ev_point <- numeric(n)
  ev_ci <- numeric(n)

  for (i in seq_len(n)) {
    coef_i <- est$coef[i]
    se_i <- est$se[i]

    if (is.na(coef_i) || is.na(se_i)) {
      rr_vec[i] <- NA_real_
      ev_point[i] <- NA_real_
      ev_ci[i] <- NA_real_
      next
    }

    # Convert ATE to approximate risk ratio
    rr <- ate_to_rr(abs(coef_i), p0, transform)
    rr_vec[i] <- rr

    # E-value for point estimate
    ev_point[i] <- evalue_rr(rr)

    # E-value for CI bound closest to null
    ci_lower_ate <- abs(coef_i) - stats::qnorm(1 - ate$alpha / 2) * se_i
    if (ci_lower_ate <= 0) {
      ev_ci[i] <- 1
    } else {
      rr_ci <- ate_to_rr(ci_lower_ate, p0, transform)
      ev_ci[i] <- evalue_rr(rr_ci)
    }
  }

  ate$estimates$rr <- rr_vec
  ate$estimates$evalue_point <- ev_point
  ate$estimates$evalue_ci <- ev_ci

  if (verbose) {
    n_robust <- sum(ev_ci > 1, na.rm = TRUE)
    cli::cli_inform(c(
      "E-value sensitivity analysis: {n} variables",
      "i" = "{n_robust}/{n} have E-value (CI) > 1 (robust to confounding)"
    ))
  }

  ate
}


#' Convert ATE (Risk Difference) to Approximate Risk Ratio
#' @keywords internal
#' @noRd
ate_to_rr <- function(ate_abs, p0, transform) {
  if (transform == "RR") {
    rr <- 1 + ate_abs / max(p0, 0.01)
  } else {
    # OR transform
    p1 <- min(p0 + ate_abs, 0.999)
    or <- (p1 / (1 - p1)) / (p0 / (1 - p0))
    rr <- sqrt(or)  # square-root bound approximation
  }
  max(rr, 1)
}


#' Compute E-value from Risk Ratio
#' @keywords internal
#' @noRd
evalue_rr <- function(rr) {
  if (is.na(rr) || rr < 1) return(1)
  rr + sqrt(rr * (rr - 1))
}
