#' VIF-Based Collinearity Screening
#'
#' Iteratively removes environmental variables with Variance Inflation Factor
#' (VIF) above the specified threshold. Optionally applies expert-guided
#' pre-filtering before VIF elimination.
#'
#' @param data A `data.frame` of environmental variables (numeric columns
#'   only). Rows are sites, columns are variables.
#' @param threshold Numeric. VIF threshold above which variables are removed.
#'   Default `10` (Zuur et al. 2010).
#' @param exclude Character vector of column names to exclude from screening
#'   (e.g., coordinates, presence). Default `NULL`.
#' @param expert_filter Character vector of variable names to remove before
#'   VIF screening (domain-knowledge redundancies). Default `NULL`.
#' @param verbose Logical. Print iteration details. Default `TRUE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`selected`}{Character vector of retained variable names.}
#'   \item{`removed`}{Character vector of removed variable names (in order).}
#'   \item{`vif_log`}{A `data.frame` with iteration-by-iteration VIF values.}
#'   \item{`data`}{Filtered `data.frame` with only retained variables.}
#' }
#'
#' @references
#' Zuur, A.F., Ieno, E.N., & Elphick, C.S. (2010). A protocol for data
#' exploration to avoid common statistical problems. *Methods in Ecology and
#' Evolution*, 1(1), 3-14.
#'
#' @export
cast_vif <- function(data,
                     threshold = 10,
                     exclude = NULL,
                     expert_filter = NULL,
                     verbose = TRUE) {
  check_suggested("car", "for VIF computation")
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data.frame.")
  }

  # Identify numeric columns only
  all_cols <- names(data)
  exclude <- exclude %||% character(0)
  num_cols <- setdiff(all_cols, exclude)
  num_cols <- num_cols[vapply(data[num_cols], is.numeric, logical(1))]

  # Stage 1: Expert pre-filter
  removed <- character(0)
  if (!is.null(expert_filter)) {
    hit <- intersect(expert_filter, num_cols)
    if (length(hit) > 0 && verbose) {
      cli::cli_inform(
        "Expert pre-filter removes {length(hit)} variable{?s}."
      )
    }
    removed <- hit
    num_cols <- setdiff(num_cols, hit)
  }

  # Complete cases for lm
  dat_clean <- stats::na.omit(data[, num_cols, drop = FALSE])
  if (nrow(dat_clean) < 10) {
    cli::cli_abort("Fewer than 10 complete cases after NA removal.")
  }

  # Stage 2: Iterative VIF elimination
  current_vars <- num_cols
  vif_log <- list()
  iteration <- 0L

  repeat {
    iteration <- iteration + 1L
    if (length(current_vars) < 3) break

    # Build dummy response for VIF (use row index)
    tmp <- dat_clean[, current_vars, drop = FALSE]
    tmp$.y <- seq_len(nrow(tmp))
    fmla <- stats::as.formula(
      paste(".y ~", paste(current_vars, collapse = " + "))
    )

    # Check for aliased (perfectly collinear) terms
    mod <- stats::lm(fmla, data = tmp)
    al <- stats::alias(mod)
    if (!is.null(al$Complete)) {
      aliased_var <- rownames(al$Complete)[1]
      if (verbose) {
        cli::cli_inform(
          "Perfect collinearity: removing {.val {aliased_var}}."
        )
      }
      removed <- c(removed, aliased_var)
      current_vars <- setdiff(current_vars, aliased_var)
      next
    }

    vif_vals <- tryCatch(
      car::vif(mod),
      error = function(e) {
        # Remove near-constant columns
        sds <- vapply(
          dat_clean[, current_vars, drop = FALSE],
          stats::sd, numeric(1), na.rm = TRUE
        )
        bad <- names(sds[sds < 1e-10])
        if (length(bad) > 0) {
          removed <<- c(removed, bad)
          current_vars <<- setdiff(current_vars, bad)
        }
        NULL
      }
    )
    if (is.null(vif_vals)) next

    max_vif <- max(vif_vals)
    max_var <- names(which.max(vif_vals))

    vif_log[[iteration]] <- data.frame(
      iteration = iteration,
      n_vars = length(current_vars),
      max_vif = round(max_vif, 2),
      removed = max_var,
      stringsAsFactors = FALSE
    )

    if (max_vif <= threshold) {
      if (verbose) {
        cli::cli_inform(
          "VIF converged: all VIF <= {threshold} ({length(current_vars)} variables)."
        )
      }
      break
    }

    if (verbose) {
      cli::cli_inform(
        "Iter {iteration}: max VIF = {round(max_vif, 1)} -> remove {.val {max_var}}"
      )
    }
    removed <- c(removed, max_var)
    current_vars <- setdiff(current_vars, max_var)

    if (length(current_vars) < 5) {
      if (verbose) cli::cli_warn("Fewer than 5 variables remaining.")
      break
    }
  }

  log_df <- if (length(vif_log) > 0) {
    do.call(rbind, vif_log)
  } else {
    data.frame(
      iteration = integer(0), n_vars = integer(0),
      max_vif = numeric(0), removed = character(0)
    )
  }

  list(
    selected = current_vars,
    removed = removed,
    vif_log = log_df,
    data = data[, current_vars, drop = FALSE]
  )
}


#' Prepare Species Data for Modeling
#'
#' Validates and prepares species occurrence data for the CAST pipeline.
#' Checks required columns, handles missing values, and splits into
#' training/test sets.
#'
#' @param data A `data.frame` with columns: `lon`, `lat`, `presence` (0/1),
#'   and environmental variables.
#' @param train_fraction Numeric. Fraction for training set. Default `0.7`.
#' @param seed Integer or `NULL`. Random seed for reproducible splitting.
#'
#' @return A list with components:
#' \describe{
#'   \item{`train`}{Training `data.frame`.}
#'   \item{`test`}{Test `data.frame`.}
#'   \item{`env_vars`}{Character vector of environmental variable names.}
#' }
#'
#' @export
cast_prepare <- function(data, train_fraction = 0.7, seed = NULL) {
  validate_species_data(data)
  env_vars <- get_env_vars(data)
  if (length(env_vars) == 0) {
    cli::cli_abort("No numeric environmental variables found in {.arg data}.")
  }

  n <- nrow(data)
  if (!is.null(seed)) set.seed(seed)
  train_idx <- sample.int(n, size = round(train_fraction * n))

  list(
    train = data[train_idx, , drop = FALSE],
    test = data[-train_idx, , drop = FALSE],
    env_vars = env_vars
  )
}
