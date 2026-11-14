# Qualification is a semantic contract, not a statistical top-k filter.
.cast_is_varying <- function(x) {
  x <- x[is.finite(x)]
  length(x) > 1L && any(x != x[1L])
}

.cast_predictor_roles <- function(vars, roles = NULL) {
  out <- stats::setNames(rep("numeric", length(vars)), vars)
  # Known library identifiers must never be interpreted as environmental axes.
  out[grepl("(^|_)id$|^(hydro_h90m_)?(basin|segment|sub_catchment)$", vars)] <- "identifier"
  out[grepl("(^|_)(direction|aspect)$", vars)] <- "cyclic"
  out[grepl("(^|_)(dry_hue|wet_hue)(_|$)", vars)] <- "nominal"
  out[grepl("(^|_)order_(hack|horton|shreve|strahler|topo)$", vars)] <- "ordinal"
  if (!is.null(roles)) {
    allowed <- c("numeric", "binary", "ordinal", "identifier", "nominal", "cyclic")
    if (!is.character(roles) || is.null(names(roles)) || anyNA(roles) ||
        anyDuplicated(names(roles)) || any(!roles %in% allowed) || any(!vars %in% names(roles)))
      cli::cli_abort("roles must provide a unique named supported role for every candidate predictor.")
    out[] <- roles[vars]
  }
  out
}

.cast_predictor_qualification <- function(data, vars, roles = NULL) {
  role <- .cast_predictor_roles(vars, roles)
  eligible <- role %in% c("numeric", "binary", "ordinal")
  reason <- ifelse(eligible, "eligible", ifelse(role == "identifier", "identifier_not_environment",
    "requires_explicit_semantic_encoding"))
  numeric <- vapply(data[, vars, drop = FALSE], function(x) is.numeric(x) || is.logical(x), logical(1))
  reason[eligible & !numeric] <- "unsupported_non_numeric"
  eligible <- eligible & numeric
  varying <- vapply(vars, function(v) numeric[[v]] && .cast_is_varying(data[[v]]), logical(1))
  reason[eligible & !varying] <- "constant_or_all_missing"
  eligible <- eligible & varying
  data.frame(variable = vars, role = unname(role), eligible = unname(eligible),
    reason = unname(reason), stringsAsFactors = FALSE)
}
