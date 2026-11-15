# Descriptive identifiability diagnostics for the saved GAM; no refitting.
.cast_gam_response_diagnostics <- function(object) {
  fit <- object$fit_full %||% object$fit
  model <- fit$models$gam$model
  if (is.null(model) || !inherits(model, "gam") || !length(model$smooth) ||
      !requireNamespace("mgcv", quietly = TRUE)) return(NULL)
  key <- .cast_digest(list(model$coefficients, model$sp, model$R, model$smooth, model$model,
    as.character(utils::packageVersion("mgcv"))))
  old <- object$diagnostics$gam_concurvity
  if (is.data.frame(old) && identical(attr(old, "fit_identity"), key)) return(old)
  value <- tryCatch({
    concurvity <- mgcv::concurvity(model, full = TRUE)
    labels <- vapply(model$smooth, `[[`, character(1), "label")
    variables <- vapply(model$smooth, function(x) paste(x$term, collapse = ":"), character(1))
    edf <- summary(model)$s.table[, "edf"]
    data.frame(term = colnames(concurvity), variable = variables[match(colnames(concurvity), labels)],
      estimate = unname(concurvity["estimate", ]), observed = unname(concurvity["observed", ]),
      worst = unname(concurvity["worst", ]), edf = unname(edf[match(colnames(concurvity), labels)]),
      interpretation = "Smooth-term substitutability; not causal identification", stringsAsFactors = FALSE)
  }, error = function(e) {
    cli::cli_warn("Saved GAM response diagnostics unavailable: {conditionMessage(e)}")
    NULL
  })
  if (!is.null(value)) attr(value, "fit_identity") <- key
  value
}
