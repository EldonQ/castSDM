library(disdat)

cat("=== PO (presence-only) ===\n")
po <- disPo("SWI")
cat(dim(po), "\n")
print(str(po))
print(head(po, 3))

cat("\n=== BG (background) ===\n")
bg <- disBg("SWI")
cat(dim(bg), "\n")
print(head(bg, 3))

cat("\n=== PA (presence-absence test) ===\n")
pa <- disPa("SWI")
cat(dim(pa), "\n")
print(names(pa)[1:10])
print(head(pa[, 1:8], 3))

cat("\n=== ENV (test env) ===\n")
env <- disEnv("SWI")
cat(dim(env), "\n")
print(head(env, 3))

cat("\n=== Species list ===\n")
sp_cols <- setdiff(names(pa), c("group", "siteid", "x", "y"))
cat("N species:", length(sp_cols), "\n")
cat(sp_cols, sep = ", ")
cat("\n")

cat("\n=== Predictors ===\n")
cat(disPredictors("SWI"), sep = ", ")
cat("\n")

cat("\n=== Unique species in PO ===\n")
cat(length(unique(po$spid)), "species\n")
cat(sort(unique(po$spid))[1:5], sep = ", ")
cat("\n")

cat("\n=== Prevalence per species ===\n")
sp_prev <- colSums(pa[, sp_cols], na.rm = TRUE)
print(sort(sp_prev, decreasing = TRUE)[1:10])
