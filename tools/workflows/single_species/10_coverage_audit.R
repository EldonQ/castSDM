# 10_coverage_audit.R -----------------------------------------------------------
# Per-layer finite-value coverage over the China land mask. The failed run
# (Presences: 2 / 537 valid background cells) proves at least one layer in the
# 1099-layer full pool has large NA holes on the plateau; any such layer kills
# every cell it touches because training rows require finite values in ALL
# predictors. This audit quantifies coverage per layer so the pool can be
# rebuilt with a coverage gate.
#
# Reads each tif once (terra::global on the full grid) and reports
# finite / land-mask ratio. ~1099 x 0.1-0.3 s; run in background.

pool <- utils::read.csv("E:/Package/castSDM_workflows/single_species/pool_layers_full.csv",
                        stringsAsFactors = FALSE)
ref <- terra::rast("E:/CHECO26/N-SDM/data/mask_china.tif")
n_land <- sum(terra::values(ref, mat = FALSE) == 1, na.rm = TRUE)
message("China land cells: ", format(n_land, big.mark = ","))

rows <- vector("list", nrow(pool))
t0 <- Sys.time()
for (i in seq_len(nrow(pool))) {
  rows[[i]] <- tryCatch({
    r <- terra::rast(pool$path[i])
    n_fin <- terra::global(r, function(x) sum(is.finite(x), na.rm = TRUE))[[1]][1]
    data.frame(layer = pool$layer[i], group = pool$group[i],
               finite = n_fin, cover = n_fin / n_land,
               stringsAsFactors = FALSE)
  }, error = function(e) data.frame(layer = pool$layer[i], group = pool$group[i],
                                    finite = NA_integer_, cover = NA_real_,
                                    stringsAsFactors = FALSE))
  if (i %% 100L == 0L) {
    message(i, "/", nrow(pool), " elapsed ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  }
}
cov <- do.call(rbind, rows)
cov <- cov[order(cov$cover), , drop = FALSE]
rownames(cov) <- NULL
utils::write.csv(cov, "E:/Package/castSDM_workflows/single_species/pool_coverage_audit.csv",
                 row.names = FALSE)

cat("\n---- coverage distribution ----\n")
print(summary(cov$cover))
for (th in c(0.5, 0.8, 0.9, 0.95, 0.98, 0.99)) {
  cat(sprintf("layers with cover >= %.2f: %d\n", th, sum(cov$cover >= th, na.rm = TRUE)))
}
cat("\n---- 40 worst layers ----\n")
print(head(cov, 40), right = FALSE)
