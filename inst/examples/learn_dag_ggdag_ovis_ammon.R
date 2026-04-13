# ==============================================================================
# Causal DAG with dagitty + ggdag (Ovis ammon / 盘羊 example)
#
# Data: same use case as run_ovis_ammon.R — prefers saved pipeline output under
#   E:/Package/castSDM_multi_species/Ovis_ammon/cast_result.rds (uses $dag).
#   If RDS is missing, falls back to bundled `ovis_ammon` + cast_dag().
#
# Structure: castSDM already learned a consensus DAG (bnlearn bootstrap HC).
#   This script **exports that graph into dagitty** for layout, d-separation
#   tests, and publication-quality ggdag figures (no redefinition of edges).
#
# Suggested packages (install once):
#   install.packages(c("dagitty", "ggdag", "ggplot2", "dplyr", "patchwork"))
#
# Run:
#   source("E:/Package/cast/inst/examples/learn_dag_ggdag_ovis_ammon.R")
# ==============================================================================

# ── Paths & output ───────────────────────────────────────────────────────────
RDS_PATH   <- "E:/Package/castSDM_multi_species/Ovis_ammon/cast_result.rds"
OUT_DIR    <- "E:/Package/castSDM_multi_species/Ovis_ammon/figures"
OUT_MAIN   <- file.path(OUT_DIR, "dagitty_ggdag_ovis_ammon_main.png")
OUT_TESTS  <- file.path(OUT_DIR, "dagitty_ggdag_ovis_ammon_local_tests.png")
FIG_WIDTH  <- 12
FIG_HEIGHT <- 7.5
FIG_DPI    <- 300

SEED <- 42L

# ── Packages ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dagitty)
  library(ggdag)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
})

if (file.exists("DESCRIPTION") &&
    grepl("castSDM", readLines("DESCRIPTION", 1))) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(castSDM)
}

# ── Load DAG + training matrix for local tests ───────────────────────────────
load_ovis_dag_and_train <- function() {
  if (file.exists(RDS_PATH)) {
    res <- readRDS(RDS_PATH)
    dag <- res$dag
    if (is.null(dag) || !inherits(dag, "cast_dag")) {
      stop("RDS found but $dag is not a cast_dag object: ", RDS_PATH)
    }
    # Training frame for dagitty::localTests (env columns only)
    if (!is.null(res$fit) && is.list(res$fit$scaling)) {
      env_vars <- res$fit$env_vars
    } else {
      env_vars <- dag$nodes
    }
    tr <- NULL
    if (!is.null(res$fit) && !is.null(res$fit$models)) {
      m0 <- res$fit$models[[1]]
      if (!is.null(m0$training_data)) {
        tr <- m0$training_data
      }
    }
    if (is.null(tr)) {
      message("[info] No training_data in RDS; rebuilding train split from package ovis_ammon.")
      data(ovis_ammon, package = "castSDM", envir = environment())
      sp <- ovis_ammon
      split <- cast_prepare(sp, train_fraction = 0.7, seed = SEED)
      tr <- split$train
      env_vars <- split$env_vars
    }
    list(dag = dag, train = tr, env_vars = env_vars)
  } else {
    message("[info] RDS not found; using bundled ovis_ammon + cast_dag().")
    data(ovis_ammon, package = "castSDM", envir = environment())
    split <- cast_prepare(ovis_ammon, train_fraction = 0.7, seed = SEED)
    dag <- cast_dag(
      split$train,
      env_vars            = split$env_vars,
      R                   = 80L,
      strength_threshold  = 0.7,
      direction_threshold = 0.6,
      seed                = SEED,
      verbose             = FALSE
    )
    list(dag = dag, train = split$train, env_vars = split$env_vars)
  }
}

#' Convert cast_dag edges to a dagitty model string (single connected DAG).
cast_dag_to_dagitty_string <- function(dag) {
  ed <- dag$edges
  if (is.null(ed) || nrow(ed) == 0) {
    stop("No edges in cast_dag; cannot build a dagitty graph.")
  }
  # dagitty identifiers: use backticks for non-standard names if needed
  qn <- function(x) {
    x <- as.character(x)
    ifelse(grepl("^[A-Za-z][A-Za-z0-9_]*$", x), x, sprintf("`%s`", x))
  }
  parts <- sprintf("%s -> %s", qn(ed$from), qn(ed$to))
  sprintf("dag {\n%s\n}", paste(parts, collapse = "\n"))
}

obj <- load_ovis_dag_and_train()
dag <- obj$dag
train <- obj$train
env_vars <- obj$env_vars

dag_txt <- cast_dag_to_dagitty_string(dag)
g <- dagitty(dag_txt)

# Environmental data matrix (complete cases, numeric only)
X <- train[, intersect(env_vars, names(train)), drop = FALSE]
X <- X[, vapply(X, is.numeric, logical(1)), drop = FALSE]
X <- as.data.frame(na.omit(X))
stopifnot(nrow(X) >= 20L, ncol(X) >= 2L)

# ── dagitty: implied independencies + local tests vs data ───────────────────
# localTests: use type "cis" (Gaussian / partial correlation) for numeric env
# matrices; "cis.chisq" can overflow when many continuous predictors are
# discretized internally.
lt <- tryCatch(
  localTests(g, data = X, type = "cis"),
  error = function(e) {
    message("[warn] localTests failed (try fewer nodes or different type): ",
            conditionMessage(e))
    NULL
  }
)

# ── ggdag: layout from dagitty (same graph as castSDM) ───────────────────────
td <- tidy_dagitty(g)

# Map bootstrap strength onto directed edges for labels (exact castSDM values)
strength_df <- dag$edges %>%
  transmute(
    name  = paste(from, to, sep = "<-"),
    from  = as.character(from),
    to    = as.character(to),
    strength,
    direction
  )

td_dat <- td$data
if (all(c("from", "to") %in% names(td_dat))) {
  td_dat <- td_dat %>%
    dplyr::left_join(
      strength_df %>% dplyr::select(from, to, strength, direction),
      by = c("from", "to")
    )
  td$data <- td_dat
}

# Edge strength summary for caption (geom_dag_label_repel on mixed node/edge rows
# breaks when `strength` is missing on node rows or ggplot3/ggdag reject extra args)
edge_strength_caption <- paste(
  sprintf(
    "%s \u2192 %s: %.2f",
    dag$edges$from,
    dag$edges$to,
    as.numeric(dag$edges$strength)
  ),
  collapse = " | "
)
if (nchar(edge_strength_caption) > 500L) {
  edge_strength_caption <- paste0(substr(edge_strength_caption, 1L, 497L), "...")
}

# ---- Main figure: ggdag high-level API (stable across ggdag / ggplot2 versions)
p_main <- ggdag(td, text_size = 3.4, use_labels = "name") +
  theme_dag(base_size = 13, base_family = "sans") +
  labs(
    title = "Environmental DAG — Ovis ammon",
    subtitle = paste0(
      "dagitty + ggdag layout on edges learned by castSDM (bnlearn bootstrap HC; ",
      "R = ", dag$boot_R, ", score = ", dag$score,
      ", strength \u2265 ", dag$strength_threshold,
      ", direction \u2265 ", dag$direction_threshold, ")"
    ),
    caption = paste0("Bootstrap edge strengths: ", edge_strength_caption)
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, colour = "gray35", size = 10),
    plot.caption  = element_text(size = 7, colour = "gray40", hjust = 0.5),
    plot.margin   = margin(12, 12, 12, 12)
  )

# ---- Optional second figure: local test z-scores (dagitty vs same data matrix)
if (!is.null(lt) && nrow(as.data.frame(lt)) > 0L) {
  ltd <- as.data.frame(lt)
  ltd$test <- gsub("^\\(|\\)$", "", rownames(ltd))
  ltd$significant <- ifelse(abs(ltd$estimate) > 1.96, "outside", "within")
  ltd$significant <- factor(ltd$significant, levels = c("within", "outside"))

  p_tests <- ggplot(ltd, aes(x = reorder(test, estimate), y = estimate, fill = significant)) +
    geom_col(width = 0.72, colour = "gray25", linewidth = 0.25) +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", colour = "gray45") +
    coord_flip() +
    scale_fill_manual(
      name   = NULL,
      values = c(within = "#64748B", outside = "#B91C1C"),
      labels = c(within = "|z| \u2264 1.96", outside = "|z| > 1.96")
    ) +
    labs(
      title = "dagitty::localTests vs training environment matrix",
      subtitle = "Approximate z-statistics for implied independencies (localTests, type = cis).",
      x = NULL,
      y = "Estimate (approx. z)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "gray35"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(OUT_TESTS, p_tests, width = FIG_WIDTH, height = 0.55 * FIG_HEIGHT, dpi = FIG_DPI)
  message("Saved: ", OUT_TESTS)
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_MAIN, p_main, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
message("Saved: ", OUT_MAIN)

# Console summary (dagitty)
cat("\n--- dagitty model (first 800 chars) ---\n")
print(g)
cat("\n--- Implied adjustment sets (example: first two nodes) ---\n")
nds <- dag$nodes
if (length(nds) >= 2L) {
  print(adjustmentSets(g, exposure = nds[1], outcome = nds[2]))
}
