# =============================================================================
# Test: All similarity methods with Brachypodium full data
# =============================================================================
#
# Computes CCS/ORS for all combinations of similarity methods:
#   1. PCC       (Pearson correlation only)
#   2. SCC       (Spearman correlation only)
#   3. PCC+MR    (Pearson + Mutual Rank)
#   4. SCC+MR    (Spearman + Mutual Rank)
#   5. MI+CLR    (Mutual Information + Context Likelihood Ratio)
#
# Uses all_pairs multicopy strategy and all orthologs as reference.
# =============================================================================

library(coexpr)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

cat("=== Loading Brachypodium data ===\n")

expr_bdis <- readRDS("data-raw/expr_bdis.rds")
expr_bsyl <- readRDS("data-raw/expr_bsyl.rds")
orthologs <- readRDS("data-raw/brachypodium_orthologs.rds")

cat("  B. distachyon:", nrow(expr_bdis), "genes x", ncol(expr_bdis), "samples\n")
cat("  B. sylvaticum:", nrow(expr_bsyl), "genes x", ncol(expr_bsyl), "samples\n")
cat("  Ortholog pairs:", nrow(orthologs), "\n")
cat("  Type distribution:\n")
print(table(orthologs$type))

# Use all_pairs strategy
orthologs_all <- handle_multicopy_orthologs(orthologs, strategy = "all_pairs")
cat("\n  After all_pairs:", nrow(orthologs_all), "pairs\n")

# Detect cores
n_cores <- max(1, parallel::detectCores() - 1)
cat("  Using", n_cores, "cores\n")

# Cache directory for similarity matrices
cache_dir <- "inst/scripts/cache"

# -----------------------------------------------------------------------------
# 2. Define similarity methods
# -----------------------------------------------------------------------------

methods <- list(
  PCC = list(
    fn = function(expr) calculate_pcc_mr(
      expr, method = "pcc", cor_method = "pearson",
      return_tri = TRUE, cache_dir = cache_dir
    ),
    label = "PCC (Pearson)"
  ),
  SCC = list(
    fn = function(expr) calculate_pcc_mr(
      expr, method = "pcc", cor_method = "spearman",
      return_tri = TRUE, cache_dir = cache_dir
    ),
    label = "SCC (Spearman)"
  ),
  PCC_MR = list(
    fn = function(expr) calculate_pcc_mr(
      expr, method = "pcc_mr", cor_method = "pearson",
      n_cores = n_cores, return_tri = TRUE, cache_dir = cache_dir
    ),
    label = "PCC+MR"
  ),
  SCC_MR = list(
    fn = function(expr) calculate_pcc_mr(
      expr, method = "pcc_mr", cor_method = "spearman",
      n_cores = n_cores, return_tri = TRUE, cache_dir = cache_dir
    ),
    label = "SCC+MR"
  ),
  MI_CLR = list(
    fn = function(expr) calculate_mi_clr(
      expr, n_bins = 10, n_cores = n_cores,
      return_tri = TRUE, cache_dir = cache_dir
    ),
    label = "MI+CLR"
  )
)

# -----------------------------------------------------------------------------
# 3. Run pipeline for each method
# -----------------------------------------------------------------------------

all_results <- list()

for (method_name in names(methods)) {
  method <- methods[[method_name]]
  cat("\n=== ", method$label, " ===\n", sep = "")

  # Compute similarity matrices
  cat("  Computing similarity for B. distachyon...\n")
  t0 <- proc.time()
  sim_bdis <- method$fn(expr_bdis)
  t1 <- proc.time()
  cat("  Done in", round((t1 - t0)[3], 1), "sec\n")

  cat("  Computing similarity for B. sylvaticum...\n")
  t0 <- proc.time()
  sim_bsyl <- method$fn(expr_bsyl)
  t1 <- proc.time()
  cat("  Done in", round((t1 - t0)[3], 1), "sec\n")

  # CCS with all orthologs as reference
  cat("  Computing CCS (all orthologs as reference)...\n")
  t0 <- proc.time()
  ccs <- calculate_ccs(
    sim_bdis, sim_bsyl, orthologs_all,
    use_only_1to1 = FALSE,
    n_cores = n_cores
  )
  t1 <- proc.time()
  cat("  CCS done in", round((t1 - t0)[3], 1), "sec\n")

  # ORS + logORS
  cat("  Computing ORS...\n")
  ors <- calculate_ors(ccs, return_log = TRUE)
  ors$method <- method$label

  all_results[[method_name]] <- ors

  # Quick summary
  cat("  Pairs:", nrow(ors), "\n")
  cat("  CCS:    median =", round(median(ors$CCS, na.rm = TRUE), 4),
      " mean =", round(mean(ors$CCS, na.rm = TRUE), 4), "\n")
  cat("  logORS: median =", round(median(ors$logORS, na.rm = TRUE), 4),
      " mean =", round(mean(ors$logORS, na.rm = TRUE), 4), "\n")

  # Clean up similarity matrices to free memory
  rm(sim_bdis, sim_bsyl)
  gc(verbose = FALSE)
}

# -----------------------------------------------------------------------------
# 4. Combined summary
# -----------------------------------------------------------------------------

cat("\n\n=== SUMMARY ACROSS ALL METHODS ===\n\n")

combined <- bind_rows(all_results)
combined$method <- factor(combined$method, levels = sapply(methods, `[[`, "label"))

summary_table <- combined |>
  group_by(method) |>
  summarize(
    n = n(),
    CCS_median = round(median(CCS, na.rm = TRUE), 4),
    CCS_mean   = round(mean(CCS, na.rm = TRUE), 4),
    CCS_sd     = round(sd(CCS, na.rm = TRUE), 4),
    CCS_na     = sum(is.na(CCS)),
    logORS_median = round(median(logORS, na.rm = TRUE), 4),
    logORS_mean   = round(mean(logORS, na.rm = TRUE), 4),
    pct_top10  = round(100 * mean(logORS > 1, na.rm = TRUE), 1),
    pct_top1   = round(100 * mean(logORS > 2, na.rm = TRUE), 1),
    pct_top01  = round(100 * mean(logORS > 3, na.rm = TRUE), 1),
    .groups = "drop"
  )

print(as.data.frame(summary_table), row.names = FALSE, right = FALSE)

# By ortholog type within each method
cat("\n\n=== BY ORTHOLOG TYPE ===\n\n")

if ("type" %in% colnames(combined)) {
  type_summary <- combined |>
    group_by(method, type) |>
    summarize(
      n = n(),
      CCS_median = round(median(CCS, na.rm = TRUE), 4),
      logORS_median = round(median(logORS, na.rm = TRUE), 4),
      pct_top10  = round(100 * mean(logORS > 1, na.rm = TRUE), 1),
      .groups = "drop"
    )

  print(as.data.frame(type_summary), row.names = FALSE, right = FALSE)
}

# Save results
output_file <- "inst/scripts/results_all_methods.rds"
saveRDS(all_results, output_file)
cat("\n\nResults saved to:", output_file, "\n")
