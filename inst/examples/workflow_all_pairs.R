# =============================================================================
# Cross-Species Co-Expression Analysis: Complete Workflow with All Ortholog Pairs
# =============================================================================
#
# This script demonstrates the full coexpr workflow using ALL ortholog pairs
# (1:1, 1:N, N:1, and N:M) from the Brachypodium example dataset.
#
# Species:
#   - Brachypodium distachyon (annual grass)
#   - Brachypodium sylvaticum (perennial grass)
#
# =============================================================================

library(coexpr)
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load Example Data
# -----------------------------------------------------------------------------

cat("=== Loading Brachypodium Example Data ===\n\n")

# Load expression matrices and ortholog relationships
expr_bdis <- readRDS(system.file("extdata", "expr_bdis_small.rds", package = "coexpr"))
expr_bsyl <- readRDS(system.file("extdata", "expr_bsyl_small.rds", package = "coexpr"))
orthologs <- readRDS(system.file("extdata", "brachypodium_orthologs_small.rds", package = "coexpr"))

cat("Expression data:\n")
cat("  B. distachyon:", nrow(expr_bdis), "genes ×", ncol(expr_bdis), "samples\n")
cat("  B. sylvaticum:", nrow(expr_bsyl), "genes ×", ncol(expr_bsyl), "samples\n")

cat("\nOrtholog pairs:", nrow(orthologs), "\n")
cat("Type distribution:\n")
print(table(orthologs$type))

# -----------------------------------------------------------------------------
# 2. Calculate Similarity Matrices
# -----------------------------------------------------------------------------

cat("\n=== Calculating Similarity Matrices (PCC+MR) ===\n\n")

sim_bdis <- calculate_pcc_mr(expr_bdis)
sim_bsyl <- calculate_pcc_mr(expr_bsyl)

cat("Similarity matrix dimensions:\n")
cat("  B. distachyon:", dim(sim_bdis), "\n")
cat("  B. sylvaticum:", dim(sim_bsyl), "\n")

# -----------------------------------------------------------------------------
# 3. Process Orthologs with "all_pairs" Strategy
# -----------------------------------------------------------------------------

cat("\n=== Processing Orthologs (all_pairs strategy) ===\n\n")

# Use all_pairs strategy to keep ALL ortholog combinations
orthologs_processed <- handle_multicopy_orthologs(
 orthologs,
 strategy = "all_pairs"
)

cat("After 'all_pairs' strategy:\n")
cat("  Total pairs:", nrow(orthologs_processed), "\n")
cat("  Type distribution:\n")
print(table(orthologs_processed$type))

# -----------------------------------------------------------------------------
# 4. Calculate Co-Expression Correlation Score (CCS)
# -----------------------------------------------------------------------------

cat("\n=== Calculating CCS ===\n\n")

# Calculate CCS using 1:1 orthologs as reference set
ccs_results <- calculate_ccs(
 sim_sp1 = sim_bdis,
 sim_sp2 = sim_bsyl,
 orthologs = orthologs_processed,
 use_only_1to1 = TRUE
)

cat("CCS results:\n")
cat("  Total pairs analyzed:", nrow(ccs_results), "\n")
cat("  CCS range:", round(range(ccs_results$ccs), 3), "\n")
cat("\nCCS distribution:\n")
print(summary(ccs_results$ccs))

# -----------------------------------------------------------------------------
# 5. Calculate Ortholog Rank Score (ORS)
# -----------------------------------------------------------------------------

cat("\n=== Calculating ORS ===\n\n")

ors_results <- calculate_ors(ccs_results, return_log = TRUE)

cat("ORS results:\n")
cat("  ORS range:", round(range(ors_results$ors_mean), 3), "\n")
cat("  logORS range:", round(range(ors_results$logors), 3), "\n")
cat("\nlogORS distribution:\n")
print(summary(ors_results$logors))

# -----------------------------------------------------------------------------
# 6. Test Significance
# -----------------------------------------------------------------------------

cat("\n=== Testing Significance ===\n\n")

ors_results <- test_ors_significance(ors_results, alpha = 0.05)

cat("Conservation levels:\n")
cat("  Below median (logORS < 0):", sum(ors_results$logors < 0), "\n")
cat("  Top 50% (logORS > 0):", sum(ors_results$logors > 0), "\n")
cat("  Top 10% (logORS > 1):", sum(ors_results$logors > 1), "\n")
cat("  Top 1% (logORS > 2):", sum(ors_results$logors > 2), "\n")
cat("  Top 0.1% (logORS > 3):", sum(ors_results$logors > 3), "\n")

# -----------------------------------------------------------------------------
# 7. Summarize Conservation by Ortholog Type
# -----------------------------------------------------------------------------

cat("\n=== Conservation by Ortholog Type ===\n\n")

conservation_summary <- summarize_conservation(ors_results, by_type = TRUE)
print(conservation_summary)

# -----------------------------------------------------------------------------
# 8. Visualizations
# -----------------------------------------------------------------------------

cat("\n=== Generating Plots ===\n\n")

# Plot 1: logORS Distribution
p1 <- ggplot(ors_results, aes(x = logors)) +
 geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8, color = "white") +
 geom_vline(xintercept = c(1, 2, 3), linetype = "dashed", color = "red") +
 labs(
   title = "Distribution of Expression Conservation (All Ortholog Pairs)",
   subtitle = paste(nrow(ors_results), "pairs including 1:1, 1:N, N:1, and N:M"),
   x = "logORS",
   y = "Count"
 ) +
 theme_minimal()

# Plot 2: CCS vs logORS colored by ortholog type
p2 <- ggplot(ors_results, aes(x = ccs, y = logors, color = type)) +
 geom_point(alpha = 0.6, size = 2) +
 geom_hline(yintercept = c(0, 1, 2), linetype = "dashed", alpha = 0.3) +
 scale_color_manual(
   values = c("1:1" = "steelblue", "1:N" = "orange",
              "N:1" = "darkgreen", "N:M" = "purple")
 ) +
 labs(
   title = "CCS vs logORS by Ortholog Type",
   x = "Co-Expression Correlation Score (CCS)",
   y = "logORS",
   color = "Ortholog Type"
 ) +
 theme_minimal() +
 theme(legend.position = "bottom")

# Plot 3: Boxplot by ortholog type
p3 <- ggplot(ors_results, aes(x = type, y = logors, fill = type)) +
 geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
 geom_hline(yintercept = c(1, 2), linetype = "dashed", color = "red", alpha = 0.5) +
 scale_fill_manual(
   values = c("1:1" = "steelblue", "1:N" = "orange",
              "N:1" = "darkgreen", "N:M" = "purple")
 ) +
 labs(
   title = "Expression Conservation by Ortholog Type",
   subtitle = "1:1 orthologs show higher conservation than multi-copy",
   x = "Ortholog Type",
   y = "logORS"
 ) +
 theme_minimal() +
 theme(legend.position = "none")

# Plot 4: Density plot by type
p4 <- ggplot(ors_results, aes(x = logors, fill = type)) +
 geom_density(alpha = 0.5) +
 geom_vline(xintercept = c(1, 2), linetype = "dashed", color = "red", alpha = 0.5) +
 scale_fill_manual(
   values = c("1:1" = "steelblue", "1:N" = "orange",
              "N:1" = "darkgreen", "N:M" = "purple")
 ) +
 labs(
   title = "logORS Density by Ortholog Type",
   x = "logORS",
   y = "Density",
   fill = "Ortholog Type"
 ) +
 theme_minimal() +
 theme(legend.position = "bottom")

# Display plots
print(p1)
print(p2)
print(p3)
print(p4)

# -----------------------------------------------------------------------------
# 9. Identify Highly Conserved Genes
# -----------------------------------------------------------------------------

cat("\n=== Highly Conserved Orthologs (Top 1%) ===\n\n")

highly_conserved <- ors_results |>
 filter(logors > 2) |>
 arrange(desc(logors))

cat("Count:", nrow(highly_conserved), "/", nrow(ors_results),
   "(", round(100 * nrow(highly_conserved) / nrow(ors_results), 1), "%)\n\n")

if (nrow(highly_conserved) > 0) {
 cat("Top 10 most conserved orthologs:\n")
 print(highly_conserved |>
         select(gene_sp1, gene_sp2, type, ccs, logors) |>
         head(10))
}

# -----------------------------------------------------------------------------
# 10. Compare Strategies
# -----------------------------------------------------------------------------

cat("\n=== Strategy Comparison ===\n\n")

# Also run with strict strategy for comparison
orthologs_strict <- handle_multicopy_orthologs(orthologs, strategy = "strict")
ccs_strict <- calculate_ccs(sim_bdis, sim_bsyl, orthologs_strict)
ors_strict <- calculate_ors(ccs_strict)

cat("Comparison:\n")
cat("  Strict (1:1 only):\n")
cat("    - Pairs:", nrow(ors_strict), "\n")
cat("    - Median logORS:", round(median(ors_strict$logors), 3), "\n")
cat("    - Top 1%:", sum(ors_strict$logors > 2), "\n\n")

cat("  All pairs:\n")
cat("    - Pairs:", nrow(ors_results), "\n")
cat("    - Median logORS:", round(median(ors_results$logors), 3), "\n")
cat("    - Top 1%:", sum(ors_results$logors > 2), "\n")

# By type breakdown
cat("\n  Median logORS by type:\n")
ors_results |>
 group_by(type) |>
 summarize(
   n = n(),
   median_logors = round(median(logors), 3),
   pct_top10 = round(100 * mean(logors > 1), 1)
 ) |>
 print()

cat("\n=== Workflow Complete ===\n")
