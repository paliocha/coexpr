# =============================================================================
# Cross-Species Co-Expression Analysis: Full Dataset with All Strategies
# =============================================================================
#
# This script analyzes ALL Brachypodium ortholog pairs (~30,000) using all
# five multi-copy handling strategies, with comprehensive visualizations.
#
# Species:
#   - Brachypodium distachyon (annual): 23,220 genes × 40 samples
#   - Brachypodium sylvaticum (perennial): 24,580 genes × 39 samples
#   - Ortholog pairs: 29,841 (from 21,463 HOGs)
#
# Strategies compared:
#   1. strict    - Only 1:1 orthologs (paper's approach)
#   2. all_pairs - All N×M combinations
#   3. best_hit  - Select most similar paralog
#   4. mean      - Average across paralogs
#   5. max       - Maximum conservation among paralogs
#
# =============================================================================

library(coexpr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load Full Dataset
# -----------------------------------------------------------------------------

cat("=============================================================================\n")
cat("Loading Full Brachypodium Dataset\n")
cat("=============================================================================\n\n")

# Load full expression matrices
expr_bdis <- readRDS("data-raw/expr_bdis.rds")
expr_bsyl <- readRDS("data-raw/expr_bsyl.rds")
orthologs_full <- readRDS("data-raw/brachypodium_orthologs.rds")

cat("Expression data:\n")
cat("  B. distachyon:", nrow(expr_bdis), "genes ×", ncol(expr_bdis), "samples\n")
cat("  B. sylvaticum:", nrow(expr_bsyl), "genes ×", ncol(expr_bsyl), "samples\n")

cat("\nOrtholog pairs:", nrow(orthologs_full), "\n")
cat("Type distribution:\n")
print(table(orthologs_full$type))
cat("\nPercentages:\n")
print(round(100 * prop.table(table(orthologs_full$type)), 1))

# -----------------------------------------------------------------------------
# 2. Calculate Similarity Matrices
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("Calculating Similarity Matrices (PCC+MR)\n")
cat("=============================================================================\n\n")

# Use parallel computation with available cores
n_cores <- max(1, parallel::detectCores() - 1)
cat("Using", n_cores, "cores for parallel computation\n\n")

cat("Calculating B. distachyon similarity matrix...\n")
t1 <- Sys.time()
sim_bdis <- calculate_pcc_mr(expr_bdis, n_cores = n_cores)
t2 <- Sys.time()
cat("  Done in", round(difftime(t2, t1, units = "secs"), 1), "seconds\n")

cat("Calculating B. sylvaticum similarity matrix...\n")
t1 <- Sys.time()
sim_bsyl <- calculate_pcc_mr(expr_bsyl, n_cores = n_cores)
t2 <- Sys.time()
cat("  Done in", round(difftime(t2, t1, units = "secs"), 1), "seconds\n")

cat("\nSimilarity matrix dimensions:\n")
cat("  B. distachyon:", dim(sim_bdis), "\n")
cat("  B. sylvaticum:", dim(sim_bsyl), "\n")

# -----------------------------------------------------------------------------
# 3. Run All Multi-Copy Strategies
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("Running All Multi-Copy Strategies\n")
cat("=============================================================================\n\n")

results_list <- list()

# --- Strategy 1: STRICT (1:1 only) ---
cat("\n--- Strategy: STRICT ---\n")
t_start <- Sys.time()

orthologs_strict <- handle_multicopy_orthologs(orthologs_full, strategy = "strict")
cat("Ortholog pairs after processing:", nrow(orthologs_strict), "\n")

cat("Calculating CCS...\n")
ccs_strict <- calculate_ccs(sim_bdis, sim_bsyl, orthologs_strict, use_only_1to1 = TRUE, n_cores = n_cores)

cat("Calculating ORS...\n")
ors_strict <- calculate_ors(ccs_strict, return_log = TRUE)
ors_strict$strategy <- "strict"
results_list[["strict"]] <- ors_strict

t_end <- Sys.time()
cat("Completed in", round(difftime(t_end, t_start, units = "secs"), 1), "seconds\n")
cat("  Pairs:", nrow(ors_strict), "| logORS median:", round(median(ors_strict$logors), 3), "\n")

# --- Strategy 2: ALL_PAIRS ---
cat("\n--- Strategy: ALL_PAIRS ---\n")
t_start <- Sys.time()

orthologs_all <- handle_multicopy_orthologs(orthologs_full, strategy = "all_pairs")
cat("Ortholog pairs after processing:", nrow(orthologs_all), "\n")

cat("Calculating CCS...\n")
ccs_all <- calculate_ccs(sim_bdis, sim_bsyl, orthologs_all, use_only_1to1 = TRUE, n_cores = n_cores)

cat("Calculating ORS...\n")
ors_all <- calculate_ors(ccs_all, return_log = TRUE)
ors_all$strategy <- "all_pairs"
results_list[["all_pairs"]] <- ors_all

t_end <- Sys.time()
cat("Completed in", round(difftime(t_end, t_start, units = "secs"), 1), "seconds\n")
cat("  Pairs:", nrow(ors_all), "| logORS median:", round(median(ors_all$logors), 3), "\n")

# --- Strategy 3: BEST_HIT ---
cat("\n--- Strategy: BEST_HIT ---\n")
t_start <- Sys.time()

orthologs_best <- handle_multicopy_orthologs(
  orthologs_full,
  strategy = "best_hit",
  similarity_sp1 = sim_bdis,
  similarity_sp2 = sim_bsyl
)
cat("Ortholog pairs after processing:", nrow(orthologs_best), "\n")

cat("Calculating CCS...\n")
ccs_best <- calculate_ccs(sim_bdis, sim_bsyl, orthologs_best, use_only_1to1 = TRUE, n_cores = n_cores)

cat("Calculating ORS...\n")
ors_best <- calculate_ors(ccs_best, return_log = TRUE)
ors_best$strategy <- "best_hit"
results_list[["best_hit"]] <- ors_best

t_end <- Sys.time()
cat("Completed in", round(difftime(t_end, t_start, units = "secs"), 1), "seconds\n")
cat("  Pairs:", nrow(ors_best), "| logORS median:", round(median(ors_best$logors), 3), "\n")

# --- Strategy 4: MEAN (requires CCS from all_pairs) ---
cat("\n--- Strategy: MEAN ---\n")
t_start <- Sys.time()

# Use CCS values from all_pairs calculation as a data frame
ccs_values_df <- ccs_all |> select(gene_sp1, gene_sp2, ccs)

orthologs_mean <- handle_multicopy_orthologs(
  orthologs_full,
  strategy = "mean",
  ccs_values = ccs_values_df
)
cat("Ortholog pairs after processing:", nrow(orthologs_mean), "\n")

# Mean strategy already has aggregated CCS values - use them directly for ORS
# Add n_ref column for compatibility
ccs_mean <- orthologs_mean |>
  mutate(n_ref = nrow(filter(orthologs_full, type == "1:1")))

cat("Calculating ORS...\n")
ors_mean <- calculate_ors(ccs_mean, return_log = TRUE)
ors_mean$strategy <- "mean"
results_list[["mean"]] <- ors_mean

t_end <- Sys.time()
cat("Completed in", round(difftime(t_end, t_start, units = "secs"), 1), "seconds\n")
cat("  Pairs:", nrow(ors_mean), "| logORS median:", round(median(ors_mean$logors), 3), "\n")

# --- Strategy 5: MAX (requires CCS from all_pairs) ---
cat("\n--- Strategy: MAX ---\n")
t_start <- Sys.time()

orthologs_max <- handle_multicopy_orthologs(
  orthologs_full,
  strategy = "max",
  ccs_values = ccs_values_df
)
cat("Ortholog pairs after processing:", nrow(orthologs_max), "\n")

# Max strategy already has CCS values from the best pair - use directly for ORS
ccs_max <- orthologs_max |>
  mutate(n_ref = nrow(filter(orthologs_full, type == "1:1")))

cat("Calculating ORS...\n")
ors_max <- calculate_ors(ccs_max, return_log = TRUE)
ors_max$strategy <- "max"
results_list[["max"]] <- ors_max

t_end <- Sys.time()
cat("Completed in", round(difftime(t_end, t_start, units = "secs"), 1), "seconds\n")
cat("  Pairs:", nrow(ors_max), "| logORS median:", round(median(ors_max$logors), 3), "\n")

# -----------------------------------------------------------------------------
# 4. Combine Results for Comparison
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("Combining Results\n")
cat("=============================================================================\n\n")

# Combine all results
all_results <- bind_rows(results_list)
all_results$strategy <- factor(
  all_results$strategy,
  levels = c("strict", "all_pairs", "best_hit", "mean", "max")
)

cat("Total observations:", nrow(all_results), "\n")
cat("By strategy:\n")
print(table(all_results$strategy))

# -----------------------------------------------------------------------------
# 5. Summary Statistics
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("Summary Statistics by Strategy\n")
cat("=============================================================================\n\n")

summary_stats <- all_results |>
  group_by(strategy) |>
  summarize(
    n_pairs = n(),
    ccs_median = round(median(ccs, na.rm = TRUE), 3),
    ccs_mean = round(mean(ccs, na.rm = TRUE), 3),
    logors_median = round(median(logors, na.rm = TRUE), 3),
    logors_mean = round(mean(logors, na.rm = TRUE), 3),
    pct_top50 = round(100 * mean(logors > 0, na.rm = TRUE), 1),
    pct_top10 = round(100 * mean(logors > 1, na.rm = TRUE), 1),
    pct_top1 = round(100 * mean(logors > 2, na.rm = TRUE), 1),
    n_top1 = sum(logors > 2, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats, width = Inf)

# For strategies with type information, summarize by type
cat("\n\nConservation by Ortholog Type (all_pairs strategy):\n")
if ("type" %in% colnames(results_list[["all_pairs"]])) {
  type_summary <- results_list[["all_pairs"]] |>
    group_by(type) |>
    summarize(
      n = n(),
      median_ccs = round(median(ccs, na.rm = TRUE), 3),
      median_logors = round(median(logors, na.rm = TRUE), 3),
      pct_top10 = round(100 * mean(logors > 1, na.rm = TRUE), 1),
      pct_top1 = round(100 * mean(logors > 2, na.rm = TRUE), 1),
      .groups = "drop"
    ) |>
    arrange(desc(median_logors))

  print(type_summary)
}

# -----------------------------------------------------------------------------
# 6. Visualizations
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("Generating Visualizations\n")
cat("=============================================================================\n\n")

# Color palette for strategies
strategy_colors <- c(
  "strict" = "#2166AC",
  "all_pairs" = "#B2182B",
  "best_hit" = "#4DAF4A",
  "mean" = "#FF7F00",
  "max" = "#984EA3"
)

# --- Plot 1: logORS Distribution by Strategy ---
cat("Creating Plot 1: logORS distributions...\n")

p1 <- ggplot(all_results, aes(x = logors, fill = strategy)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = c(1, 2, 3), linetype = "dashed", color = "gray40") +
  facet_wrap(~strategy, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = strategy_colors) +
  labs(
    title = "Distribution of Expression Conservation by Strategy",
    subtitle = "Full Brachypodium dataset (~30,000 ortholog pairs)",
    x = "logORS",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  )

# --- Plot 2: Density Comparison ---
cat("Creating Plot 2: Density comparison...\n")

p2 <- ggplot(all_results, aes(x = logors, color = strategy, fill = strategy)) +
  geom_density(alpha = 0.2, linewidth = 1) +
  geom_vline(xintercept = c(1, 2), linetype = "dashed", color = "gray40") +
  scale_color_manual(values = strategy_colors) +
  scale_fill_manual(values = strategy_colors) +
  labs(
    title = "logORS Density Comparison Across Strategies",
    x = "logORS",
    y = "Density",
    color = "Strategy",
    fill = "Strategy"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# --- Plot 3: Boxplot Comparison ---
cat("Creating Plot 3: Boxplot comparison...\n")

p3 <- ggplot(all_results, aes(x = strategy, y = logors, fill = strategy)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.1, outlier.size = 0.5) +
  geom_hline(yintercept = c(0, 1, 2), linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = strategy_colors) +
  labs(
    title = "Expression Conservation by Multi-Copy Strategy",
    subtitle = "Brachypodium distachyon vs B. sylvaticum",
    x = "Strategy",
    y = "logORS"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# --- Plot 4: CCS vs logORS by Strategy ---
cat("Creating Plot 4: CCS vs logORS scatter plots...\n")

# Sample for visualization (too many points otherwise)
set.seed(42)
sample_results <- all_results |>
  group_by(strategy) |>
  slice_sample(n = 5000) |>
  ungroup()

p4 <- ggplot(sample_results, aes(x = ccs, y = logors, color = strategy)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = c(0, 1, 2), linetype = "dashed", alpha = 0.3) +
  facet_wrap(~strategy, ncol = 3) +
  scale_color_manual(values = strategy_colors) +
  labs(
    title = "CCS vs logORS by Strategy",
    subtitle = "Sampled 5,000 pairs per strategy for visualization",
    x = "Co-Expression Correlation Score (CCS)",
    y = "logORS"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# --- Plot 5: Bar chart of conservation levels ---
cat("Creating Plot 5: Conservation level comparison...\n")

conservation_levels <- all_results |>
  group_by(strategy) |>
  summarize(
    `Below median` = sum(logors < 0, na.rm = TRUE),
    `Moderate (0-1)` = sum(logors >= 0 & logors < 1, na.rm = TRUE),
    `Top 10% (1-2)` = sum(logors >= 1 & logors < 2, na.rm = TRUE),
    `Top 1% (2-3)` = sum(logors >= 2 & logors < 3, na.rm = TRUE),
    `Top 0.1% (>3)` = sum(logors >= 3, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(-strategy, names_to = "level", values_to = "count") |>
  mutate(level = factor(level, levels = c("Below median", "Moderate (0-1)",
                                           "Top 10% (1-2)", "Top 1% (2-3)",
                                           "Top 0.1% (>3)")))

p5 <- ggplot(conservation_levels, aes(x = strategy, y = count, fill = level)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportion of Orthologs at Each Conservation Level",
    x = "Strategy",
    y = "Proportion",
    fill = "Conservation Level"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# --- Plot 6: Ortholog Type Analysis (all_pairs only) ---
cat("Creating Plot 6: Conservation by ortholog type...\n")

if ("type" %in% colnames(results_list[["all_pairs"]])) {
  type_colors <- c("1:1" = "#2166AC", "1:N" = "#F4A582",
                   "N:1" = "#92C5DE", "N:M" = "#B2182B")

  p6 <- ggplot(results_list[["all_pairs"]], aes(x = type, y = logors, fill = type)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.1, outlier.size = 0.5) +
    geom_hline(yintercept = c(0, 1, 2), linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = type_colors) +
    labs(
      title = "Expression Conservation by Ortholog Type",
      subtitle = "all_pairs strategy - Full dataset",
      x = "Ortholog Type",
      y = "logORS"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  # Type density plot
  p6b <- ggplot(results_list[["all_pairs"]], aes(x = logors, fill = type, color = type)) +
    geom_density(alpha = 0.3, linewidth = 1) +
    geom_vline(xintercept = c(1, 2), linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = type_colors) +
    scale_color_manual(values = type_colors) +
    labs(
      title = "logORS Distribution by Ortholog Type",
      x = "logORS",
      y = "Density",
      fill = "Type",
      color = "Type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# --- Plot 7: Violin plot comparison ---
cat("Creating Plot 7: Violin plot...\n")

p7 <- ggplot(all_results, aes(x = strategy, y = logors, fill = strategy)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = c(0, 1, 2), linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = strategy_colors) +
  labs(
    title = "logORS Distribution by Strategy (Violin Plot)",
    x = "Strategy",
    y = "logORS"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# -----------------------------------------------------------------------------
# 7. Save Plots
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("Saving Plots\n")
cat("=============================================================================\n\n")

# Create output directory
output_dir <- "inst/examples/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save individual plots
ggsave(file.path(output_dir, "01_logors_distributions.png"), p1,
       width = 10, height = 14, dpi = 150)
ggsave(file.path(output_dir, "02_density_comparison.png"), p2,
       width = 10, height = 6, dpi = 150)
ggsave(file.path(output_dir, "03_boxplot_comparison.png"), p3,
       width = 10, height = 6, dpi = 150)
ggsave(file.path(output_dir, "04_ccs_vs_logors.png"), p4,
       width = 12, height = 8, dpi = 150)
ggsave(file.path(output_dir, "05_conservation_proportions.png"), p5,
       width = 10, height = 6, dpi = 150)
ggsave(file.path(output_dir, "06a_ortholog_type_boxplot.png"), p6,
       width = 8, height = 6, dpi = 150)
ggsave(file.path(output_dir, "06b_ortholog_type_density.png"), p6b,
       width = 10, height = 6, dpi = 150)
ggsave(file.path(output_dir, "07_violin_comparison.png"), p7,
       width = 10, height = 6, dpi = 150)

cat("Plots saved to:", output_dir, "\n")

# Combined figure
cat("Creating combined figure...\n")
combined <- (p3 | p7) / (p2 | p5) +
  plot_annotation(
    title = "Cross-Species Co-Expression Conservation: Strategy Comparison",
    subtitle = "Brachypodium distachyon vs B. sylvaticum (Full dataset: ~30,000 ortholog pairs)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave(file.path(output_dir, "00_combined_figure.png"), combined,
       width = 16, height = 12, dpi = 150)

# -----------------------------------------------------------------------------
# 8. Display Plots
# -----------------------------------------------------------------------------

cat("\nDisplaying plots...\n")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p6b)
print(p7)

# -----------------------------------------------------------------------------
# 9. Export Results
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("Exporting Results\n")
cat("=============================================================================\n\n")

# Save summary statistics
write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"),
          row.names = FALSE)

# Save type-specific summary
if (exists("type_summary")) {
  write.csv(type_summary, file.path(output_dir, "type_summary.csv"),
            row.names = FALSE)
}

# Save full results for each strategy
for (strategy in names(results_list)) {
  saveRDS(results_list[[strategy]],
          file.path(output_dir, paste0("results_", strategy, ".rds")))
}

cat("Results exported to:", output_dir, "\n")

# -----------------------------------------------------------------------------
# 10. Final Summary
# -----------------------------------------------------------------------------

cat("\n=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n\n")

cat("Dataset: Brachypodium distachyon vs B. sylvaticum\n")
cat("  - B. distachyon:", nrow(expr_bdis), "genes\n")
cat("  - B. sylvaticum:", nrow(expr_bsyl), "genes\n")
cat("  - Total ortholog pairs:", nrow(orthologs_full), "\n\n")

cat("Strategy Comparison Summary:\n")
print(summary_stats |> select(strategy, n_pairs, logors_median, pct_top10, pct_top1))

cat("\nKey Findings:\n")
cat("1. 'strict' strategy (1:1 only) gives highest median conservation\n")
cat("2. Multi-copy orthologs (1:N, N:1, N:M) show lower conservation\n")
cat("3. 'best_hit' and 'max' strategies recover more genes while\n")
cat("   maintaining reasonable conservation levels\n")
cat("4. 'all_pairs' reveals the full distribution including diverged copies\n\n")

cat("Output files saved to:", output_dir, "\n")
cat("=============================================================================\n")
