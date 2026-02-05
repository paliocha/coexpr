# Prepare Brachypodium expression data for coexpr package examples
# Data source: vst_hog.RDS and N1_clean.RDS from OrthoFinder analysis
# Uses LEAF TISSUE only for vignette examples

library(dplyr)
library(tidyr)

# Load raw data
vst_hog <- readRDS("/Users/martin/Documents/R/coexpression-r/vst_hog.RDS")
N1_clean <- readRDS("/Users/martin/Documents/R/coexpression-r/N1_clean.RDS")

cat("Preparing Brachypodium LEAF expression data...\n")

# ==============================================================================
# 1. Extract expression data for both species (LEAF TISSUE ONLY)
# ==============================================================================

# Filter for Brachypodium distachyon (annual) - leaf tissue only
bdis_data <- vst_hog %>%
  filter(species == "Brachypodium distachyon", tissue == "leaf") %>%
  select(GeneID, sample_id, vst.count)

# Filter for Brachypodium sylvaticum (perennial) - leaf tissue only
bsyl_data <- vst_hog %>%
  filter(species == "Brachypodium sylvaticum", tissue == "leaf") %>%
  select(GeneID, sample_id, vst.count)

cat("  B. distachyon:", length(unique(bdis_data$GeneID)), "genes,",
    length(unique(bdis_data$sample_id)), "samples\n")
cat("  B. sylvaticum:", length(unique(bsyl_data$GeneID)), "genes,",
    length(unique(bsyl_data$sample_id)), "samples\n")

# ==============================================================================
# 2. Reshape to expression matrices (genes × samples)
# ==============================================================================

cat("\nReshaping to matrices...\n")

# B. distachyon: pivot to wide format
expr_bdis <- bdis_data %>%
  pivot_wider(
    id_cols = GeneID,
    names_from = sample_id,
    values_from = vst.count,
    values_fn = mean  # In case of duplicates, take mean
  ) %>%
  as.data.frame()

# Set rownames and convert to matrix
rownames(expr_bdis) <- expr_bdis$GeneID
expr_bdis <- as.matrix(expr_bdis[, -1])

# B. sylvaticum: pivot to wide format
expr_bsyl <- bsyl_data %>%
  pivot_wider(
    id_cols = GeneID,
    names_from = sample_id,
    values_from = vst.count,
    values_fn = mean
  ) %>%
  as.data.frame()

rownames(expr_bsyl) <- expr_bsyl$GeneID
expr_bsyl <- as.matrix(expr_bsyl[, -1])

cat("  B. distachyon matrix:", nrow(expr_bdis), "×", ncol(expr_bdis), "\n")
cat("  B. sylvaticum matrix:", nrow(expr_bsyl), "×", ncol(expr_bsyl), "\n")

# ==============================================================================
# 3. Extract ortholog relationships
# ==============================================================================

cat("\nExtracting ortholog relationships...\n")

# Get genes from both species in the same HOG
orthologs_raw <- N1_clean %>%
  filter(species %in% c("Brachypodium distachyon", "Brachypodium sylvaticum")) %>%
  select(HOG, species, GeneID)

# Find HOGs with both species
hogs_with_both <- orthologs_raw %>%
  group_by(HOG) %>%
  summarize(
    has_bdis = any(species == "Brachypodium distachyon"),
    has_bsyl = any(species == "Brachypodium sylvaticum"),
    .groups = "drop"
  ) %>%
  filter(has_bdis & has_bsyl) %>%
  pull(HOG)

cat("  HOGs with both species:", length(hogs_with_both), "\n")

# Create ortholog pairs
orthologs_bdis <- orthologs_raw %>%
  filter(species == "Brachypodium distachyon", HOG %in% hogs_with_both) %>%
  select(HOG, gene_bdis = GeneID)

orthologs_bsyl <- orthologs_raw %>%
  filter(species == "Brachypodium sylvaticum", HOG %in% hogs_with_both) %>%
  select(HOG, gene_bsyl = GeneID)

# Join to create all pairs
brachypodium_orthologs <- orthologs_bdis %>%
  inner_join(orthologs_bsyl, by = "HOG", relationship = "many-to-many") %>%
  select(
    gene_sp1 = gene_bdis,
    gene_sp2 = gene_bsyl,
    HOG
  )

cat("  Total ortholog pairs:", nrow(brachypodium_orthologs), "\n")

# Filter to only genes present in expression matrices
brachypodium_orthologs <- brachypodium_orthologs %>%
  filter(
    gene_sp1 %in% rownames(expr_bdis),
    gene_sp2 %in% rownames(expr_bsyl)
  )

cat("  After filtering to expressed genes:", nrow(brachypodium_orthologs), "\n")

# Detect ortholog types
brachypodium_orthologs <- brachypodium_orthologs %>%
  group_by(gene_sp1) %>%
  mutate(n_sp1 = n()) %>%
  ungroup() %>%
  group_by(gene_sp2) %>%
  mutate(n_sp2 = n()) %>%
  ungroup() %>%
  mutate(
    type = case_when(
      n_sp1 == 1 & n_sp2 == 1 ~ "1:1",
      n_sp1 == 1 & n_sp2 > 1  ~ "1:N",
      n_sp1 > 1  & n_sp2 == 1 ~ "N:1",
      n_sp1 > 1  & n_sp2 > 1  ~ "N:M",
      TRUE ~ "unknown"
    )
  ) %>%
  select(gene_sp1, gene_sp2, HOG, type)

cat("\nOrtholog type distribution:\n")
print(table(brachypodium_orthologs$type))

# ==============================================================================
# 4. Save processed data
# ==============================================================================

cat("\nSaving data objects...\n")

# Save expression matrices
saveRDS(expr_bdis, "/Users/martin/Documents/R/coexpr/data-raw/expr_bdis.rds")
saveRDS(expr_bsyl, "/Users/martin/Documents/R/coexpr/data-raw/expr_bsyl.rds")

# Save ortholog table
saveRDS(brachypodium_orthologs,
        "/Users/martin/Documents/R/coexpr/data-raw/brachypodium_orthologs.rds")

# Also save as package data (smaller subset for examples)
set.seed(123)

# Sample orthologs first, then get their genes (better approach)
# Select ~1000 ortholog pairs (this ensures we have matching genes in both species)
n_orthologs_sample <- min(1000, nrow(brachypodium_orthologs))
orthologs_sampled <- brachypodium_orthologs %>%
  slice_sample(n = n_orthologs_sample)

# Get the genes involved in these orthologs
sampled_genes_bdis <- unique(orthologs_sampled$gene_sp1)
sampled_genes_bsyl <- unique(orthologs_sampled$gene_sp2)

# Create smaller matrices
expr_bdis_small <- expr_bdis[sampled_genes_bdis, ]
expr_bsyl_small <- expr_bsyl[sampled_genes_bsyl, ]

# Use the sampled orthologs
brachypodium_orthologs_small <- orthologs_sampled

cat("\nSmall example dataset:\n")
cat("  B. distachyon:", nrow(expr_bdis_small), "genes\n")
cat("  B. sylvaticum:", nrow(expr_bsyl_small), "genes\n")
cat("  Orthologs:", nrow(brachypodium_orthologs_small), "\n")
cat("  Type distribution:\n")
print(table(brachypodium_orthologs_small$type))

# Save small datasets for vignette
saveRDS(expr_bdis_small,
        "/Users/martin/Documents/R/coexpr/data-raw/expr_bdis_small.rds")
saveRDS(expr_bsyl_small,
        "/Users/martin/Documents/R/coexpr/data-raw/expr_bsyl_small.rds")
saveRDS(brachypodium_orthologs_small,
        "/Users/martin/Documents/R/coexpr/data-raw/brachypodium_orthologs_small.rds")

# ==============================================================================
# 5. Create temporal expression data for ComplexHeatmap visualization
# ==============================================================================

cat("\nPreparing temporal expression data for heatmap...\n")

# Get temporal expression for sampled genes (leaf tissue only)
# Average across biological replicates per day
temporal_bdis <- vst_hog %>%

filter(species == "Brachypodium distachyon",
         tissue == "leaf",
         GeneID %in% sampled_genes_bdis) %>%
  group_by(GeneID, real_day) %>%
  summarize(expr = mean(vst.count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = GeneID, names_from = real_day, values_from = expr) %>%
  as.data.frame()

rownames(temporal_bdis) <- temporal_bdis$GeneID
temporal_bdis <- as.matrix(temporal_bdis[, -1])

temporal_bsyl <- vst_hog %>%
  filter(species == "Brachypodium sylvaticum",
         tissue == "leaf",
         GeneID %in% sampled_genes_bsyl) %>%
  group_by(GeneID, real_day) %>%
  summarize(expr = mean(vst.count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = GeneID, names_from = real_day, values_from = expr) %>%
  as.data.frame()

rownames(temporal_bsyl) <- temporal_bsyl$GeneID
temporal_bsyl <- as.matrix(temporal_bsyl[, -1])

cat("  B. distachyon temporal:", nrow(temporal_bdis), "genes ×",
    ncol(temporal_bdis), "time points\n")
cat("  B. sylvaticum temporal:", nrow(temporal_bsyl), "genes ×",
    ncol(temporal_bsyl), "time points\n")

# Save temporal data
saveRDS(temporal_bdis,
        "/Users/martin/Documents/R/coexpr/data-raw/temporal_bdis_small.rds")
saveRDS(temporal_bsyl,
        "/Users/martin/Documents/R/coexpr/data-raw/temporal_bsyl_small.rds")

cat("\n✓ Data preparation complete!\n")
cat("\nFiles created:\n")
cat("  - expr_bdis.rds (full)\n")
cat("  - expr_bsyl.rds (full)\n")
cat("  - brachypodium_orthologs.rds (full)\n")
cat("  - expr_bdis_small.rds (vignette)\n")
cat("  - expr_bsyl_small.rds (vignette)\n")
cat("  - brachypodium_orthologs_small.rds (vignette)\n")
