# Prepare LOLD species-pair expression data for coexpr analysis
# Data source: vst_hog.RDS and N1_clean.RDS from ../coexpr-lold/
# Uses LEAF TISSUE only, following the pattern of prepare_brachypodium_data.R
#
# Species pairs:
#   1. Festuca pratensis (FPRA) vs Vulpia bromoides (VBRO)
#   2. Hordeum vulgare (HVUL) vs Hordeum jubatum (HJUB)
#   3. Briza maxima (BMAX) vs Briza media (BMED)

library(dplyr)
library(tidyr)

# Load raw data
vst_hog <- readRDS("/Users/martin/Documents/R/coexpr-lold/vst_hog.RDS")
N1_clean <- readRDS("/Users/martin/Documents/R/coexpr-lold/N1_clean.RDS")

cat("Available species in vst_hog:\n")
print(sort(unique(vst_hog$species)))
cat("\nAvailable species in N1_clean:\n")
print(sort(unique(N1_clean$species)))

# Helper: extract expression matrix for a species (leaf tissue only)
extract_expr <- function(vst, species_name, label) {
  sp_data <- vst |>
    filter(species == species_name, tissue == "leaf") |>
    select(GeneID, sample_id, vst.count)

  cat(sprintf("  %s: %d genes, %d samples\n",
              label, length(unique(sp_data$GeneID)),
              length(unique(sp_data$sample_id))))

  expr <- sp_data |>
    pivot_wider(
      id_cols = GeneID,
      names_from = sample_id,
      values_from = vst.count,
      values_fn = mean
    ) |>
    as.data.frame()

  rownames(expr) <- expr$GeneID
  as.matrix(expr[, -1])
}

# Helper: build ortholog table from N1_clean for a species pair
build_orthologs <- function(n1, sp1_name, sp2_name, sp1_label, sp2_label,
                            expr_sp1, expr_sp2) {
  orthologs_raw <- n1 |>
    filter(species %in% c(sp1_name, sp2_name)) |>
    select(HOG, species, GeneID)

  hogs_with_both <- orthologs_raw |>
    group_by(HOG) |>
    summarize(
      has_sp1 = any(species == sp1_name),
      has_sp2 = any(species == sp2_name),
      .groups = "drop"
    ) |>
    filter(has_sp1 & has_sp2) |>
    pull(HOG)

  cat(sprintf("  HOGs with both species: %d\n", length(hogs_with_both)))

  orth_sp1 <- orthologs_raw |>
    filter(species == sp1_name, HOG %in% hogs_with_both) |>
    select(HOG, gene_sp1 = GeneID)

  orth_sp2 <- orthologs_raw |>
    filter(species == sp2_name, HOG %in% hogs_with_both) |>
    select(HOG, gene_sp2 = GeneID)

  orth <- orth_sp1 |>
    inner_join(orth_sp2, by = "HOG", relationship = "many-to-many") |>
    select(gene_sp1, gene_sp2, HOG)

  cat(sprintf("  Total ortholog pairs: %d\n", nrow(orth)))

  # Filter to genes present in expression matrices
  orth <- orth |>
    filter(
      gene_sp1 %in% rownames(expr_sp1),
      gene_sp2 %in% rownames(expr_sp2)
    )

  cat(sprintf("  After filtering to expressed genes: %d\n", nrow(orth)))

  # Detect ortholog types
  orth <- orth |>
    group_by(gene_sp1) |>
    mutate(n_sp1 = n()) |>
    ungroup() |>
    group_by(gene_sp2) |>
    mutate(n_sp2 = n()) |>
    ungroup() |>
    mutate(
      type = case_when(
        n_sp1 == 1 & n_sp2 == 1 ~ "1:1",
        n_sp1 == 1 & n_sp2 > 1  ~ "1:N",
        n_sp1 > 1  & n_sp2 == 1 ~ "N:1",
        n_sp1 > 1  & n_sp2 > 1  ~ "N:M",
        TRUE ~ "unknown"
      )
    ) |>
    select(gene_sp1, gene_sp2, HOG, type)

  cat(sprintf("\n  Ortholog type distribution (%s vs %s):\n", sp1_label, sp2_label))
  print(table(orth$type))

  orth
}

output_dir <- "/Users/martin/Documents/R/coexpr/data-raw"

# ==============================================================================
# 1. Festuca pratensis vs Vulpia bromoides
# ==============================================================================

cat("\n==============================================================================\n")
cat("1. Festuca pratensis vs Vulpia bromoides (LEAF)\n")
cat("==============================================================================\n\n")

expr_fpra <- extract_expr(vst_hog, "Festuca pratensis", "F. pratensis")
expr_vbro <- extract_expr(vst_hog, "Vulpia bromoides", "V. bromoides")

cat(sprintf("\n  F. pratensis matrix: %d x %d\n", nrow(expr_fpra), ncol(expr_fpra)))
cat(sprintf("  V. bromoides matrix: %d x %d\n", nrow(expr_vbro), ncol(expr_vbro)))

festuca_vulpia_orthologs <- build_orthologs(
  N1_clean, "Festuca pratensis", "Vulpia bromoides",
  "F. pratensis", "V. bromoides",
  expr_fpra, expr_vbro
)

saveRDS(expr_fpra, file.path(output_dir, "expr_fpra.rds"))
saveRDS(expr_vbro, file.path(output_dir, "expr_vbro.rds"))
saveRDS(festuca_vulpia_orthologs,
        file.path(output_dir, "festuca_vulpia_orthologs.rds"))

# ==============================================================================
# 2. Hordeum vulgare vs Hordeum jubatum
# ==============================================================================

cat("\n==============================================================================\n")
cat("2. Hordeum vulgare vs Hordeum jubatum (LEAF)\n")
cat("==============================================================================\n\n")

expr_hvul <- extract_expr(vst_hog, "Hordeum vulgare", "H. vulgare")
expr_hjub <- extract_expr(vst_hog, "Hordeum jubatum", "H. jubatum")

cat(sprintf("\n  H. vulgare matrix: %d x %d\n", nrow(expr_hvul), ncol(expr_hvul)))
cat(sprintf("  H. jubatum matrix: %d x %d\n", nrow(expr_hjub), ncol(expr_hjub)))

hordeum_orthologs <- build_orthologs(
  N1_clean, "Hordeum vulgare", "Hordeum jubatum",
  "H. vulgare", "H. jubatum",
  expr_hvul, expr_hjub
)

saveRDS(expr_hvul, file.path(output_dir, "expr_hvul.rds"))
saveRDS(expr_hjub, file.path(output_dir, "expr_hjub.rds"))
saveRDS(hordeum_orthologs,
        file.path(output_dir, "hordeum_orthologs.rds"))

# ==============================================================================
# 3. Briza maxima vs Briza media
# ==============================================================================

cat("\n==============================================================================\n")
cat("3. Briza maxima vs Briza media (LEAF)\n")
cat("==============================================================================\n\n")

expr_bmax <- extract_expr(vst_hog, "Briza maxima", "B. maxima")
expr_bmed <- extract_expr(vst_hog, "Briza media", "B. media")

cat(sprintf("\n  B. maxima matrix: %d x %d\n", nrow(expr_bmax), ncol(expr_bmax)))
cat(sprintf("  B. media matrix: %d x %d\n", nrow(expr_bmed), ncol(expr_bmed)))

briza_orthologs <- build_orthologs(
  N1_clean, "Briza maxima", "Briza media",
  "B. maxima", "B. media",
  expr_bmax, expr_bmed
)

saveRDS(expr_bmax, file.path(output_dir, "expr_bmax.rds"))
saveRDS(expr_bmed, file.path(output_dir, "expr_bmed.rds"))
saveRDS(briza_orthologs,
        file.path(output_dir, "briza_orthologs.rds"))

# ==============================================================================
# Summary
# ==============================================================================

cat("\n==============================================================================\n")
cat("Data preparation complete!\n")
cat("==============================================================================\n\n")

cat("Files created in", output_dir, ":\n")
cat("  - expr_fpra.rds, expr_vbro.rds, festuca_vulpia_orthologs.rds\n")
cat("  - expr_hvul.rds, expr_hjub.rds, hordeum_orthologs.rds\n")
cat("  - expr_bmax.rds, expr_bmed.rds, briza_orthologs.rds\n")
