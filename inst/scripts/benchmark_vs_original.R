#!/usr/bin/env Rscript
#
# Benchmark script for coexpr package
#
# This script downloads plant expression data and ortholog information from
# public databases, then runs the coexpr pipeline to validate against the
# methodology described in Grønvold & Hvidsten.
#
# Usage: Rscript benchmark_vs_original.R
#
# The script will:
# 1. Download expression data for 3-5 plant species (or use simulated fallback)
# 2. Download ortholog relationships from Ensembl Plants
# 3. Calculate PCC+MR similarity matrices
# 4. Calculate CCS and ORS for species pairs
# 5. Generate validation plots
#
# Output: benchmark_results/ directory with plots and summary statistics
#

library(coexpr)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Configuration
OUTPUT_DIR <- "benchmark_results"
CACHE_DIR <- file.path(OUTPUT_DIR, "cache")
USE_SIMULATED_FALLBACK <- TRUE  # Set to FALSE to require real data

# Species to analyze (subset of the paper's 5 species)
# Using common names and their Ensembl identifiers
SPECIES <- list(
  At = list(
    name = "Arabidopsis thaliana",
    ensembl = "athaliana_eg_gene",
    short = "At"
  ),
  Os = list(
    name = "Oryza sativa",
    ensembl = "osativa_eg_gene",
    short = "Os"
  ),
  Sl = list(
    name = "Solanum lycopersicum",
    ensembl = "slycopersicum_eg_gene",
    short = "Sl"
  )
)

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)


#' Download expression data from EBI Expression Atlas
#'
#' @param species_ensembl Ensembl species identifier
#' @param max_samples Maximum samples to download
#'
#' @return Expression matrix (genes x samples) or NULL if failed
download_expression_atlas <- function(species_ensembl, max_samples = 200) {
  message(sprintf("Attempting to download expression data for %s...", species_ensembl))

  # Check for httr package
  if (!requireNamespace("httr", quietly = TRUE)) {
    message("Package 'httr' not available, skipping download")
    return(NULL)
  }

  # Expression Atlas API endpoint
  # This is a simplified approach - real implementation would query specific experiments
  base_url <- "https://www.ebi.ac.uk/gxa/experiments-content"

  tryCatch({
    # For demonstration, we'll indicate this would query the API
    # In practice, you'd need to:
    # 1. Search for experiments for this species
    # 2. Download normalized expression matrices
    # 3. Combine across experiments

    message("  Note: Full Expression Atlas integration requires API setup")
    message("  Using simulated data instead")
    return(NULL)
  }, error = function(e) {
    message(sprintf("  Failed to download: %s", e$message))
    return(NULL)
  })
}


#' Download orthologs from Ensembl Plants via biomaRt
#'
#' @param sp1_ensembl Ensembl identifier for species 1
#' @param sp2_ensembl Ensembl identifier for species 2
#'
#' @return Data frame with ortholog pairs or NULL if failed
download_orthologs_biomart <- function(sp1_ensembl, sp2_ensembl) {
  message(sprintf("Attempting to download orthologs between %s and %s...",
                  sp1_ensembl, sp2_ensembl))

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("Package 'biomaRt' not available, skipping download")
    return(NULL)
  }

  tryCatch({
    # Connect to Ensembl Plants
    mart <- biomaRt::useEnsembl(
      biomart = "plants_mart",
      host = "https://plants.ensembl.org",
      dataset = sp1_ensembl
    )

    # Get ortholog attribute name for species 2
    sp2_short <- gsub("_eg_gene$", "", sp2_ensembl)
    homolog_attr <- sprintf("%s_homolog_ensembl_gene", sp2_short)

    # Query orthologs
    orthologs <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", homolog_attr, "homolog_type"),
      mart = mart
    )

    if (nrow(orthologs) == 0) {
      message("  No orthologs found")
      return(NULL)
    }

    # Rename columns
    colnames(orthologs) <- c("gene_sp1", "gene_sp2", "ortholog_type")

    # Map ortholog types
    orthologs <- orthologs |>
      mutate(type = case_when(
        ortholog_type == "ortholog_one2one" ~ "1:1",
        ortholog_type == "ortholog_one2many" ~ "1:N",
        ortholog_type == "ortholog_many2many" ~ "N:M",
        TRUE ~ "other"
      )) |>
      filter(!is.na(gene_sp2), gene_sp2 != "")

    message(sprintf("  Found %d ortholog pairs", nrow(orthologs)))
    return(orthologs)
  }, error = function(e) {
    message(sprintf("  Failed to download: %s", e$message))
    return(NULL)
  })
}


#' Generate simulated expression data for benchmarking
#'
#' Creates expression matrices with realistic co-expression patterns
#'
#' @param n_genes Number of genes
#' @param n_samples Number of samples
#' @param n_modules Number of co-expression modules
#' @param seed Random seed
#'
#' @return Expression matrix (genes x samples)
simulate_expression_data <- function(n_genes = 1000, n_samples = 100,
                                     n_modules = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate base expression with co-expression modules
  genes_per_module <- n_genes %/% n_modules
  expr <- matrix(0, nrow = n_genes, ncol = n_samples)

  for (m in seq_len(n_modules)) {
    # Module-specific pattern
    module_pattern <- rnorm(n_samples)

    # Genes in this module
    start_gene <- (m - 1) * genes_per_module + 1
    end_gene <- min(m * genes_per_module, n_genes)

    for (g in start_gene:end_gene) {
      # Gene expression = module pattern + noise
      noise_level <- runif(1, 0.3, 0.8)
      expr[g, ] <- module_pattern * (1 - noise_level) + rnorm(n_samples) * noise_level
    }
  }

  # Handle remaining genes (not in modules)
  remaining <- (n_modules * genes_per_module + 1):n_genes
  if (length(remaining) > 0 && remaining[1] <= n_genes) {
    expr[remaining, ] <- matrix(rnorm(length(remaining) * n_samples),
                                nrow = length(remaining))
  }

  rownames(expr) <- paste0("Gene", seq_len(n_genes))
  colnames(expr) <- paste0("Sample", seq_len(n_samples))

  return(expr)
}


#' Generate simulated orthologs with realistic patterns
#'
#' @param genes_sp1 Gene IDs for species 1
#' @param genes_sp2 Gene IDs for species 2
#' @param prop_1to1 Proportion of 1:1 orthologs
#' @param seed Random seed
#'
#' @return Data frame with ortholog pairs
simulate_orthologs <- function(genes_sp1, genes_sp2, prop_1to1 = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n1 <- length(genes_sp1)
  n2 <- length(genes_sp2)
  n_orthologs <- min(n1, n2)

  # Create 1:1 orthologs
  n_1to1 <- round(n_orthologs * prop_1to1)

  orthologs <- data.frame(
    gene_sp1 = genes_sp1[seq_len(n_1to1)],
    gene_sp2 = genes_sp2[seq_len(n_1to1)],
    type = "1:1",
    stringsAsFactors = FALSE
  )

  # Add some 1:N orthologs
  n_1toN <- round(n_orthologs * 0.2)
  if (n_1toN > 0) {
    for (i in seq_len(n_1toN)) {
      idx1 <- n_1to1 + i
      if (idx1 > n1) break

      # Map to 2-3 genes in species 2
      n_targets <- sample(2:3, 1)
      idx2_start <- n_1to1 + (i - 1) * 3 + 1
      idx2_end <- min(idx2_start + n_targets - 1, n2)

      if (idx2_start > n2) break

      for (idx2 in idx2_start:idx2_end) {
        orthologs <- rbind(orthologs, data.frame(
          gene_sp1 = genes_sp1[idx1],
          gene_sp2 = genes_sp2[idx2],
          type = "1:N",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(orthologs)
}


#' Run the full benchmark pipeline
#'
#' @param species_data List with expression matrices per species
#' @param orthologs_data List with ortholog tables per species pair
#'
#' @return List with CCS and ORS results
run_benchmark <- function(species_data, orthologs_data) {
  results <- list()

  # Calculate similarity matrices for each species
  message("\n=== Calculating similarity matrices ===")
  sim_matrices <- list()

  for (sp_name in names(species_data)) {
    message(sprintf("\nProcessing %s (%d genes x %d samples)...",
                    sp_name,
                    nrow(species_data[[sp_name]]),
                    ncol(species_data[[sp_name]])))

    cache_file <- file.path(CACHE_DIR, sprintf("sim_%s.rds", sp_name))

    if (file.exists(cache_file)) {
      message("  Loading cached similarity matrix...")
      sim_matrices[[sp_name]] <- readRDS(cache_file)
    } else {
      sim_matrices[[sp_name]] <- calculate_pcc_mr(
        species_data[[sp_name]],
        method = "pcc_mr",
        n_cores = 1
      )
      saveRDS(sim_matrices[[sp_name]], cache_file)
    }

    message(sprintf("  Similarity matrix: %d x %d genes",
                    sim_matrices[[sp_name]]$n,
                    sim_matrices[[sp_name]]$n))
  }

  # Calculate CCS for each species pair
  message("\n=== Calculating CCS for species pairs ===")

  for (pair_name in names(orthologs_data)) {
    species_pair <- strsplit(pair_name, "_")[[1]]
    sp1 <- species_pair[1]
    sp2 <- species_pair[2]

    message(sprintf("\nProcessing %s vs %s...", sp1, sp2))

    orthologs <- orthologs_data[[pair_name]]
    message(sprintf("  %d ortholog pairs (%d 1:1)",
                    nrow(orthologs),
                    sum(orthologs$type == "1:1")))

    # Calculate CCS
    ccs_results <- calculate_ccs(
      sim_matrices[[sp1]],
      sim_matrices[[sp2]],
      orthologs,
      use_only_1to1 = TRUE,
      handle_self_diagonal = "mean"
    )

    # Calculate ORS
    message("Calculating ORS...")
    ors_results <- calculate_ors(ccs_results)

    results[[pair_name]] <- list(
      ccs = ccs_results,
      ors = ors_results
    )

    # Summary statistics
    message(sprintf("  CCS: median = %.3f, mean = %.3f",
                    median(ccs_results$ccs, na.rm = TRUE),
                    mean(ccs_results$ccs, na.rm = TRUE)))
    message(sprintf("  logORS: median = %.3f, mean = %.3f",
                    median(ors_results$logORS, na.rm = TRUE),
                    mean(ors_results$logORS, na.rm = TRUE)))
  }

  return(results)
}


#' Generate validation plots
#'
#' @param results Results from run_benchmark
#' @param output_dir Directory for plot files
generate_plots <- function(results, output_dir) {
  message("\n=== Generating validation plots ===")

  # 1. logORS density plot
  message("  Creating logORS density plot...")

  ors_combined <- bind_rows(
    lapply(names(results), function(pair_name) {
      results[[pair_name]]$ors |>
        mutate(species_pair = pair_name)
    })
  )

  p1 <- ggplot(ors_combined, aes(x = logORS, fill = species_pair)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "logORS Distribution by Species Pair",
      subtitle = "Higher values indicate more conserved co-expression",
      x = "logORS",
      y = "Density",
      fill = "Species Pair"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "logORS_density.pdf"), p1,
         width = 8, height = 6)

  # 2. CCS density plot
  message("  Creating CCS density plot...")

  ccs_combined <- bind_rows(
    lapply(names(results), function(pair_name) {
      results[[pair_name]]$ccs |>
        mutate(species_pair = pair_name)
    })
  )

  p2 <- ggplot(ccs_combined, aes(x = ccs, fill = species_pair)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "CCS Distribution by Species Pair",
      x = "Co-expression Correlation Score (CCS)",
      y = "Density",
      fill = "Species Pair"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "CCS_density.pdf"), p2,
         width = 8, height = 6)

  # 3. logORS boxplot by ortholog type
  message("  Creating logORS boxplot by ortholog type...")

  ors_with_type <- bind_rows(
    lapply(names(results), function(pair_name) {
      results[[pair_name]]$ors |>
        mutate(species_pair = pair_name)
    })
  )

  if ("type" %in% colnames(ors_with_type)) {
    p3 <- ggplot(ors_with_type, aes(x = type, y = logORS, fill = type)) +
      geom_boxplot() +
      facet_wrap(~species_pair) +
      labs(
        title = "logORS by Ortholog Type",
        subtitle = "1:1 orthologs expected to have higher ORS",
        x = "Ortholog Type",
        y = "logORS"
      ) +
      theme_minimal() +
      theme(legend.position = "none")

    ggsave(file.path(output_dir, "logORS_by_type.pdf"), p3,
           width = 10, height = 6)
  }

  # 4. Summary statistics table
  message("  Creating summary statistics...")

  summary_stats <- bind_rows(
    lapply(names(results), function(pair_name) {
      ors <- results[[pair_name]]$ors
      ccs <- results[[pair_name]]$ccs

      data.frame(
        species_pair = pair_name,
        n_pairs = nrow(ors),
        median_CCS = median(ccs$ccs, na.rm = TRUE),
        mean_CCS = mean(ccs$ccs, na.rm = TRUE),
        median_logORS = median(ors$logORS, na.rm = TRUE),
        mean_logORS = mean(ors$logORS, na.rm = TRUE),
        prop_ORS_gt_1 = mean(ors$logORS > 1, na.rm = TRUE),
        prop_ORS_gt_2 = mean(ors$logORS > 2, na.rm = TRUE)
      )
    })
  )

  write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"),
            row.names = FALSE)

  message("  Plots saved to: ", output_dir)
  return(summary_stats)
}


#' Validate results against expected patterns
#'
#' @param results Results from run_benchmark
#'
#' @return List with validation checks
validate_results <- function(results) {
  message("\n=== Validation Checks ===")
  checks <- list()

  for (pair_name in names(results)) {
    ors <- results[[pair_name]]$ors

    # Check 1: Median logORS should be around 0 (random expectation)
    median_logORS <- median(ors$logORS, na.rm = TRUE)
    checks[[paste0(pair_name, "_median_near_zero")]] <- abs(median_logORS) < 0.5

    message(sprintf("  %s: Median logORS = %.3f (expected ~0)",
                    pair_name, median_logORS))

    # Check 2: Some proportion should have high ORS (top 10% = logORS > 1)
    prop_high_ors <- mean(ors$logORS > 1, na.rm = TRUE)
    checks[[paste0(pair_name, "_has_conserved_pairs")]] <- prop_high_ors > 0.05

    message(sprintf("  %s: %.1f%% have logORS > 1 (expected >5%%)",
                    pair_name, prop_high_ors * 100))

    # Check 3: If type column exists, 1:1 should have higher ORS than 1:N
    if ("type" %in% colnames(ors)) {
      ors_1to1 <- ors$logORS[ors$type == "1:1"]
      ors_1toN <- ors$logORS[ors$type == "1:N"]

      if (length(ors_1toN) > 10) {
        diff <- median(ors_1to1, na.rm = TRUE) - median(ors_1toN, na.rm = TRUE)
        checks[[paste0(pair_name, "_1to1_higher")]] <- diff > 0

        message(sprintf("  %s: 1:1 vs 1:N median difference = %.3f (expected >0)",
                        pair_name, diff))
      }
    }
  }

  # Overall pass/fail
  n_passed <- sum(unlist(checks))
  n_total <- length(checks)
  message(sprintf("\n  Validation: %d/%d checks passed", n_passed, n_total))

  return(checks)
}


# Main execution
main <- function() {
  message("=== coexpr Benchmark Script ===")
  message(sprintf("Output directory: %s", OUTPUT_DIR))

  # Try to download real data, fall back to simulated
  species_data <- list()
  orthologs_data <- list()

  # Attempt real data download
  for (sp_name in names(SPECIES)) {
    sp_info <- SPECIES[[sp_name]]
    expr <- download_expression_atlas(sp_info$ensembl)

    if (is.null(expr) && USE_SIMULATED_FALLBACK) {
      message(sprintf("  Using simulated data for %s", sp_info$name))
      expr <- simulate_expression_data(
        n_genes = 500,
        n_samples = 80,
        n_modules = 8,
        seed = which(names(SPECIES) == sp_name) * 100
      )
      # Add species prefix to gene names
      rownames(expr) <- paste0(sp_name, "_", rownames(expr))
    }

    if (!is.null(expr)) {
      species_data[[sp_name]] <- expr
    }
  }

  # Check we have at least 2 species
  if (length(species_data) < 2) {
    stop("Need at least 2 species for benchmark")
  }

  # Get orthologs for each species pair
  sp_names <- names(species_data)
  for (i in seq_len(length(sp_names) - 1)) {
    for (j in (i + 1):length(sp_names)) {
      sp1 <- sp_names[i]
      sp2 <- sp_names[j]
      pair_name <- paste0(sp1, "_", sp2)

      orthologs <- download_orthologs_biomart(
        SPECIES[[sp1]]$ensembl,
        SPECIES[[sp2]]$ensembl
      )

      if (is.null(orthologs) && USE_SIMULATED_FALLBACK) {
        message(sprintf("  Using simulated orthologs for %s", pair_name))
        orthologs <- simulate_orthologs(
          rownames(species_data[[sp1]]),
          rownames(species_data[[sp2]]),
          prop_1to1 = 0.7,
          seed = i * 1000 + j
        )
      }

      if (!is.null(orthologs)) {
        orthologs_data[[pair_name]] <- orthologs
      }
    }
  }

  # Run benchmark
  results <- run_benchmark(species_data, orthologs_data)

  # Generate plots
  summary_stats <- generate_plots(results, OUTPUT_DIR)

  # Validate results
  validation <- validate_results(results)

  # Print summary
  message("\n=== Summary Statistics ===")
  print(summary_stats)

  message("\n=== Benchmark Complete ===")
  message(sprintf("Results saved to: %s", OUTPUT_DIR))

  # Return results for interactive use
  invisible(list(
    results = results,
    summary = summary_stats,
    validation = validation
  ))
}

# Run if executed as script
if (!interactive()) {
  main()
}
