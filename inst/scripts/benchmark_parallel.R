#!/usr/bin/env Rscript
#' Benchmark script for parallelization performance
#'
#' This script benchmarks the parallel vs sequential implementations of:
#' 1. calculate_pcc_mr() - similarity matrix calculation
#' 2. calculate_ccs() - co-expression correlation score calculation
#'
#' Usage:
#'   Rscript inst/scripts/benchmark_parallel.R
#'
#' Or source interactively:
#'   source("inst/scripts/benchmark_parallel.R")

library(coexpr)
library(microbenchmark)

# --- Configuration ---
N_GENES <- 5000    # Number of genes for similarity benchmark
N_SAMPLES <- 40    # Number of samples
N_ORTHOLOGS <- 10000  # Number of ortholog pairs for CCS benchmark
N_CORES <- parallel::detectCores() - 1
N_REPS <- 3        # Number of benchmark repetitions

message("=== Parallelization Benchmark ===")
message(sprintf("Configuration:"))
message(sprintf("  - Genes: %d", N_GENES))
message(sprintf("  - Samples: %d", N_SAMPLES))
message(sprintf("  - Ortholog pairs: %d", N_ORTHOLOGS))
message(sprintf("  - Cores available: %d", N_CORES))
message(sprintf("  - Benchmark repetitions: %d", N_REPS))
message("")

# --- Generate test data ---
message("Generating test data...")
set.seed(42)

# Expression matrix
expr <- matrix(rnorm(N_GENES * N_SAMPLES), nrow = N_GENES, ncol = N_SAMPLES)
rownames(expr) <- paste0("Gene", seq_len(N_GENES))
colnames(expr) <- paste0("Sample", seq_len(N_SAMPLES))

# Pre-calculate similarity matrices for CCS benchmark
message("Pre-calculating similarity matrices for CCS benchmark...")
sim <- calculate_pcc_mr(expr)

# Generate fake ortholog pairs (gene pairs for CCS calculation)
orthologs <- data.frame(
  gene_sp1 = sample(rownames(sim), N_ORTHOLOGS, replace = TRUE),
  gene_sp2 = sample(rownames(sim), N_ORTHOLOGS, replace = TRUE),
  type = "1:1",
  stringsAsFactors = FALSE
)

message("")

# --- Benchmark 1: Similarity Matrix Calculation ---
message("=== Benchmark 1: Similarity Matrix (calculate_pcc_mr) ===")
message(sprintf("Matrix size: %d x %d genes", N_GENES, N_GENES))

if (N_GENES >= 5000) {
  bm_sim <- microbenchmark(
    sequential = calculate_pcc_mr(expr, n_cores = 1),
    parallel = calculate_pcc_mr(expr, n_cores = N_CORES),
    times = N_REPS
  )

  print(bm_sim)

  # Calculate speedup
  median_seq <- median(bm_sim$time[bm_sim$expr == "sequential"])
  median_par <- median(bm_sim$time[bm_sim$expr == "parallel"])
  speedup_sim <- median_seq / median_par

  message(sprintf("\nSpeedup (parallel vs sequential): %.2fx", speedup_sim))
  message("")
} else {
  message("Skipping parallel benchmark for similarity (N_GENES < 5000)")
  message("Parallel mode only activates for matrices with >5000 genes")
  message("")
}

# --- Benchmark 2: CCS Calculation ---
message("=== Benchmark 2: CCS Calculation (calculate_ccs) ===")
message(sprintf("Ortholog pairs: %d", N_ORTHOLOGS))

bm_ccs <- microbenchmark(
  sequential = calculate_ccs(sim, sim, orthologs, n_cores = 1),
  parallel = calculate_ccs(sim, sim, orthologs, n_cores = N_CORES),
  times = N_REPS
)

print(bm_ccs)

# Calculate speedup
median_seq_ccs <- median(bm_ccs$time[bm_ccs$expr == "sequential"])
median_par_ccs <- median(bm_ccs$time[bm_ccs$expr == "parallel"])
speedup_ccs <- median_seq_ccs / median_par_ccs

message(sprintf("\nSpeedup (parallel vs sequential): %.2fx", speedup_ccs))
message("")

# --- Summary ---
message("=== Summary ===")
message(sprintf("CCS calculation speedup: %.2fx with %d cores", speedup_ccs, N_CORES))
if (exists("speedup_sim")) {
  message(sprintf("Similarity calculation speedup: %.2fx with %d cores", speedup_sim, N_CORES))
}
message("")
message("Note: Actual speedup depends on:")
message("  - Number of CPU cores")
message("  - Data size (larger = better parallel efficiency)")
message("  - Memory bandwidth")
message("  - Overhead from data transfer to workers")
