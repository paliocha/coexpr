test_that("calculate_ccs works with basic input", {
  # Create simple test similarity matrices
  set.seed(123)
  n_genes <- 20

  # Similarity matrices for two species
  sim_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp1 <- (sim_sp1 + t(sim_sp1)) / 2  # Make symmetric
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  sim_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp2 <- (sim_sp2 + t(sim_sp2)) / 2
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Create 1:1 orthologs
  orthologs <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_genes),
    gene_sp2 = paste0("Gene_sp2_", 1:n_genes),
    type = "1:1"
  )

  # Calculate CCS
  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)

  # Check output structure
  expect_true(is.data.frame(ccs_results))
  expect_true("CCS" %in% colnames(ccs_results))
  expect_true("n_ref" %in% colnames(ccs_results))
  expect_equal(nrow(ccs_results), nrow(orthologs))

  # Check CCS values are in valid range
  expect_true(all(ccs_results$CCS >= -1 & ccs_results$CCS <= 1))

  # Check reference count is correct
  expect_true(all(ccs_results$n_ref == n_genes))
})

test_that("calculate_ccs filters reference orthologs correctly", {
  set.seed(456)
  n_genes <- 15

  # Create similarity matrices
  sim_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp1 <- (sim_sp1 + t(sim_sp1)) / 2
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  sim_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp2 <- (sim_sp2 + t(sim_sp2)) / 2
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Create mixed orthologs (1:1 and 1:N)
  orthologs <- data.frame(
    gene_sp1 = c(paste0("Gene_sp1_", 1:10), "Gene_sp1_11", "Gene_sp1_11"),
    gene_sp2 = c(paste0("Gene_sp2_", 1:10), "Gene_sp2_11", "Gene_sp2_12"),
    type = c(rep("1:1", 10), "1:N", "1:N")
  )

  # Test with use_only_1to1 = TRUE (default)
  ccs_1to1 <- calculate_ccs(sim_sp1, sim_sp2, orthologs, use_only_1to1 = TRUE)

  # Should use only 10 1:1 orthologs as reference
  expect_true(all(ccs_1to1$n_ref == 10))

  # Test with use_only_1to1 = FALSE
  ccs_all <- calculate_ccs(sim_sp1, sim_sp2, orthologs, use_only_1to1 = FALSE)

  # Should use all 12 orthologs as reference
  expect_true(all(ccs_all$n_ref == 12))
})

test_that("calculate_ccs validates inputs", {
  # Test non-matrix similarity
  expect_error(
    calculate_ccs(data.frame(a = 1:10), matrix(1:100, 10, 10),
                  data.frame(gene_sp1 = "A", gene_sp2 = "B")),
    "Similarity matrices must be matrices"
  )

  # Test missing required columns
  sim <- matrix(runif(100), 10, 10)
  rownames(sim) <- colnames(sim) <- paste0("G", 1:10)

  expect_error(
    calculate_ccs(sim, sim, data.frame(wrong_col1 = "A", wrong_col2 = "B")),
    "orthologs must have columns: gene_sp1, gene_sp2"
  )

  # Test insufficient reference orthologs
  orthologs_few <- data.frame(
    gene_sp1 = paste0("G", 1:5),
    gene_sp2 = paste0("G", 1:5),
    type = "1:1"
  )

  expect_error(
    calculate_ccs(sim, sim, orthologs_few),
    "Need at least 10 reference orthologs"
  )
})

test_that("calculate_ccs handles missing genes", {
  set.seed(789)
  n_genes <- 15

  # Create similarity matrices
  sim_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp1 <- (sim_sp1 + t(sim_sp1)) / 2
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  sim_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp2 <- (sim_sp2 + t(sim_sp2)) / 2
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Create orthologs with some genes not in similarity matrices
  orthologs <- data.frame(
    gene_sp1 = c(paste0("Gene_sp1_", 1:10), "Gene_sp1_MISSING1", "Gene_sp1_MISSING2"),
    gene_sp2 = c(paste0("Gene_sp2_", 1:10), "Gene_sp2_MISSING1", "Gene_sp2_MISSING2"),
    type = "1:1"
  )

  # Calculate CCS - should filter out missing genes
  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)

  # Should only keep genes present in both matrices
  expect_equal(nrow(ccs_results), 10)
  expect_true(all(ccs_results$gene_sp1 %in% rownames(sim_sp1)))
  expect_true(all(ccs_results$gene_sp2 %in% rownames(sim_sp2)))
})

test_that("calculate_ccs with auto-detected 1:1 orthologs", {
  set.seed(111)
  n_genes <- 15

  # Create similarity matrices
  sim_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp1 <- (sim_sp1 + t(sim_sp1)) / 2
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  sim_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp2 <- (sim_sp2 + t(sim_sp2)) / 2
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Orthologs WITHOUT type column
  orthologs_no_type <- data.frame(
    gene_sp1 = c(paste0("Gene_sp1_", 1:10), "Gene_sp1_11", "Gene_sp1_11"),
    gene_sp2 = c(paste0("Gene_sp2_", 1:10), "Gene_sp2_11", "Gene_sp2_12")
  )

  # Should auto-detect 1:1 orthologs
  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs_no_type, use_only_1to1 = TRUE)

  # Should have detected 10 1:1 orthologs as reference
  expect_true(all(ccs_results$n_ref == 10))
  expect_equal(nrow(ccs_results), 12)  # All 12 pairs get CCS, but only 10 used as ref
})

test_that("calculate_ccs parallel produces same results as sequential", {
  skip_on_cran()  # Parallel tests can be slow on CRAN

  set.seed(999)
  n_genes <- 50

  # Similarity matrices for two species
  sim_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp1 <- (sim_sp1 + t(sim_sp1)) / 2
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  sim_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp2 <- (sim_sp2 + t(sim_sp2)) / 2
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Create orthologs
  orthologs <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_genes),
    gene_sp2 = paste0("Gene_sp2_", 1:n_genes),
    type = "1:1"
  )

  # Calculate CCS sequentially
  ccs_seq <- calculate_ccs(sim_sp1, sim_sp2, orthologs, n_cores = 1)

  # Calculate CCS in parallel (with 2 cores for testing)
  ccs_par <- calculate_ccs(sim_sp1, sim_sp2, orthologs, n_cores = 2)

  # Results should be identical (or very close due to floating point)
  expect_equal(ccs_seq$CCS, ccs_par$CCS, tolerance = 1e-10)
  expect_equal(ccs_seq$n_ref, ccs_par$n_ref)
  expect_equal(nrow(ccs_seq), nrow(ccs_par))
})


# Tests for self-diagonal handling

test_that("handle_self_diagonal methods produce valid CCS values", {
  set.seed(1001)
  n_genes <- 25

  # Create similarity matrices
  sim_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp1 <- (sim_sp1 + t(sim_sp1)) / 2
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  sim_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp2 <- (sim_sp2 + t(sim_sp2)) / 2
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Create orthologs
  orthologs <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_genes),
    gene_sp2 = paste0("Gene_sp2_", 1:n_genes),
    type = "1:1"
  )

  # Test all three methods
  ccs_mean <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "mean")
  ccs_na <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "na")
  ccs_none <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "none")

  # All should return valid CCS values in range [-1, 1]
  expect_true(all(ccs_mean$CCS >= -1 & ccs_mean$CCS <= 1))
  expect_true(all(ccs_na$CCS >= -1 & ccs_na$CCS <= 1))
  expect_true(all(ccs_none$CCS >= -1 & ccs_none$CCS <= 1))
})

test_that("handle_self_diagonal 'mean' and 'na' give lower CCS than 'none' for reference genes", {
  # This test verifies that diagonal handling reduces CCS for genes in the reference set
  # (where the self-correlation = 1 artificially inflates the correlation)

  set.seed(1002)
  n_genes <- 20

  # Create similarity matrices with low off-diagonal correlations
  # This makes the diagonal inflation more pronounced
  sim_sp1 <- matrix(runif(n_genes * n_genes, 0.1, 0.3), nrow = n_genes)
  sim_sp1 <- (sim_sp1 + t(sim_sp1)) / 2
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  sim_sp2 <- matrix(runif(n_genes * n_genes, 0.1, 0.3), nrow = n_genes)
  sim_sp2 <- (sim_sp2 + t(sim_sp2)) / 2
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Create orthologs (all 1:1, so all genes are in reference set)
  orthologs <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_genes),
    gene_sp2 = paste0("Gene_sp2_", 1:n_genes),
    type = "1:1"
  )

  ccs_mean <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "mean")
  ccs_na <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "na")
  ccs_none <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "none")

  # Mean CCS with 'mean' or 'na' handling should be lower than 'none'
  # because the self-correlation (1.0, 1.0) pair is removed/replaced
  expect_lt(mean(ccs_mean$CCS), mean(ccs_none$CCS) + 0.01)  # Allow small tolerance
  expect_lt(mean(ccs_na$CCS), mean(ccs_none$CCS) + 0.01)
})

test_that("handle_self_diagonal does not affect non-reference genes", {
  set.seed(1003)
  n_genes <- 20

  # Create similarity matrices
  sim_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp1 <- (sim_sp1 + t(sim_sp1)) / 2
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  sim_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  sim_sp2 <- (sim_sp2 + t(sim_sp2)) / 2
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Create orthologs: genes 1-15 are 1:1 (reference), genes 16-20 are 1:N (not reference)
  orthologs <- data.frame(
    gene_sp1 = c(paste0("Gene_sp1_", 1:15),
                 "Gene_sp1_16", "Gene_sp1_16",  # 1:N
                 "Gene_sp1_17", "Gene_sp1_17",
                 "Gene_sp1_18"),
    gene_sp2 = c(paste0("Gene_sp2_", 1:15),
                 "Gene_sp2_16", "Gene_sp2_17",
                 "Gene_sp2_18", "Gene_sp2_19",
                 "Gene_sp2_20"),
    type = c(rep("1:1", 15), rep("1:N", 5))
  )

  ccs_mean <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "mean")
  ccs_none <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "none")

  # For non-reference genes (1:N type), CCS should be similar regardless of diagonal handling
  # because they are not in the reference set, so no diagonal to handle
  non_ref_idx <- which(ccs_mean$type == "1:N")
  ref_idx <- which(ccs_mean$type == "1:1")

  # For 1:N genes, results should be very close (only floating point differences)
  # because diagonal handling only affects genes in the reference set
  # Note: The genes 16-17 from sp1 ARE in the matrix, but not in 1:1 reference set
  expect_equal(ccs_mean$CCS[non_ref_idx], ccs_none$CCS[non_ref_idx], tolerance = 1e-10)
})

test_that("handle_self_diagonal works with tri_similarity objects", {
  set.seed(1004)
  n_genes <- 20

  # Create similarity matrices
  mat_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  mat_sp1 <- (mat_sp1 + t(mat_sp1)) / 2
  diag(mat_sp1) <- 1
  rownames(mat_sp1) <- colnames(mat_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  mat_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  mat_sp2 <- (mat_sp2 + t(mat_sp2)) / 2
  diag(mat_sp2) <- 1
  rownames(mat_sp2) <- colnames(mat_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Convert to tri_similarity
  sim_sp1 <- as.TriSimilarity(mat_sp1)
  sim_sp2 <- as.TriSimilarity(mat_sp2)

  orthologs <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_genes),
    gene_sp2 = paste0("Gene_sp2_", 1:n_genes),
    type = "1:1"
  )

  # All methods should work with tri_similarity
  ccs_mean <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "mean")
  ccs_na <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "na")
  ccs_none <- calculate_ccs(sim_sp1, sim_sp2, orthologs, handle_self_diagonal = "none")

  expect_true(all(is.finite(ccs_mean$CCS)))
  expect_true(all(is.finite(ccs_na$CCS)))
  expect_true(all(is.finite(ccs_none$CCS)))
})

test_that("handle_self_diagonal results match between matrix and tri_similarity", {
  set.seed(1005)
  n_genes <- 15

  # Create similarity matrices
  mat_sp1 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  mat_sp1 <- (mat_sp1 + t(mat_sp1)) / 2
  diag(mat_sp1) <- 1
  rownames(mat_sp1) <- colnames(mat_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  mat_sp2 <- matrix(runif(n_genes * n_genes), nrow = n_genes)
  mat_sp2 <- (mat_sp2 + t(mat_sp2)) / 2
  diag(mat_sp2) <- 1
  rownames(mat_sp2) <- colnames(mat_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Convert to tri_similarity
  tri_sp1 <- as.TriSimilarity(mat_sp1)
  tri_sp2 <- as.TriSimilarity(mat_sp2)

  orthologs <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_genes),
    gene_sp2 = paste0("Gene_sp2_", 1:n_genes),
    type = "1:1"
  )

  # Compare results for all methods
  for (method in c("mean", "na", "none")) {
    ccs_mat <- calculate_ccs(mat_sp1, mat_sp2, orthologs, handle_self_diagonal = method)
    ccs_tri <- calculate_ccs(tri_sp1, tri_sp2, orthologs, handle_self_diagonal = method)

    expect_equal(ccs_mat$CCS, ccs_tri$CCS, tolerance = 1e-10,
                 info = sprintf("Method '%s' results differ", method))
  }
})
