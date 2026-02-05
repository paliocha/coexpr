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
  expect_true("ccs" %in% colnames(ccs_results))
  expect_true("n_ref" %in% colnames(ccs_results))
  expect_equal(nrow(ccs_results), nrow(orthologs))

  # Check CCS values are in valid range
  expect_true(all(ccs_results$ccs >= -1 & ccs_results$ccs <= 1))

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
  expect_equal(ccs_seq$ccs, ccs_par$ccs, tolerance = 1e-10)
  expect_equal(ccs_seq$n_ref, ccs_par$n_ref)
  expect_equal(nrow(ccs_seq), nrow(ccs_par))
})
