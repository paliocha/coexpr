test_that("calculate_pcc_mr works with basic input", {
  # Create simple test data
  set.seed(123)
  expr <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:100)

  # Calculate similarity with PCC only (return_tri = FALSE for legacy tests)
  sim <- calculate_pcc_mr(expr, method = "pcc", return_tri = FALSE)

  # Check properties
  expect_true(is.matrix(sim))
  expect_equal(nrow(sim), 100)
  expect_equal(ncol(sim), 100)
  expect_true(all(sim >= -1 & sim <= 1))  # Correlation bounds
  expect_equal(as.vector(diag(sim)), rep(1, 100))  # Self-correlation = 1 (ignore names)
})

test_that("calculate_pcc_mr handles mutual rank correctly", {
  set.seed(123)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  # Calculate with mutual rank (return_tri = FALSE for legacy tests)
  sim_mr <- calculate_pcc_mr(expr, method = "pcc_mr", return_tri = FALSE)

  # Check properties
  expect_true(is.matrix(sim_mr))
  expect_true(all(sim_mr >= 0 & sim_mr <= 1))  # Should be [0,1] after MR normalization
  expect_equal(as.vector(diag(sim_mr)), rep(1, 50))  # Diagonal = 1 (ignore names)
  expect_true(isSymmetric(sim_mr))  # Should be symmetric
})

test_that("calculate_pcc_mr validates input correctly", {
  # Test non-matrix input
  expect_error(
    calculate_pcc_mr(data.frame(a = 1:10)),
    "must be a matrix"
  )

  # Test matrix without rownames
  expr <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  expect_warning(
    calculate_pcc_mr(expr, return_tri = FALSE),
    NA  # Should not warn if rownames missing (will be auto-generated)
  )
})

test_that("calculate_pcc_mr parallel produces same results as sequential", {
  skip_on_cran()  # Parallel tests can be slow on CRAN

  # Use a large enough matrix to trigger parallel mode (>5000 genes)
  # For testing, we'll use a smaller matrix and test the chunking logic
  set.seed(42)
  n_genes <- 6000  # Just above the 5000 threshold
  n_samples <- 10

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene", 1:n_genes)

  # Calculate similarity sequentially
  sim_seq <- calculate_pcc_mr(expr, method = "pcc", n_cores = 1, return_tri = FALSE)

  # Calculate similarity in parallel (with 2 cores for testing)
  # Suppress package version warnings from future/furrr workers
  sim_par <- suppressWarnings(calculate_pcc_mr(expr, method = "pcc", n_cores = 2, return_tri = FALSE))

  # Results should be identical (or very close due to floating point)
  expect_equal(sim_seq, sim_par, tolerance = 1e-10)
})


# MI+CLR tests

# MI+CLR tests with C++ implementation

test_that("calculate_mi_clr returns valid TriSimilarity", {
  skip_on_cran()

  set.seed(123)
  n_genes <- 30
  n_samples <- 50  # Need enough samples for MI estimation

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene", 1:n_genes)

  sim <- calculate_mi_clr(expr, n_bins = 5)

  # Should return TriSimilarity (S4)
  expect_s4_class(sim, "TriSimilarity")
  expect_equal(sim@n, n_genes)
  expect_equal(sim@genes, paste0("Gene", 1:n_genes))

  # CLR values are non-negative
  expect_true(all(sim@data >= 0))
})

test_that("calculate_mi_clr returns matrix when return_tri = FALSE", {
  skip_on_cran()

  set.seed(456)
  n_genes <- 20
  n_samples <- 40

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("G", 1:n_genes)

  sim <- calculate_mi_clr(expr, n_bins = 5, return_tri = FALSE)

  expect_true(is.matrix(sim))
  expect_equal(dim(sim), c(n_genes, n_genes))
  expect_equal(rownames(sim), paste0("G", 1:n_genes))
  expect_equal(colnames(sim), paste0("G", 1:n_genes))
})

test_that("calculate_mi_clr validates input", {
  # Non-matrix input
  expect_error(
    calculate_mi_clr(data.frame(a = 1:10)),
    "must be a matrix"
  )
})

test_that("calculate_mi_clr handles zero variance genes", {
  skip_on_cran()

  set.seed(111)
  n_genes <- 20
  n_samples <- 40

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene", 1:n_genes)

  # Make one gene have zero variance
  expr[5, ] <- 0

  # Should warn about zero variance
  expect_warning(
    sim <- calculate_mi_clr(expr, n_bins = 5),
    "zero variance"
  )

  # Should still return valid result
  expect_s4_class(sim, "TriSimilarity")
})

test_that("calculate_mi_clr warns about few samples", {
  set.seed(222)
  n_genes <- 10
  n_samples <- 15  # Less than 3 * n_bins for n_bins = 10

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene", 1:n_genes)

  # Should warn about insufficient samples
  expect_warning(
    calculate_mi_clr(expr, n_bins = 10),
    "works best with at least"
  )
})

test_that("calculate_mi_clr CLR transformation produces symmetric matrix", {
  skip_on_cran()

  set.seed(333)
  n_genes <- 25
  n_samples <- 50

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene", 1:n_genes)

  sim <- calculate_mi_clr(expr, n_bins = 5, return_tri = FALSE)

  # CLR should be symmetric
  expect_true(isSymmetric(sim, tol = 1e-10))
})

test_that("calculate_mi_clr works with CCS calculation", {
  skip_on_cran()

  set.seed(444)
  n_genes <- 20
  n_samples <- 40

  # Create expression for two species
  expr_sp1 <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  expr_sp2 <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Calculate MI+CLR similarity
  sim_sp1 <- calculate_mi_clr(expr_sp1, n_bins = 5)
  sim_sp2 <- calculate_mi_clr(expr_sp2, n_bins = 5)

  # Create orthologs
  orthologs <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_genes),
    gene_sp2 = paste0("Gene_sp2_", 1:n_genes),
    type = "1:1"
  )

  # CCS should work with MI+CLR similarity matrices
  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)

  expect_true(is.data.frame(ccs_results))
  expect_equal(nrow(ccs_results), n_genes)
  expect_true(all(ccs_results$CCS >= -1 & ccs_results$CCS <= 1))
})

test_that("calculate_mi_clr parallel produces consistent results", {
  skip_on_cran()

  set.seed(555)
  n_genes <- 30
  n_samples <- 50

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene", 1:n_genes)

  # Sequential
  sim_seq <- calculate_mi_clr(expr, n_bins = 5, n_cores = 1, return_tri = FALSE)

  # Parallel (if OpenMP available)
  if (has_openmp()) {
    sim_par <- calculate_mi_clr(expr, n_bins = 5, n_cores = 2, return_tri = FALSE)
    expect_equal(sim_seq, sim_par, tolerance = 1e-10)
  }
})


# cor_method tests

test_that("Spearman correlation produces valid output", {
  set.seed(700)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  sim <- calculate_pcc_mr(expr, method = "pcc", cor_method = "spearman",
                          return_tri = FALSE)

  expect_true(is.matrix(sim))
  expect_equal(dim(sim), c(50, 50))
  expect_true(all(sim >= -1 & sim <= 1))
  expect_true(isSymmetric(sim, tol = 1e-10))
  expect_equal(as.vector(diag(sim)), rep(1, 50))
})

test_that("Spearman differs from Pearson on non-normal data", {
  set.seed(701)
  # Create skewed data where Spearman and Pearson should differ
  expr <- matrix(rexp(50 * 20, rate = 0.5), nrow = 50, ncol = 20)
  rownames(expr) <- paste0("Gene", 1:50)

  sim_pearson <- calculate_pcc_mr(expr, method = "pcc", cor_method = "pearson",
                                  return_tri = FALSE)
  sim_spearman <- calculate_pcc_mr(expr, method = "pcc", cor_method = "spearman",
                                   return_tri = FALSE)

  # They should not be identical

  expect_false(identical(sim_pearson, sim_spearman))
  # But both should be valid correlation matrices
  expect_true(all(sim_pearson >= -1 & sim_pearson <= 1))
  expect_true(all(sim_spearman >= -1 & sim_spearman <= 1))
})

test_that("Spearman + MR produces values in [0,1]", {
  set.seed(702)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  sim <- calculate_pcc_mr(expr, method = "pcc_mr", cor_method = "spearman",
                          return_tri = FALSE)

  expect_true(all(sim >= 0 & sim <= 1))
  expect_equal(as.vector(diag(sim)), rep(1, 50))
  expect_true(isSymmetric(sim))
})

test_that("Kendall correlation produces valid output", {
  set.seed(703)
  # Small matrix since Kendall is slow
  expr <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:20)

  sim <- calculate_pcc_mr(expr, method = "pcc", cor_method = "kendall",
                          return_tri = FALSE)

  expect_true(is.matrix(sim))
  expect_equal(dim(sim), c(20, 20))
  expect_true(all(sim >= -1 & sim <= 1))
  expect_true(isSymmetric(sim, tol = 1e-10))
  expect_equal(as.vector(diag(sim)), rep(1, 20))
})

test_that("Kendall warns for large matrices", {
  set.seed(704)
  expr <- matrix(rnorm(501 * 5), nrow = 501, ncol = 5)
  rownames(expr) <- paste0("Gene", 1:501)

  expect_warning(
    calculate_pcc_mr(expr, method = "pcc", cor_method = "kendall",
                     return_tri = FALSE),
    "Kendall correlation is slow"
  )
})

test_that("default cor_method produces identical results to explicit pearson", {
  set.seed(705)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  sim_default <- calculate_pcc_mr(expr, method = "pcc", return_tri = FALSE)
  sim_pearson <- calculate_pcc_mr(expr, method = "pcc", cor_method = "pearson",
                                  return_tri = FALSE)

  expect_identical(sim_default, sim_pearson)
})
