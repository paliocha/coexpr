test_that("calculate_pcc_mr works with basic input", {
  # Create simple test data
  set.seed(123)
  expr <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:100)

  # Calculate similarity with PCC only
  sim <- calculate_pcc_mr(expr, method = "pcc")

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

  # Calculate with mutual rank
  sim_mr <- calculate_pcc_mr(expr, method = "pcc_mr")

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
    calculate_pcc_mr(expr),
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
  sim_seq <- calculate_pcc_mr(expr, method = "pcc", n_cores = 1)

  # Calculate similarity in parallel (with 2 cores for testing)
  # Suppress package version warnings from future/furrr workers
  sim_par <- suppressWarnings(calculate_pcc_mr(expr, method = "pcc", n_cores = 2))

  # Results should be identical (or very close due to floating point)
  expect_equal(sim_seq, sim_par, tolerance = 1e-10)
})
