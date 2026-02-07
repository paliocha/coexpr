# Tests for C++ mutual rank implementation

test_that("has_openmp returns boolean", {
  result <- has_openmp()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("get_max_threads returns positive integer", {
  result <- get_max_threads()
  expect_type(result, "integer")
  expect_gte(result, 1)
})

test_that("compute_ranks_desc_cpp computes descending ranks correctly", {
  # Simple case: values in ascending order should get descending ranks
  x <- c(1, 2, 3, 4, 5)
  ranks <- as.vector(compute_ranks_desc_cpp(x))
  expect_equal(ranks, c(5, 4, 3, 2, 1))

  # Values in descending order should get ascending ranks
  x <- c(5, 4, 3, 2, 1)
  ranks <- as.vector(compute_ranks_desc_cpp(x))
  expect_equal(ranks, c(1, 2, 3, 4, 5))

  # Ties should get average ranks
  x <- c(1, 2, 2, 3)
  ranks <- as.vector(compute_ranks_desc_cpp(x))
  expect_equal(ranks, c(4, 2.5, 2.5, 1))
})

test_that("mutual_rank_transform_cached_cpp produces valid output", {
  # Create small symmetric correlation matrix
  set.seed(123)
  n <- 50
  pcc <- matrix(rnorm(n * n), n, n)
  pcc <- (pcc + t(pcc)) / 2  # Make symmetric
  diag(pcc) <- 1  # Diagonal = 1

  result <- mutual_rank_transform_cached_cpp(pcc, n_cores = 1)

  # Check dimensions
  expect_equal(dim(result), c(n, n))

  # Check symmetry
  expect_equal(result, t(result))


  # Check diagonal is 1
  expect_equal(diag(result), rep(1, n))

  # Check values in [0, 1]
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("mutual_rank_transform_cpp (streaming) produces valid output", {
  set.seed(123)
  n <- 50
  pcc <- matrix(rnorm(n * n), n, n)
  pcc <- (pcc + t(pcc)) / 2
  diag(pcc) <- 1

  result <- mutual_rank_transform_cpp(pcc, n_cores = 1)

  expect_equal(dim(result), c(n, n))
  expect_equal(result, t(result))
  expect_equal(diag(result), rep(1, n))
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("cached and streaming methods produce identical results", {
  set.seed(456)
  n <- 100
  pcc <- matrix(rnorm(n * n), n, n)
  pcc <- (pcc + t(pcc)) / 2
  diag(pcc) <- 1

  result_cached <- mutual_rank_transform_cached_cpp(pcc, n_cores = 1)
  result_streaming <- mutual_rank_transform_cpp(pcc, n_cores = 1)

  expect_equal(result_cached, result_streaming, tolerance = 1e-10)
})

test_that("calculate_pcc_mr with mr_method parameter works", {
  set.seed(789)
  n_genes <- 100
  n_samples <- 20
  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(expr) <- paste0("Gene", 1:n_genes)

  # Test cached method (return_tri = FALSE for legacy tests)
  sim_cached <- calculate_pcc_mr(expr, mr_method = "cached", return_tri = FALSE)
  expect_equal(dim(sim_cached), c(n_genes, n_genes))
  expect_true(all(diag(sim_cached) == 1))

  # Test streaming method
  sim_streaming <- calculate_pcc_mr(expr, mr_method = "streaming", return_tri = FALSE)
  expect_equal(dim(sim_streaming), c(n_genes, n_genes))
  expect_true(all(diag(sim_streaming) == 1))

  # Both should produce same results
  expect_equal(sim_cached, sim_streaming, tolerance = 1e-10)
})

test_that("mutual rank handles edge cases", {
  # 2x2 matrix (minimum size)
  pcc <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  result <- mutual_rank_transform_cached_cpp(pcc, n_cores = 1)

  expect_equal(dim(result), c(2, 2))
  expect_equal(diag(result), c(1, 1))
  expect_equal(result, t(result))

  # Perfect correlation (all same value except diagonal)
  n <- 10
  pcc <- matrix(0.9, n, n)
  diag(pcc) <- 1
  result <- mutual_rank_transform_cached_cpp(pcc, n_cores = 1)

  expect_equal(dim(result), c(n, n))
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})
