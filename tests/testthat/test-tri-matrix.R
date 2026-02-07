# Tests for S4 TriSimilarity class

test_that("as.TriSimilarity converts full matrix correctly", {
  # Create a symmetric test matrix
  set.seed(123)
  n <- 10
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2  # Make symmetric
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("Gene", 1:n)

  # Convert to triangular
  tri <- as.TriSimilarity(mat)

  expect_s4_class(tri, "TriSimilarity")
  expect_equal(tri@n, n)
  expect_equal(tri@genes, paste0("Gene", 1:n))
  expect_equal(tri@diag_value, 1)
  expect_equal(length(tri@data), n * (n - 1) / 2)
})

test_that("as.matrix.TriSimilarity reconstructs full matrix", {
  set.seed(456)
  n <- 15
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("G", 1:n)

  # Round-trip conversion
  tri <- as.TriSimilarity(mat)
  reconstructed <- as.matrix(tri)

  expect_equal(reconstructed, mat, tolerance = 1e-10)
})

test_that("TriSimilarity subsetting works correctly", {
  set.seed(789)
  n <- 8
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("Gene", 1:n)

  tri <- as.TriSimilarity(mat)

  # Single element access by index
  expect_equal(tri[1, 1], 1)  # Diagonal
  expect_equal(tri[1, 2], mat[1, 2])
  expect_equal(tri[2, 1], mat[2, 1])  # Lower triangle (should mirror)
  expect_equal(tri[3, 5], mat[3, 5])

  # Single element access by name
  expect_equal(tri["Gene1", "Gene1"], 1)
  expect_equal(tri["Gene2", "Gene5"], mat["Gene2", "Gene5"])

  # Multiple elements
  sub <- tri[1:3, 1:3]
  expect_equal(sub, mat[1:3, 1:3], tolerance = 1e-10)
})

test_that("extractColumn.TriSimilarity matches full matrix column", {
  set.seed(111)
  n <- 20
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("Gene", 1:n)

  tri <- as.TriSimilarity(mat)

  # Test column extraction for various positions
  for (j in c(1, 5, 10, 15, 20)) {
    col_tri <- extractColumn(tri, j)
    col_mat <- mat[, j]
    expect_equal(col_tri, col_mat, tolerance = 1e-10,
                 info = sprintf("Column %d extraction failed", j))
  }

  # Test by name
  col_tri <- extractColumn(tri, "Gene7")
  col_mat <- mat[, "Gene7"]
  expect_equal(col_tri, col_mat, tolerance = 1e-10)
})

test_that("extractColumns.TriSimilarity returns correct matrix", {
  set.seed(222)
  n <- 12
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("G", 1:n)

  tri <- as.TriSimilarity(mat)

  # Extract multiple columns
  cols <- extractColumns(tri, c("G2", "G5", "G10"))
  expected <- mat[, c("G2", "G5", "G10")]

  expect_equal(cols, expected, tolerance = 1e-10)
})

test_that("extractRows.TriSimilarity returns correct matrix", {
  set.seed(333)
  n <- 10
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("G", 1:n)

  tri <- as.TriSimilarity(mat)

  # Extract multiple rows
  rows <- extractRows(tri, c(1, 3, 7))
  expected <- mat[c(1, 3, 7), ]

  expect_equal(rows, expected, tolerance = 1e-10)
})

test_that("TriSimilarity memory usage is approximately 50%", {
  # Create a larger matrix to get meaningful size comparison
  set.seed(444)
  n <- 100
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("Gene", 1:n)

  tri <- as.TriSimilarity(mat)

  # The triangular data should be about 50% of the full matrix
  data_ratio <- length(tri@data) / (n * n)
  expect_lt(data_ratio, 0.51)
  expect_gt(data_ratio, 0.49)
})

test_that("dim and dimnames work for TriSimilarity", {
  set.seed(555)
  n <- 5
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- letters[1:n]

  tri <- as.TriSimilarity(mat)

  expect_equal(dim(tri), c(n, n))
  expect_equal(dimnames(tri), list(letters[1:n], letters[1:n]))
})

test_that("is.TriSimilarity identifies objects correctly", {
  set.seed(666)
  mat <- matrix(runif(25), 5, 5)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1

  tri <- as.TriSimilarity(mat)

  expect_true(is.TriSimilarity(tri))
  expect_false(is.TriSimilarity(mat))
  expect_false(is.TriSimilarity(list(data = 1:10)))
})

test_that("save and load TriSimilarity works", {
  skip_on_cran()

  set.seed(777)
  n <- 10
  mat <- matrix(runif(n * n), n, n)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("G", 1:n)

  tri <- as.TriSimilarity(mat)

  # Save and load
  tmp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_file), add = TRUE)

  saveTriSimilarity(tri, tmp_file)
  loaded <- loadTriSimilarity(tmp_file)

  expect_s4_class(loaded, "TriSimilarity")
  expect_equal(loaded@data, tri@data)
  expect_equal(loaded@genes, tri@genes)
  expect_equal(loaded@n, tri@n)
  expect_equal(loaded@diag_value, tri@diag_value)
})

test_that("TriSimilarity constructor validates input length", {
  # Correct: 5 genes need 5*4/2 = 10 values
  expect_error(
    TriSimilarity(1:10, paste0("G", 1:5)),
    NA
  )

  # Wrong: 5 genes but only 5 values
  expect_error(
    TriSimilarity(1:5, paste0("G", 1:5)),
    "data length"
  )
})

test_that("show.TriSimilarity displays info without error", {
  set.seed(888)
  mat <- matrix(runif(25), 5, 5)
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- paste0("G", 1:5)

  tri <- as.TriSimilarity(mat)

  # Should show without error
  expect_output(show(tri), "TriSimilarity")
  expect_output(show(tri), "5 x 5 genes")
})

test_that("CCS calculation produces same results with TriSimilarity vs matrix", {
  skip_on_cran()

  set.seed(999)
  n_genes <- 30

  # Create two similarity matrices
  mat_sp1 <- matrix(runif(n_genes * n_genes), n_genes, n_genes)
  mat_sp1 <- (mat_sp1 + t(mat_sp1)) / 2
  diag(mat_sp1) <- 1
  rownames(mat_sp1) <- colnames(mat_sp1) <- paste0("Gene_sp1_", 1:n_genes)

  mat_sp2 <- matrix(runif(n_genes * n_genes), n_genes, n_genes)
  mat_sp2 <- (mat_sp2 + t(mat_sp2)) / 2
  diag(mat_sp2) <- 1
  rownames(mat_sp2) <- colnames(mat_sp2) <- paste0("Gene_sp2_", 1:n_genes)

  # Convert to TriSimilarity (S4)
  tri_sp1 <- as.TriSimilarity(mat_sp1)
  tri_sp2 <- as.TriSimilarity(mat_sp2)

  # Create orthologs
  orthologs <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_genes),
    gene_sp2 = paste0("Gene_sp2_", 1:n_genes),
    type = "1:1"
  )

  # Calculate CCS with both formats
  ccs_mat <- calculate_ccs(mat_sp1, mat_sp2, orthologs, handle_self_diagonal = "none")
  ccs_tri <- calculate_ccs(tri_sp1, tri_sp2, orthologs, handle_self_diagonal = "none")

  # Results should be identical
  expect_equal(ccs_mat$CCS, ccs_tri$CCS, tolerance = 1e-10)
  expect_equal(ccs_mat$n_ref, ccs_tri$n_ref)
})

test_that("calculate_pcc_mr returns TriSimilarity by default", {
  skip_on_cran()

  set.seed(1111)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  # Default should return TriSimilarity (S4)
  sim <- calculate_pcc_mr(expr, method = "pcc")
  expect_s4_class(sim, "TriSimilarity")

  # With return_tri = FALSE should return matrix
  sim_mat <- calculate_pcc_mr(expr, method = "pcc", return_tri = FALSE)
  expect_true(is.matrix(sim_mat))

  # Values should be equivalent
  expect_equal(as.matrix(sim), sim_mat, tolerance = 1e-10)
})
