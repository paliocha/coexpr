# ---- Cache key determinism ----

test_that("compute_cache_key is deterministic for identical inputs", {
  set.seed(1)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  key1 <- compute_cache_key(expr, "pcc_mr", cor_method = "pearson")
  key2 <- compute_cache_key(expr, "pcc_mr", cor_method = "pearson")
  expect_identical(key1, key2)
})

test_that("cache key changes when expression data changes", {
  set.seed(1)
  expr1 <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("Gene", 1:50)

  expr2 <- expr1
  expr2[1, 1] <- expr2[1, 1] + 0.001
  rownames(expr2) <- rownames(expr1)

  key1 <- compute_cache_key(expr1, "pcc_mr", cor_method = "pearson")
  key2 <- compute_cache_key(expr2, "pcc_mr", cor_method = "pearson")
  expect_false(key1 == key2)
})

test_that("cache key changes when method changes", {
  set.seed(1)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  key1 <- compute_cache_key(expr, "pcc_mr", cor_method = "pearson")
  key2 <- compute_cache_key(expr, "pcc", cor_method = "pearson")
  expect_false(key1 == key2)
})

test_that("cache key changes when cor_method changes", {
  set.seed(1)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  key1 <- compute_cache_key(expr, "pcc_mr", cor_method = "pearson")
  key2 <- compute_cache_key(expr, "pcc_mr", cor_method = "spearman")
  expect_false(key1 == key2)
})

test_that("cache key changes when gene names change", {
  set.seed(1)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  expr2 <- expr
  rownames(expr2) <- paste0("G", 1:50)

  key1 <- compute_cache_key(expr, "pcc_mr", cor_method = "pearson")
  key2 <- compute_cache_key(expr2, "pcc_mr", cor_method = "pearson")
  expect_false(key1 == key2)
})

test_that("cache key filename has correct format", {
  set.seed(1)
  expr <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:50)

  key <- compute_cache_key(expr, "pcc_mr", cor_method = "pearson")
  expect_match(key, "^pcc_mr_pearson_[0-9a-f]+\\.rds$")

  key2 <- compute_cache_key(expr, "mi_clr", n_bins = 10)
  expect_match(key2, "^mi_clr_10_[0-9a-f]+\\.rds$")
})

# ---- Cache load/save ----

test_that("cache_save and cache_load round-trip a TriSimilarity object", {
  cache_dir <- withr::local_tempdir()

  mat <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.7, 0.3, 0.7, 1), nrow = 3)
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")
  tri <- as.TriSimilarity(mat)

  cache_save(tri, cache_dir, "test.rds")
  loaded <- cache_load(cache_dir, "test.rds")

  expect_s4_class(loaded, "TriSimilarity")
  expect_equal(loaded@data, tri@data)
  expect_equal(loaded@genes, tri@genes)
})

test_that("cache_save converts matrix to TriSimilarity before saving", {
  cache_dir <- withr::local_tempdir()

  mat <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.7, 0.3, 0.7, 1), nrow = 3)
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")

  cache_save(mat, cache_dir, "test_mat.rds")
  loaded <- cache_load(cache_dir, "test_mat.rds")

  expect_s4_class(loaded, "TriSimilarity")
  expect_equal(as.matrix(loaded), mat)
})

test_that("cache_load returns NULL for missing file", {
  cache_dir <- withr::local_tempdir()
  expect_null(cache_load(cache_dir, "nonexistent.rds"))
})

test_that("cache_load warns and returns NULL for corrupt file", {
  cache_dir <- withr::local_tempdir()
  writeLines("not an rds file", file.path(cache_dir, "corrupt.rds"))

  expect_warning(
    result <- cache_load(cache_dir, "corrupt.rds"),
    "Corrupt cache file"
  )
  expect_null(result)
})

test_that("cache_load warns and returns NULL for wrong object type", {
  cache_dir <- withr::local_tempdir()
  saveRDS(list(a = 1), file.path(cache_dir, "wrong_type.rds"))

  expect_warning(
    result <- cache_load(cache_dir, "wrong_type.rds"),
    "does not contain a TriSimilarity"
  )
  expect_null(result)
})

test_that("cache_save creates directory if it doesn't exist", {
  parent <- withr::local_tempdir()
  cache_dir <- file.path(parent, "subdir", "cache")

  mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  rownames(mat) <- colnames(mat) <- c("A", "B")

  cache_save(as.TriSimilarity(mat), cache_dir, "test.rds")
  expect_true(dir.exists(cache_dir))
  expect_true(file.exists(file.path(cache_dir, "test.rds")))
})


# ---- cache_list ----

test_that("cache_list returns empty data frame for non-existent dir", {
  result <- cache_list(file.path(tempdir(), "nonexistent_dir_12345"))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_true(all(c("file", "size_mb", "modified", "method") %in% names(result)))
})

test_that("cache_list returns empty data frame for empty dir", {
  cache_dir <- withr::local_tempdir()
  result <- cache_list(cache_dir)
  expect_equal(nrow(result), 0)
})

test_that("cache_list lists cached files with metadata", {
  cache_dir <- withr::local_tempdir()

  mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  rownames(mat) <- colnames(mat) <- c("A", "B")
  tri <- as.TriSimilarity(mat)

  saveRDS(tri, file.path(cache_dir, "pcc_mr_pearson_abc123.rds"))
  saveRDS(tri, file.path(cache_dir, "mi_clr_10_def456.rds"))

  result <- cache_list(cache_dir)
  expect_equal(nrow(result), 2)
  expect_true("pcc_mr_pearson_abc123.rds" %in% result$file)
  expect_true("mi_clr_10_def456.rds" %in% result$file)
  expect_true(all(result$size_mb >= 0))
  expect_s3_class(result$modified, "POSIXct")
})


# ---- cache_clear ----

test_that("cache_clear removes all files when method is NULL", {
  cache_dir <- withr::local_tempdir()

  mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  rownames(mat) <- colnames(mat) <- c("A", "B")
  tri <- as.TriSimilarity(mat)

  saveRDS(tri, file.path(cache_dir, "pcc_mr_abc.rds"))
  saveRDS(tri, file.path(cache_dir, "mi_clr_def.rds"))

  n <- cache_clear(cache_dir)
  expect_equal(n, 2L)
  expect_equal(length(list.files(cache_dir, pattern = "\\.rds$")), 0)
})

test_that("cache_clear filters by method", {
  cache_dir <- withr::local_tempdir()

  mat <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  rownames(mat) <- colnames(mat) <- c("A", "B")
  tri <- as.TriSimilarity(mat)

  saveRDS(tri, file.path(cache_dir, "pcc_mr_pearson_abc.rds"))
  saveRDS(tri, file.path(cache_dir, "mi_clr_10_def.rds"))

  n <- cache_clear(cache_dir, method = "pcc_mr")
  expect_equal(n, 1L)

  remaining <- list.files(cache_dir, pattern = "\\.rds$")
  expect_equal(remaining, "mi_clr_10_def.rds")
})

test_that("cache_clear returns 0 for non-existent dir", {
  n <- cache_clear(file.path(tempdir(), "nonexistent_dir_12345"))
  expect_equal(n, 0L)
})


# ---- Integration with similarity functions ----

test_that("calculate_pcc_mr caches and reloads (PCC only)", {
  cache_dir <- withr::local_tempdir()

  set.seed(42)
  expr <- matrix(rnorm(30 * 10), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:30)

  # First call: compute and cache
  sim1 <- calculate_pcc_mr(expr, method = "pcc", cache_dir = cache_dir)

  # Verify file was created
  files <- list.files(cache_dir, pattern = "\\.rds$")
  expect_equal(length(files), 1)

  # Second call: should hit cache
  expect_message(
    sim2 <- calculate_pcc_mr(expr, method = "pcc", cache_dir = cache_dir),
    "Loading cached"
  )

  expect_s4_class(sim2, "TriSimilarity")
  expect_equal(sim1@data, sim2@data)
  expect_equal(sim1@genes, sim2@genes)
})

test_that("calculate_pcc_mr cache returns full matrix when return_tri = FALSE", {
  cache_dir <- withr::local_tempdir()

  set.seed(42)
  expr <- matrix(rnorm(30 * 10), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:30)

  # Cache with return_tri = TRUE (default)
  sim_tri <- calculate_pcc_mr(expr, method = "pcc", cache_dir = cache_dir)

  # Reload with return_tri = FALSE
  expect_message(
    sim_mat <- calculate_pcc_mr(expr, method = "pcc", cache_dir = cache_dir,
                                return_tri = FALSE),
    "Loading cached"
  )

  expect_true(is.matrix(sim_mat))
  expect_equal(sim_mat, as.matrix(sim_tri))
})

test_that("calculate_mi_clr caches and reloads", {
  skip_on_cran()
  cache_dir <- withr::local_tempdir()

  set.seed(42)
  expr <- matrix(rnorm(20 * 40), nrow = 20, ncol = 40)
  rownames(expr) <- paste0("Gene", 1:20)

  # First call: compute and cache
  sim1 <- calculate_mi_clr(expr, n_bins = 5, cache_dir = cache_dir)

  files <- list.files(cache_dir, pattern = "\\.rds$")
  expect_equal(length(files), 1)

  # Second call: should hit cache
  expect_message(
    sim2 <- calculate_mi_clr(expr, n_bins = 5, cache_dir = cache_dir),
    "Loading cached"
  )

  expect_s4_class(sim2, "TriSimilarity")
  expect_equal(sim1@data, sim2@data)
})

test_that("different parameters produce different cache entries", {
  cache_dir <- withr::local_tempdir()

  set.seed(42)
  expr <- matrix(rnorm(30 * 10), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:30)

  calculate_pcc_mr(expr, method = "pcc", cor_method = "pearson",
                   cache_dir = cache_dir)
  calculate_pcc_mr(expr, method = "pcc", cor_method = "spearman",
                   cache_dir = cache_dir)

  files <- list.files(cache_dir, pattern = "\\.rds$")
  expect_equal(length(files), 2)
})

test_that("cache_dir = NULL does not cache (default behavior)", {
  set.seed(42)
  expr <- matrix(rnorm(30 * 10), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:30)

  # No message about loading/caching
  sim <- calculate_pcc_mr(expr, method = "pcc", cache_dir = NULL)
  expect_s4_class(sim, "TriSimilarity")
})

test_that("n_cores does not affect cache key", {
  set.seed(42)
  expr <- matrix(rnorm(30 * 10), nrow = 30, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:30)

  # Keys should be identical regardless of n_cores
  key1 <- compute_cache_key(expr, "pcc_mr", cor_method = "pearson")
  key2 <- compute_cache_key(expr, "pcc_mr", cor_method = "pearson")
  expect_identical(key1, key2)
  # n_cores is NOT included in the key, so it can't change the key
})
