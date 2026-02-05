test_that("detect_ortholog_types works correctly", {
  # Create test orthologs
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3", "A4", "A4"),
    gene_sp2 = c("B1", "B2", "B3", "B4", "B5", "B5")
  )

  # Detect types
  orthologs_typed <- coexpr:::detect_ortholog_types(orthologs)

  # Check types
  expect_equal(orthologs_typed$type[1], "1:1")  # A1-B1
  expect_equal(orthologs_typed$type[2], "1:1")  # A2-B2
  expect_equal(orthologs_typed$type[3], "N:1")  # A3-B3, A3-B4 (B3,B4 both map to A3)
  expect_equal(orthologs_typed$type[6], "N:M")  # A4-B5 (both A4 and B5 appear twice)
})

test_that("handle_multicopy_orthologs strict strategy works", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3"),
    gene_sp2 = c("B1", "B2", "B3", "B4"),
    type = c("1:1", "1:1", "1:N", "1:N")
  )

  result <- handle_multicopy_orthologs(orthologs, strategy = "strict")

  # Should only keep 1:1
  expect_equal(nrow(result), 2)
  expect_true(all(result$type == "1:1"))
})

test_that("handle_multicopy_orthologs all_pairs strategy works", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3"),
    gene_sp2 = c("B1", "B2", "B3", "B4")
  )

  result <- handle_multicopy_orthologs(orthologs, strategy = "all_pairs")

  # Should keep all pairs
  expect_equal(nrow(result), 4)
})
