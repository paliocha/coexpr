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


# --- Tests for mean strategy ---

# Type naming convention from detect_ortholog_types:
#   "N:1" = n_sp1 > 1 (sp1 gene has N sp2 partners), n_sp2 == 1
#   "1:N" = n_sp1 == 1, n_sp2 > 1 (sp2 gene has N sp1 partners)

test_that("mean strategy averages CCS for N:1 orthologs (sp1 gene has multiple sp2 partners)", {
  # A2 maps to B2 and B3 → detect_ortholog_types says N:1 (A2 has n_sp1=2)
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B3")
  )

  ccs_values <- data.frame(
    gene_sp1 = c("A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B3"),
    CCS = c(0.8, 0.6, 0.4)
  )

  result <- handle_multicopy_orthologs(
    orthologs, strategy = "mean", ccs_values = ccs_values
  )

  # A1-B1 is 1:1 => kept as-is with CCS = 0.8
  a1_row <- result[result$gene_sp1 == "A1", ]
  expect_equal(nrow(a1_row), 1)
  expect_equal(a1_row$CCS, 0.8)

  # A2 is N:1 => grouped by gene_sp1, averaged to (0.6 + 0.4) / 2 = 0.5
  a2_row <- result[result$gene_sp1 == "A2", ]
  expect_equal(nrow(a2_row), 1)
  expect_equal(a2_row$CCS, 0.5)

  # gene_sp2 should be a composite for the aggregated row
  expect_true(grepl(";", a2_row$gene_sp2))
  expect_true(grepl("B2", a2_row$gene_sp2))
  expect_true(grepl("B3", a2_row$gene_sp2))
})

test_that("mean strategy averages CCS for 1:N orthologs (sp2 gene has multiple sp1 partners)", {
  # A1,A2 both map to B1 → detect_ortholog_types says 1:N (B1 has n_sp2=2)
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3"),
    gene_sp2 = c("B1", "B1", "B2")
  )

  ccs_values <- data.frame(
    gene_sp1 = c("A1", "A2", "A3"),
    gene_sp2 = c("B1", "B1", "B2"),
    CCS = c(0.7, 0.3, 0.9)
  )

  result <- handle_multicopy_orthologs(
    orthologs, strategy = "mean", ccs_values = ccs_values
  )

  # A3-B2 is 1:1 => kept as-is
  a3_row <- result[result$gene_sp2 == "B2", ]
  expect_equal(nrow(a3_row), 1)
  expect_equal(a3_row$CCS, 0.9)

  # A1,A2 -> B1 is 1:N => grouped by gene_sp2 (B1), averaged to (0.7 + 0.3) / 2 = 0.5
  b1_row <- result[result$gene_sp2 == "B1", ]
  expect_equal(nrow(b1_row), 1)
  expect_equal(b1_row$CCS, 0.5)

  # gene_sp1 should be composite
  expect_true(grepl(";", b1_row$gene_sp1))
  expect_true(grepl("A1", b1_row$gene_sp1))
  expect_true(grepl("A2", b1_row$gene_sp1))
})

test_that("mean strategy handles N:M orthologs", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B1", "B2")
  )

  ccs_values <- data.frame(
    gene_sp1 = c("A1", "A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B1", "B2"),
    CCS = c(0.8, 0.6, 0.4, 0.2)
  )

  result <- handle_multicopy_orthologs(
    orthologs, strategy = "mean", ccs_values = ccs_values
  )

  # All pairs are N:M, should be aggregated by gene_sp1
  # A1: mean(0.8, 0.6) = 0.7;  A2: mean(0.4, 0.2) = 0.3
  expect_equal(nrow(result), 2)
  expect_true(all(grepl("aggregated", result$type)))

  a1_row <- result[result$gene_sp1 == "A1", ]
  expect_equal(a1_row$CCS, 0.7)

  a2_row <- result[result$gene_sp1 == "A2", ]
  expect_equal(a2_row$CCS, 0.3)
})

test_that("mean strategy requires ccs_values", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2"),
    gene_sp2 = c("B1", "B2")
  )

  expect_error(
    handle_multicopy_orthologs(orthologs, strategy = "mean"),
    "requires pre-calculated ccs_values"
  )
})


# --- Tests for max strategy ---

test_that("max strategy selects best sp2 copy for N:1 (sp1 gene has multiple sp2 partners)", {
  # A2 maps to B2 and B3 → N:1
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B3")
  )

  ccs_values <- data.frame(
    gene_sp1 = c("A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B3"),
    CCS = c(0.8, 0.3, 0.9)
  )

  result <- handle_multicopy_orthologs(
    orthologs, strategy = "max", ccs_values = ccs_values
  )

  # A1-B1 (1:1) kept as-is
  expect_true("B1" %in% result$gene_sp2)

  # A2 has N:1 => should pick B3 (CCS = 0.9 > 0.3)
  a2_row <- result[result$gene_sp1 == "A2", ]
  expect_equal(nrow(a2_row), 1)
  expect_equal(a2_row$gene_sp2, "B3")
  expect_equal(a2_row$CCS, 0.9)
})

test_that("max strategy selects best sp1 copy for 1:N (sp2 gene has multiple sp1 partners)", {
  # A1,A2 both map to B1 → 1:N
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3"),
    gene_sp2 = c("B1", "B1", "B2")
  )

  ccs_values <- data.frame(
    gene_sp1 = c("A1", "A2", "A3"),
    gene_sp2 = c("B1", "B1", "B2"),
    CCS = c(0.7, 0.3, 0.9)
  )

  result <- handle_multicopy_orthologs(
    orthologs, strategy = "max", ccs_values = ccs_values
  )

  # A3-B2 (1:1) kept
  expect_true("A3" %in% result$gene_sp1)

  # 1:N to B1 => should pick A1 (CCS = 0.7 > 0.3)
  b1_row <- result[result$gene_sp2 == "B1", ]
  expect_equal(nrow(b1_row), 1)
  expect_equal(b1_row$gene_sp1, "A1")
  expect_equal(b1_row$CCS, 0.7)
})

test_that("max strategy handles N:M with two-pass selection", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B1", "B2")
  )

  ccs_values <- data.frame(
    gene_sp1 = c("A1", "A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B1", "B2"),
    CCS = c(0.9, 0.3, 0.4, 0.8)
  )

  result <- handle_multicopy_orthologs(
    orthologs, strategy = "max", ccs_values = ccs_values
  )

  # Two-pass: first pick best sp2 per sp1 (A1->B1=0.9, A2->B2=0.8)
  # Then pick best sp1 per sp2 (B1->A1=0.9, B2->A2=0.8)
  # Result: A1-B1 and A2-B2
  expect_equal(nrow(result), 2)
  expect_true("A1" %in% result$gene_sp1)
  expect_true("A2" %in% result$gene_sp1)
})

test_that("max strategy requires ccs_values", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2"),
    gene_sp2 = c("B1", "B2")
  )

  expect_error(
    handle_multicopy_orthologs(orthologs, strategy = "max"),
    "requires pre-calculated ccs_values"
  )
})


# --- Tests for best_hit strategy ---

test_that("best_hit strategy requires similarity matrices", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2"),
    gene_sp2 = c("B1", "B2")
  )

  expect_error(
    handle_multicopy_orthologs(orthologs, strategy = "best_hit"),
    "requires similarity_sp1 and similarity_sp2"
  )
})

test_that("best_hit strategy selects by average similarity", {
  # Create small similarity matrices
  sim_sp1 <- matrix(c(
    1.0, 0.8, 0.2,
    0.8, 1.0, 0.3,
    0.2, 0.3, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3")

  sim_sp2 <- matrix(c(
    1.0, 0.9, 0.1,
    0.9, 1.0, 0.4,
    0.1, 0.4, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3")

  # A2 maps to B2 and B3 → N:1 (A2 has n_sp1=2)
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B3")
  )

  result <- handle_multicopy_orthologs(
    orthologs,
    strategy = "best_hit",
    similarity_sp1 = sim_sp1,
    similarity_sp2 = sim_sp2
  )

  # A1-B1 (1:1) should be kept
  expect_true("A1" %in% result$gene_sp1)

  # A2->B2 or A2->B3: B2 has higher avg similarity ((1+0.9+0.4)/3=0.767 vs (0.1+0.4+1)/3=0.5)
  a2_row <- result[result$gene_sp1 == "A2", ]
  expect_equal(nrow(a2_row), 1)
  expect_equal(a2_row$gene_sp2, "B2")
})


# --- Edge cases ---

test_that("strategies handle all-1:1 input gracefully", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3"),
    gene_sp2 = c("B1", "B2", "B3")
  )

  ccs_values <- data.frame(
    gene_sp1 = c("A1", "A2", "A3"),
    gene_sp2 = c("B1", "B2", "B3"),
    CCS = c(0.5, 0.7, 0.3)
  )

  # All strategies should return all 3 rows for pure 1:1 input
  for (strat in c("strict", "all_pairs")) {
    result <- handle_multicopy_orthologs(orthologs, strategy = strat)
    expect_equal(nrow(result), 3, info = paste("Strategy:", strat))
  }

  for (strat in c("mean", "max")) {
    result <- handle_multicopy_orthologs(
      orthologs, strategy = strat, ccs_values = ccs_values
    )
    expect_equal(nrow(result), 3, info = paste("Strategy:", strat))
  }
})

test_that("type column is auto-detected when missing", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B3")
  )

  ccs_values <- data.frame(
    gene_sp1 = c("A1", "A2", "A2"),
    gene_sp2 = c("B1", "B2", "B3"),
    CCS = c(0.5, 0.7, 0.3)
  )

  # Should not error — type detection happens internally
  result_mean <- handle_multicopy_orthologs(
    orthologs, strategy = "mean", ccs_values = ccs_values
  )
  expect_true(nrow(result_mean) > 0)

  result_max <- handle_multicopy_orthologs(
    orthologs, strategy = "max", ccs_values = ccs_values
  )
  expect_true(nrow(result_max) > 0)
})
