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

# --- Tests for collapse_orthologs ---

# Helper: build a symmetric similarity matrix from named values
make_sim_matrix <- function(gene_names, off_diag_vals = NULL) {
  n <- length(gene_names)
  mat <- diag(n)
  rownames(mat) <- colnames(mat) <- gene_names
  if (!is.null(off_diag_vals)) {
    # off_diag_vals is a list of (i, j, val) triples
    for (entry in off_diag_vals) {
      mat[entry[1], entry[2]] <- entry[3]
      mat[entry[2], entry[1]] <- entry[3]
    }
  }
  mat
}

test_that("collapse_orthologs basic N:1 collapse selects best homeolog", {
  # 3 natural 1:1 pairs: A1-B1, A2-B2, A3-B3
  # 1 N:1 group: A4 maps to B4a and B4b (sp1 gene A4 has 2 sp2 partners)
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A4", "A4"),
    gene_sp2 = c("B1", "B2", "B3", "B4a", "B4b")
  )

  # Species 1: varied off-diag values so reference vectors have non-zero variance
  sim_sp1 <- matrix(c(
    1.0, 0.8, 0.3, 0.7,
    0.8, 1.0, 0.4, 0.6,
    0.3, 0.4, 1.0, 0.5,
    0.7, 0.6, 0.5, 1.0
  ), nrow = 4, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3", "A4")

  # Species 2: B4a mimics the ref pattern of A4, B4b does not
  sim_sp2 <- matrix(c(
    1.0, 0.8, 0.3, 0.7, 0.1,
    0.8, 1.0, 0.4, 0.6, 0.2,
    0.3, 0.4, 1.0, 0.5, 0.3,
    0.7, 0.6, 0.5, 1.0, 0.2,
    0.1, 0.2, 0.3, 0.2, 1.0
  ), nrow = 5, byrow = TRUE)
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3", "B4a", "B4b")

  result <- collapse_orthologs(
    orthologs, sim_sp1, sim_sp2,
    multicopy_sp = "sp2", max_copy_number = 2L
  )

  # Should have 4 rows: 3 natural 1:1 + 1 collapsed
  expect_equal(nrow(result), 4)

  # B4a should be selected (better match to ref pattern)
  collapsed_row <- result[!is.na(result$original_type), ]
  expect_equal(nrow(collapsed_row), 1)
  expect_equal(collapsed_row$gene_sp1, "A4")
  expect_equal(collapsed_row$gene_sp2, "B4a")
})

test_that("collapse_orthologs output has correct columns and types", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3"),
    gene_sp2 = c("B1", "B2", "B3a", "B3b")
  )

  sim_sp1 <- matrix(c(
    1.0, 0.7, 0.3,
    0.7, 1.0, 0.5,
    0.3, 0.5, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3")

  sim_sp2 <- matrix(c(
    1.0, 0.7, 0.4, 0.2,
    0.7, 1.0, 0.5, 0.3,
    0.4, 0.5, 1.0, 0.6,
    0.2, 0.3, 0.6, 1.0
  ), nrow = 4, byrow = TRUE)
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3a", "B3b")

  result <- collapse_orthologs(orthologs, sim_sp1, sim_sp2, multicopy_sp = "sp2")

  expect_true(all(c("gene_sp1", "gene_sp2", "type", "original_type",
                     "homeolog_score", "n_candidates") %in% colnames(result)))
  expect_type(result$gene_sp1, "character")
  expect_type(result$gene_sp2, "character")
  expect_type(result$type, "character")
  expect_type(result$original_type, "character")
  expect_type(result$homeolog_score, "double")
  expect_type(result$n_candidates, "integer")
})

test_that("collapsed pairs get type='1:1' and original_type records source", {
  # Need at least 2 natural 1:1 for non-degenerate reference vectors
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3"),
    gene_sp2 = c("B1", "B2", "B3a", "B3b")
  )

  sim_sp1 <- matrix(c(
    1.0, 0.7, 0.4,
    0.7, 1.0, 0.5,
    0.4, 0.5, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3")

  sim_sp2 <- matrix(c(
    1.0, 0.7, 0.6, 0.2,
    0.7, 1.0, 0.5, 0.3,
    0.6, 0.5, 1.0, 0.4,
    0.2, 0.3, 0.4, 1.0
  ), nrow = 4, byrow = TRUE)
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3a", "B3b")

  result <- collapse_orthologs(orthologs, sim_sp1, sim_sp2, multicopy_sp = "sp2")

  # All rows should have type = "1:1"
  expect_true(all(result$type == "1:1"))

  # Natural 1:1 should have NA original_type
  natural <- result[result$gene_sp1 == "A1", ]
  expect_true(is.na(natural$original_type))
  expect_true(is.na(natural$homeolog_score))
  expect_equal(natural$n_candidates, 1L)

  # Collapsed pair should have original_type = "N:1"
  collapsed <- result[result$gene_sp1 == "A3", ]
  expect_equal(nrow(collapsed), 1)
  expect_equal(collapsed$original_type, "N:1")
  expect_false(is.na(collapsed$homeolog_score))
  expect_equal(collapsed$n_candidates, 2L)
})

test_that("max_copy_number filters out larger groups", {
  # A3 maps to B3a, B3b (2 copies) — should be collapsed with max_copy_number=2
  # A4 maps to B4a, B4b, B4c (3 copies) — should be skipped
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3", "A4", "A4", "A4"),
    gene_sp2 = c("B1", "B2", "B3a", "B3b", "B4a", "B4b", "B4c")
  )

  sim_sp1 <- matrix(c(
    1.0, 0.7, 0.4, 0.6,
    0.7, 1.0, 0.5, 0.3,
    0.4, 0.5, 1.0, 0.8,
    0.6, 0.3, 0.8, 1.0
  ), nrow = 4, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3", "A4")

  sim_sp2 <- matrix(0.4, 7, 7)
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3a", "B3b",
                                                "B4a", "B4b", "B4c")
  # Add varied values for ref genes
  sim_sp2["B1", "B2"] <- sim_sp2["B2", "B1"] <- 0.7
  sim_sp2["B3a", "B1"] <- sim_sp2["B1", "B3a"] <- 0.6
  sim_sp2["B3a", "B2"] <- sim_sp2["B2", "B3a"] <- 0.5
  sim_sp2["B3b", "B1"] <- sim_sp2["B1", "B3b"] <- 0.2
  sim_sp2["B3b", "B2"] <- sim_sp2["B2", "B3b"] <- 0.3

  result <- collapse_orthologs(
    orthologs, sim_sp1, sim_sp2,
    multicopy_sp = "sp2", max_copy_number = 2L
  )

  # 2 natural 1:1 + 1 collapsed (A3 group) = 3; A4 group excluded
  expect_equal(nrow(result), 3)
  expect_false("A4" %in% result$gene_sp1)
  expect_true("A3" %in% result$gene_sp1)
})

test_that("max_copy_number = NULL collapses all group sizes", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3", "A3"),
    gene_sp2 = c("B1", "B2", "B3a", "B3b", "B3c")
  )

  sim_sp1 <- matrix(c(
    1.0, 0.7, 0.4,
    0.7, 1.0, 0.5,
    0.4, 0.5, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3")

  sim_sp2 <- matrix(0.4, 5, 5)
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3a", "B3b", "B3c")
  sim_sp2["B1", "B2"] <- sim_sp2["B2", "B1"] <- 0.7
  sim_sp2["B3a", "B1"] <- sim_sp2["B1", "B3a"] <- 0.6
  sim_sp2["B3a", "B2"] <- sim_sp2["B2", "B3a"] <- 0.5
  sim_sp2["B3b", "B1"] <- sim_sp2["B1", "B3b"] <- 0.3
  sim_sp2["B3b", "B2"] <- sim_sp2["B2", "B3b"] <- 0.2
  sim_sp2["B3c", "B1"] <- sim_sp2["B1", "B3c"] <- 0.1
  sim_sp2["B3c", "B2"] <- sim_sp2["B2", "B3c"] <- 0.3

  result <- collapse_orthologs(
    orthologs, sim_sp1, sim_sp2,
    multicopy_sp = "sp2", max_copy_number = NULL
  )

  # 2 natural 1:1 + 1 collapsed (3-copy group included) = 3
  expect_equal(nrow(result), 3)
  expect_true("A3" %in% result$gene_sp1)
})

test_that("collapse_orthologs works with TriSimilarity input", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3"),
    gene_sp2 = c("B1", "B2", "B3a", "B3b")
  )

  sim_sp1_mat <- matrix(c(
    1.0, 0.7, 0.4,
    0.7, 1.0, 0.5,
    0.4, 0.5, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp1_mat) <- colnames(sim_sp1_mat) <- c("A1", "A2", "A3")

  sim_sp2_mat <- matrix(c(
    1.0, 0.7, 0.9, 0.1,
    0.7, 1.0, 0.6, 0.2,
    0.9, 0.6, 1.0, 0.3,
    0.1, 0.2, 0.3, 1.0
  ), nrow = 4, byrow = TRUE)
  rownames(sim_sp2_mat) <- colnames(sim_sp2_mat) <- c("B1", "B2", "B3a", "B3b")

  # Convert to TriSimilarity
  tri_sp1 <- as.TriSimilarity(sim_sp1_mat)
  tri_sp2 <- as.TriSimilarity(sim_sp2_mat)

  result_tri <- collapse_orthologs(
    orthologs, tri_sp1, tri_sp2, multicopy_sp = "sp2"
  )

  result_mat <- collapse_orthologs(
    orthologs, sim_sp1_mat, sim_sp2_mat, multicopy_sp = "sp2"
  )

  # Both should give same result
  expect_equal(nrow(result_tri), nrow(result_mat))
  expect_equal(result_tri$gene_sp2, result_mat$gene_sp2)
  expect_equal(result_tri$homeolog_score, result_mat$homeolog_score,
               tolerance = 1e-10)
})

test_that("collapse_orthologs works with full matrix input", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3"),
    gene_sp2 = c("B1", "B2", "B3a", "B3b")
  )

  sim_sp1 <- matrix(c(
    1.0, 0.7, 0.4,
    0.7, 1.0, 0.5,
    0.4, 0.5, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3")

  sim_sp2 <- matrix(c(
    1.0, 0.7, 0.5, 0.3,
    0.7, 1.0, 0.6, 0.4,
    0.5, 0.6, 1.0, 0.5,
    0.3, 0.4, 0.5, 1.0
  ), nrow = 4, byrow = TRUE)
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3a", "B3b")

  result <- collapse_orthologs(orthologs, sim_sp1, sim_sp2, multicopy_sp = "sp2")

  expect_equal(nrow(result), 3)
  expect_true(all(result$type == "1:1"))
})

test_that("multicopy_sp = 'sp1' collapses 1:N groups correctly", {
  # 1:N: gene_sp2 appears multiple times (multiple sp1 genes map to same sp2 gene)
  # A1a, A1b both map to B1; A2-B2 and A3-B3 are natural 1:1
  orthologs <- data.frame(
    gene_sp1 = c("A1a", "A1b", "A2", "A3"),
    gene_sp2 = c("B1", "B1", "B2", "B3")
  )

  sim_sp1 <- matrix(c(
    1.0, 0.5, 0.9, 0.3,
    0.5, 1.0, 0.1, 0.4,
    0.9, 0.1, 1.0, 0.6,
    0.3, 0.4, 0.6, 1.0
  ), nrow = 4, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1a", "A1b", "A2", "A3")

  sim_sp2 <- matrix(c(
    1.0, 0.8, 0.3,
    0.8, 1.0, 0.5,
    0.3, 0.5, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3")

  result <- collapse_orthologs(
    orthologs, sim_sp1, sim_sp2, multicopy_sp = "sp1"
  )

  # 2 natural 1:1 (A2-B2, A3-B3) + 1 collapsed = 3
  expect_equal(nrow(result), 3)

  # A1a should be selected (stronger co-expression correlation with refs)
  collapsed <- result[!is.na(result$original_type), ]
  expect_equal(nrow(collapsed), 1)
  expect_equal(collapsed$gene_sp1, "A1a")
  expect_equal(collapsed$gene_sp2, "B1")
  expect_equal(collapsed$original_type, "1:N")
})

test_that("multicopy_sp = 'both' collapses both sides", {
  # N:1 group: A3 maps to B3a, B3b
  # 1:N group: A4a, A4b both map to B4
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3", "A4a", "A4b"),
    gene_sp2 = c("B1", "B2", "B3a", "B3b", "B4", "B4")
  )

  sim_sp1 <- matrix(c(
    1.0, 0.7, 0.4, 0.6, 0.3,
    0.7, 1.0, 0.5, 0.8, 0.2,
    0.4, 0.5, 1.0, 0.6, 0.4,
    0.6, 0.8, 0.6, 1.0, 0.5,
    0.3, 0.2, 0.4, 0.5, 1.0
  ), nrow = 5, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3", "A4a", "A4b")

  sim_sp2 <- matrix(c(
    1.0, 0.7, 0.6, 0.2, 0.8,
    0.7, 1.0, 0.5, 0.3, 0.6,
    0.6, 0.5, 1.0, 0.4, 0.5,
    0.2, 0.3, 0.4, 1.0, 0.3,
    0.8, 0.6, 0.5, 0.3, 1.0
  ), nrow = 5, byrow = TRUE)
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3a", "B3b", "B4")

  result <- collapse_orthologs(
    orthologs, sim_sp1, sim_sp2, multicopy_sp = "both"
  )

  # 2 natural 1:1 + 1 N:1 collapsed + 1 1:N collapsed = 4
  expect_equal(nrow(result), 4)
  expect_true(all(result$type == "1:1"))

  # Check both types were collapsed
  collapsed_types <- result$original_type[!is.na(result$original_type)]
  expect_true("N:1" %in% collapsed_types)
  expect_true("1:N" %in% collapsed_types)
})

test_that("gene not in similarity matrix warns and skips", {
  # Need 2+ natural 1:1. A1-B1, A2-B2 are 1:1; A3 maps to B3a, B3b (N:1)
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A3"),
    gene_sp2 = c("B1", "B2", "B3a", "B3b")
  )

  sim_sp1 <- matrix(c(
    1.0, 0.7, 0.4,
    0.7, 1.0, 0.5,
    0.4, 0.5, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3")

  # B3b not in similarity matrix
  sim_sp2 <- matrix(c(
    1.0, 0.7, 0.5,
    0.7, 1.0, 0.6,
    0.5, 0.6, 1.0
  ), nrow = 3, byrow = TRUE)
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3a")

  expect_warning(
    result <- collapse_orthologs(
      orthologs, sim_sp1, sim_sp2, multicopy_sp = "sp2"
    ),
    "not found in similarity_sp2"
  )

  # Should still get a result with B3a selected
  collapsed <- result[result$gene_sp1 == "A3", ]
  expect_equal(nrow(collapsed), 1)
  expect_equal(collapsed$gene_sp2, "B3a")
})

test_that("no multi-copy groups returns 1:1 with message", {
  orthologs <- data.frame(
    gene_sp1 = c("A1", "A2", "A3"),
    gene_sp2 = c("B1", "B2", "B3")
  )

  sim_sp1 <- matrix(0.5, 3, 3)
  diag(sim_sp1) <- 1
  rownames(sim_sp1) <- colnames(sim_sp1) <- c("A1", "A2", "A3")

  sim_sp2 <- matrix(0.5, 3, 3)
  diag(sim_sp2) <- 1
  rownames(sim_sp2) <- colnames(sim_sp2) <- c("B1", "B2", "B3")

  expect_message(
    result <- collapse_orthologs(orthologs, sim_sp1, sim_sp2),
    "No multi-copy groups found"
  )

  expect_equal(nrow(result), 3)
  expect_true(all(result$type == "1:1"))
  expect_true(all(is.na(result$original_type)))
})

test_that("collapse_orthologs result integrates with calculate_ccs", {
  # Build realistic similarity matrices with enough genes
  set.seed(42)
  n_sp1 <- 15
  n_sp2 <- 17

  # Create random correlation matrices
  raw1 <- matrix(stats::runif(n_sp1^2, 0.1, 0.9), n_sp1, n_sp1)
  sim_sp1 <- (raw1 + t(raw1)) / 2
  diag(sim_sp1) <- 1
  gene_names_sp1 <- paste0("A", seq_len(n_sp1))
  rownames(sim_sp1) <- colnames(sim_sp1) <- gene_names_sp1

  raw2 <- matrix(stats::runif(n_sp2^2, 0.1, 0.9), n_sp2, n_sp2)
  sim_sp2 <- (raw2 + t(raw2)) / 2
  diag(sim_sp2) <- 1
  gene_names_sp2 <- paste0("B", seq_len(n_sp2))
  rownames(sim_sp2) <- colnames(sim_sp2) <- gene_names_sp2

  # Create orthologs: 10 natural 1:1 + 2 N:1 groups (2 copies each)
  orthologs <- data.frame(
    gene_sp1 = c(paste0("A", 1:10), "A11", "A11", "A12", "A12"),
    gene_sp2 = c(paste0("B", 1:10), "B11", "B12", "B13", "B14")
  )

  # Collapse
  expanded <- collapse_orthologs(
    orthologs, sim_sp1, sim_sp2,
    multicopy_sp = "sp2", max_copy_number = 2L
  )

  # Should have 12 rows (10 + 2 collapsed groups)
  expect_equal(nrow(expanded), 12)

  # Use as reference in calculate_ccs — should not error
  ccs_result <- calculate_ccs(sim_sp1, sim_sp2, expanded, use_only_1to1 = TRUE)

  expect_true(is.data.frame(ccs_result))
  expect_true("CCS" %in% colnames(ccs_result))
  expect_true(all(!is.na(ccs_result$CCS)))
  # CCS should expand: 12 pairs evaluated
  expect_equal(nrow(ccs_result), 12)
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
