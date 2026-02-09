test_that("calculate_ors computes ORS correctly", {
  # Create test CCS results
  set.seed(234)
  n_orthologs <- 50

  ccs_results <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_orthologs),
    gene_sp2 = paste0("Gene_sp2_", 1:n_orthologs),
    CCS = rnorm(n_orthologs, mean = 0.5, sd = 0.2),
    n_ref = 50,
    type = "1:1"
  )

  # Calculate ORS
  ors_results <- calculate_ors(ccs_results, return_log = FALSE)

  # Check output structure
  expect_true(is.data.frame(ors_results))
  expect_true("ORS" %in% colnames(ors_results))
  expect_false("ORS_sp1_to_sp2" %in% colnames(ors_results))
  expect_false("ORS_sp2_to_sp1" %in% colnames(ors_results))
  expect_equal(nrow(ors_results), n_orthologs)

  # Check ORS values are in [0, 1]
  expect_true(all(ors_results$ORS >= 0 & ors_results$ORS <= 1))

  # Check ORS uses global ranking (proper distribution from 1/n to 1)
  # ORS should span [1/n, 1] for n unique CCS values
  expected_global <- rank(ccs_results$CCS, ties.method = "average") / n_orthologs
  expect_equal(ors_results$ORS, expected_global)

  # Check that ORS gives proper distribution (not all 1s)
  expect_true(min(ors_results$ORS) < 0.5)
  expect_true(max(ors_results$ORS) > 0.5)
})

test_that("calculate_ors logORS transformation works", {
  set.seed(345)
  n_orthologs <- 100

  # Create CCS results with varying scores
  ccs_results <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_orthologs),
    gene_sp2 = paste0("Gene_sp2_", 1:n_orthologs),
    CCS = runif(n_orthologs, -0.5, 1),
    n_ref = 50
  )

  # Test with return_log = TRUE
  ors_with_log <- calculate_ors(ccs_results, return_log = TRUE)
  expect_true("logORS" %in% colnames(ors_with_log))

  # Test with return_log = FALSE
  ors_without_log <- calculate_ors(ccs_results, return_log = FALSE)
  expect_false("logORS" %in% colnames(ors_without_log))

  # Check logORS transformation is correct
  # logORS = -log10(1 + 1e-4 - ORS)
  expected_logORS <- -log10(1 + 1e-4 - ors_with_log$ORS)
  expect_equal(ors_with_log$logORS, expected_logORS)

  # Check that high ORS gives positive logORS
  # Top 10% should have logORS > 1
  top_10pct_idx <- which(ors_with_log$ORS > 0.9)
  if (length(top_10pct_idx) > 0) {
    expect_true(all(ors_with_log$logORS[top_10pct_idx] > 1))
  }
})

test_that("calculate_ors directional scoring works with directional = TRUE", {
  # Create asymmetric scenario with enough orthologs (need >= 10)
  ccs_results <- data.frame(
    gene_sp1 = c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A10", "A10"),
    gene_sp2 = c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12"),
    CCS = c(0.9, 0.8, 0.75, 0.7, 0.65, 0.62, 0.58, 0.55, 0.52, 0.6, 0.5, 0.4),
    n_ref = 10
  )

  ors_results <- calculate_ors(ccs_results, return_log = FALSE, directional = TRUE)

  # Directional columns should be present
  expect_true("ORS_sp1_to_sp2" %in% colnames(ors_results))
  expect_true("ORS_sp2_to_sp1" %in% colnames(ors_results))

  # For A10 (3 copies in sp2), ORS_sp1_to_sp2 should rank among its 3 copies
  a10_orthologs <- ors_results[ors_results$gene_sp1 == "A10", ]
  expect_equal(nrow(a10_orthologs), 3)

  # Highest CCS (A10-B10: 0.6) should have highest ORS_sp1_to_sp2
  expect_true(a10_orthologs$ORS_sp1_to_sp2[a10_orthologs$gene_sp2 == "B10"] == 1.0)
  expect_true(a10_orthologs$ORS_sp1_to_sp2[a10_orthologs$gene_sp2 == "B11"] == 2/3)
  expect_true(a10_orthologs$ORS_sp1_to_sp2[a10_orthologs$gene_sp2 == "B12"] == 1/3)
})

test_that("test_ors_significance with BH correction (default)", {
  set.seed(456)
  n_orthologs <- 100

  ccs_results <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_orthologs),
    gene_sp2 = paste0("Gene_sp2_", 1:n_orthologs),
    CCS = rnorm(n_orthologs, mean = 0.5, sd = 0.3),
    n_ref = 50
  )

  ors_results <- calculate_ors(ccs_results)

  # Default is BH correction
  ors_sig <- test_ors_significance(ors_results, alpha = 0.05)

  # Check output columns
  expect_true("pvalue" %in% colnames(ors_sig))
  expect_true("padj" %in% colnames(ors_sig))
  expect_true("significant" %in% colnames(ors_sig))

  # Raw pvalue = 1 - ORS
  expected_pvalue <- 1 - ors_sig$ORS
  expect_equal(ors_sig$pvalue, expected_pvalue)

  # padj should be BH-adjusted
  expected_padj <- p.adjust(ors_sig$pvalue, method = "BH")
  expect_equal(ors_sig$padj, expected_padj)

  # Significance is based on padj, not raw pvalue
  expect_equal(ors_sig$significant, ors_sig$padj < 0.05)

  # BH-adjusted p-values are >= raw p-values (more conservative)
  expect_true(all(ors_sig$padj >= ors_sig$pvalue))
})

test_that("test_ors_significance with p_adjust_method = 'none' uses raw p-values", {
  set.seed(456)
  n_orthologs <- 100

  ccs_results <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_orthologs),
    gene_sp2 = paste0("Gene_sp2_", 1:n_orthologs),
    CCS = rnorm(n_orthologs, mean = 0.5, sd = 0.3),
    n_ref = 50
  )

  ors_results <- calculate_ors(ccs_results)
  ors_sig <- test_ors_significance(ors_results, alpha = 0.05,
                                   p_adjust_method = "none")

  # No padj column when method = "none"
  expect_false("padj" %in% colnames(ors_sig))
  expect_true("pvalue" %in% colnames(ors_sig))

  # Significance based on raw pvalue
  expect_equal(ors_sig$significant, ors_sig$pvalue < 0.05)
})

test_that("BH correction is more conservative than no correction", {
  set.seed(789)
  n_orthologs <- 100

  ccs_results <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_orthologs),
    gene_sp2 = paste0("Gene_sp2_", 1:n_orthologs),
    CCS = rnorm(n_orthologs, mean = 0.5, sd = 0.3),
    n_ref = 50
  )

  ors_results <- calculate_ors(ccs_results)

  sig_none <- test_ors_significance(ors_results, alpha = 0.05,
                                    p_adjust_method = "none")
  sig_bh <- test_ors_significance(ors_results, alpha = 0.05,
                                  p_adjust_method = "BH")

  n_none <- sum(sig_none$significant)
  n_bh <- sum(sig_bh$significant)

  # BH should find fewer or equal significant hits
  expect_true(n_bh <= n_none)
})

test_that("test_ors_significance supports different correction methods", {
  set.seed(456)
  n_orthologs <- 50

  ccs_results <- data.frame(
    gene_sp1 = paste0("Gene_sp1_", 1:n_orthologs),
    gene_sp2 = paste0("Gene_sp2_", 1:n_orthologs),
    CCS = rnorm(n_orthologs, mean = 0.5, sd = 0.3),
    n_ref = 50
  )

  ors_results <- calculate_ors(ccs_results)

  # Bonferroni
  sig_bonf <- test_ors_significance(ors_results, p_adjust_method = "bonferroni")
  expect_true("padj" %in% colnames(sig_bonf))
  expect_equal(sig_bonf$padj,
               p.adjust(sig_bonf$pvalue, method = "bonferroni"))

  # holm
  sig_holm <- test_ors_significance(ors_results, p_adjust_method = "holm")
  expect_true("padj" %in% colnames(sig_holm))
  expect_equal(sig_holm$padj,
               p.adjust(sig_holm$pvalue, method = "holm"))
})

test_that("calculate_ors validates inputs", {
  # Test missing 'ccs' column
  bad_results <- data.frame(
    gene_sp1 = c("A", "B"),
    gene_sp2 = c("C", "D"),
    wrong_col = c(0.5, 0.6)
  )

  expect_error(
    calculate_ors(bad_results),
    "ccs_results must contain 'CCS' column"
  )

  # Test insufficient ortholog pairs
  few_results <- data.frame(
    gene_sp1 = paste0("A", 1:5),
    gene_sp2 = paste0("B", 1:5),
    CCS = runif(5)
  )

  expect_error(
    calculate_ors(few_results),
    "Need at least 10 ortholog pairs"
  )
})

test_that("test_ors_significance validates inputs", {
  # Test missing 'ORS' column
  bad_results <- data.frame(
    gene_sp1 = c("A", "B"),
    gene_sp2 = c("C", "D")
  )

  expect_error(
    test_ors_significance(bad_results),
    "ors_results must contain 'ORS'"
  )
})

test_that("calculate_ors handles edge cases", {
  # Test with perfect CCS = 1 for all
  perfect_ccs <- data.frame(
    gene_sp1 = paste0("Gene", 1:20),
    gene_sp2 = paste0("Gene", 1:20),
    CCS = rep(1, 20),
    n_ref = 20
  )

  ors_perfect <- calculate_ors(perfect_ccs, return_log = FALSE)

  # All should have equal ORS (ties)
  expect_true(all(ors_perfect$ORS > 0))

  # Test with all negative CCS
  negative_ccs <- data.frame(
    gene_sp1 = paste0("Gene", 1:20),
    gene_sp2 = paste0("Gene", 1:20),
    CCS = runif(20, -1, 0),
    n_ref = 20
  )

  ors_negative <- calculate_ors(negative_ccs, return_log = TRUE)

  # Should still compute ORS correctly
  expect_true(all(!is.na(ors_negative$ORS)))
  expect_true(all(!is.na(ors_negative$logORS)))
})
