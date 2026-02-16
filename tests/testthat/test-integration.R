# Integration tests for the full coexpr pipeline
#
# Tests the end-to-end flow:
#   expression data → similarity matrices → CCS → ORS → significance → multicopy

# ---------------------------------------------------------------------------
# Shared test data: two "species" with known co-expression structure
# ---------------------------------------------------------------------------

# Create reproducible expression matrices with planted co-expression
set.seed(42)
n_samples <- 30

# Species 1: 30 genes, 30 samples
# Genes 1-10 form a co-expressed cluster
base_signal <- matrix(rnorm(10 * n_samples), nrow = 10, ncol = n_samples)
noise <- function(n, sd = 0.5) matrix(rnorm(n * n_samples, sd = sd), nrow = n, ncol = n_samples)

# Cluster 1 (genes 1-10): correlated via shared signal
shared1 <- matrix(rep(rnorm(n_samples), 10), nrow = 10, byrow = TRUE)
expr_sp1_cluster1 <- shared1 + noise(10, sd = 0.3)

# Cluster 2 (genes 11-20): different shared signal
shared2 <- matrix(rep(rnorm(n_samples), 10), nrow = 10, byrow = TRUE)
expr_sp1_cluster2 <- shared2 + noise(10, sd = 0.3)

# Independent genes (21-30)
expr_sp1_indep <- noise(10, sd = 1.0)

expr_sp1 <- rbind(expr_sp1_cluster1, expr_sp1_cluster2, expr_sp1_indep)
rownames(expr_sp1) <- paste0("SP1_G", 1:30)
colnames(expr_sp1) <- paste0("S", 1:n_samples)

# Species 2: 35 genes (some paralogs), same cluster structure for orthologs
shared1b <- matrix(rep(rnorm(n_samples), 10), nrow = 10, byrow = TRUE)
expr_sp2_cluster1 <- shared1b + noise(10, sd = 0.3)

shared2b <- matrix(rep(rnorm(n_samples), 10), nrow = 10, byrow = TRUE)
expr_sp2_cluster2 <- shared2b + noise(10, sd = 0.3)

expr_sp2_indep <- noise(15, sd = 1.0)

expr_sp2 <- rbind(expr_sp2_cluster1, expr_sp2_cluster2, expr_sp2_indep)
rownames(expr_sp2) <- paste0("SP2_G", 1:35)
colnames(expr_sp2) <- paste0("S", 1:n_samples)

# Make some ortholog pairs have truly conserved co-expression by sharing signal
# Genes in cluster 1 of both species share the same underlying structure
conserved_signal <- matrix(rep(rnorm(n_samples), 10), nrow = 10, byrow = TRUE)
expr_sp1[1:10, ] <- conserved_signal + noise(10, sd = 0.3)
expr_sp2[1:10, ] <- conserved_signal + noise(10, sd = 0.3)

# Orthologs: 20 1:1 pairs + some multi-copy
orthologs <- data.frame(
  gene_sp1 = c(
    paste0("SP1_G", 1:20),          # 1:1 orthologs
    "SP1_G21", "SP1_G21",           # 1:N (SP1_G21 -> SP2_G21, SP2_G31)
    "SP1_G22", "SP1_G23"            # N:1 (SP1_G22, SP1_G23 -> SP2_G22)
  ),
  gene_sp2 = c(
    paste0("SP2_G", 1:20),          # 1:1 orthologs
    "SP2_G21", "SP2_G31",           # 1:N copies
    "SP2_G22", "SP2_G22"            # N:1 copies
  ),
  stringsAsFactors = FALSE
)


# ===========================================================================
# Test 1: Full pipeline with PCC similarity
# ===========================================================================
test_that("full pipeline works: PCC → CCS → ORS → significance", {
  # Step 1: Similarity matrices (PCC only, no MR, for speed)
  sim_sp1 <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = TRUE)
  sim_sp2 <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = TRUE)

  expect_s4_class(sim_sp1, "TriSimilarity")
  expect_s4_class(sim_sp2, "TriSimilarity")
  expect_equal(sim_sp1@n, 30L)
  expect_equal(sim_sp2@n, 35L)

  # Step 2: CCS
  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)

  expect_true("CCS" %in% colnames(ccs_results))
  expect_true("n_ref" %in% colnames(ccs_results))
  expect_true(all(is.finite(ccs_results$CCS)))
  expect_true(all(ccs_results$CCS >= -1 & ccs_results$CCS <= 1))

  # Conserved orthologs (cluster 1) should have higher CCS than random
  cluster1_ccs <- ccs_results$CCS[ccs_results$gene_sp1 %in% paste0("SP1_G", 1:10)]
  other_ccs <- ccs_results$CCS[!ccs_results$gene_sp1 %in% paste0("SP1_G", 1:10)]
  expect_gt(mean(cluster1_ccs), mean(other_ccs))

  # Step 3: ORS
  ors_results <- calculate_ors(ccs_results)

  expect_true("ORS" %in% colnames(ors_results))
  expect_true("logORS" %in% colnames(ors_results))
  expect_true(all(ors_results$ORS >= 0 & ors_results$ORS <= 1))

  # Step 4: Significance testing
  sig_results <- test_ors_significance(ors_results)

  expect_true("pvalue" %in% colnames(sig_results))
  expect_true("padj" %in% colnames(sig_results))
  expect_true(all(sig_results$pvalue >= 0 & sig_results$pvalue <= 1))
  expect_true(all(sig_results$padj >= 0 & sig_results$padj <= 1))
})


# ===========================================================================
# Test 2: Full pipeline with PCC+MR similarity
# ===========================================================================
test_that("full pipeline works: PCC+MR → CCS → ORS", {
  sim_sp1 <- calculate_pcc_mr(expr_sp1, method = "pcc_mr", return_tri = TRUE)
  sim_sp2 <- calculate_pcc_mr(expr_sp2, method = "pcc_mr", return_tri = TRUE)

  # MR values should be in [0, 1]
  expect_true(all(sim_sp1@data >= 0 & sim_sp1@data <= 1))
  expect_true(all(sim_sp2@data >= 0 & sim_sp2@data <= 1))

  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)
  ors_results <- calculate_ors(ccs_results)

  # Basic sanity
  expect_equal(nrow(ors_results), nrow(ccs_results))
  expect_true(all(is.finite(ors_results$logORS)))
})


# ===========================================================================
# Test 3: Full pipeline with MI+CLR similarity
# ===========================================================================
test_that("full pipeline works: MI+CLR → CCS → ORS", {
  sim_sp1 <- calculate_mi_clr(expr_sp1, n_bins = 5, return_tri = TRUE)
  sim_sp2 <- calculate_mi_clr(expr_sp2, n_bins = 5, return_tri = TRUE)

  expect_s4_class(sim_sp1, "TriSimilarity")
  expect_true(all(sim_sp1@data >= 0))  # CLR values are non-negative

  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)
  ors_results <- calculate_ors(ccs_results)

  expect_equal(nrow(ors_results), nrow(ccs_results))
  expect_true(all(is.finite(ors_results$logORS)))
})


# ===========================================================================
# Test 4: Pipeline with full matrix (non-TriSimilarity) input
# ===========================================================================
test_that("pipeline works with full matrix similarity input", {
  sim_sp1 <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = FALSE)
  sim_sp2 <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = FALSE)

  expect_true(is.matrix(sim_sp1))
  expect_true(is.matrix(sim_sp2))

  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)
  ors_results <- calculate_ors(ccs_results)

  expect_equal(nrow(ors_results), nrow(ccs_results))
  expect_true(all(is.finite(ors_results$logORS)))
})


# ===========================================================================
# Test 5: CCS results are consistent between TriSimilarity and full matrix
# ===========================================================================
test_that("TriSimilarity and full matrix give identical CCS", {
  sim_sp1_tri <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = TRUE)
  sim_sp2_tri <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = TRUE)
  sim_sp1_mat <- as.matrix(sim_sp1_tri)
  sim_sp2_mat <- as.matrix(sim_sp2_tri)

  ccs_tri <- calculate_ccs(sim_sp1_tri, sim_sp2_tri, orthologs)
  ccs_mat <- calculate_ccs(sim_sp1_mat, sim_sp2_mat, orthologs)

  expect_equal(ccs_tri$CCS, ccs_mat$CCS, tolerance = 1e-10)
})


# ===========================================================================
# Test 6: Multicopy integration — strict vs best_hit vs collapse
# ===========================================================================
test_that("multicopy strategies integrate with CCS/ORS pipeline", {
  sim_sp1 <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = TRUE)
  sim_sp2 <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = TRUE)

  # Strict: only 1:1
  strict <- handle_multicopy_orthologs(orthologs, strategy = "strict",
                                       similarity_sp1 = sim_sp1,
                                       similarity_sp2 = sim_sp2)
  expect_true(all(strict$type == "1:1"))

  ccs_strict <- calculate_ccs(sim_sp1, sim_sp2, strict)
  ors_strict <- calculate_ors(ccs_strict)
  expect_equal(nrow(ors_strict), nrow(strict))

  # Best hit: resolves multi-copy to 1:1
  best <- handle_multicopy_orthologs(orthologs, strategy = "best_hit",
                                     similarity_sp1 = sim_sp1,
                                     similarity_sp2 = sim_sp2)
  expect_gte(nrow(best), nrow(strict))

  ccs_best <- calculate_ccs(sim_sp1, sim_sp2, best)
  ors_best <- calculate_ors(ccs_best)
  expect_equal(nrow(ors_best), nrow(best))

  # All pairs: keeps everything
  all_p <- handle_multicopy_orthologs(orthologs, strategy = "all_pairs",
                                      similarity_sp1 = sim_sp1,
                                      similarity_sp2 = sim_sp2)
  expect_equal(nrow(all_p), nrow(orthologs))

  ccs_all <- calculate_ccs(sim_sp1, sim_sp2, all_p)
  ors_all <- calculate_ors(ccs_all)
  expect_equal(nrow(ors_all), nrow(all_p))
})


# ===========================================================================
# Test 7: collapse_orthologs integration
# ===========================================================================
test_that("collapse_orthologs integrates with CCS/ORS pipeline", {
  sim_sp1 <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = TRUE)
  sim_sp2 <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = TRUE)

  collapsed <- collapse_orthologs(
    orthologs = orthologs,
    similarity_sp1 = sim_sp1,
    similarity_sp2 = sim_sp2,
    multicopy_sp = "both"
  )

  expect_true("homeolog_score" %in% colnames(collapsed))
  expect_true("n_candidates" %in% colnames(collapsed))

  # Should resolve multi-copy groups
  expect_lte(nrow(collapsed), nrow(orthologs))

  # Feed into CCS/ORS
  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, collapsed)
  ors_results <- calculate_ors(ccs_results)

  expect_equal(nrow(ors_results), nrow(collapsed))
  expect_true(all(is.finite(ors_results$logORS)))
})


# ===========================================================================
# Test 8: Caching round-trip preserves results
# ===========================================================================
test_that("cached similarity matrices produce identical CCS", {
  cache_dir <- tempfile("cache_test_")

  sim_sp1 <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = TRUE,
                               cache_dir = cache_dir)
  sim_sp2 <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = TRUE,
                               cache_dir = cache_dir)

  # Reload from cache
  sim_sp1_cached <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = TRUE,
                                      cache_dir = cache_dir)
  sim_sp2_cached <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = TRUE,
                                      cache_dir = cache_dir)

  # CCS should be identical
  ccs_orig <- calculate_ccs(sim_sp1, sim_sp2, orthologs)
  ccs_cached <- calculate_ccs(sim_sp1_cached, sim_sp2_cached, orthologs)

  expect_equal(ccs_orig$CCS, ccs_cached$CCS, tolerance = 1e-14)

  unlink(cache_dir, recursive = TRUE)
})


# ===========================================================================
# Test 9: Self-diagonal handling options produce valid results
# ===========================================================================
test_that("all self-diagonal handling options work end-to-end", {
  sim_sp1 <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = TRUE)
  sim_sp2 <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = TRUE)

  for (diag_method in c("mean", "na", "none")) {
    ccs <- calculate_ccs(sim_sp1, sim_sp2, orthologs,
                         handle_self_diagonal = diag_method)
    expect_true(all(is.finite(ccs$CCS)),
                info = paste("handle_self_diagonal =", diag_method))
    expect_true(all(ccs$CCS >= -1 & ccs$CCS <= 1),
                info = paste("handle_self_diagonal =", diag_method))
  }
})


# ===========================================================================
# Test 10: summarize_conservation works on real pipeline output
# ===========================================================================
test_that("summarize_conservation works on full pipeline output", {
  sim_sp1 <- calculate_pcc_mr(expr_sp1, method = "pcc", return_tri = TRUE)
  sim_sp2 <- calculate_pcc_mr(expr_sp2, method = "pcc", return_tri = TRUE)

  all_p <- handle_multicopy_orthologs(orthologs, strategy = "all_pairs",
                                      similarity_sp1 = sim_sp1,
                                      similarity_sp2 = sim_sp2)
  ccs_results <- calculate_ccs(sim_sp1, sim_sp2, all_p)
  ors_results <- calculate_ors(ccs_results)

  # Summarize by type
  summary <- summarize_conservation(ors_results, by_type = TRUE)
  expect_true("type" %in% colnames(summary))
  expect_true("median_logORS" %in% colnames(summary))
  expect_true("pct_top10" %in% colnames(summary))

  # Overall summary
  summary_all <- summarize_conservation(ors_results, by_type = FALSE)
  expect_equal(nrow(summary_all), 1)
})
