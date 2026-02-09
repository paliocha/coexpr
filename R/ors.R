#' Calculate ortholog rank score (ORS)
#'
#' Calculates ORS which measures the statistical significance of co-expression
#' conservation. Higher ORS indicates stronger evidence that orthologous genes
#' have conserved co-expression partners.
#'
#' @param ccs_results Data frame from `calculate_ccs()` with CCS values
#' @param return_log Logical. If TRUE (default), returns logORS which is more
#'   suitable for visualization. If FALSE, returns raw ORS proportion.
#' @param directional Logical. If TRUE, also computes directional ORS columns
#'   (`ORS_sp1_to_sp2`, `ORS_sp2_to_sp1`) for multi-copy ortholog analysis.
#'   Default is FALSE.
#'
#' @return Input data frame with added columns:
#'   - `ORS`: Global rank-based ORS (recommended for analysis)
#'   - `logORS` (if return_log=TRUE): -log10 transformed ORS
#'   - `ORS_sp1_to_sp2`, `ORS_sp2_to_sp1` (if directional=TRUE):
#'     within-group ORS for multi-copy analysis
#'
#' @details
#' ORS is calculated as the proportion of ortholog pairs with CCS less than or
#' equal to the CCS of the focal ortholog pair: `ORS_ij = P(CCS <= CCS_ij)`
#'
#' For multi-copy orthologs (1:N, N:1, N:M), directional ORS is calculated
#' separately for each species direction by grouping on the gene from that
#' species. For 1:1 orthologs, both directions use global ranking.
#'
#' The mean bidirectional ORS is then log-transformed:
#' `logORS_ij = -log10(1 + 10^-4 - (ORS^sp1_ij + ORS^sp2_ji)/2)`
#'
#' Interpretation of logORS:
#' \itemize{
#'   \item >1: Ortholog is in top 10
#'   \item >2: Top 1
#'   \item >3: Top 0.1
#'   \item <0: Below median (diverged)
#' }
#'
#' @references
#' Grønvold L, Hvidsten TR. Cross species comparative transcriptomics using
#' co-expression networks.
#'
#' @examples
#' \dontrun{
#' ccs_results <- calculate_ccs(sim_at, sim_os, orthologs)
#' ors_results <- calculate_ors(ccs_results)
#'
#' # Filter highly conserved orthologs
#' conserved <- ors_results |> dplyr::filter(logORS > 2)
#' }
#'
#' @export
calculate_ors <- function(ccs_results, return_log = TRUE, directional = FALSE) {

  if (!"CCS" %in% colnames(ccs_results)) {
    stop("ccs_results must contain 'CCS' column from calculate_ccs()")
  }

  if (nrow(ccs_results) < 10) {
    stop("Need at least 10 ortholog pairs to calculate ORS")
  }

  n_total <- nrow(ccs_results)

  # Calculate ORS: proportion of ortholog pairs with CCS <= CCS of focal pair
  #
  # The ORS measures how well an ortholog pair's CCS compares to all other pairs.
  # We use GLOBAL ranking across all ortholog pairs to get meaningful statistics.
  #
  # For multi-copy orthologs, we also compute within-group rankings:
  # - ORS_sp1_to_sp2: rank within all pairs sharing the same gene_sp1
  # - ORS_sp2_to_sp1: rank within all pairs sharing the same gene_sp2
  # These indicate which copy is best when a gene has multiple orthologs.
  #
  # The primary metric (ORS) uses GLOBAL ranking for interpretability.

  # Global ORS: rank against ALL ortholog pairs
  # na.last = "keep" preserves NA for pairs with missing CCS
  global_ranks <- rank(ccs_results$CCS, ties.method = "average", na.last = "keep")
  ors_global <- global_ranks / sum(!is.na(ccs_results$CCS))

  ors_results <- ccs_results |>
    dplyr::mutate(ORS = ors_global)

  # Optionally compute directional ORS for multi-copy analysis
  if (directional) {
    # Species 1 -> 2: For each gene in sp1, rank among all its ortholog pairs
    ors_results <- ors_results |>
      dplyr::group_by(.data$gene_sp1) |>
      dplyr::mutate(
        ORS_sp1_to_sp2 = rank(.data$CCS, ties.method = "average") / dplyr::n()
      ) |>
      dplyr::ungroup()

    # Species 2 -> 1: For each gene in sp2, rank among all its ortholog pairs
    ors_results <- ors_results |>
      dplyr::group_by(.data$gene_sp2) |>
      dplyr::mutate(
        ORS_sp2_to_sp1 = rank(.data$CCS, ties.method = "average") / dplyr::n()
      ) |>
      dplyr::ungroup()
  }

  # Log transformation for better visualization
  if (return_log) {
    ors_results <- ors_results |>
      dplyr::mutate(
        logORS = -log10(1 + 1e-4 - .data$ORS)
      )
  }

  return(ors_results)
}


#' Test significance of ORS
#'
#' Tests whether an ortholog pair has significantly conserved co-expression
#' compared to random expectation.
#'
#' @param ors_results Data frame from `calculate_ors()`
#' @param alpha Significance level (default 0.05)
#' @param p_adjust_method Method for multiple testing correction, passed to
#'   [stats::p.adjust()]. Default `"BH"` (Benjamini-Hochberg FDR). Use
#'   `"none"` to skip adjustment and test on raw p-values.
#'
#' @return Input data frame with added columns:
#'   - `pvalue`: raw empirical p-value (1 - ORS)
#'   - `padj`: adjusted p-value (only when `p_adjust_method != "none"`)
#'   - `significant`: logical, TRUE if adjusted p-value < alpha
#'     (or raw p-value when `p_adjust_method = "none"`)
#'
#' @details
#' The null hypothesis is that the ortholog has no more conserved co-expression
#' than a random gene pair. Under the null, ORS is uniform on \[0, 1\], so the
#' raw p-value is simply `1 - ORS`.
#'
#' Because genome-scale analyses test thousands of ortholog pairs
#' simultaneously, multiple testing correction is applied by default using
#' the Benjamini-Hochberg procedure to control the false discovery rate.
#'
#' @export
test_ors_significance <- function(ors_results, alpha = 0.05,
                                  p_adjust_method = "BH") {

  if (!"ORS" %in% colnames(ors_results)) {
    stop("ors_results must contain 'ORS' from calculate_ors()")
  }

  p_adjust_method <- match.arg(p_adjust_method,
                               choices = stats::p.adjust.methods)

  # Empirical p-value: under null ORS ~ Uniform[0,1], so p = 1 - ORS
  ors_results <- ors_results |>
    dplyr::mutate(pvalue = 1 - .data$ORS)

  if (p_adjust_method == "none") {
    ors_results <- ors_results |>
      dplyr::mutate(significant = .data$pvalue < alpha)
  } else {
    ors_results <- ors_results |>
      dplyr::mutate(
        padj = stats::p.adjust(.data$pvalue, method = p_adjust_method),
        significant = .data$padj < alpha
      )
  }

  n_sig <- sum(ors_results$significant, na.rm = TRUE)
  method_label <- if (p_adjust_method == "none") "unadjusted" else p_adjust_method
  message(sprintf("%d / %d orthologs are significant (alpha = %.3f, method = %s)",
                  n_sig, nrow(ors_results), alpha, method_label))

  return(ors_results)
}
