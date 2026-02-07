#' Calculate ortholog rank score (ORS)
#'
#' Calculates ORS which measures the statistical significance of co-expression
#' conservation. Higher ORS indicates stronger evidence that orthologous genes
#' have conserved co-expression partners.
#'
#' @param ccs_results Data frame from `calculate_ccs()` with CCS values
#' @param return_log Logical. If TRUE (default), returns logORS which is more
#'   suitable for visualization. If FALSE, returns raw ORS proportion.
#'
#' @return Input data frame with added columns:
#'   - `ORS_sp1_to_sp2`: ORS from species 1 perspective (for multi-copy)
#'   - `ORS_sp2_to_sp1`: ORS from species 2 perspective (for multi-copy)
#'   - `ORS`: Mean of bidirectional ORS (recommended for analysis)
#'   - `logORS` (if return_log=TRUE): -log10 transformed ORS
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
calculate_ors <- function(ccs_results, return_log = TRUE) {

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
  global_ranks <- rank(ccs_results$CCS, ties.method = "average")
  ors_global <- global_ranks / n_total

  # Also compute directional ORS for multi-copy analysis
  # Species 1 -> 2: For each gene in sp1, rank among all its ortholog pairs
  ors_sp1_to_sp2 <- ccs_results |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::mutate(
      ORS_sp1_to_sp2 = rank(.data$CCS, ties.method = "average") / dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::pull(.data$ORS_sp1_to_sp2)

  # Species 2 -> 1: For each gene in sp2, rank among all its ortholog pairs
  ors_sp2_to_sp1 <- ccs_results |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::mutate(
      ORS_sp2_to_sp1 = rank(.data$CCS, ties.method = "average") / dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::pull(.data$ORS_sp2_to_sp1)

  # Combine results
  ors_results <- ccs_results |>
    dplyr::mutate(
      ORS_sp1_to_sp2 = ors_sp1_to_sp2,
      ORS_sp2_to_sp1 = ors_sp2_to_sp1,
      # Use global ORS as the primary metric for interpretability
      ORS = ors_global
    )

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
#'
#' @return Input data frame with added column `significant` (logical)
#'
#' @details
#' The null hypothesis is that the ortholog has no more conserved co-expression
#' than a random gene pair. This corresponds to ORS = 0.5 (median rank).
#'
#' The alternative hypothesis is that ORS > 0.5 (better than median).
#'
#' @export
test_ors_significance <- function(ors_results, alpha = 0.05) {

  if (!"ORS" %in% colnames(ors_results)) {
    stop("ors_results must contain 'ORS' from calculate_ors()")
  }

  # Calculate empirical P-value: P(ORS >= observed | null: ORS = 0.5)
  # Under null, ORS should be uniform[0,1], so P-value = 1 - ORS
  ors_results <- ors_results |>
    dplyr::mutate(
      pvalue = 1 - .data$ORS,
      significant = .data$pvalue < alpha
    )

  n_sig <- sum(ors_results$significant, na.rm = TRUE)
  message(sprintf("%d / %d orthologs are significant (alpha = %.3f)",
                  n_sig, nrow(ors_results), alpha))

  return(ors_results)
}
