#' Handle multi-copy orthologs
#'
#' Resolves how to handle 1:N and N:M ortholog relationships for CCS calculation.
#' The paper (Grønvold & Hvidsten) uses only 1:1 orthologs, but this function
#' provides several strategies for incorporating multi-copy orthologs.
#'
#' @param orthologs Data frame with columns `gene_sp1` and `gene_sp2`.
#'   Optionally include `type` column with values "1:1", "1:N", "N:1", "N:M".
#' @param strategy Character. Strategy for handling multi-copy orthologs:
#'   - `"strict"` (default): Keep only 1:1 orthologs (matches paper)
#'   - `"best_hit"`: For 1:N, select best copy; for N:M, select best reciprocal pair
#'   - `"mean"`: Average CCS across all copies (requires CCS pre-calculated)
#'   - `"max"`: Use maximum CCS (most conserved copy, requires CCS pre-calculated)
#'   - `"all_pairs"`: Keep all combinations (increases dataset size)
#' @param ccs_values Optional. Pre-calculated CCS values (required for "mean" and "max")
#' @param similarity_sp1 Optional. Similarity matrix for species 1 (for "best_hit")
#' @param similarity_sp2 Optional. Similarity matrix for species 2 (for "best_hit")
#'
#' @return Processed orthologs data frame suitable for CCS calculation
#'
#' @details
#' **Paper's findings**: Duplicated genes show progressively lower expression
#' conservation (1:2 < 1:3 < 1:4), EXCEPT for very recent WGD duplicates
#' (~13 Mya in Glycine max) which show no divergence.
#'
#' **Strategy recommendations**:
#' - Use `"strict"` for main analysis (most conservative, publication-ready)
#' - Use `"best_hit"` for exploratory analysis when need more gene coverage
#' - Use `"all_pairs"` to study expression divergence among paralogs
#' - Avoid `"mean"` unless you have biological reason to believe partial sub-functionalization
#'
#' @references
#' Grønvold L, Hvidsten TR. Cross species comparative transcriptomics using
#' co-expression networks.
#'
#' Ohno S (1970) Evolution by Gene Duplication.
#'
#' @examples
#' \dontrun{
#' # Load orthologs (may contain 1:N and N:M relationships)
#' orthologs_all <- read.csv("orthologs.csv")
#'
#' # Strategy 1: Strict (only 1:1, matches paper)
#' orthologs_1to1 <- handle_multicopy_orthologs(
#'   orthologs_all,
#'   strategy = "strict"
#' )
#'
#' # Strategy 2: Best hit (when need more coverage)
#' orthologs_best <- handle_multicopy_orthologs(
#'   orthologs_all,
#'   strategy = "best_hit",
#'   similarity_sp1 = sim_at,
#'   similarity_sp2 = sim_os
#' )
#'
#' # Strategy 3: All pairs (study paralog divergence)
#' orthologs_all_pairs <- handle_multicopy_orthologs(
#'   orthologs_all,
#'   strategy = "all_pairs"
#' )
#' }
#'
#' @export
handle_multicopy_orthologs <- function(orthologs,
                                       strategy = c("strict", "best_hit", "mean",
                                                    "max", "all_pairs"),
                                       ccs_values = NULL,
                                       similarity_sp1 = NULL,
                                       similarity_sp2 = NULL) {

  strategy <- match.arg(strategy)

  # Validate input
  required_cols <- c("gene_sp1", "gene_sp2")
  if (!all(required_cols %in% colnames(orthologs))) {
    stop(sprintf("orthologs must have columns: %s",
                 paste(required_cols, collapse = ", ")))
  }

  # Detect ortholog types if not provided
  if (!"type" %in% colnames(orthologs)) {
    orthologs <- detect_ortholog_types(orthologs)
  }

  # Count ortholog types
  type_counts <- table(orthologs$type)
  message(sprintf("Ortholog composition: %s",
                  paste(names(type_counts), "=", type_counts, collapse = ", ")))

  # Apply strategy
  result <- switch(strategy,
    "strict" = filter_1to1_only(orthologs),
    "best_hit" = select_best_hits(orthologs, similarity_sp1, similarity_sp2),
    "mean" = aggregate_by_mean(orthologs, ccs_values),
    "max" = aggregate_by_max(orthologs, ccs_values),
    "all_pairs" = keep_all_pairs(orthologs)
  )

  message(sprintf("Strategy '%s': Retained %d / %d ortholog pairs",
                  strategy, nrow(result), nrow(orthologs)))

  return(result)
}


#' Detect ortholog relationship types
#'
#' @param orthologs Data frame with gene_sp1 and gene_sp2
#' @return Data frame with added 'type' column
#' @keywords internal
#' @noRd
detect_ortholog_types <- function(orthologs) {

  # Count occurrences of each gene
  counts_sp1 <- orthologs |>
    dplyr::count(.data$gene_sp1, name = "n_sp1")

  counts_sp2 <- orthologs |>
    dplyr::count(.data$gene_sp2, name = "n_sp2")

  # Classify ortholog types
  orthologs_typed <- orthologs |>
    dplyr::left_join(counts_sp1, by = "gene_sp1") |>
    dplyr::left_join(counts_sp2, by = "gene_sp2") |>
    dplyr::mutate(
      type = dplyr::case_when(
        .data$n_sp1 == 1 & .data$n_sp2 == 1 ~ "1:1",
        .data$n_sp1 == 1 & .data$n_sp2 > 1  ~ "1:N",
        .data$n_sp1 > 1  & .data$n_sp2 == 1 ~ "N:1",
        .data$n_sp1 > 1  & .data$n_sp2 > 1  ~ "N:M",
        TRUE ~ "unknown"
      )
    ) |>
    dplyr::select(-"n_sp1", -"n_sp2")

  return(orthologs_typed)
}


#' Filter to 1:1 orthologs only (strict strategy)
#'
#' @keywords internal
#' @noRd
filter_1to1_only <- function(orthologs) {
  orthologs |>
    dplyr::filter(.data$type == "1:1")
}


#' Select best hit for multi-copy orthologs
#'
#' @keywords internal
#' @noRd
select_best_hits <- function(orthologs, similarity_sp1, similarity_sp2) {

  if (is.null(similarity_sp1) || is.null(similarity_sp2)) {
    stop("Strategy 'best_hit' requires similarity_sp1 and similarity_sp2")
  }

  # For 1:1, keep as is
  orthologs_1to1 <- orthologs |>
    dplyr::filter(.data$type == "1:1")

  # For 1:N (sp1 has one copy, sp2 has multiple)
  # Select sp2 copy with highest average similarity to sp1
  orthologs_1toN <- orthologs |>
    dplyr::filter(.data$type == "1:N") |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::slice_max(order_by = calculate_avg_similarity(
      .data$gene_sp2, similarity_sp2
    ), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # For N:1 (sp1 has multiple, sp2 has one)
  orthologs_Nto1 <- orthologs |>
    dplyr::filter(.data$type == "N:1") |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::slice_max(order_by = calculate_avg_similarity(
      .data$gene_sp1, similarity_sp1
    ), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # For N:M, select best reciprocal pair based on combined similarity
  orthologs_NtoM <- orthologs |>
    dplyr::filter(.data$type == "N:M") |>
    dplyr::mutate(
      score = calculate_avg_similarity(.data$gene_sp1, similarity_sp1) +
              calculate_avg_similarity(.data$gene_sp2, similarity_sp2)
    ) |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::slice_max(order_by = .data$score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::slice_max(order_by = .data$score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-.data$score)

  # Combine all
  dplyr::bind_rows(orthologs_1to1, orthologs_1toN, orthologs_Nto1, orthologs_NtoM)
}


#' Calculate average similarity score for a gene
#'
#' @keywords internal
#' @noRd
calculate_avg_similarity <- function(gene, similarity_matrix) {
  # Vectorize for use with dplyr
  vapply(gene, function(g) {
    if (g %in% rownames(similarity_matrix)) {
      mean(similarity_matrix[g, ], na.rm = TRUE)
    } else {
      0
    }
  }, numeric(1))
}


#' Aggregate multi-copy orthologs by mean CCS
#'
#' @keywords internal
#' @noRd
aggregate_by_mean <- function(orthologs, ccs_values) {
  if (is.null(ccs_values)) {
    stop("Strategy 'mean' requires pre-calculated ccs_values")
  }

  # Merge with CCS values
  orthologs_with_ccs <- orthologs |>
    dplyr::left_join(ccs_values, by = c("gene_sp1", "gene_sp2"))

  # For multi-copy, average CCS
  orthologs_with_ccs |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::summarize(
      gene_sp2 = .data$gene_sp2[1],  # Keep first for simplicity
      ccs = mean(.data$ccs, na.rm = TRUE),
      type = "aggregated",
      .groups = "drop"
    )
}


#' Aggregate multi-copy orthologs by max CCS
#'
#' @keywords internal
#' @noRd
aggregate_by_max <- function(orthologs, ccs_values) {
  if (is.null(ccs_values)) {
    stop("Strategy 'max' requires pre-calculated ccs_values")
  }

  # Merge with CCS values
  orthologs_with_ccs <- orthologs |>
    dplyr::left_join(ccs_values, by = c("gene_sp1", "gene_sp2"))

  # For each gene in sp1, keep copy with max CCS
  orthologs_with_ccs |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::slice_max(order_by = .data$ccs, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}


#' Keep all ortholog pairs
#'
#' @keywords internal
#' @noRd
keep_all_pairs <- function(orthologs) {
  # Already have all pairs, just return
  orthologs
}
