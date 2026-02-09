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
    "all_pairs" = orthologs
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

  # Type naming convention from detect_ortholog_types:
  #   "N:1" = n_sp1 > 1 (sp1 gene has N sp2 partners), n_sp2 == 1
  #   "1:N" = n_sp1 == 1, n_sp2 > 1 (sp2 gene has N sp1 partners)

  # For 1:1, keep as is
  orthologs_1to1 <- orthologs |>
    dplyr::filter(.data$type == "1:1")

  # N:1 — sp1 gene has N sp2 partners: pick best sp2 by average similarity
  orthologs_Nto1 <- orthologs |>
    dplyr::filter(.data$type == "N:1") |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::slice_max(order_by = calculate_avg_similarity(
      .data$gene_sp2, similarity_sp2
    ), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # 1:N — sp2 gene has N sp1 partners: pick best sp1 by average similarity
  orthologs_1toN <- orthologs |>
    dplyr::filter(.data$type == "1:N") |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::slice_max(order_by = calculate_avg_similarity(
      .data$gene_sp1, similarity_sp1
    ), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # N:M — two-pass: best sp2 per sp1, then best sp1 per sp2
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
    dplyr::select(-"score")

  # Combine all
  dplyr::bind_rows(orthologs_1to1, orthologs_Nto1, orthologs_1toN, orthologs_NtoM)
}


#' Calculate average similarity score for a gene
#'
#' @keywords internal
#' @noRd
calculate_avg_similarity <- function(gene, similarity_matrix) {
  # Get gene names depending on object type
  sim_genes <- if (is(similarity_matrix, "TriSimilarity")) {
    similarity_matrix@genes
  } else {
    rownames(similarity_matrix)
  }

  # Vectorize for use with dplyr
  vapply(gene, function(g) {
    if (g %in% sim_genes) {
      if (is(similarity_matrix, "TriSimilarity")) {
        mean(extractColumn(similarity_matrix, g), na.rm = TRUE)
      } else {
        mean(similarity_matrix[g, ], na.rm = TRUE)
      }
    } else {
      0
    }
  }, numeric(1))
}


#' Aggregate multi-copy orthologs by mean CCS
#'
#' For 1:N orthologs, averages CCS across sp2 copies (one row per sp1 gene).
#' For N:1, averages across sp1 copies (one row per sp2 gene).
#' For N:M, averages across all sp2 copies per sp1 gene, then across
#' sp1 copies per sp2 gene (two-pass reduction).
#'
#' @keywords internal
#' @noRd
aggregate_by_mean <- function(orthologs, ccs_values) {
  if (is.null(ccs_values)) {
    stop("Strategy 'mean' requires pre-calculated ccs_values")
  }

  # Ensure type detection
  if (!"type" %in% colnames(orthologs)) {
    orthologs <- detect_ortholog_types(orthologs)
  }

  # Merge with CCS values
  orthologs_with_ccs <- orthologs |>
    dplyr::left_join(ccs_values, by = c("gene_sp1", "gene_sp2"))

  # 1:1 — keep as-is
  one_to_one <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "1:1")

  # N:1 — sp1 gene has N sp2 partners: group by gene_sp1, average across sp2 copies
  n_to_one <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "N:1") |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::summarize(
      gene_sp2 = paste(.data$gene_sp2, collapse = ";"),
      CCS = mean(.data$CCS, na.rm = TRUE),
      type = "aggregated_N:1",
      .groups = "drop"
    )

  # 1:N — sp2 gene has N sp1 partners: group by gene_sp2, average across sp1 copies
  one_to_n <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "1:N") |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::summarize(
      gene_sp1 = paste(.data$gene_sp1, collapse = ";"),
      CCS = mean(.data$CCS, na.rm = TRUE),
      type = "aggregated_1:N",
      .groups = "drop"
    )

  # N:M — average across sp2 copies per sp1 gene
  n_to_m <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "N:M")

  if (nrow(n_to_m) > 0) {
    n_to_m <- n_to_m |>
      dplyr::group_by(.data$gene_sp1) |>
      dplyr::summarize(
        gene_sp2 = paste(.data$gene_sp2, collapse = ";"),
        CCS = mean(.data$CCS, na.rm = TRUE),
        type = "aggregated_N:M",
        .groups = "drop"
      )
  }

  dplyr::bind_rows(one_to_one, n_to_one, one_to_n, n_to_m)
}


#' Aggregate multi-copy orthologs by max CCS
#'
#' For 1:N orthologs, selects the sp2 copy with the highest CCS.
#' For N:1, selects the sp1 copy with the highest CCS.
#' For N:M, two-pass selection: best sp2 per sp1, then best sp1 per sp2.
#'
#' @keywords internal
#' @noRd
aggregate_by_max <- function(orthologs, ccs_values) {
  if (is.null(ccs_values)) {
    stop("Strategy 'max' requires pre-calculated ccs_values")
  }

  # Ensure type detection
  if (!"type" %in% colnames(orthologs)) {
    orthologs <- detect_ortholog_types(orthologs)
  }

  # Merge with CCS values
  orthologs_with_ccs <- orthologs |>
    dplyr::left_join(ccs_values, by = c("gene_sp1", "gene_sp2"))

  # 1:1 — keep as-is
  one_to_one <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "1:1")

  # N:1 — sp1 gene has N sp2 partners: keep best sp2 copy per sp1 gene
  n_to_one <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "N:1") |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::slice_max(order_by = .data$CCS, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # 1:N — sp2 gene has N sp1 partners: keep best sp1 copy per sp2 gene
  one_to_n <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "1:N") |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::slice_max(order_by = .data$CCS, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # N:M — two-pass: best sp2 per sp1, then best sp1 per sp2
  n_to_m <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "N:M")

  if (nrow(n_to_m) > 0) {
    n_to_m <- n_to_m |>
      dplyr::group_by(.data$gene_sp1) |>
      dplyr::slice_max(order_by = .data$CCS, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::group_by(.data$gene_sp2) |>
      dplyr::slice_max(order_by = .data$CCS, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  dplyr::bind_rows(one_to_one, one_to_n, n_to_one, n_to_m)
}
