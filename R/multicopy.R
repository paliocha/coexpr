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


#' Collapse multi-copy orthologs to expand CCS reference set
#'
#' For each multi-copy ortholog group, runs a preliminary CCS-like scoring pass
#' using the existing 1:1 reference set, then selects the best representative
#' from each group. This expands the reference set while maintaining a proper
#' bijective (1:1) mapping.
#'
#' This is useful when very few genes have clean 1:1 ortholog relationships,
#' as occurs with polyploids (WGD-derived homeologs) or de novo transcriptome
#' assemblies (fragmented isoforms).
#'
#' @param orthologs Data frame with columns `gene_sp1`, `gene_sp2`, and
#'   optionally `type` (from `detect_ortholog_types()`).
#' @param similarity_sp1 Similarity matrix for species 1. Can be a regular
#'   matrix or a TriSimilarity object.
#' @param similarity_sp2 Similarity matrix for species 2. Can be a regular
#'   matrix or a TriSimilarity object.
#' @param multicopy_sp Character. Which species has the multi-copy issue:
#'   - `"sp2"` (default): Collapses N:1 groups (sp2 has duplicates, e.g.
#'     diploid-to-polyploid comparison).
#'   - `"sp1"`: Collapses 1:N groups (sp1 has duplicates).
#'   - `"both"`: Collapses both sides sequentially (e.g. two de novo
#'     transcriptomes). Handles N:M groups via greedy two-pass selection.
#' @param max_copy_number Integer or NULL. Only collapse groups with at most
#'   this many copies (default 2L, matching WGD). Set to NULL to collapse all
#'   group sizes.
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{gene_sp1, gene_sp2}{The selected ortholog pair.}
#'     \item{type}{Set to `"1:1"` for all rows, so collapsed pairs integrate
#'       seamlessly with [calculate_ccs()] `use_only_1to1 = TRUE`.}
#'     \item{original_type}{Pre-collapse type: `NA` for natural 1:1 orthologs;
#'       `"N:1"`, `"1:N"`, or `"N:M"` for collapsed pairs.}
#'     \item{homeolog_score}{Preliminary CCS score used for selection (`NA` for
#'       natural 1:1 pairs).}
#'     \item{n_candidates}{Number of candidates in the group (1 for natural
#'       1:1 pairs).}
#'   }
#'
#' @details
#' The scoring function for each candidate pair correlates their co-expression
#' vectors (restricted to the 1:1 reference genes) across the two species — the
#' same logic as CCS, but using only the initial 1:1 set as reference.
#'
#' **Type convention** (from `detect_ortholog_types()`):
#' - `"N:1"`: `gene_sp1` appears multiple times (sp1 gene has N sp2 partners).
#'   Multi-copy is on sp2 side. Group by `gene_sp1` to pick best `gene_sp2`.
#' - `"1:N"`: `gene_sp2` appears multiple times (sp2 gene has N sp1 partners).
#'   Multi-copy is on sp1 side. Group by `gene_sp2` to pick best `gene_sp1`.
#' - `"N:M"`: Both sides have duplicates.
#'
#' @examples
#' \dontrun{
#' # Collapse N:1 homeologs (diploid vs polyploid)
#' expanded_ref <- collapse_orthologs(
#'   orthologs, sim_sp1, sim_sp2,
#'   multicopy_sp = "sp2", max_copy_number = 2L
#' )
#'
#' # Use expanded reference for CCS
#' ccs <- calculate_ccs(sim_sp1, sim_sp2, expanded_ref)
#' }
#'
#' @export
collapse_orthologs <- function(orthologs,
                               similarity_sp1,
                               similarity_sp2,
                               multicopy_sp = c("sp2", "sp1", "both"),
                               max_copy_number = 2L) {

  multicopy_sp <- match.arg(multicopy_sp)

  # Validate inputs
  required_cols <- c("gene_sp1", "gene_sp2")
  if (!all(required_cols %in% colnames(orthologs))) {
    stop(sprintf("orthologs must have columns: %s",
                 paste(required_cols, collapse = ", ")))
  }

  is_valid_sim <- function(x) is.matrix(x) || is(x, "TriSimilarity")
  if (!is_valid_sim(similarity_sp1) || !is_valid_sim(similarity_sp2)) {
    stop("similarity_sp1 and similarity_sp2 must be matrices or TriSimilarity objects")
  }

  # Detect ortholog types if not provided
  if (!"type" %in% colnames(orthologs)) {
    orthologs <- detect_ortholog_types(orthologs)
  }

  # Extract 1:1 orthologs as initial reference set
  ref_pairs <- orthologs |>
    dplyr::filter(.data$type == "1:1")

  if (nrow(ref_pairs) == 0) {
    stop("No 1:1 orthologs found. Cannot build reference set for scoring.")
  }

  # Get gene names from similarity matrices
  get_sim_genes <- function(sim) {
    if (is(sim, "TriSimilarity")) sim@genes else rownames(sim)
  }

  genes_sp1 <- get_sim_genes(similarity_sp1)
  genes_sp2 <- get_sim_genes(similarity_sp2)

  # Filter reference pairs to those present in both similarity matrices
  ref_pairs <- ref_pairs |>
    dplyr::filter(.data$gene_sp1 %in% genes_sp1,
                  .data$gene_sp2 %in% genes_sp2)

  if (nrow(ref_pairs) == 0) {
    stop("No 1:1 reference orthologs found in similarity matrices.")
  }

  ref_genes_sp1 <- ref_pairs$gene_sp1
  ref_genes_sp2 <- ref_pairs$gene_sp2

  # Helper to extract a column from either matrix type or TriSimilarity
  get_sim_column <- function(sim, gene) {
    if (is(sim, "TriSimilarity")) {
      extractColumn(sim, gene)
    } else {
      sim[, gene]
    }
  }

  # Score a candidate pair using preliminary CCS (correlation of co-expression

  # vectors restricted to the 1:1 reference genes)
  score_pair <- function(gene_sp1, gene_sp2) {
    vec1 <- get_sim_column(similarity_sp1, gene_sp1)[ref_genes_sp1]
    vec2 <- get_sim_column(similarity_sp2, gene_sp2)[ref_genes_sp2]
    stats::cor(vec1, vec2, use = "pairwise.complete.obs")
  }

  # Build 1:1 result rows
  result_1to1 <- ref_pairs |>
    dplyr::select("gene_sp1", "gene_sp2") |>
    dplyr::mutate(
      type = "1:1",
      original_type = NA_character_,
      homeolog_score = NA_real_,
      n_candidates = 1L
    )

  # Helper to collapse one side of multi-copy groups
  # group_col: column to group by (the singleton side)
  # candidate_col: column with duplicates (candidates to pick from)
  # target_type: ortholog type to filter for (e.g. "N:1", "1:N")
  collapse_groups <- function(orth, group_col, candidate_col, target_type) {
    multi <- orth |>
      dplyr::filter(.data$type == target_type)

    if (nrow(multi) == 0) {
      return(data.frame(
        gene_sp1 = character(0), gene_sp2 = character(0),
        type = character(0), original_type = character(0),
        homeolog_score = numeric(0), n_candidates = integer(0),
        stringsAsFactors = FALSE
      ))
    }

    # Count candidates per group for max_copy_number filtering
    group_counts <- multi |>
      dplyr::count(!!rlang::sym(group_col), name = "n_cand")

    multi <- multi |>
      dplyr::left_join(group_counts, by = group_col)

    # Apply max_copy_number filter
    if (!is.null(max_copy_number)) {
      skipped <- multi |>
        dplyr::filter(.data$n_cand > max_copy_number)
      if (nrow(skipped) > 0) {
        n_groups_skipped <- length(unique(skipped[[group_col]]))
        message(sprintf(
          "Skipping %d %s groups with >%d copies",
          n_groups_skipped, target_type, max_copy_number
        ))
      }
      multi <- multi |>
        dplyr::filter(.data$n_cand <= max_copy_number)
    }

    if (nrow(multi) == 0) {
      return(data.frame(
        gene_sp1 = character(0), gene_sp2 = character(0),
        type = character(0), original_type = character(0),
        homeolog_score = numeric(0), n_candidates = integer(0),
        stringsAsFactors = FALSE
      ))
    }

    # Score each candidate
    scores <- vapply(seq_len(nrow(multi)), function(i) {
      g1 <- multi$gene_sp1[i]
      g2 <- multi$gene_sp2[i]

      # Check both genes exist in similarity matrices
      if (!(g1 %in% genes_sp1)) {
        warning(sprintf("Gene '%s' not found in similarity_sp1, skipping", g1))
        return(NA_real_)
      }
      if (!(g2 %in% genes_sp2)) {
        warning(sprintf("Gene '%s' not found in similarity_sp2, skipping", g2))
        return(NA_real_)
      }

      score_pair(g1, g2)
    }, numeric(1))

    multi$homeolog_score <- scores

    # Select best candidate per group
    selected <- multi |>
      dplyr::group_by(!!rlang::sym(group_col)) |>
      dplyr::filter(!is.na(.data$homeolog_score)) |>
      dplyr::slice_max(order_by = .data$homeolog_score, n = 1,
                       with_ties = FALSE) |>
      dplyr::ungroup()

    # Warn about groups where all candidates had NA scores
    groups_with_scores <- unique(selected[[group_col]])
    all_groups <- unique(multi[[group_col]])
    na_groups <- setdiff(all_groups, groups_with_scores)
    if (length(na_groups) > 0) {
      warning(sprintf(
        "%d %s group(s) skipped: all candidates had NA scores",
        length(na_groups), target_type
      ))
    }

    selected |>
      dplyr::transmute(
        gene_sp1 = .data$gene_sp1,
        gene_sp2 = .data$gene_sp2,
        type = "1:1",
        original_type = target_type,
        homeolog_score = .data$homeolog_score,
        n_candidates = .data$n_cand
      )
  }

  # Collapse depending on multicopy_sp
  collapsed_parts <- list()

  if (multicopy_sp %in% c("sp2", "both")) {
    # N:1: gene_sp1 appears multiple times, group by gene_sp1, pick best gene_sp2
    collapsed_parts$n_to_1 <- collapse_groups(
      orthologs, "gene_sp1", "gene_sp2", "N:1"
    )
  }

  if (multicopy_sp %in% c("sp1", "both")) {
    # 1:N: gene_sp2 appears multiple times, group by gene_sp2, pick best gene_sp1
    collapsed_parts$one_to_n <- collapse_groups(
      orthologs, "gene_sp2", "gene_sp1", "1:N"
    )
  }

  if (multicopy_sp == "both") {
    # N:M: greedy two-pass — best sp2 per sp1, then deduplicate sp2
    nm_pairs <- orthologs |>
      dplyr::filter(.data$type == "N:M")

    if (nrow(nm_pairs) > 0) {
      # Count candidates and apply max_copy_number
      # For N:M, count unique sp2 per sp1 as proxy for group size
      if (!is.null(max_copy_number)) {
        sp2_per_sp1 <- nm_pairs |>
          dplyr::count(.data$gene_sp1, name = "n_sp2")
        sp1_per_sp2 <- nm_pairs |>
          dplyr::count(.data$gene_sp2, name = "n_sp1")

        nm_pairs <- nm_pairs |>
          dplyr::left_join(sp2_per_sp1, by = "gene_sp1") |>
          dplyr::left_join(sp1_per_sp2, by = "gene_sp2") |>
          dplyr::filter(.data$n_sp2 <= max_copy_number,
                        .data$n_sp1 <= max_copy_number) |>
          dplyr::select(-"n_sp2", -"n_sp1")
      }

      if (nrow(nm_pairs) > 0) {
        # Score all N:M candidates
        scores <- vapply(seq_len(nrow(nm_pairs)), function(i) {
          g1 <- nm_pairs$gene_sp1[i]
          g2 <- nm_pairs$gene_sp2[i]
          if (!(g1 %in% genes_sp1) || !(g2 %in% genes_sp2)) {
            return(NA_real_)
          }
          score_pair(g1, g2)
        }, numeric(1))

        nm_pairs$homeolog_score <- scores

        # Pass 1: best sp2 per sp1
        nm_selected <- nm_pairs |>
          dplyr::filter(!is.na(.data$homeolog_score)) |>
          dplyr::group_by(.data$gene_sp1) |>
          dplyr::slice_max(order_by = .data$homeolog_score, n = 1,
                           with_ties = FALSE) |>
          dplyr::ungroup()

        # Pass 2: deduplicate sp2 (best sp1 per sp2)
        nm_selected <- nm_selected |>
          dplyr::group_by(.data$gene_sp2) |>
          dplyr::slice_max(order_by = .data$homeolog_score, n = 1,
                           with_ties = FALSE) |>
          dplyr::ungroup()

        collapsed_parts$n_to_m <- nm_selected |>
          dplyr::transmute(
            gene_sp1 = .data$gene_sp1,
            gene_sp2 = .data$gene_sp2,
            type = "1:1",
            original_type = "N:M",
            homeolog_score = .data$homeolog_score,
            n_candidates = nrow(nm_pairs)
          )
      }
    }
  }

  # Combine all parts
  result <- dplyr::bind_rows(c(list(result_1to1), collapsed_parts))

  n_collapsed <- sum(!is.na(result$original_type))
  if (n_collapsed == 0) {
    message("No multi-copy groups found to collapse. Returning 1:1 pairs only.")
  } else {
    message(sprintf(
      "Collapsed %d multi-copy groups: %d total pairs (from %d original 1:1)",
      n_collapsed, nrow(result), nrow(result_1to1)
    ))
  }

  result
}
