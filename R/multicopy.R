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


#' Select best hit for multi-copy orthologs using CCS-based scoring
#'
#' Scores each candidate pair by correlating co-expression vectors against the
#' 1:1 reference set (preliminary CCS). Selects the highest-scoring candidate
#' per multi-copy group.
#'
#' @keywords internal
#' @noRd
select_best_hits <- function(orthologs, similarity_sp1, similarity_sp2) {

  if (is.null(similarity_sp1) || is.null(similarity_sp2)) {
    stop("Strategy 'best_hit' requires similarity_sp1 and similarity_sp2")
  }

  # Build 1:1 reference set for CCS-based scoring
  ref_pairs <- orthologs |>
    dplyr::filter(.data$type == "1:1")

  get_sim_genes <- function(sim) {
    if (is(sim, "TriSimilarity")) sim@genes else rownames(sim)
  }
  genes_sp1 <- get_sim_genes(similarity_sp1)
  genes_sp2 <- get_sim_genes(similarity_sp2)

  ref_pairs <- ref_pairs |>
    dplyr::filter(.data$gene_sp1 %in% genes_sp1,
                  .data$gene_sp2 %in% genes_sp2)

  if (nrow(ref_pairs) < 3) {
    stop("Strategy 'best_hit' requires at least 3 reference 1:1 orthologs ",
         "in both similarity matrices for CCS-based scoring")
  }

  # Score all non-1:1 candidates via preliminary CCS
  score_fn <- make_pair_scorer(similarity_sp1, similarity_sp2,
                               ref_pairs$gene_sp1, ref_pairs$gene_sp2)

  multi <- orthologs |>
    dplyr::filter(.data$type != "1:1")

  if (nrow(multi) == 0) {
    return(ref_pairs)
  }

  multi$score <- vapply(seq_len(nrow(multi)), function(i) {
    score_fn(multi$gene_sp1[i], multi$gene_sp2[i])
  }, numeric(1))

  # N:1 — sp1 gene has N sp2 partners: pick best sp2 per sp1 gene
  orthologs_Nto1 <- multi |>
    dplyr::filter(.data$type == "N:1") |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::slice_max(order_by = .data$score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-"score")

  # 1:N — sp2 gene has N sp1 partners: pick best sp1 per sp2 gene
  orthologs_1toN <- multi |>
    dplyr::filter(.data$type == "1:N") |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::slice_max(order_by = .data$score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(-"score")

  # N:M — optimal 1:1 assignment (Hungarian algorithm or greedy fallback)
  orthologs_NtoM <- multi |>
    dplyr::filter(.data$type == "N:M")

  if (nrow(orthologs_NtoM) > 0) {
    orthologs_NtoM <- resolve_nm_optimal(orthologs_NtoM, score_col = "score") |>
      dplyr::select(-"score")
  } else {
    orthologs_NtoM <- orthologs_NtoM |>
      dplyr::select(-"score")
  }

  dplyr::bind_rows(ref_pairs, orthologs_Nto1, orthologs_1toN, orthologs_NtoM)
}


#' Create a pair scoring function using preliminary CCS
#'
#' Returns a closure that scores an (sp1, sp2) gene pair by correlating
#' their co-expression vectors restricted to the reference ortholog genes.
#'
#' @param similarity_sp1,similarity_sp2 Similarity matrices or TriSimilarity
#' @param ref_genes_sp1,ref_genes_sp2 Character vectors of reference gene names
#' @return Function(gene_sp1, gene_sp2) -> numeric CCS score
#' @keywords internal
#' @noRd
make_pair_scorer <- function(similarity_sp1, similarity_sp2,
                             ref_genes_sp1, ref_genes_sp2) {

  get_sim_genes <- function(sim) {
    if (is(sim, "TriSimilarity")) sim@genes else rownames(sim)
  }
  genes_sp1 <- get_sim_genes(similarity_sp1)
  genes_sp2 <- get_sim_genes(similarity_sp2)

  get_col <- function(sim, gene) {
    if (is(sim, "TriSimilarity")) extractColumn(sim, gene) else sim[, gene]
  }

  function(gene_sp1, gene_sp2) {
    if (!(gene_sp1 %in% genes_sp1) || !(gene_sp2 %in% genes_sp2)) {
      return(NA_real_)
    }
    vec1 <- get_col(similarity_sp1, gene_sp1)[ref_genes_sp1]
    vec2 <- get_col(similarity_sp2, gene_sp2)[ref_genes_sp2]
    stats::cor(vec1, vec2, use = "pairwise.complete.obs")
  }
}


#' Resolve N:M ortholog groups to optimal 1:1 assignments
#'
#' Given a data frame of N:M pairs with scores, finds the 1:1 assignment
#' that maximizes total score. Uses the Hungarian algorithm (via
#' `clue::solve_LSAP`) when available, otherwise falls back to greedy
#' two-pass selection.
#'
#' @param pairs Data frame with columns `gene_sp1`, `gene_sp2`, and a score
#'   column named by `score_col`.
#' @param score_col Character. Name of the score column (default `"score"`).
#' @return Data frame: one row per assigned pair (subset of input rows).
#' @keywords internal
#' @noRd
resolve_nm_optimal <- function(pairs, score_col = "score") {
  if (nrow(pairs) == 0) return(pairs)

  # Remove rows with NA scores
  pairs <- pairs[!is.na(pairs[[score_col]]), ]
  if (nrow(pairs) == 0) return(pairs)

  sp1_genes <- unique(pairs$gene_sp1)
  sp2_genes <- unique(pairs$gene_sp2)

  # Trivial cases: 1 gene on either side

  if (length(sp1_genes) == 1 || length(sp2_genes) == 1) {
    # Just pick best per group
    if (length(sp1_genes) == 1) {
      return(pairs[which.max(pairs[[score_col]]), , drop = FALSE])
    } else {
      return(pairs[which.max(pairs[[score_col]]), , drop = FALSE])
    }
  }

  # Try Hungarian algorithm
  if (requireNamespace("clue", quietly = TRUE)) {
    # Build score matrix (sp1 rows x sp2 cols)
    score_mat <- matrix(-Inf, nrow = length(sp1_genes), ncol = length(sp2_genes))
    rownames(score_mat) <- sp1_genes
    colnames(score_mat) <- sp2_genes

    score_mat[cbind(pairs$gene_sp1, pairs$gene_sp2)] <- pairs[[score_col]]

    # solve_LSAP minimizes cost with non-negative entries, so transform:
    # cost = max_score - score (large cost for low scores, zero for best)
    # Replace -Inf with a sentinel before computing max
    nr <- nrow(score_mat)
    nc <- ncol(score_mat)
    finite_scores <- score_mat[is.finite(score_mat)]
    if (length(finite_scores) == 0) return(pairs[0, , drop = FALSE])
    max_score <- max(finite_scores)

    cost_mat <- score_mat
    cost_mat[!is.finite(cost_mat)] <- max_score - 1e6  # will become +1e6 cost
    cost_mat <- max_score - cost_mat  # now 0 = best, large = worst

    # Pad to square if needed (LSAP requires square matrix)
    if (nr != nc) {
      n <- max(nr, nc)
      padded <- matrix(max(cost_mat) + 1, n, n)  # dummy rows/cols = high cost
      padded[seq_len(nr), seq_len(nc)] <- cost_mat
      assignment <- clue::solve_LSAP(padded)
      assignment <- as.integer(assignment)
    } else {
      assignment <- clue::solve_LSAP(cost_mat)
      assignment <- as.integer(assignment)
    }

    # Extract valid assignments (within original dimensions and with real scores)
    valid_i <- seq_len(min(nr, length(assignment)))
    valid_j <- assignment[valid_i]
    keep <- valid_j <= nc & score_mat[cbind(valid_i, valid_j)] > -Inf
    selected <- data.frame(
      gene_sp1 = sp1_genes[valid_i[keep]],
      gene_sp2 = sp2_genes[valid_j[keep]],
      stringsAsFactors = FALSE
    )

    # Rejoin with original pairs to preserve all columns
    result <- dplyr::semi_join(pairs, selected, by = c("gene_sp1", "gene_sp2"))
    return(result)
  }

  # Fallback: greedy two-pass
  result <- pairs |>
    dplyr::group_by(.data$gene_sp1) |>
    dplyr::slice_max(order_by = .data[[score_col]], n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::slice_max(order_by = .data[[score_col]], n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  result
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
    dplyr::summarise(
      gene_sp2 = paste(.data$gene_sp2, collapse = ";"),
      CCS = mean(.data$CCS, na.rm = TRUE),
      type = "aggregated_N:1",
      .groups = "drop"
    )

  # 1:N — sp2 gene has N sp1 partners: group by gene_sp2, average across sp1 copies
  one_to_n <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "1:N") |>
    dplyr::group_by(.data$gene_sp2) |>
    dplyr::summarise(
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
      dplyr::summarise(
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

  # N:M — optimal 1:1 assignment (Hungarian algorithm or greedy fallback)
  n_to_m <- orthologs_with_ccs |>
    dplyr::filter(.data$type == "N:M")

  if (nrow(n_to_m) > 0) {
    n_to_m <- n_to_m |>
      dplyr::mutate(score = .data$CCS) |>
      resolve_nm_optimal(score_col = "score") |>
      dplyr::select(-"score")
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
#'     \item{selection_confidence}{(only if `n_bootstrap > 0`) Proportion of
#'       bootstrap replicates selecting the same candidate. Values near 1
#'       indicate a robust choice; low values suggest the selection is sensitive
#'       to reference composition. `NA` for natural 1:1 pairs.}
#'     \item{score_ci_low, score_ci_high}{(only if `n_bootstrap > 0`)
#'       2.5th and 97.5th percentile of the homeolog score across bootstrap
#'       replicates. `NA` for natural 1:1 pairs.}
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
#' @param n_bootstrap Integer. Number of bootstrap iterations to assess
#'   selection confidence. If > 0, subsamples 80% of the reference set
#'   `n_bootstrap` times and re-scores candidates. Adds a
#'   `selection_confidence` column (proportion of bootstraps that select the
#'   same candidate) and `score_ci_low`/`score_ci_high` (2.5th/97.5th
#'   percentile of homeolog scores). Default 0 (no bootstrap).
#'
#' @export
collapse_orthologs <- function(orthologs,
                               similarity_sp1,
                               similarity_sp2,
                               multicopy_sp = c("sp2", "sp1", "both"),
                               max_copy_number = 2L,
                               n_bootstrap = 0L) {

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

  # CCS-based scoring: correlate co-expression vectors against 1:1 reference
  score_pair <- make_pair_scorer(similarity_sp1, similarity_sp2,
                                 ref_genes_sp1, ref_genes_sp2)

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

        # Count candidates per sp1 gene (for n_candidates metadata)
        nm_group_counts <- nm_pairs |>
          dplyr::count(.data$gene_sp1, name = "n_cand")
        nm_pairs <- nm_pairs |>
          dplyr::left_join(nm_group_counts, by = "gene_sp1")

        # Optimal 1:1 assignment (Hungarian algorithm or greedy fallback)
        nm_selected <- nm_pairs |>
          dplyr::mutate(score = .data$homeolog_score) |>
          resolve_nm_optimal(score_col = "score") |>
          dplyr::select(-"score")

        collapsed_parts$n_to_m <- nm_selected |>
          dplyr::transmute(
            gene_sp1 = .data$gene_sp1,
            gene_sp2 = .data$gene_sp2,
            type = "1:1",
            original_type = "N:M",
            homeolog_score = .data$homeolog_score,
            n_candidates = .data$n_cand
          )
      }
    }
  }

  # Combine all parts
  result <- dplyr::bind_rows(c(list(result_1to1), collapsed_parts))

  # Bootstrap confidence assessment (C++ accelerated)
  n_bootstrap <- as.integer(n_bootstrap)
  if (n_bootstrap > 0 && any(!is.na(result$original_type))) {
    collapsed_idx <- which(!is.na(result$original_type))

    # Pre-extract full reference submatrices once (n_ref x n_genes)
    get_col <- function(sim, gene) {
      if (is(sim, "TriSimilarity")) extractColumn(sim, gene) else sim[, gene]
    }

    boot_confidence <- rep(NA_real_, nrow(result))
    boot_score_low <- rep(NA_real_, nrow(result))
    boot_score_high <- rep(NA_real_, nrow(result))

    base_seed <- sample.int(.Machine$integer.max, 1)

    for (ci in collapsed_idx) {
      row <- result[ci, ]
      otype <- row$original_type

      if (otype == "N:1") {
        group_col <- "gene_sp1"
        cand_col <- "gene_sp2"
      } else if (otype == "1:N") {
        group_col <- "gene_sp2"
        cand_col <- "gene_sp1"
      } else {
        group_col <- "gene_sp1"
        cand_col <- "gene_sp2"
      }

      group_val <- row[[group_col]]
      selected_cand <- row[[cand_col]]

      cands <- orthologs[orthologs$type == otype &
                           orthologs[[group_col]] == group_val, ]
      if (nrow(cands) <= 1) {
        boot_confidence[ci] <- 1.0
        next
      }

      # Build reference submatrices for this group: n_ref x n_candidates
      # Each column = candidate's co-expression vector against references
      n_cand <- nrow(cands)
      mat_sp1 <- vapply(cands$gene_sp1, \(g) get_col(similarity_sp1, g)[ref_genes_sp1],
                         numeric(length(ref_genes_sp1)))
      mat_sp2 <- vapply(cands$gene_sp2, \(g) get_col(similarity_sp2, g)[ref_genes_sp2],
                         numeric(length(ref_genes_sp2)))

      selected_idx <- which(cands[[cand_col]] == selected_cand) - 1L  # 0-based

      # Call C++ bootstrap
      boot_result <- bootstrap_group_scores_cpp(
        mat_sp1, mat_sp2, selected_idx[1],
        n_bootstrap, as.integer(base_seed + ci)
      )

      boot_confidence[ci] <- boot_result$n_same / n_bootstrap
      boot_scores <- boot_result$scores[is.finite(boot_result$scores)]
      if (length(boot_scores) >= 2) {
        boot_score_low[ci] <- stats::quantile(boot_scores, 0.025, names = FALSE)
        boot_score_high[ci] <- stats::quantile(boot_scores, 0.975, names = FALSE)
      }
    }

    result$selection_confidence <- boot_confidence
    result$score_ci_low <- boot_score_low
    result$score_ci_high <- boot_score_high
  }

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


#' Iteratively expand ortholog reference set
#'
#' Starts with 1:1 orthologs as the reference set, then iteratively adds
#' the best multi-copy ortholog representatives. Each iteration re-scores
#' candidates using the expanded reference from the previous round, so
#' later additions benefit from a richer reference.
#'
#' This is useful when few genes have clean 1:1 relationships (e.g.,
#' polyploid comparisons) and a single-pass collapse is insufficient.
#'
#' @param orthologs Data frame with columns `gene_sp1`, `gene_sp2`, and
#'   optionally `type`.
#' @param similarity_sp1,similarity_sp2 Similarity matrices or TriSimilarity
#'   objects for each species.
#' @param multicopy_sp Character. Which species has multi-copy genes:
#'   `"sp2"` (default), `"sp1"`, or `"both"`.
#' @param max_copy_number Integer or NULL. Only consider groups with at most
#'   this many copies (default 2L). NULL = no limit.
#' @param max_iterations Integer. Maximum number of expansion rounds
#'   (default 5). Iteration stops early if no new pairs are added.
#' @param score_threshold Numeric. Minimum preliminary CCS (homeolog score)
#'   required to accept a candidate (default 0.1). Higher values are more
#'   conservative.
#'
#' @return Data frame with the same columns as [collapse_orthologs()], plus
#'   an `iteration` column indicating when each pair was added (0 = original
#'   1:1, 1 = first expansion, etc.).
#'
#' @details
#' Algorithm:
#' 1. Extract 1:1 orthologs as initial reference (iteration 0)
#' 2. Score all multi-copy candidates against current reference
#' 3. For each group, select the best candidate above `score_threshold`
#' 4. Add selected pairs to the reference set
#' 5. Repeat from step 2 with remaining unassigned groups
#' 6. Stop when no new pairs are added or `max_iterations` is reached
#'
#' @examples
#' \dontrun{
#' expanded <- expand_reference_iteratively(
#'   orthologs, sim_sp1, sim_sp2,
#'   multicopy_sp = "sp2", max_copy_number = 2L,
#'   max_iterations = 5, score_threshold = 0.1
#' )
#'
#' ccs <- calculate_ccs(sim_sp1, sim_sp2, expanded)
#' }
#'
#' @export
expand_reference_iteratively <- function(orthologs,
                                         similarity_sp1,
                                         similarity_sp2,
                                         multicopy_sp = c("sp2", "sp1", "both"),
                                         max_copy_number = 2L,
                                         max_iterations = 5L,
                                         score_threshold = 0.1) {

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

  if (!"type" %in% colnames(orthologs)) {
    orthologs <- detect_ortholog_types(orthologs)
  }

  get_sim_genes <- function(sim) {
    if (is(sim, "TriSimilarity")) sim@genes else rownames(sim)
  }
  genes_sp1 <- get_sim_genes(similarity_sp1)
  genes_sp2 <- get_sim_genes(similarity_sp2)

  # Initialize reference with 1:1 orthologs
  ref_pairs <- orthologs |>
    dplyr::filter(.data$type == "1:1",
                  .data$gene_sp1 %in% genes_sp1,
                  .data$gene_sp2 %in% genes_sp2)

  if (nrow(ref_pairs) == 0) {
    stop("No 1:1 orthologs found in similarity matrices.")
  }

  result <- ref_pairs |>
    dplyr::select("gene_sp1", "gene_sp2") |>
    dplyr::mutate(
      type = "1:1",
      original_type = NA_character_,
      homeolog_score = NA_real_,
      n_candidates = 1L,
      iteration = 0L
    )

  # Determine which ortholog types to expand
  target_types <- switch(multicopy_sp,
    "sp2"  = "N:1",
    "sp1"  = "1:N",
    "both" = c("N:1", "1:N")
  )

  # Build pool of multi-copy candidates
  candidates <- orthologs |>
    dplyr::filter(.data$type %in% target_types,
                  .data$gene_sp1 %in% genes_sp1,
                  .data$gene_sp2 %in% genes_sp2)

  # Apply max_copy_number filter
  if (!is.null(max_copy_number) && nrow(candidates) > 0) {
    for (tt in target_types) {
      group_col <- if (tt == "N:1") "gene_sp1" else "gene_sp2"
      group_counts <- candidates |>
        dplyr::filter(.data$type == tt) |>
        dplyr::count(!!rlang::sym(group_col), name = "n_cand")

      too_large <- group_counts |>
        dplyr::filter(.data$n_cand > max_copy_number) |>
        dplyr::pull(!!rlang::sym(group_col))

      if (length(too_large) > 0) {
        candidates <- candidates |>
          dplyr::filter(!(!!rlang::sym(group_col) %in% too_large &
                            .data$type == tt))
      }
    }
  }

  if (nrow(candidates) == 0) {
    message("No eligible multi-copy groups. Returning 1:1 pairs only.")
    return(result)
  }

  # Compute n_candidates per group (vectorized via dplyr)
  # N:1 groups: count per gene_sp1; 1:N groups: count per gene_sp2
  n1_counts <- candidates |>
    dplyr::filter(.data$type == "N:1") |>
    dplyr::add_count(.data$gene_sp1, name = "n_cand")
  on_counts <- candidates |>
    dplyr::filter(.data$type != "N:1") |>
    dplyr::add_count(.data$gene_sp2, name = "n_cand")
  candidates <- dplyr::bind_rows(n1_counts, on_counts)

  # Track which groups have been assigned
  assigned_groups <- character(0)

  for (iter in seq_len(max_iterations)) {
    ref_genes_sp1 <- result$gene_sp1
    ref_genes_sp2 <- result$gene_sp2

    # Score remaining candidates against current reference
    score_fn <- make_pair_scorer(similarity_sp1, similarity_sp2,
                                 ref_genes_sp1, ref_genes_sp2)

    remaining <- candidates |>
      dplyr::filter(!paste(.data$gene_sp1, .data$gene_sp2) %in% assigned_groups)

    if (nrow(remaining) == 0) break

    # Remove candidates whose group is already assigned
    for (tt in target_types) {
      group_col <- if (tt == "N:1") "gene_sp1" else "gene_sp2"
      assigned_in_result <- result[[group_col]]
      remaining <- remaining |>
        dplyr::filter(!(!!rlang::sym(group_col) %in% assigned_in_result &
                          .data$type == tt))
    }

    if (nrow(remaining) == 0) break

    # Score each candidate
    remaining$homeolog_score <- vapply(seq_len(nrow(remaining)), function(i) {
      score_fn(remaining$gene_sp1[i], remaining$gene_sp2[i])
    }, numeric(1))

    # Select best per group, above threshold
    new_pairs <- data.frame()
    for (tt in target_types) {
      group_col <- if (tt == "N:1") "gene_sp1" else "gene_sp2"
      type_remaining <- remaining |>
        dplyr::filter(.data$type == tt, !is.na(.data$homeolog_score),
                      .data$homeolog_score >= score_threshold)

      if (nrow(type_remaining) > 0) {
        best <- type_remaining |>
          dplyr::group_by(!!rlang::sym(group_col)) |>
          dplyr::slice_max(order_by = .data$homeolog_score, n = 1,
                           with_ties = FALSE) |>
          dplyr::ungroup()
        new_pairs <- dplyr::bind_rows(new_pairs, best)
      }
    }

    if (nrow(new_pairs) == 0) {
      message(sprintf("Iteration %d: no candidates above threshold. Stopping.", iter))
      break
    }

    # Add to result
    new_rows <- new_pairs |>
      dplyr::transmute(
        gene_sp1 = .data$gene_sp1,
        gene_sp2 = .data$gene_sp2,
        type = "1:1",
        original_type = .data$type,
        homeolog_score = .data$homeolog_score,
        n_candidates = as.integer(.data$n_cand),
        iteration = iter
      )

    result <- dplyr::bind_rows(result, new_rows)

    message(sprintf("Iteration %d: added %d pairs (total: %d)",
                    iter, nrow(new_rows), nrow(result)))

    # Mark these groups as assigned
    for (tt in target_types) {
      group_col <- if (tt == "N:1") "gene_sp1" else "gene_sp2"
      assigned_groups <- c(assigned_groups,
                           new_pairs[[group_col]][new_pairs$type == tt])
    }
  }

  n_expanded <- sum(!is.na(result$original_type))
  n_iters <- max(result$iteration)
  message(sprintf(
    "Reference expansion complete: %d pairs (%d original 1:1 + %d expanded) over %d iteration(s)",
    nrow(result), sum(result$iteration == 0), n_expanded, n_iters
  ))

  result
}


#' Analyze expression divergence among paralogs
#'
#' Characterizes expression divergence patterns in multi-copy ortholog groups
#' using CCS and ORS values from an `all_pairs` analysis. For each group of
#' paralogs (1:N, N:1, N:M), identifies the primary (most conserved) copy,
#' quantifies divergence between copies, and classifies divergence levels.
#'
#' @param ors_results Data frame from `calculate_ors()` with at minimum columns
#'   `gene_sp1`, `gene_sp2`, `CCS`, `ORS`, and `logORS`. Should typically be
#'   generated using `all_pairs` strategy and `directional = TRUE` in
#'   `calculate_ors()`. Must include a `type` column with ortholog types.
#' @param ccs_threshold Numeric. CCS threshold for classifying divergence.
#'   Pairs with CCS above this are "conserved". Default 0.3.
#' @param logors_threshold Numeric. logORS threshold for "highly conserved"
#'   classification. Default 1 (top 10%).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{per_group}{Data frame with one row per multi-copy group:
#'       \itemize{
#'         \item `group_gene`: The gene with multiple orthologs
#'         \item `type`: Ortholog type (1:N, N:1, or N:M)
#'         \item `copy_number`: Number of copies in the group
#'         \item `primary_gene`: Best-conserved copy (highest CCS)
#'         \item `primary_ccs`: CCS of the primary copy
#'         \item `secondary_ccs`: CCS of the second-best copy (NA if only 1)
#'         \item `delta_ccs`: Difference between primary and secondary CCS
#'         \item `mean_ccs`: Mean CCS across all copies
#'         \item `classification`: One of "conserved", "partially_diverged",
#'           "fully_diverged"
#'       }
#'     }
#'     \item{by_copy_number}{Data frame summarizing divergence by copy number
#'       (1:2, 1:3, 1:4, ...):
#'       \itemize{
#'         \item `copy_number`: Number of copies
#'         \item `n_groups`: Number of ortholog groups
#'         \item `median_primary_ccs`: Median CCS of primary copy
#'         \item `median_delta_ccs`: Median CCS gap between primary and secondary
#'         \item `median_logORS`: Median logORS of primary copy
#'         \item `pct_conserved`: Percent of groups classified as conserved
#'       }
#'     }
#'   }
#'
#' @details
#' The paper's key finding is that 1:2, 1:3, 1:4 orthologs show progressively
#' lower ORS (expression divergence), except for recent WGD duplicates (~13 Mya
#' in *Glycine max*) which show no divergence.
#'
#' **Classification rules**:
#' \itemize{
#'   \item **conserved**: Primary copy CCS >= `ccs_threshold` AND
#'     delta_CCS < 0.1 (copies are similarly conserved)
#'   \item **partially_diverged**: Primary copy CCS >= `ccs_threshold` AND
#'     delta_CCS >= 0.1 (one copy retains function, others diverge)
#'   \item **fully_diverged**: Primary copy CCS < `ccs_threshold`
#'     (all copies have diverged)
#' }
#'
#' @examples
#' \dontrun{
#' # Run all_pairs CCS and directional ORS
#' orthologs_all <- handle_multicopy_orthologs(orthologs, strategy = "all_pairs")
#' ccs <- calculate_ccs(sim_sp1, sim_sp2, orthologs_all)
#' ors <- calculate_ors(ccs, directional = TRUE)
#'
#' # Analyze paralog divergence
#' div <- analyze_paralog_divergence(ors)
#' div$per_group     # Per-group details
#' div$by_copy_number  # Summary by copy number
#' }
#'
#' @export
analyze_paralog_divergence <- function(ors_results,
                                       ccs_threshold = 0.3,
                                       logors_threshold = 1) {

  # Validate input
  required_cols <- c("gene_sp1", "gene_sp2", "CCS", "ORS", "logORS", "type")
  missing <- setdiff(required_cols, colnames(ors_results))
  if (length(missing) > 0) {
    stop(sprintf("ors_results missing required columns: %s. ",
                 paste(missing, collapse = ", ")),
         "Run calculate_ors(calculate_ccs(...), return_log = TRUE) with type column.")
  }

  # Identify multi-copy groups
  multicopy <- ors_results |>
    dplyr::filter(.data$type %in% c("1:N", "N:1", "N:M"))

  if (nrow(multicopy) == 0) {
    stop("No multi-copy orthologs found in ors_results. ",
         "Use strategy = 'all_pairs' in handle_multicopy_orthologs().")
  }

  # Type convention (from detect_ortholog_types):
  # N:1: n_sp1 > 1 = multiple sp1 genes map to one sp2 gene -> group by gene_sp2
  # 1:N: n_sp2 > 1 = one sp1 gene maps to multiple sp2 genes -> group by gene_sp1

  # Helper to process one type using split-apply to avoid .data pronoun issues
  process_type <- function(data, type, group_col, copy_col) {
    if (nrow(data) == 0) return(NULL)

    groups <- split(data, data[[group_col]])
    rows <- lapply(names(groups), function(gname) {
      g <- groups[[gname]]
      best_idx <- which.max(g$CCS)
      sorted_ccs <- sort(g$CCS, decreasing = TRUE, na.last = TRUE)
      data.frame(
        group_gene = gname,
        type = type,
        copy_number = nrow(g),
        primary_gene = g[[copy_col]][best_idx],
        primary_ccs = g$CCS[best_idx],
        primary_logORS = g$logORS[best_idx],
        secondary_ccs = if (nrow(g) > 1) sorted_ccs[2] else NA_real_,
        mean_ccs = mean(g$CCS, na.rm = TRUE),
        min_ccs = min(g$CCS, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
    result_df <- do.call(rbind, rows)
    result_df$delta_ccs <- result_df$primary_ccs -
      ifelse(is.na(result_df$secondary_ccs), result_df$primary_ccs,
             result_df$secondary_ccs)
    result_df$classification <- dplyr::case_when(
      result_df$primary_ccs < ccs_threshold ~ "fully_diverged",
      result_df$delta_ccs >= 0.1 ~ "partially_diverged",
      TRUE ~ "conserved"
    )
    result_df
  }

  per_group_list <- list()

  # N:1: multiple sp1 genes map to one sp2 gene -> group by gene_sp2
  n1_data <- multicopy |> dplyr::filter(.data$type == "N:1")
  per_group_list$n1 <- process_type(n1_data, "N:1", "gene_sp2", "gene_sp1")

  # 1:N: one sp1 gene maps to multiple sp2 genes -> group by gene_sp1
  onetomany_data <- multicopy |> dplyr::filter(.data$type == "1:N")
  per_group_list$onetomany <- process_type(onetomany_data, "1:N", "gene_sp1",
                                           "gene_sp2")

  # N:M: group by gene_sp1 to see how each sp1 gene's copies diverge
  nm_data <- multicopy |> dplyr::filter(.data$type == "N:M")
  if (nrow(nm_data) > 0) {
    per_group_list$nm <- process_type(nm_data, "N:M", "gene_sp1", "gene_sp2")
  }

  per_group <- dplyr::bind_rows(per_group_list) |>
    dplyr::select("group_gene", "type", "copy_number", "primary_gene",
                  "primary_ccs", "secondary_ccs", "delta_ccs", "mean_ccs",
                  "min_ccs", "primary_logORS", "classification")

  # Summarize by copy number
  by_copy_number <- per_group |>
    dplyr::group_by(.data$copy_number) |>
    dplyr::summarise(
      n_groups = dplyr::n(),
      median_primary_ccs = stats::median(.data$primary_ccs, na.rm = TRUE),
      median_delta_ccs = stats::median(.data$delta_ccs, na.rm = TRUE),
      median_logORS = stats::median(.data$primary_logORS, na.rm = TRUE),
      pct_conserved = sum(.data$classification == "conserved") /
        dplyr::n() * 100,
      pct_partially_diverged = sum(.data$classification == "partially_diverged") /
        dplyr::n() * 100,
      pct_fully_diverged = sum(.data$classification == "fully_diverged") /
        dplyr::n() * 100,
      .groups = "drop"
    )

  message(sprintf("Analyzed %d multi-copy groups (%d total paralog pairs)",
                  nrow(per_group), nrow(multicopy)))
  message(sprintf("  Copy numbers: %s",
                  paste(sort(unique(per_group$copy_number)), collapse = ", ")))
  message(sprintf("  Classification: %d conserved, %d partially diverged, %d fully diverged",
                  sum(per_group$classification == "conserved"),
                  sum(per_group$classification == "partially_diverged"),
                  sum(per_group$classification == "fully_diverged")))

  list(
    per_group = per_group,
    by_copy_number = by_copy_number
  )
}


#' Diagnose reference set quality
#'
#' Assesses whether the 1:1 ortholog reference set is adequate for CCS
#' calculation. Reports size, effective size (genes present in both similarity
#' matrices), and optionally CCS distribution statistics to detect bias.
#'
#' @param orthologs Data frame with columns `gene_sp1` and `gene_sp2`.
#'   Optionally `type` column.
#' @param similarity_sp1 Similarity matrix or TriSimilarity for species 1.
#' @param similarity_sp2 Similarity matrix or TriSimilarity for species 2.
#' @param ccs_results Optional. Data frame from `calculate_ccs()` to assess
#'   CCS distribution within the reference set.
#'
#' @return A list with:
#'   \describe{
#'     \item{n_total}{Total number of ortholog pairs}
#'     \item{n_1to1}{Number of 1:1 ortholog pairs}
#'     \item{pct_1to1}{Percent of orthologs that are 1:1}
#'     \item{n_effective}{Number of 1:1 pairs present in both similarity matrices}
#'     \item{pct_effective}{Percent of 1:1 pairs that are effective}
#'     \item{missing_sp1}{Gene IDs in reference but missing from similarity_sp1}
#'     \item{missing_sp2}{Gene IDs in reference but missing from similarity_sp2}
#'     \item{ccs_stats}{If `ccs_results` provided: list with median, mean, sd,
#'       skewness, and IQR of reference CCS values}
#'     \item{warnings}{Character vector of diagnostic warnings}
#'   }
#'
#' @details
#' A good reference set should:
#' \itemize{
#'   \item Have at least 50 effective 1:1 orthologs (>100 preferred)
#'   \item Have nearly all 1:1 orthologs present in both similarity matrices
#'   \item Show a roughly symmetric CCS distribution (not heavily skewed)
#' }
#'
#' Warning thresholds:
#' \itemize{
#'   \item `< 10 effective references`: Error-level (CCS unreliable)
#'   \item `< 20 effective references`: Strong warning
#'   \item `< 50 effective references`: Warning
#'   \item `> 20% missing from similarity matrices`: Warning
#' }
#'
#' @examples
#' \dontrun{
#' diag <- diagnose_reference(orthologs, sim_sp1, sim_sp2)
#' diag$warnings  # Check for issues
#'
#' # With CCS results for distribution analysis
#' ccs <- calculate_ccs(sim_sp1, sim_sp2, orthologs)
#' diag <- diagnose_reference(orthologs, sim_sp1, sim_sp2, ccs_results = ccs)
#' diag$ccs_stats  # Distribution statistics
#' }
#'
#' @export
diagnose_reference <- function(orthologs, similarity_sp1, similarity_sp2,
                               ccs_results = NULL) {

  if (!is.data.frame(orthologs) ||
      !all(c("gene_sp1", "gene_sp2") %in% colnames(orthologs))) {
    stop("orthologs must be a data frame with columns gene_sp1 and gene_sp2")
  }

  # Add types if missing
  if (!"type" %in% colnames(orthologs)) {
    orthologs <- detect_ortholog_types(orthologs)
  }

  n_total <- nrow(orthologs)
  ref <- orthologs[orthologs$type == "1:1", ]
  n_1to1 <- nrow(ref)

  # Get gene names from similarity matrices
  get_sim_genes <- function(sim) {
    if (methods::is(sim, "TriSimilarity")) sim@genes else rownames(sim)
  }
  genes_sp1 <- get_sim_genes(similarity_sp1)
  genes_sp2 <- get_sim_genes(similarity_sp2)

  # Check which reference genes are present
  in_sp1 <- ref$gene_sp1 %in% genes_sp1
  in_sp2 <- ref$gene_sp2 %in% genes_sp2
  in_both <- in_sp1 & in_sp2
  n_effective <- sum(in_both)

  missing_sp1 <- ref$gene_sp1[!in_sp1]
  missing_sp2 <- ref$gene_sp2[!in_sp2]

  # Build warnings
  warnings <- character(0)
  if (n_effective < 10) {
    warnings <- c(warnings,
      sprintf("CRITICAL: Only %d effective reference orthologs. CCS values will be unreliable.", n_effective))
  } else if (n_effective < 20) {
    warnings <- c(warnings,
      sprintf("WARNING: Only %d effective reference orthologs. CCS estimates may be noisy.", n_effective))
  } else if (n_effective < 50) {
    warnings <- c(warnings,
      sprintf("Note: %d effective reference orthologs. Consider >50 for robust CCS.", n_effective))
  }

  pct_missing <- if (n_1to1 > 0) (1 - n_effective / n_1to1) * 100 else 0
  if (pct_missing > 20) {
    warnings <- c(warnings,
      sprintf("WARNING: %.1f%% of 1:1 reference genes missing from similarity matrices.", pct_missing))
  }

  # CCS distribution analysis
  ccs_stats <- NULL
  if (!is.null(ccs_results)) {
    if (!"CCS" %in% colnames(ccs_results)) {
      warnings <- c(warnings, "ccs_results provided but missing CCS column.")
    } else {
      # Filter to 1:1 reference pairs
      if ("type" %in% colnames(ccs_results)) {
        ref_ccs <- ccs_results$CCS[ccs_results$type == "1:1"]
      } else {
        ref_ccs <- ccs_results$CCS
      }
      ref_ccs <- ref_ccs[!is.na(ref_ccs)]

      if (length(ref_ccs) >= 3) {
        ccs_mean <- mean(ref_ccs)
        ccs_sd <- stats::sd(ref_ccs)
        ccs_median <- stats::median(ref_ccs)
        ccs_iqr <- stats::IQR(ref_ccs)
        # Skewness (Pearson's moment coefficient)
        ccs_skew <- if (ccs_sd > 0) {
          mean(((ref_ccs - ccs_mean) / ccs_sd)^3)
        } else {
          0
        }

        ccs_stats <- list(
          n = length(ref_ccs),
          median = ccs_median,
          mean = ccs_mean,
          sd = ccs_sd,
          skewness = ccs_skew,
          iqr = ccs_iqr,
          q25 = stats::quantile(ref_ccs, 0.25, names = FALSE),
          q75 = stats::quantile(ref_ccs, 0.75, names = FALSE)
        )

        if (abs(ccs_skew) > 1) {
          warnings <- c(warnings,
            sprintf("CCS distribution is skewed (skewness = %.2f). Reference may be biased.", ccs_skew))
        }
      }
    }
  }

  # Print summary
  message(sprintf("Reference set diagnostics:"))
  message(sprintf("  Total orthologs: %d", n_total))
  message(sprintf("  1:1 orthologs: %d (%.1f%%)", n_1to1,
                  if (n_total > 0) n_1to1 / n_total * 100 else 0))
  message(sprintf("  Effective (in both matrices): %d (%.1f%% of 1:1)",
                  n_effective,
                  if (n_1to1 > 0) n_effective / n_1to1 * 100 else 0))
  if (length(missing_sp1) > 0) {
    message(sprintf("  Missing from sp1 matrix: %d genes", length(missing_sp1)))
  }
  if (length(missing_sp2) > 0) {
    message(sprintf("  Missing from sp2 matrix: %d genes", length(missing_sp2)))
  }
  if (!is.null(ccs_stats)) {
    message(sprintf("  CCS: median=%.3f, mean=%.3f, sd=%.3f, skew=%.2f",
                    ccs_stats$median, ccs_stats$mean, ccs_stats$sd,
                    ccs_stats$skewness))
  }
  if (length(warnings) > 0) {
    message(paste(sprintf("  * %s", warnings), collapse = "\n"))
  }

  list(
    n_total = n_total,
    n_1to1 = n_1to1,
    pct_1to1 = if (n_total > 0) n_1to1 / n_total * 100 else 0,
    n_effective = n_effective,
    pct_effective = if (n_1to1 > 0) n_effective / n_1to1 * 100 else 0,
    missing_sp1 = missing_sp1,
    missing_sp2 = missing_sp2,
    ccs_stats = ccs_stats,
    warnings = warnings
  )
}
