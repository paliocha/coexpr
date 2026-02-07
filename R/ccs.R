#' Calculate co-expression correlation score (CCS)
#'
#' Calculates CCS between orthologous genes by correlating their co-expression
#' patterns with a set of reference 1:1 orthologs. This measures how well
#' orthologous genes retain the same co-expression partners across species.
#'
#' @param sim_sp1 Similarity matrix for species 1. Can be a regular matrix or
#'   a TriSimilarity object.
#' @param sim_sp2 Similarity matrix for species 2. Can be a regular matrix or
#'   a TriSimilarity object.
#' @param orthologs Data frame with ortholog pairs. Must have columns:
#'   - `gene_sp1`: Gene IDs from species 1
#'   - `gene_sp2`: Gene IDs from species 2
#'   - `type` (optional): Ortholog type ("1:1", "1:N", "N:1", "N:M")
#' @param use_only_1to1 Logical. If TRUE (default), uses only 1:1 orthologs
#'   as reference (recommended, matches paper). If FALSE, uses all orthologs.
#' @param n_cores Integer. Number of cores for parallel computation. Default is
#'   1 (sequential). Set to a higher value to enable parallel processing via
#'   the future/furrr framework. Use `parallel::detectCores() - 1` for maximum
#'   parallelization while leaving one core free.
#' @param handle_self_diagonal Character. How to handle self-correlation values
#'   when a gene being evaluated is also in the reference set:
#'   - `"mean"` (default): Replace diagonal with mean of non-diagonal values.
#'     This matches the original CoExCorr implementation.
#'   - `"na"`: Set diagonal to NA and use pairwise complete observations.
#'   - `"none"`: Keep diagonal value as-is (may inflate CCS).
#'
#' @return Data frame with columns:
#'   - `gene_sp1`, `gene_sp2`: Ortholog pair
#'   - `CCS`: Co-expression correlation score
#'   - `n_ref`: Number of reference orthologs used
#'
#' @details
#' For ortholog pair (i, j), CCS is calculated as:
#'   CCS_ij = cor(S^sp1_r,i, S^sp2_r,j)
#'
#' Where r is the set of 1:1 reference orthologs. Only the rows/columns
#' corresponding to reference orthologs are correlated.
#'
#' The reference set should contain unduplicated genes performing ancestral
#' functions with conserved expression patterns.
#'
#' ## Self-Diagonal Handling
#'
#' When calculating CCS for an ortholog pair (i, j), if gene i is also in the
#' reference set, the co-expression vector for species 1 will include
#' `sim[i,i] = 1.0` (self-correlation). Similarly for gene j in species 2.
#' This can inflate CCS values because self-correlations are always 1.0 in
#' both species.
#'
#' The `handle_self_diagonal` parameter controls how to address this:
#' - `"mean"`: Replaces the diagonal value with the mean of the other values
#'   in that column of the reference subset. This is the most conservative.
#' - `"na"`: Sets diagonal to NA and uses `cor(..., use = "pairwise.complete.obs")`.
#' - `"none"`: Keeps the diagonal value, which may slightly inflate CCS.
#'
#' ## Parallelization
#'
#' CCS calculation is embarrassingly parallel - each ortholog pair can be
#' computed independently. For large datasets (>10,000 pairs), parallel
#' processing can provide 4-8x speedup on multi-core machines.
#'
#' @references
#' Grønvold L, Hvidsten TR. Cross species comparative transcriptomics using
#' co-expression networks.
#'
#' @examples
#' \dontrun{
#' # Calculate similarity matrices
#' sim_at <- calculate_pcc_mr(expr_at)
#' sim_os <- calculate_pcc_mr(expr_os)
#'
#' # Define orthologs
#' orthologs <- data.frame(
#'   gene_sp1 = c("AT1G01010", "AT1G01020"),
#'   gene_sp2 = c("Os01g0100100", "Os01g0100200"),
#'   type = c("1:1", "1:1")
#' )
#'
#' # Calculate CCS (sequential)
#' ccs_results <- calculate_ccs(sim_at, sim_os, orthologs)
#'
#' # Calculate CCS (parallel with 4 cores)
#' ccs_results <- calculate_ccs(sim_at, sim_os, orthologs, n_cores = 4)
#'
#' # Without self-diagonal handling (may give slightly higher values)
#' ccs_no_fix <- calculate_ccs(sim_at, sim_os, orthologs,
#'                              handle_self_diagonal = "none")
#' }
#'
#' @export
calculate_ccs <- function(sim_sp1, sim_sp2, orthologs,
                          use_only_1to1 = TRUE,
                          n_cores = 1,
                          handle_self_diagonal = c("mean", "na", "none")) {

  handle_self_diagonal <- match.arg(handle_self_diagonal)

  # Validate inputs - accept matrices or TriSimilarity objects
  is_valid_sim <- function(x) {
    is.matrix(x) || is(x, "TriSimilarity")
  }

  if (!is_valid_sim(sim_sp1) || !is_valid_sim(sim_sp2)) {
    stop("Similarity matrices must be matrices or TriSimilarity objects")
  }

  # Get gene names from similarity matrices
  get_genes <- function(sim) {
    if (is(sim, "TriSimilarity")) {
      sim@genes
    } else {
      rownames(sim)
    }
  }

  genes_sp1 <- get_genes(sim_sp1)
  genes_sp2 <- get_genes(sim_sp2)

  required_cols <- c("gene_sp1", "gene_sp2")
  if (!all(required_cols %in% colnames(orthologs))) {
    stop(sprintf("orthologs must have columns: %s",
                 paste(required_cols, collapse = ", ")))
  }

  # Define reference orthologs (1:1 only)
  if (use_only_1to1) {
    if ("type" %in% colnames(orthologs)) {
      ref_orthologs <- orthologs |>
        dplyr::filter(.data$type == "1:1")
    } else {
      # If no type column, detect 1:1 by uniqueness
      ref_orthologs <- orthologs |>
        dplyr::group_by(.data$gene_sp1) |>
        dplyr::filter(dplyr::n() == 1) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$gene_sp2) |>
        dplyr::filter(dplyr::n() == 1) |>
        dplyr::ungroup()
    }

    message(sprintf("Using %d 1:1 reference orthologs (from %d total)",
                    nrow(ref_orthologs), nrow(orthologs)))
  } else {
    ref_orthologs <- orthologs
    message(sprintf("Using all %d orthologs as reference", nrow(ref_orthologs)))
  }

  if (nrow(ref_orthologs) < 10) {
    stop("Need at least 10 reference orthologs for reliable CCS calculation")
  }

  # Filter reference genes that exist in both similarity matrices
  ref_sp1 <- intersect(ref_orthologs$gene_sp1, genes_sp1)
  ref_sp2 <- intersect(ref_orthologs$gene_sp2, genes_sp2)

  ref_orthologs_filt <- ref_orthologs |>
    dplyr::filter(.data$gene_sp1 %in% ref_sp1, .data$gene_sp2 %in% ref_sp2)

  if (nrow(ref_orthologs_filt) < 10) {
    stop(sprintf(
      "Only %d reference orthologs found in similarity matrices (need >=10)",
      nrow(ref_orthologs_filt)
    ))
  }

  # Filter orthologs to those present in similarity matrices
  orthologs_filt <- orthologs |>
    dplyr::filter(.data$gene_sp1 %in% genes_sp1,
                  .data$gene_sp2 %in% genes_sp2)

  n_pairs <- nrow(orthologs_filt)
  n_ref <- nrow(ref_orthologs_filt)

  message(sprintf("Calculating CCS for %d ortholog pairs using %d references",
                  n_pairs, n_ref))

  # Pre-extract reference gene vectors
  ref_genes_sp1 <- ref_orthologs_filt$gene_sp1
  ref_genes_sp2 <- ref_orthologs_filt$gene_sp2

  # Helper function to extract column from similarity matrix
  extract_sim_column <- function(sim, gene) {
    if (is(sim, "TriSimilarity")) {
      extractColumn(sim, gene)
    } else {
      sim[, gene]
    }
  }

  # Helper function to handle self-diagonal
  handle_diagonal <- function(coexpr, gene, ref_genes, method) {
    if (method == "none") {
      return(coexpr)
    }

    # Check if gene is in reference set
    self_idx <- which(ref_genes == gene)

    if (length(self_idx) == 0) {
      # Gene not in reference set, no diagonal to handle
      return(coexpr)
    }

    if (method == "mean") {
      # Replace diagonal with mean of non-diagonal values
      non_diag_vals <- coexpr[-self_idx]
      coexpr[self_idx] <- mean(non_diag_vals, na.rm = TRUE)
    } else if (method == "na") {
      # Set diagonal to NA
      coexpr[self_idx] <- NA
    }

    coexpr
  }

  # Determine correlation method based on diagonal handling
  cor_use <- if (handle_self_diagonal == "na") {
    "pairwise.complete.obs"
  } else {
    "everything"
  }

  # Calculate CCS for each ortholog pair
  if (n_cores > 1) {
    # Use fork-based parallelism on Unix/macOS (shares memory, no copying)
    if (.Platform$OS.type == "unix") {
      message(sprintf("Using parallel computation with %d cores (fork-based)", n_cores))

      # mclapply uses forking which shares memory - no need to copy large matrices
      ccs_values <- parallel::mclapply(
        seq_len(n_pairs),
        function(i) {
          gene_sp1 <- orthologs_filt$gene_sp1[i]
          gene_sp2 <- orthologs_filt$gene_sp2[i]

          coexpr_sp1 <- extract_sim_column(sim_sp1, gene_sp1)[ref_genes_sp1]
          coexpr_sp2 <- extract_sim_column(sim_sp2, gene_sp2)[ref_genes_sp2]

          # Handle self-diagonal
          coexpr_sp1 <- handle_diagonal(coexpr_sp1, gene_sp1, ref_genes_sp1, handle_self_diagonal)
          coexpr_sp2 <- handle_diagonal(coexpr_sp2, gene_sp2, ref_genes_sp2, handle_self_diagonal)

          stats::cor(coexpr_sp1, coexpr_sp2, method = "pearson", use = cor_use)
        },
        mc.cores = n_cores
      )
      ccs_values <- unlist(ccs_values)

    } else {
      # Windows: use furrr with chunked processing to reduce memory
      message(sprintf("Using parallel computation with %d cores (multisession)", n_cores))

      # Pre-extract the relevant columns for all ortholog pairs
      if (is(sim_sp1, "TriSimilarity")) {
        sim_sp1_ref <- extractRows(sim_sp1, ref_genes_sp1)
        sim_sp1_cols <- sim_sp1_ref[, orthologs_filt$gene_sp1, drop = FALSE]
      } else {
        sim_sp1_cols <- sim_sp1[ref_genes_sp1, orthologs_filt$gene_sp1, drop = FALSE]
      }

      if (is(sim_sp2, "TriSimilarity")) {
        sim_sp2_ref <- extractRows(sim_sp2, ref_genes_sp2)
        sim_sp2_cols <- sim_sp2_ref[, orthologs_filt$gene_sp2, drop = FALSE]
      } else {
        sim_sp2_cols <- sim_sp2[ref_genes_sp2, orthologs_filt$gene_sp2, drop = FALSE]
      }

      # Set up parallel backend with increased global size limit
      old_limit <- getOption("future.globals.maxSize")
      options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB limit
      on.exit(options(future.globals.maxSize = old_limit), add = TRUE)

      oplan <- future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(oplan), add = TRUE)

      # Get gene names for diagonal handling
      pair_genes_sp1 <- orthologs_filt$gene_sp1
      pair_genes_sp2 <- orthologs_filt$gene_sp2

      # Process in parallel using furrr with pre-extracted data
      ccs_values <- furrr::future_map_dbl(
        seq_len(n_pairs),
        function(i) {
          coexpr_sp1 <- sim_sp1_cols[, i]
          coexpr_sp2 <- sim_sp2_cols[, i]

          # Handle self-diagonal
          coexpr_sp1 <- handle_diagonal(coexpr_sp1, pair_genes_sp1[i], ref_genes_sp1, handle_self_diagonal)
          coexpr_sp2 <- handle_diagonal(coexpr_sp2, pair_genes_sp2[i], ref_genes_sp2, handle_self_diagonal)

          stats::cor(coexpr_sp1, coexpr_sp2, method = "pearson", use = cor_use)
        },
        .options = furrr::furrr_options(seed = TRUE)
      )
    }

    ccs_results <- orthologs_filt |>
      dplyr::mutate(
        CCS = ccs_values,
        n_ref = n_ref
      )

  } else {
    # Sequential computation
    ccs_results <- orthologs_filt |>
      dplyr::rowwise() |>
      dplyr::mutate(
        CCS = calculate_ccs_pair_internal(
          .data$gene_sp1, .data$gene_sp2,
          sim_sp1, sim_sp2,
          ref_genes_sp1, ref_genes_sp2,
          handle_self_diagonal, cor_use
        ),
        n_ref = n_ref
      ) |>
      dplyr::ungroup()
  }

  return(ccs_results)
}


#' Calculate CCS for a single ortholog pair (internal)
#'
#' @param gene_sp1 Gene ID in species 1
#' @param gene_sp2 Gene ID in species 2
#' @param sim_sp1 Similarity matrix for species 1
#' @param sim_sp2 Similarity matrix for species 2
#' @param ref_genes_sp1 Reference gene IDs for species 1
#' @param ref_genes_sp2 Reference gene IDs for species 2
#' @param handle_self_diagonal Method for diagonal handling
#' @param cor_use Argument for cor() 'use' parameter
#'
#' @return CCS value
#' @keywords internal
#' @noRd
calculate_ccs_pair_internal <- function(gene_sp1, gene_sp2, sim_sp1, sim_sp2,
                                        ref_genes_sp1, ref_genes_sp2,
                                        handle_self_diagonal, cor_use) {

  # Helper function to extract column from similarity matrix
  extract_sim_column <- function(sim, gene) {
    if (is(sim, "TriSimilarity")) {
      extractColumn(sim, gene)
    } else {
      sim[, gene]
    }
  }

  # Helper function to handle self-diagonal
  handle_diagonal <- function(coexpr, gene, ref_genes, method) {
    if (method == "none") {
      return(coexpr)
    }

    self_idx <- which(ref_genes == gene)

    if (length(self_idx) == 0) {
      return(coexpr)
    }

    if (method == "mean") {
      non_diag_vals <- coexpr[-self_idx]
      coexpr[self_idx] <- mean(non_diag_vals, na.rm = TRUE)
    } else if (method == "na") {
      coexpr[self_idx] <- NA
    }

    coexpr
  }

  # Extract co-expression vectors with reference genes
  coexpr_sp1 <- extract_sim_column(sim_sp1, gene_sp1)[ref_genes_sp1]
  coexpr_sp2 <- extract_sim_column(sim_sp2, gene_sp2)[ref_genes_sp2]

  # Handle self-diagonal
  coexpr_sp1 <- handle_diagonal(coexpr_sp1, gene_sp1, ref_genes_sp1, handle_self_diagonal)
  coexpr_sp2 <- handle_diagonal(coexpr_sp2, gene_sp2, ref_genes_sp2, handle_self_diagonal)

  # Calculate Pearson correlation
  ccs <- cor(coexpr_sp1, coexpr_sp2, method = "pearson", use = cor_use)

  return(ccs)
}


#' Calculate CCS for a single ortholog pair (legacy wrapper)
#'
#' @param gene_sp1 Gene ID in species 1
#' @param gene_sp2 Gene ID in species 2
#' @param sim_sp1 Similarity matrix for species 1
#' @param sim_sp2 Similarity matrix for species 2
#' @param ref_orthologs Reference ortholog pairs
#'
#' @return CCS value
#' @keywords internal
#' @noRd
calculate_ccs_pair <- function(gene_sp1, gene_sp2, sim_sp1, sim_sp2, ref_orthologs) {
  calculate_ccs_pair_internal(
    gene_sp1, gene_sp2, sim_sp1, sim_sp2,
    ref_orthologs$gene_sp1, ref_orthologs$gene_sp2,
    "mean", "everything"
  )
}
