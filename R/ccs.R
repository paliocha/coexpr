#' Calculate co-expression correlation score (CCS)
#'
#' Calculates CCS between orthologous genes by correlating their co-expression
#' patterns with a set of reference 1:1 orthologs. This measures how well
#' orthologous genes retain the same co-expression partners across species.
#'
#' @param sim_sp1 Similarity matrix for species 1 (genes × genes)
#' @param sim_sp2 Similarity matrix for species 2 (genes × genes)
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
#'
#' @return Data frame with columns:
#'   - `gene_sp1`, `gene_sp2`: Ortholog pair
#'   - `ccs`: Co-expression correlation score
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
#' }
#'
#' @export
calculate_ccs <- function(sim_sp1, sim_sp2, orthologs,
                          use_only_1to1 = TRUE,
                          n_cores = 1) {

  # Validate inputs
  if (!is.matrix(sim_sp1) || !is.matrix(sim_sp2)) {
    stop("Similarity matrices must be matrices")
  }

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
  ref_sp1 <- intersect(ref_orthologs$gene_sp1, rownames(sim_sp1))
  ref_sp2 <- intersect(ref_orthologs$gene_sp2, rownames(sim_sp2))

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
    dplyr::filter(.data$gene_sp1 %in% rownames(sim_sp1),
                  .data$gene_sp2 %in% rownames(sim_sp2))

  n_pairs <- nrow(orthologs_filt)
  n_ref <- nrow(ref_orthologs_filt)

  message(sprintf("Calculating CCS for %d ortholog pairs using %d references",
                  n_pairs, n_ref))

  # Calculate CCS for each ortholog pair
  if (n_cores > 1) {
    # Pre-extract reference gene vectors for efficiency
    ref_genes_sp1 <- ref_orthologs_filt$gene_sp1
    ref_genes_sp2 <- ref_orthologs_filt$gene_sp2

    # Use fork-based parallelism on Unix/macOS (shares memory, no copying)
    if (.Platform$OS.type == "unix") {
      message(sprintf("Using parallel computation with %d cores (fork-based)", n_cores))

      # mclapply uses forking which shares memory - no need to copy large matrices
      ccs_values <- parallel::mclapply(
        seq_len(n_pairs),
        function(i) {
          gene_sp1 <- orthologs_filt$gene_sp1[i]
          gene_sp2 <- orthologs_filt$gene_sp2[i]
          coexpr_sp1 <- sim_sp1[ref_genes_sp1, gene_sp1]
          coexpr_sp2 <- sim_sp2[ref_genes_sp2, gene_sp2]
          stats::cor(coexpr_sp1, coexpr_sp2, method = "pearson")
        },
        mc.cores = n_cores
      )
      ccs_values <- unlist(ccs_values)

    } else {
      # Windows: use furrr with chunked processing to reduce memory
      message(sprintf("Using parallel computation with %d cores (multisession)", n_cores))

      # Pre-extract the relevant rows from similarity matrices to reduce memory
      sim_sp1_ref <- sim_sp1[ref_genes_sp1, orthologs_filt$gene_sp1, drop = FALSE]
      sim_sp2_ref <- sim_sp2[ref_genes_sp2, orthologs_filt$gene_sp2, drop = FALSE]

      # Set up parallel backend with increased global size limit
      old_limit <- getOption("future.globals.maxSize")
      options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB limit
      on.exit(options(future.globals.maxSize = old_limit), add = TRUE)

      oplan <- future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(oplan), add = TRUE)

      # Process in parallel using furrr with pre-extracted data
      ccs_values <- furrr::future_map_dbl(
        seq_len(n_pairs),
        function(i) {
          coexpr_sp1 <- sim_sp1_ref[, i]
          coexpr_sp2 <- sim_sp2_ref[, i]
          stats::cor(coexpr_sp1, coexpr_sp2, method = "pearson")
        },
        .options = furrr::furrr_options(seed = TRUE)
      )
    }

    ccs_results <- orthologs_filt |>
      dplyr::mutate(
        ccs = ccs_values,
        n_ref = n_ref
      )

  } else {
    # Sequential computation (original method)
    ccs_results <- orthologs_filt |>
      dplyr::rowwise() |>
      dplyr::mutate(
        ccs = calculate_ccs_pair(
          .data$gene_sp1, .data$gene_sp2,
          sim_sp1, sim_sp2,
          ref_orthologs_filt
        ),
        n_ref = n_ref
      ) |>
      dplyr::ungroup()
  }

  return(ccs_results)
}


#' Calculate CCS for a single ortholog pair
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

  # Extract co-expression vectors with reference genes
  coexpr_sp1 <- sim_sp1[ref_orthologs$gene_sp1, gene_sp1]
  coexpr_sp2 <- sim_sp2[ref_orthologs$gene_sp2, gene_sp2]

  # Calculate Pearson correlation
  ccs <- cor(coexpr_sp1, coexpr_sp2, method = "pearson")

  return(ccs)
}
