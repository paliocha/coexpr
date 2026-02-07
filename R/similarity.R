#' Calculate PCC+MR similarity matrix
#'
#' Calculates co-expression similarity using Pearson correlation followed by
#' mutual rank normalization (Obayashi et al. 2009). This is the best-performing
#' method identified in Grønvold & Hvidsten for cross-species comparisons.
#'
#' @param expr_matrix Numeric matrix of gene expression (genes × samples).
#'   Rows are genes, columns are samples. Should be log-transformed.
#' @param method Character. Either "pcc_mr" (default, recommended) for Pearson
#'   correlation with mutual rank, or "pcc" for Pearson only.
#' @param n_cores Integer. Number of cores for parallel mutual rank computation
#'   via OpenMP. Default is 1 (sequential).
#' @param chunk_size Deprecated. No longer used (kept for API compatibility).
#' @param mr_method Character. Method for mutual rank computation:
#'   - "cached" (default): Precomputes all row ranks, faster but uses more memory
#'     (O(n^2) for rank matrix). Recommended for most cases.
#'   - "streaming": Computes ranks on-demand, slower but more memory-efficient.
#'     Use for very large matrices that cause memory issues with "cached".
#' @param return_tri Logical. If TRUE (default), returns a TriSimilarity object
#'   that stores only the upper triangle for ~50% memory savings. If FALSE,
#'   returns a full symmetric matrix. Note: CCS calculation works with both formats.
#'
#' @return If return_tri = TRUE (default), a TriSimilarity object.
#'   If return_tri = FALSE, a symmetric similarity matrix (genes × genes).
#'   Values are in the range 0 to 1. Higher values indicate stronger co-expression.
#'
#' @details
#' The PCC+MR method:
#' 1. Calculates Pearson correlation: S^PCC_ij = cor(E_i, E_j)
#' 2. Transforms using log mutual rank: S^(PCC+MR)_ij = 1 - log(sqrt(R_ij * R_ji)) / log(n)
#'
#' Where R_ij is the rank of S^PCC_ij in row i (ordered high to low).
#' The log transformation emphasizes strong correlations and the mutual rank
#' reduces bias from large co-expression clusters.
#'
#' ## Triangular Storage
#'
#' By default, the function returns a TriSimilarity object that stores only
#' the upper triangle of the symmetric matrix. This reduces memory usage by
#' approximately 50%. Use `as.matrix()` to convert to a full matrix if needed.
#'
#' ## Mutual Rank Methods
#'
#' The mutual rank transformation is implemented in C++ for efficiency:
#' - **cached**: O(n^2 log n) time, O(n^2) extra memory. Precomputes all row
#'   ranks before computing mutual ranks. Much faster for typical datasets.
#' - **streaming**: O(n^3) time, O(n) extra memory per thread. Computes ranks
#'   on-demand for each row pair. Use only if "cached" causes memory issues.
#'
#' ## Parallelization
#'
#' The mutual rank step uses OpenMP for parallel computation when available.
#' The correlation step uses Rfast::cora which is a fast C++ implementation.
#'
#' @references
#' Obayashi T, Kinoshita K (2009) Rank of correlation coefficient as a comparable
#' measure for biological significance of gene coexpression. DNA Research 16:249-260.
#'
#' Grønvold L, Hvidsten TR. Cross species comparative transcriptomics using
#' co-expression networks.
#'
#' @examples
#' \dontrun{
#' # Load expression data
#' expr <- matrix(rnorm(1000 * 50), nrow = 1000, ncol = 50)
#' rownames(expr) <- paste0("Gene", 1:1000)
#'
#' # Calculate similarity (returns TriSimilarity by default)
#' sim <- calculate_pcc_mr(expr)
#' print(sim)  # Shows memory savings
#'
#' # Convert to full matrix if needed
#' sim_full <- as.matrix(sim)
#'
#' # Or return full matrix directly (uses more memory)
#' sim_full <- calculate_pcc_mr(expr, return_tri = FALSE)
#'
#' # For very large matrices with memory constraints, use streaming
#' sim <- calculate_pcc_mr(large_expr, mr_method = "streaming")
#'
#' # Parallel with 4 cores
#' sim <- calculate_pcc_mr(expr, n_cores = 4)
#' }
#'
#' @export
calculate_pcc_mr <- function(expr_matrix, method = c("pcc_mr", "pcc"),
                             n_cores = 1, chunk_size = NULL,
                             mr_method = c("cached", "streaming"),
                             return_tri = TRUE) {
  method <- match.arg(method)
  mr_method <- match.arg(mr_method)

  # Validate input
  if (!is.matrix(expr_matrix)) {
    stop("expr_matrix must be a matrix")
  }
  if (any(is.na(expr_matrix))) {
    warning("expr_matrix contains NA values, these will propagate to similarity matrix")
  }

  n_genes <- nrow(expr_matrix)

  # Calculate Pearson correlation using Rfast (very fast C++ implementation)
  # Note: We avoid furrr/future here as it spawns separate processes that
  # copy the expression matrix, causing high memory usage. Rfast::cora is
  # already highly optimized and fast enough for most matrices.
  message("Calculating Pearson correlation...")
  sim_pcc <- Rfast::cora(t(expr_matrix))

  # Clip to [-1, 1] due to numerical precision issues
  sim_pcc[sim_pcc > 1] <- 1
  sim_pcc[sim_pcc < -1] <- -1

  rownames(sim_pcc) <- rownames(expr_matrix)
  colnames(sim_pcc) <- rownames(expr_matrix)

  if (method == "pcc") {
    if (return_tri) {
      return(as.TriSimilarity(sim_pcc))
    } else {
      return(sim_pcc)
    }
  }

  # Apply mutual rank normalization using C++ implementation
  if (has_openmp()) {
    effective_cores <- min(n_cores, get_max_threads())
  } else {
    effective_cores <- 1
  }

  if (mr_method == "cached") {
    # Cached method: precompute all row ranks, faster but uses O(n^2) extra memory
    if (effective_cores > 1) {
      message(sprintf("Applying mutual rank normalization (C++ cached, %d cores)...", effective_cores))
    } else {
      message("Applying mutual rank normalization (C++ cached)...")
    }

    if (return_tri) {
      # Use triangular output function for memory efficiency
      result <- mutual_rank_transform_tri_cpp(sim_pcc, effective_cores)
      return(TriSimilarity(
        data = result$data,
        genes = rownames(expr_matrix),
        diag_value = result$diag_value
      ))
    } else {
      sim_pcc_mr <- mutual_rank_transform_cached_cpp(sim_pcc, effective_cores)
    }
  } else {
    # Streaming method: compute ranks on-demand, slower but memory-efficient
    if (effective_cores > 1) {
      message(sprintf("Applying mutual rank normalization (C++ streaming, %d cores)...", effective_cores))
    } else {
      message("Applying mutual rank normalization (C++ streaming)...")
    }
    sim_pcc_mr <- mutual_rank_transform_cpp(sim_pcc, effective_cores)

    if (return_tri) {
      rownames(sim_pcc_mr) <- rownames(expr_matrix)
      colnames(sim_pcc_mr) <- rownames(expr_matrix)
      return(as.TriSimilarity(sim_pcc_mr))
    }
  }

  rownames(sim_pcc_mr) <- rownames(expr_matrix)
  colnames(sim_pcc_mr) <- rownames(expr_matrix)

  return(sim_pcc_mr)
}


#' Calculate MI+CLR similarity matrix
#'
#' Calculates co-expression similarity using mutual information (MI) followed by
#' context likelihood ratio (CLR) normalization. This is an alternative to PCC+MR
#' that was tested in Grønvold & Hvidsten. CLR helps address hub gene bias by
#' normalizing by the background distribution of MI values.
#'
#' The computation is implemented in C++ with OpenMP parallelization for
#' efficient processing of large gene expression matrices.
#'
#' @param expr_matrix Numeric matrix of gene expression (genes × samples).
#'   Rows are genes, columns are samples. Should be log-transformed.
#' @param n_bins Integer. Number of bins for discretization. Default is 10.
#'   More bins capture finer relationships but require more samples.
#' @param n_cores Integer. Number of cores for parallel computation via OpenMP.
#'   Default is 1 (sequential).
#' @param return_tri Logical. If TRUE (default), returns a TriSimilarity object.
#'   If FALSE, returns a full symmetric matrix.
#'
#' @return If return_tri = TRUE (default), a TriSimilarity object (S4 class).
#'   If return_tri = FALSE, a symmetric similarity matrix (genes × genes).
#'   Values are non-negative (CLR transformation can produce values > 1).
#'
#' @details
#' The MI+CLR method:
#' 1. Discretizes expression values into equal-frequency bins
#' 2. Calculates mutual information: MI(X,Y) = H(X) + H(Y) - H(X,Y)
#' 3. Applies CLR transformation to normalize by background:
#'    CLR_ij = sqrt(max(0, Z_row)^2 + max(0, Z_col)^2)
#'    where Z_row and Z_col are z-scores of MI_ij within its row and column
#'
#' The CLR transformation (Faith et al. 2007) reduces the bias from hub genes
#' that have high MI with many other genes.
#'
#' ## C++ Implementation
#'
#' This function uses a fast C++ implementation with OpenMP parallelization:
#' - Discretization: O(n * m log m) where n = genes, m = samples
#' - MI computation: O(n^2 * m) for all pairwise comparisons
#' - CLR transformation: O(n^2)
#'
#' For a 1000-gene matrix with 100 samples, expect ~10-20x speedup compared
#' to pure R implementation when using 4 cores.
#'
#' ## Sample Requirements
#'
#' MI estimation requires sufficient samples for reliable discretization.
#' A rule of thumb is at least 3 * n_bins samples. For 10 bins, this means
#' at least 30 samples.
#'
#' @references
#' Faith JJ et al. (2007) Large-scale mapping and validation of Escherichia coli
#' transcriptional regulation from a compendium of expression profiles.
#' PLoS Biology 5:e8.
#'
#' Grønvold L, Hvidsten TR. Cross species comparative transcriptomics using
#' co-expression networks.
#'
#' @examples
#' \dontrun{
#' # Calculate MI+CLR similarity
#' sim <- calculate_mi_clr(expr, n_bins = 10)
#'
#' # Use with CCS calculation
#' ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)
#'
#' # Parallel computation with 4 cores
#' sim <- calculate_mi_clr(expr, n_bins = 10, n_cores = 4)
#' }
#'
#' @export
calculate_mi_clr <- function(expr_matrix,
                             n_bins = 10,
                             n_cores = 1,
                             return_tri = TRUE) {

  # Validate input
  if (!is.matrix(expr_matrix)) {
    stop("expr_matrix must be a matrix")
  }

  if (ncol(expr_matrix) < 3 * n_bins) {
    warning(sprintf(
      "MI estimation works best with at least %d samples (3 * n_bins). You have %d samples. Consider reducing n_bins.",
      3 * n_bins, ncol(expr_matrix)
    ))
  }

  n_genes <- nrow(expr_matrix)
  gene_names <- rownames(expr_matrix)

  if (is.null(gene_names)) {
    gene_names <- paste0("Gene", seq_len(n_genes))
    rownames(expr_matrix) <- gene_names
  }

  # Handle zero-variance genes - warn but let C++ handle it
  gene_vars <- apply(expr_matrix, 1, stats::var, na.rm = TRUE)
  zero_var_genes <- which(gene_vars == 0 | is.na(gene_vars))

  if (length(zero_var_genes) > 0) {
    warning(sprintf(
      "%d genes have zero variance and will have MI = 0 with all other genes",
      length(zero_var_genes)
    ))
  }

  # Determine effective core count
  if (has_openmp()) {
    effective_cores <- min(n_cores, get_max_threads())
  } else {
    effective_cores <- 1
  }

  if (effective_cores > 1) {
    message(sprintf("Computing MI+CLR with %d cores...", effective_cores))
  }

  # Call C++ implementation
  if (return_tri) {
    result <- compute_mi_clr_tri_cpp(expr_matrix, as.integer(n_bins), effective_cores)
    return(TriSimilarity(
      data = result$data,
      genes = gene_names,
      diag_value = result$diag_value
    ))
  } else {
    clr_matrix <- compute_mi_clr_cpp(expr_matrix, as.integer(n_bins), effective_cores)
    rownames(clr_matrix) <- gene_names
    colnames(clr_matrix) <- gene_names
    return(clr_matrix)
  }
}


#' Validate expression matrix format
#'
#' @param expr_matrix Expression matrix to validate
#' @param name Character name for error messages
#'
#' @return TRUE if valid, throws error otherwise
#' @keywords internal
#' @noRd
validate_expr_matrix <- function(expr_matrix, name = "expr_matrix") {
  if (!is.matrix(expr_matrix)) {
    stop(sprintf("%s must be a matrix, got %s", name, class(expr_matrix)[1]))
  }

  if (is.null(rownames(expr_matrix))) {
    stop(sprintf("%s must have rownames (gene IDs)", name))
  }

  if (any(duplicated(rownames(expr_matrix)))) {
    stop(sprintf("%s has duplicated rownames", name))
  }

  if (ncol(expr_matrix) < 3) {
    stop(sprintf("%s must have at least 3 samples", name))
  }

  TRUE
}
