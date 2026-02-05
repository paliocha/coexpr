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
#'
#' @return Symmetric similarity matrix (genes × genes) with values in the range
#'   0 to 1. Higher values indicate stronger co-expression.
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
#' # Calculate similarity (default cached method, fastest)
#' sim <- calculate_pcc_mr(expr)
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
                             mr_method = c("cached", "streaming")) {
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
    return(sim_pcc)
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
    sim_pcc_mr <- mutual_rank_transform_cached_cpp(sim_pcc, effective_cores)
  } else {
    # Streaming method: compute ranks on-demand, slower but memory-efficient
    if (effective_cores > 1) {
      message(sprintf("Applying mutual rank normalization (C++ streaming, %d cores)...", effective_cores))
    } else {
      message("Applying mutual rank normalization (C++ streaming)...")
    }
    sim_pcc_mr <- mutual_rank_transform_cpp(sim_pcc, effective_cores)
  }

  rownames(sim_pcc_mr) <- rownames(expr_matrix)
  colnames(sim_pcc_mr) <- rownames(expr_matrix)

  return(sim_pcc_mr)
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
