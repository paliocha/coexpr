# ---- Internal similarity matrix helpers ----
# Shared across ccs.R, multicopy.R, and other files that work with
# similarity matrices (regular matrices or TriSimilarity objects).

# Get gene names from a similarity matrix or TriSimilarity object
# @param sim A matrix or TriSimilarity object
# @return Character vector of gene names
# @noRd
sim_genes <- function(sim) {
  if (is(sim, "TriSimilarity")) sim@genes else rownames(sim)
}

# Check whether an object is a valid similarity matrix
# @param x Object to test
# @return Logical
# @noRd
is_valid_sim <- function(x) {
  is.matrix(x) || is(x, "TriSimilarity")
}

# Extract a column from a similarity matrix or TriSimilarity object
# @param sim A matrix or TriSimilarity object
# @param gene Character. Gene name to extract
# @return Named numeric vector
# @noRd
sim_column <- function(sim, gene) {
  if (is(sim, "TriSimilarity")) extractColumn(sim, gene) else sim[, gene]
}

# Validate that a data frame has the required ortholog columns
# @param orthologs Data frame to validate
# @param extra Character vector of additional required columns
# @noRd
check_ortholog_cols <- function(orthologs, extra = character(0)) {
  required <- c("gene_sp1", "gene_sp2", extra)
  missing <- setdiff(required, colnames(orthologs))
  if (length(missing) > 0) {
    stop(sprintf("orthologs must have columns: %s",
                 paste(required, collapse = ", ")))
  }
  invisible(TRUE)
}


#' Summarize ortholog conservation
#'
#' Provides summary statistics for ORS results to quickly assess overall
#' conservation patterns.
#'
#' @param ors_results Data frame from `calculate_ors()`
#' @param by_type Logical. If TRUE and 'type' column exists, summarize by
#'   ortholog type (1:1, 1:N, etc.)
#'
#' @return Data frame with summary statistics
#'
#' @export
summarize_conservation <- function(ors_results, by_type = TRUE) {

  if (!"logORS" %in% colnames(ors_results)) {
    stop("ors_results must contain 'logORS' from calculate_ors()")
  }

  if (by_type && "type" %in% colnames(ors_results)) {
    summary_df <- ors_results |>
      dplyr::group_by(.data$type) |>
      dplyr::summarise(
        n = dplyr::n(),
        median_logORS = median(.data$logORS, na.rm = TRUE),
        mean_logORS = mean(.data$logORS, na.rm = TRUE),
        sd_logORS = sd(.data$logORS, na.rm = TRUE),
        pct_top10 = sum(.data$logORS > 1, na.rm = TRUE) / dplyr::n() * 100,
        pct_top1 = sum(.data$logORS > 2, na.rm = TRUE) / dplyr::n() * 100,
        .groups = "drop"
      )
  } else {
    summary_df <- ors_results |>
      dplyr::summarise(
        n = dplyr::n(),
        median_logORS = median(.data$logORS, na.rm = TRUE),
        mean_logORS = mean(.data$logORS, na.rm = TRUE),
        sd_logORS = sd(.data$logORS, na.rm = TRUE),
        pct_top10 = sum(.data$logORS > 1, na.rm = TRUE) / dplyr::n() * 100,
        pct_top1 = sum(.data$logORS > 2, na.rm = TRUE) / dplyr::n() * 100
      )
  }

  summary_df
}
