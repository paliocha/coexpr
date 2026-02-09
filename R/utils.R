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
      dplyr::summarize(
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
      dplyr::summarize(
        n = dplyr::n(),
        median_logORS = median(.data$logORS, na.rm = TRUE),
        mean_logORS = mean(.data$logORS, na.rm = TRUE),
        sd_logORS = sd(.data$logORS, na.rm = TRUE),
        pct_top10 = sum(.data$logORS > 1, na.rm = TRUE) / dplyr::n() * 100,
        pct_top1 = sum(.data$logORS > 2, na.rm = TRUE) / dplyr::n() * 100
      )
  }

  return(summary_df)
}


#' Validate ortholog table format
#'
#' @param orthologs Ortholog data frame to validate
#'
#' @return TRUE if valid, throws error otherwise
#'
#' @keywords internal
validate_orthologs <- function(orthologs) {

  if (!is.data.frame(orthologs)) {
    stop("orthologs must be a data frame")
  }

  required_cols <- c("gene_sp1", "gene_sp2")
  missing_cols <- setdiff(required_cols, colnames(orthologs))
  if (length(missing_cols) > 0) {
    stop(sprintf("orthologs missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  if (nrow(orthologs) == 0) {
    stop("orthologs is empty")
  }

  if (any(is.na(orthologs$gene_sp1)) || any(is.na(orthologs$gene_sp2))) {
    stop("orthologs contains NA values in gene columns")
  }

  TRUE
}
