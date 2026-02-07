#' coexpr: Cross-Species Co-Expression Network Analysis
#'
#' Implementation of cross-species comparative transcriptomics using co-expression
#' networks based on Grønvold & Hvidsten. Provides functions to calculate
#' co-expression correlation scores (CCS) and ortholog rank scores (ORS) to
#' identify conserved gene regulation patterns across species.
#'
#' @section Main functions:
#' - [calculate_pcc_mr()]: Calculate PCC+MR similarity matrix
#' - [calculate_ccs()]: Calculate co-expression correlation score
#' - [calculate_ors()]: Calculate ortholog rank score
#' - [handle_multicopy_orthologs()]: Handle 1:N and N:M ortholog relationships
#'
#' @section Typical workflow:
#' 1. Calculate similarity matrices for each species using [calculate_pcc_mr()]
#' 2. Process orthologs with [handle_multicopy_orthologs()]
#' 3. Calculate CCS with [calculate_ccs()]
#' 4. Calculate ORS with [calculate_ors()]
#' 5. Test significance with [test_ors_significance()]
#'
#' @section Parallelization:
#' Both [calculate_pcc_mr()] and [calculate_ccs()] support parallel computation
#' via the `n_cores` parameter. This uses the future/furrr framework for
#' cross-platform parallel processing.
#'
#' Example:
#' ```
#' # Use 4 cores for CCS calculation
#' ccs_results <- calculate_ccs(sim1, sim2, orthologs, n_cores = 4)
#'
#' # Use all available cores minus one
#' n <- parallel::detectCores() - 1
#' sim <- calculate_pcc_mr(expr, n_cores = n)
#' ```
#'
#' CCS calculation sees the largest speedup (4-8x on typical hardware) since
#' it is embarrassingly parallel. Similarity matrix calculation benefits from
#' parallelization primarily for large datasets (>20,000 genes).
#'
#' @docType package
#' @name coexpr-package
#' @aliases coexpr
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @importFrom stats cor median sd
#' @importFrom future plan multisession
#' @importFrom furrr future_map future_map_dbl furrr_options
## usethis namespace: end
NULL


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom dplyr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


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

  if (!"logors" %in% colnames(ors_results)) {
    stop("ors_results must contain 'logors' from calculate_ors()")
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
#' @export
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
