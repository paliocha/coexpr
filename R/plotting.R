#' Plot logORS distribution
#'
#' Histogram of logORS values with significance threshold lines at 1
#' (top 10%), 2 (top 1%), and 3 (top 0.1%).
#'
#' @param ors_results Data frame from [calculate_ors()] with `logORS` column.
#' @param bins Integer. Number of histogram bins. Default 50.
#' @param title Character. Plot title. Default "Distribution of logORS".
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' ors_results <- calculate_ors(ccs_results)
#' plot_ors_distribution(ors_results)
#' }
#'
#' @export
plot_ors_distribution <- function(ors_results, bins = 50,
                                  title = "Distribution of logORS") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install with install.packages('ggplot2').")
  }
  if (!"logORS" %in% colnames(ors_results)) {
    stop("ors_results must contain 'logORS' column from calculate_ors()")
  }

  ggplot2::ggplot(ors_results, ggplot2::aes(x = .data$logORS)) +
    ggplot2::geom_histogram(bins = bins, fill = "steelblue", alpha = 0.8,
                            color = "white") +
    ggplot2::geom_vline(xintercept = c(1, 2, 3), linetype = "dashed",
                        color = "red", linewidth = 0.7) +
    ggplot2::annotate("text", x = c(1, 2, 3), y = Inf,
                      label = c("Top 10%", "Top 1%", "Top 0.1%"),
                      vjust = 2, hjust = -0.1, color = "red", size = 3) +
    ggplot2::labs(
      x = "logORS",
      y = "Count",
      title = title,
      subtitle = sprintf("n = %d ortholog pairs", nrow(ors_results))
    ) +
    ggplot2::theme_minimal()
}


#' Plot logORS by ortholog type
#'
#' Boxplot comparing logORS distributions across ortholog types (1:1, 1:N,
#' N:1, N:M). Requires a `type` column in the input data.
#'
#' @param ors_results Data frame from [calculate_ors()] with `logORS` and
#'   `type` columns.
#' @param title Character. Plot title. Default "Conservation by ortholog type".
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' ors_results <- calculate_ors(ccs_results)
#' plot_ors_by_type(ors_results)
#' }
#'
#' @export
plot_ors_by_type <- function(ors_results,
                             title = "Conservation by ortholog type") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install with install.packages('ggplot2').")
  }
  if (!"logORS" %in% colnames(ors_results)) {
    stop("ors_results must contain 'logORS' column from calculate_ors()")
  }
  if (!"type" %in% colnames(ors_results)) {
    stop("ors_results must contain 'type' column (use handle_multicopy_orthologs() first)")
  }

  ggplot2::ggplot(ors_results,
                  ggplot2::aes(x = .data$type, y = .data$logORS,
                               fill = .data$type)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    ggplot2::geom_hline(yintercept = c(1, 2), linetype = "dashed",
                        color = "red", alpha = 0.5) +
    ggplot2::labs(
      x = "Ortholog type",
      y = "logORS",
      title = title,
      fill = "Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}


#' Plot CCS vs logORS
#'
#' Scatter plot of CCS against logORS with conservation zones colored.
#'
#' @param ors_results Data frame from [calculate_ors()] with `CCS` and
#'   `logORS` columns.
#' @param title Character. Plot title. Default "CCS vs logORS".
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' ors_results <- calculate_ors(ccs_results)
#' plot_ccs_vs_ors(ors_results)
#' }
#'
#' @export
plot_ccs_vs_ors <- function(ors_results,
                            title = "CCS vs logORS") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install with install.packages('ggplot2').")
  }
  if (!all(c("CCS", "logORS") %in% colnames(ors_results))) {
    stop("ors_results must contain 'CCS' and 'logORS' columns")
  }

  ors_results$conservation <- cut(
    ors_results$logORS,
    breaks = c(-Inf, 0, 1, 2, Inf),
    labels = c("Diverged (<0)", "Moderate (0-1)",
               "Top 10% (1-2)", "Top 1% (>2)")
  )

  ggplot2::ggplot(ors_results,
                  ggplot2::aes(x = .data$CCS, y = .data$logORS,
                               color = .data$conservation)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_hline(yintercept = c(0, 1, 2, 3), linetype = "dashed",
                        alpha = 0.3) +
    ggplot2::scale_color_manual(
      values = c("Diverged (<0)" = "grey60",
                 "Moderate (0-1)" = "steelblue",
                 "Top 10% (1-2)" = "orange",
                 "Top 1% (>2)" = "red")
    ) +
    ggplot2::labs(
      x = "CCS (co-expression correlation score)",
      y = "logORS (ortholog rank score)",
      title = title,
      color = "Conservation"
    ) +
    ggplot2::theme_minimal()
}
