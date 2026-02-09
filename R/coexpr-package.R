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
#' @useDynLib coexpr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods new is setClass setGeneric setMethod setValidity callNextMethod show
#' @importFrom rlang .data
#' @importFrom stats cor median sd setNames
#' @importFrom utils object.size head
#' @importFrom future plan multisession
#' @importFrom furrr future_map future_map_dbl furrr_options
## usethis namespace: end
NULL
