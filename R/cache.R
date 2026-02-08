#' @importFrom utils object.size
NULL

# Internal: Compute cache key for a similarity matrix computation
#
# Hashes expression matrix content + method parameters to create a
# deterministic filename for caching.
#
# @param expr_matrix Numeric matrix (genes x samples)
# @param method_label Character. Method identifier (e.g., "pcc_mr", "mi_clr")
# @param ... Additional method parameters to include in the key
#
# @return Character string suitable as a filename (without directory)
# @keywords internal
# @noRd
compute_cache_key <- function(expr_matrix, method_label, ...) {
  params <- list(...)


  # Build content to hash: matrix data + dimensions + names + params
  key_content <- list(
    data = as.vector(expr_matrix),
    dim = dim(expr_matrix),
    rownames = rownames(expr_matrix),
    colnames = colnames(expr_matrix),
    method_label = method_label,
    params = params
  )

  hash <- digest::digest(key_content, algo = "xxhash64")


  # Build human-readable prefix from params

  param_str <- paste(vapply(params, as.character, character(1)), collapse = "_")
  if (nzchar(param_str)) {
    sprintf("%s_%s_%s.rds", method_label, param_str, hash)
  } else {
    sprintf("%s_%s.rds", method_label, hash)
  }
}


# Internal: Load a cached similarity matrix
#
# @param cache_dir Character. Directory containing cache files
# @param cache_key Character. Filename to load
#
# @return A TriSimilarity object, or NULL if not found or corrupt
# @keywords internal
# @noRd
cache_load <- function(cache_dir, cache_key) {
  path <- file.path(cache_dir, cache_key)

  if (!file.exists(path)) {
    return(NULL)
  }

  obj <- tryCatch(
    readRDS(path),
    error = function(e) {
      warning(sprintf("Corrupt cache file '%s', will recompute: %s",
                      cache_key, conditionMessage(e)), call. = FALSE)
      return(NULL)
    }
  )

  if (is.null(obj)) return(NULL)

  if (!methods::is(obj, "TriSimilarity")) {
    warning(sprintf("Cache file '%s' does not contain a TriSimilarity object, will recompute",
                    cache_key), call. = FALSE)
    return(NULL)
  }

  obj
}


# Internal: Save a similarity matrix to cache
#
# Converts to TriSimilarity if needed before saving.
#
# @param obj A TriSimilarity object or matrix to cache
# @param cache_dir Character. Directory to write to (created if needed)
# @param cache_key Character. Filename to save as
#
# @return Invisibly returns the file path
# @keywords internal
# @noRd
cache_save <- function(obj, cache_dir, cache_key) {
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Always store as TriSimilarity for compact storage
  if (is.matrix(obj)) {
    obj <- as.TriSimilarity(obj)
  }

  path <- file.path(cache_dir, cache_key)
  saveRDS(obj, path)
  invisible(path)
}


#' List cached similarity matrices
#'
#' Lists all cached similarity matrix files in a cache directory, with
#' file size and modification time.
#'
#' @param cache_dir Character. Path to the cache directory.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{file}{Character. Cache filename.}
#'     \item{size_mb}{Numeric. File size in megabytes.}
#'     \item{modified}{POSIXct. Last modification time.}
#'     \item{method}{Character. Method label extracted from filename.}
#'   }
#'   Returns a zero-row data frame if the directory doesn't exist or is empty.
#'
#' @examples
#' \dontrun{
#' cache_list("my_cache")
#' }
#'
#' @export
cache_list <- function(cache_dir) {
  empty <- data.frame(
    file = character(0),
    size_mb = numeric(0),
    modified = as.POSIXct(character(0)),
    method = character(0),
    stringsAsFactors = FALSE
  )

  if (!dir.exists(cache_dir)) {
    return(empty)
  }

  files <- list.files(cache_dir, pattern = "\\.rds$", full.names = FALSE)

  if (length(files) == 0) {
    return(empty)
  }

  info <- file.info(file.path(cache_dir, files))

  # Extract method from filename (first part before _ or .)
  methods <- sub("^([^_]+(?:_[^_]+)?)_.*\\.rds$", "\\1", files)

  data.frame(
    file = files,
    size_mb = round(info$size / 1e6, 2),
    modified = info$mtime,
    method = methods,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' Clear cached similarity matrices
#'
#' Removes cached similarity matrix files from a cache directory.
#'
#' @param cache_dir Character. Path to the cache directory.
#' @param method Character or NULL. If specified, only remove files matching
#'   this method (e.g., "pcc_mr", "mi_clr"). If NULL (default), removes all
#'   cache files.
#'
#' @return Invisibly returns the number of files removed.
#'
#' @examples
#' \dontrun{
#' # Remove all cached files
#' cache_clear("my_cache")
#'
#' # Remove only PCC+MR caches
#' cache_clear("my_cache", method = "pcc_mr")
#' }
#'
#' @export
cache_clear <- function(cache_dir, method = NULL) {
  if (!dir.exists(cache_dir)) {
    return(invisible(0L))
  }

  if (is.null(method)) {
    files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  } else {
    pattern <- sprintf("^%s_.*\\.rds$", method)
    files <- list.files(cache_dir, pattern = pattern, full.names = TRUE)
  }

  if (length(files) == 0) {
    return(invisible(0L))
  }

  removed <- file.remove(files)
  invisible(sum(removed))
}
