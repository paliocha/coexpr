#' @include RcppExports.R
NULL

#' TriSimilarity Class
#'
#' S4 class for efficient storage of symmetric similarity matrices.
#' Stores only the upper triangle (excluding diagonal), reducing memory
#' by approximately 50%.
#'
#' @slot data Numeric vector of upper triangle values in column-major order.
#' @slot genes Character vector of gene names.
#' @slot n Integer. Number of genes.
#' @slot diag_value Numeric. Value for diagonal elements (typically 1.0).
#'
#' @details
#' The upper triangle is stored in column-major order, matching R's
#' \code{upper.tri()} ordering. For a matrix with n genes, the data
#' slot contains n*(n-1)/2 elements.
#'
#' For position (i,j) with i < j, the index in the data vector is:
#' \code{(j-1)*(j-2)/2 + i}
#'
#' @seealso \code{\link{TriSimilarity}} for the constructor function.
#'
#' @exportClass TriSimilarity
#' @name TriSimilarity-class
#' @rdname TriSimilarity-class
setClass("TriSimilarity",
  slots = c(
    data = "numeric",
    genes = "character",
    n = "integer",
    diag_value = "numeric"
  ),
  prototype = list(
    data = numeric(0),
    genes = character(0),
    n = 0L,
    diag_value = 1.0
  )
)


# Validity check for TriSimilarity
# (internal, no roxygen docs needed)
setValidity("TriSimilarity", function(object) {
  errors <- character()

  # Check data length
  expected_len <- object@n * (object@n - 1L) / 2L
  if (length(object@data) != expected_len) {
    errors <- c(errors, sprintf(
      "data length (%d) doesn't match expected (%d) for %d genes",
      length(object@data), expected_len, object@n
    ))
  }

  # Check genes length
  if (length(object@genes) != object@n) {
    errors <- c(errors, sprintf(
      "genes length (%d) doesn't match n (%d)",
      length(object@genes), object@n
    ))
  }

  # Check diag_value length
  if (length(object@diag_value) != 1) {
    errors <- c(errors, "diag_value must be length 1")
  }

  if (length(errors) == 0) TRUE else errors
})


#' Create a TriSimilarity object
#'
#' Constructor for triangular similarity matrices.
#'
#' @param data Numeric vector of upper triangle values (column-major order,
#'   excluding diagonal). Length should be n*(n-1)/2 where n is number of genes.
#' @param genes Character vector of gene names.
#' @param diag_value Numeric. Value for diagonal elements (default 1.0).
#'
#' @return A TriSimilarity object.
#'
#' @examples
#' \dontrun{
#' # Create from components
#' tri <- TriSimilarity(
#'   data = c(0.5, 0.3, 0.7),  # Upper triangle of 3x3
#'   genes = c("A", "B", "C"),
#'   diag_value = 1.0
#' )
#' }
#'
#' @export
TriSimilarity <- function(data, genes, diag_value = 1.0) {
  n <- length(genes)

  new("TriSimilarity",
    data = as.numeric(data),
    genes = as.character(genes),
    n = as.integer(n),
    diag_value = as.numeric(diag_value)
  )
}


#' Convert a full matrix to TriSimilarity
#'
#' Extracts the upper triangle of a symmetric matrix and stores it efficiently.
#'
#' @param mat Numeric matrix. Should be symmetric with gene names as
#'   row/column names.
#'
#' @return A TriSimilarity object.
#'
#' @examples
#' \dontrun{
#' # Create symmetric matrix
#' mat <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.7, 0.3, 0.7, 1), nrow = 3)
#' rownames(mat) <- colnames(mat) <- c("A", "B", "C")
#'
#' # Convert to triangular
#' tri <- as.TriSimilarity(mat)
#' }
#'
#' @export
as.TriSimilarity <- function(mat) {
  if (!is.matrix(mat)) {
    stop("mat must be a matrix")
  }

  if (nrow(mat) != ncol(mat)) {
    stop("mat must be square")
  }

  n <- nrow(mat)
  genes <- rownames(mat)

  if (is.null(genes)) {
    genes <- paste0("Gene", seq_len(n))
  }

  # Assume constant diagonal
  diag_val <- mat[1, 1]

  # Extract upper triangle (excluding diagonal) in column-major order
  upper_idx <- which(upper.tri(mat))
  data <- mat[upper_idx]

  TriSimilarity(data, genes, diag_val)
}


#' @describeIn TriSimilarity-class Convert to full matrix
#' @param x TriSimilarity object
#' @param ... Ignored (for S4 method consistency)
#' @export
setMethod("as.matrix", "TriSimilarity", function(x, ...) {
  n <- x@n
  mat <- matrix(x@diag_value, nrow = n, ncol = n)
  rownames(mat) <- colnames(mat) <- x@genes

  # Fill upper triangle
  upper_idx <- which(upper.tri(mat))
  mat[upper_idx] <- x@data

  # Mirror to lower triangle
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

  mat
})


#' @describeIn TriSimilarity-class Get dimensions
#' @export
setMethod("dim", "TriSimilarity", function(x) {
  c(x@n, x@n)
})


#' @describeIn TriSimilarity-class Get dimension names
#' @export
setMethod("dimnames", "TriSimilarity", function(x) {
  list(x@genes, x@genes)
})


#' @describeIn TriSimilarity-class Show method
#' @param object TriSimilarity object
#' @export
setMethod("show", "TriSimilarity", function(object) {
  mem_used <- object.size(object@data)
  mem_full <- 8 * object@n * object@n

  cat("TriSimilarity: ", object@n, " x ", object@n, " genes\n", sep = "")
  cat(sprintf("Storage: %.2f MB (%.1f%% of full matrix)\n",
              as.numeric(mem_used) / 1e6,
              100 * as.numeric(mem_used) / mem_full))
  cat("Diagonal value:", object@diag_value, "\n")

  if (object@n <= 5) {
    cat("\nFull matrix:\n")
    print(as.matrix(object))
  } else {
    cat("First 5 genes:", paste(head(object@genes, 5), collapse = ", "), "\n")
  }

  invisible(object)
})


#' Subset a TriSimilarity object
#'
#' @param x TriSimilarity object
#' @param i Row indices or gene names
#' @param j Column indices or gene names
#' @param ... Ignored
#' @param drop Logical. If TRUE and result is single value, return scalar.
#'
#' @return Numeric value(s) or matrix subset.
#'
#' @export
setMethod("[", c("TriSimilarity", "ANY", "ANY"),
  function(x, i, j, ..., drop = TRUE) {
    # Handle missing j (row extraction: x[i, ])
    if (missing(j)) {
      j <- seq_len(x@n)
    }
    # Handle missing i (column extraction: x[, j])
    if (missing(i)) {
      i <- seq_len(x@n)
    }
    # Convert character to indices
    if (is.character(i)) i <- match(i, x@genes)
    if (is.character(j)) j <- match(j, x@genes)

    if (any(is.na(i)) || any(is.na(j))) {
      stop("Gene name not found in similarity matrix")
    }

    n <- x@n

    # Single element access
    if (length(i) == 1 && length(j) == 1) {
      if (i == j) {
        return(x@diag_value)
      }

      # Ensure ii < jj for upper triangle
      ii <- min(i, j)
      jj <- max(i, j)

      # Index in upper triangle vector (column-major)
      idx <- (jj - 1) * (jj - 2) / 2 + ii
      return(x@data[idx])
    }

    # Vectorized multi-element access
    rows_exp <- rep(i, times = length(j))
    cols_exp <- rep(j, each = length(i))

    is_diag <- rows_exp == cols_exp
    ii <- pmin(rows_exp, cols_exp)
    jj <- pmax(rows_exp, cols_exp)

    idx <- (jj - 1L) * (jj - 2L) / 2L + ii
    values <- x@data[idx]
    values[is_diag] <- x@diag_value

    result <- matrix(values, nrow = length(i), ncol = length(j))
    rownames(result) <- x@genes[i]
    colnames(result) <- x@genes[j]
    result
  }
)

#' @rdname sub-TriSimilarity-ANY-ANY-method
#' @export
setMethod("[", c("TriSimilarity", "ANY", "missing"),
  function(x, i, j, ..., drop = TRUE) {
    # Extract row: x[i, ]
    j <- seq_len(x@n)
    callNextMethod(x, i, j, ..., drop = drop)
  }
)

#' @rdname sub-TriSimilarity-ANY-ANY-method
#' @export
setMethod("[", c("TriSimilarity", "missing", "ANY"),
  function(x, i, j, ..., drop = TRUE) {
    # Extract column: x[, j]
    i <- seq_len(x@n)
    callNextMethod(x, i, j, ..., drop = drop)
  }
)


#' Extract a column from TriSimilarity
#'
#' Efficiently extracts a full column (co-expression vector) from a
#' triangular similarity matrix.
#'
#' @param x TriSimilarity object
#' @param gene Gene name or index to extract column for
#'
#' @return Named numeric vector
#'
#' @export
setGeneric("extractColumn", function(x, gene) standardGeneric("extractColumn"))

#' @rdname extractColumn
#' @export
setMethod("extractColumn", "TriSimilarity", function(x, gene) {
  # Convert character to index
  if (is.character(gene)) {
    j <- match(gene, x@genes)
    if (is.na(j)) {
      stop(sprintf("Gene '%s' not found in similarity matrix", gene))
    }
  } else {
    j <- gene
  }

  n <- x@n
  col <- numeric(n)

  # Diagonal
  col[j] <- x@diag_value

  # Values where row < j (in upper triangle): index = (j-1)*(j-2)/2 + row
  if (j > 1L) {
    rows_above <- seq_len(j - 1L)
    col[rows_above] <- x@data[(j - 1L) * (j - 2L) / 2L + rows_above]
  }

  # Values where row > j (mirror from upper triangle): index = (row-1)*(row-2)/2 + j
  if (j < n) {
    rows_below <- seq.int(j + 1L, n)
    col[rows_below] <- x@data[(rows_below - 1L) * (rows_below - 2L) / 2L + j]
  }

  names(col) <- x@genes
  col
})


#' Extract multiple columns from TriSimilarity
#'
#' @param x TriSimilarity object
#' @param genes Gene names or indices to extract
#'
#' @return Numeric matrix with rows = all genes, columns = requested genes
#'
#' @export
setGeneric("extractColumns", function(x, genes) standardGeneric("extractColumns"))

#' @rdname extractColumns
#' @export
setMethod("extractColumns", "TriSimilarity", function(x, genes) {
  if (is.character(genes)) {
    indices <- match(genes, x@genes)
    if (any(is.na(indices))) {
      missing <- genes[is.na(indices)]
      stop(sprintf("Genes not found: %s", paste(missing, collapse = ", ")))
    }
  } else {
    indices <- genes
  }

  result <- matrix(NA_real_, nrow = x@n, ncol = length(indices))

  for (ci in seq_along(indices)) {
    result[, ci] <- extractColumn(x, indices[ci])
  }

  rownames(result) <- x@genes
  colnames(result) <- x@genes[indices]
  result
})


#' Extract rows from TriSimilarity
#'
#' @param x TriSimilarity object
#' @param genes Gene names or indices to extract
#'
#' @return Numeric matrix with rows = requested genes, columns = all genes
#'
#' @export
setGeneric("extractRows", function(x, genes) standardGeneric("extractRows"))

#' @rdname extractRows
#' @export
setMethod("extractRows", "TriSimilarity", function(x, genes) {
  t(extractColumns(x, genes))
})


#' Check if object is TriSimilarity
#'
#' @param x Object to check
#'
#' @return Logical
#'
#' @export
is.TriSimilarity <- function(x) {
  inherits(x, "TriSimilarity")
}


#' Save TriSimilarity to file
#'
#' @param x TriSimilarity object
#' @param file File path (should end in .rds)
#'
#' @return Invisibly returns file path
#'
#' @export
saveTriSimilarity <- function(x, file) {
  if (!is(x, "TriSimilarity")) {
    stop("x must be a TriSimilarity object")
  }
  saveRDS(x, file)
  invisible(file)
}


#' Load TriSimilarity from file
#'
#' @param file File path
#'
#' @return TriSimilarity object
#'
#' @export
loadTriSimilarity <- function(file) {
  x <- readRDS(file)
  if (!is(x, "TriSimilarity")) {
    stop("File does not contain a TriSimilarity object")
  }
  x
}
