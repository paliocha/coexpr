// mutual_rank.cpp
// Memory-efficient streaming mutual rank transformation for co-expression analysis
//
// This implementation computes the mutual rank normalization in a streaming fashion,
// avoiding the need to store multiple n x n intermediate matrices in memory.
// The algorithm is O(n^3) in computation but only O(n^2) in memory (for input + output).
//
// Author: Martin Paliocha <martin.paliocha@nmbu.no>
// Based on: Obayashi & Kinoshita (2009) DNA Research 16:249-260

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//' Compute average ranks for a vector (descending order)
//'
//' For each element, computes its rank when sorted in descending order.
//' Ties are handled using average ranks.
//'
//' @param x Input vector
//' @return Vector of ranks (1-based, average method for ties)
// [[Rcpp::export]]
arma::vec compute_ranks_desc_cpp(const arma::vec& x) {
    const uword n = x.n_elem;
    arma::vec ranks(n);

    // Create index vector for sorting
    std::vector<uword> indices(n);
    for (uword i = 0; i < n; ++i) {
        indices[i] = i;
    }

    // Sort indices by x values in descending order
    std::sort(indices.begin(), indices.end(),
              [&x](uword a, uword b) { return x(a) > x(b); });

    // Assign ranks with tie handling (average method)
    uword i = 0;
    while (i < n) {
        uword j = i;

        // Find all elements with the same value (ties)
        while (j < n - 1 && x(indices[j]) == x(indices[j + 1])) {
            ++j;
        }

        // Average rank for tied elements
        // Ranks are 1-based: position i gets rank i+1
        double avg_rank = (static_cast<double>(i + 1) + static_cast<double>(j + 1)) / 2.0;

        // Assign the average rank to all tied elements
        for (uword k = i; k <= j; ++k) {
            ranks(indices[k]) = avg_rank;
        }

        i = j + 1;
    }

    return ranks;
}


//' Streaming mutual rank transformation (memory-efficient, parallelized)
//'
//' Transforms a Pearson correlation matrix into a PCC+MR similarity matrix
//' using the mutual rank normalization from Obayashi & Kinoshita (2009).
//'
//' This implementation uses a streaming approach that computes ranks on-demand
//' for each row pair, avoiding storage of intermediate rank matrices.
//' The algorithm is fully parallelizable via OpenMP.
//'
//' Memory usage: O(n^2) for input + O(n^2) for output + O(n) per thread
//' Time complexity: O(n^3) due to on-demand rank computation
//'
//' @param sim_pcc Symmetric Pearson correlation matrix (n x n)
//' @param n_cores Number of OpenMP threads to use (default: 1)
//' @return Mutual rank normalized similarity matrix (n x n) with values in range 0 to 1
//'
//' @details
//' The transformation is:
//' \deqn{S^{PCC+MR}_{ij} = 1 - \frac{\log(\sqrt{R_{ij} \cdot R_{ji}})}{\log(n)}}
//'
//' Where R_ij is the rank of gene j in gene i's correlation row (descending order).
//'
//' Since the input PCC matrix is symmetric, we exploit that colRanks(i,j) = rowRanks(j,i).
//' This means we only need to compute row ranks for each row.
//'
//' @export
// [[Rcpp::export]]
arma::mat mutual_rank_transform_cpp(const arma::mat& sim_pcc, int n_cores = 1) {
    const uword n = sim_pcc.n_rows;

    if (n != sim_pcc.n_cols) {
        stop("sim_pcc must be a square matrix");
    }

    if (n < 2) {
        stop("Matrix must have at least 2 rows/columns");
    }

    const double log_n = std::log(static_cast<double>(n));

    // Allocate result matrix
    arma::mat result(n, n, arma::fill::zeros);

    // Set number of threads
#ifdef _OPENMP
    if (n_cores > 1) {
        omp_set_num_threads(n_cores);
    }
#endif

    // Progress tracking (only from thread 0)
    const uword report_interval = std::max(static_cast<uword>(1), n / 20);

    // Parallel loop over rows
    // We process the upper triangle (including diagonal) and mirror to lower
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (uword i = 0; i < n; ++i) {

        // Progress reporting (only from thread 0, non-parallel section)
#ifdef _OPENMP
        if (omp_get_thread_num() == 0 && i % report_interval == 0 && i > 0) {
            Rcpp::Rcout << "  Progress: " << (100 * i / n) << "%" << std::endl;
        }
#else
        if (i % report_interval == 0 && i > 0) {
            Rcpp::Rcout << "  Progress: " << (100 * i / n) << "%" << std::endl;
        }
#endif

        // Compute ranks for row i (descending order of correlations)
        arma::vec row_i = sim_pcc.row(i).t();
        arma::vec ranks_i = compute_ranks_desc_cpp(row_i);

        // Process columns j >= i (upper triangle)
        for (uword j = i; j < n; ++j) {
            if (i == j) {
                // Diagonal: self-similarity is 1
                result(i, j) = 1.0;
            } else {
                // Compute ranks for row j
                // Due to symmetry of PCC: rank of i in column j = rank of i in row j
                arma::vec row_j = sim_pcc.row(j).t();
                arma::vec ranks_j = compute_ranks_desc_cpp(row_j);

                // R_ij = rank of j in row i
                // R_ji = rank of i in row j (= rank of i in column j due to symmetry)
                double R_ij = ranks_i(j);
                double R_ji = ranks_j(i);

                // Mutual rank transformation
                double mutual_rank = std::sqrt(R_ij * R_ji);
                double val = 1.0 - std::log(mutual_rank) / log_n;

                // Ensure value is in valid range [0, 1]
                if (val < 0.0) val = 0.0;
                if (val > 1.0) val = 1.0;

                // Store in both upper and lower triangle (symmetric)
                result(i, j) = val;
                result(j, i) = val;
            }
        }
    }

    return result;
}


//' Cached mutual rank transformation (speed-optimized, parallelized)
//'
//' Transforms a Pearson correlation matrix into a PCC+MR similarity matrix
//' using the mutual rank normalization from Obayashi & Kinoshita (2009).
//'
//' This implementation precomputes all row ranks once, then uses them for
//' the mutual rank computation. This is faster than the streaming approach
//' but uses more memory (O(n^2) for the rank matrix).
//'
//' Memory usage: O(n^2) for input + O(n^2) for ranks + O(n^2) for output = 3 matrices
//' Time complexity: O(n^2 log n) for ranking + O(n^2) for transformation
//'
//' @param sim_pcc Symmetric Pearson correlation matrix (n x n)
//' @param n_cores Number of OpenMP threads to use (default: 1)
//' @return Mutual rank normalized similarity matrix (n x n) with values in range 0 to 1
//'
//' @details
//' Since the input PCC matrix is symmetric, column ranks equal transposed row ranks.
//' This allows us to precompute row ranks once and reuse them efficiently.
//'
//' @export
// [[Rcpp::export]]
arma::mat mutual_rank_transform_cached_cpp(const arma::mat& sim_pcc, int n_cores = 1) {
    const uword n = sim_pcc.n_rows;

    if (n != sim_pcc.n_cols) {
        stop("sim_pcc must be a square matrix");
    }

    if (n < 2) {
        stop("Matrix must have at least 2 rows/columns");
    }

    const double log_n = std::log(static_cast<double>(n));

    // Allocate rank matrix and result matrix
    arma::mat row_ranks(n, n);
    arma::mat result(n, n, arma::fill::zeros);

    // Set number of threads
#ifdef _OPENMP
    if (n_cores > 1) {
        omp_set_num_threads(n_cores);
    }
#endif

    // Step 1: Precompute all row ranks in parallel
    Rcpp::Rcout << "  Computing row ranks..." << std::endl;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (uword i = 0; i < n; ++i) {
        arma::vec row_i = sim_pcc.row(i).t();
        row_ranks.row(i) = compute_ranks_desc_cpp(row_i).t();
    }

    // Step 2: Compute mutual rank and transform in parallel
    Rcpp::Rcout << "  Computing mutual ranks..." << std::endl;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (uword i = 0; i < n; ++i) {
        for (uword j = i; j < n; ++j) {
            if (i == j) {
                result(i, j) = 1.0;
            } else {
                // R_ij = rank of j in row i
                // R_ji = rank of i in row j (= col rank due to symmetry)
                double R_ij = row_ranks(i, j);
                double R_ji = row_ranks(j, i);

                double mutual_rank = std::sqrt(R_ij * R_ji);
                double val = 1.0 - std::log(mutual_rank) / log_n;

                // Clamp to [0, 1]
                if (val < 0.0) val = 0.0;
                if (val > 1.0) val = 1.0;

                result(i, j) = val;
                result(j, i) = val;
            }
        }
    }

    return result;
}


//' Check if OpenMP is available
//'
//' @return TRUE if OpenMP support is compiled in, FALSE otherwise
//' @export
// [[Rcpp::export]]
bool has_openmp() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}


//' Get maximum number of OpenMP threads
//'
//' @return Maximum number of threads available, or 1 if OpenMP not available
//' @export
// [[Rcpp::export]]
int get_max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}
