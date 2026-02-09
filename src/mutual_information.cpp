// mutual_information.cpp
// Fast mutual information and CLR computation for co-expression analysis
//
// This implementation computes MI using empirical estimator and applies
// CLR (Context Likelihood Ratio) transformation for normalization.
// Parallelized via OpenMP for efficient computation of large gene matrices.
//
// Author: Martin Paliocha <martin.paliocha@nmbu.no>
// Based on: Faith et al. (2007) PLoS Biology 5:e8

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


//' Discretize a vector into equal-frequency bins
//'
//' @param x Numeric vector to discretize
//' @param n_bins Number of bins
//' @return Integer vector of bin assignments (0-indexed)
// [[Rcpp::export]]
arma::ivec discretize_equalfreq_cpp(const arma::vec& x, int n_bins) {
    const uword n = x.n_elem;
    arma::ivec result(n);

    // Get sorted indices
    arma::uvec sorted_idx = arma::sort_index(x);

    // Assign to bins based on rank
    for (uword i = 0; i < n; ++i) {
        int bin = static_cast<int>((static_cast<double>(i) * n_bins) / n);
        if (bin >= n_bins) bin = n_bins - 1;
        result(sorted_idx(i)) = bin;
    }

    return result;
}


//' Discretize expression matrix (all genes)
//'
//' @param expr Expression matrix (genes x samples)
//' @param n_bins Number of bins
//' @param n_cores Number of OpenMP threads (default 1)
//' @return Integer matrix of bin assignments (genes x samples)
// [[Rcpp::export]]
arma::imat discretize_matrix_cpp(const arma::mat& expr, int n_bins, int n_cores = 1) {
    const uword n_genes = expr.n_rows;
    const uword n_samples = expr.n_cols;

    arma::imat result(n_genes, n_samples);

#ifdef _OPENMP
    if (n_cores > 1) {
        omp_set_num_threads(n_cores);
    }
#endif

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (uword i = 0; i < n_genes; ++i) {
        arma::vec row = expr.row(i).t();

        // Check for zero variance
        double var = arma::var(row);
        if (var < 1e-10 || !std::isfinite(var)) {
            // All same bin for zero variance genes
            result.row(i).fill(0);
        } else {
            result.row(i) = discretize_equalfreq_cpp(row, n_bins).t();
        }
    }

    return result;
}


// Compute entropy of a discrete variable (internal helper)
inline double compute_entropy(const arma::ivec& x, int n_bins) {
    const uword n = x.n_elem;

    // Count occurrences
    std::vector<int> counts(n_bins, 0);
    for (uword i = 0; i < n; ++i) {
        counts[x(i)]++;
    }

    // Compute entropy: -sum(p * log(p))
    double entropy = 0.0;
    for (int i = 0; i < n_bins; ++i) {
        if (counts[i] > 0) {
            double p = static_cast<double>(counts[i]) / n;
            entropy -= p * std::log(p);
        }
    }

    return entropy;
}


// Compute joint entropy of two discrete variables (internal helper)
inline double compute_joint_entropy(const arma::ivec& x, const arma::ivec& y, int n_bins) {
    const uword n = x.n_elem;

    // Count joint occurrences using flat index
    int n_bins_sq = n_bins * n_bins;
    std::vector<int> joint_counts(n_bins_sq, 0);

    for (uword i = 0; i < n; ++i) {
        int idx = x(i) * n_bins + y(i);
        joint_counts[idx]++;
    }

    // Compute joint entropy
    double entropy = 0.0;
    for (int i = 0; i < n_bins_sq; ++i) {
        if (joint_counts[i] > 0) {
            double p = static_cast<double>(joint_counts[i]) / n;
            entropy -= p * std::log(p);
        }
    }

    return entropy;
}


// Compute mutual information between two discrete variables (internal helper)
// MI(X,Y) = H(X) + H(Y) - H(X,Y)
inline double compute_mi(const arma::ivec& x, const arma::ivec& y, int n_bins) {
    double h_x = compute_entropy(x, n_bins);
    double h_y = compute_entropy(y, n_bins);
    double h_xy = compute_joint_entropy(x, y, n_bins);

    double mi = h_x + h_y - h_xy;

    // MI should be non-negative (can be slightly negative due to numerical issues)
    if (mi < 0.0) mi = 0.0;

    return mi;
}


//' Compute mutual information matrix
//'
//' Computes pairwise MI for all gene pairs from discretized expression.
//'
//' @param expr_disc Discretized expression matrix (genes x samples, integer)
//' @param n_bins Number of bins used in discretization
//' @param n_cores Number of OpenMP threads
//' @return Symmetric MI matrix (genes x genes)
//'
//' @details
//' This function computes MI(i,j) = H(i) + H(j) - H(i,j) for all gene pairs.
//' The diagonal contains self-information (entropy).
//' Parallelized via OpenMP for efficiency.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_mi_matrix_cpp(const arma::imat& expr_disc, int n_bins, int n_cores = 1) {
    const uword n = expr_disc.n_rows;

    arma::mat mi_matrix(n, n, arma::fill::zeros);

    // Precompute entropies for all genes
    arma::vec entropies(n);

#ifdef _OPENMP
    if (n_cores > 1) {
        omp_set_num_threads(n_cores);
    }
#endif

    Rcpp::Rcout << "  Computing gene entropies..." << std::endl;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (uword i = 0; i < n; ++i) {
        arma::ivec row_i = arma::conv_to<arma::ivec>::from(expr_disc.row(i).t());
        entropies(i) = compute_entropy(row_i, n_bins);
    }

    // Set diagonal to entropy (self-information)
    mi_matrix.diag() = entropies;

    Rcpp::Rcout << "  Computing pairwise mutual information..." << std::endl;

    // Compute upper triangle of MI matrix
    const uword report_interval = std::max(static_cast<uword>(1), n / 20);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (uword i = 0; i < n; ++i) {
#ifdef _OPENMP
        if (omp_get_thread_num() == 0 && i % report_interval == 0 && i > 0) {
            Rcpp::Rcout << "  Progress: " << (100 * i / n) << "%" << std::endl;
        }
#else
        if (i % report_interval == 0 && i > 0) {
            Rcpp::Rcout << "  Progress: " << (100 * i / n) << "%" << std::endl;
        }
#endif

        arma::ivec row_i = arma::conv_to<arma::ivec>::from(expr_disc.row(i).t());

        for (uword j = i + 1; j < n; ++j) {
            arma::ivec row_j = arma::conv_to<arma::ivec>::from(expr_disc.row(j).t());

            // MI = H(i) + H(j) - H(i,j)
            double h_ij = compute_joint_entropy(row_i, row_j, n_bins);
            double mi = entropies(i) + entropies(j) - h_ij;

            if (mi < 0.0) mi = 0.0;

            mi_matrix(i, j) = mi;
            mi_matrix(j, i) = mi;
        }
    }

    return mi_matrix;
}


//' Apply CLR (Context Likelihood Ratio) transformation
//'
//' Normalizes MI matrix using z-scores to reduce hub gene bias.
//'
//' @param mi_matrix Symmetric mutual information matrix
//' @param n_cores Number of OpenMP threads (currently unused, for API consistency)
//' @return CLR-transformed matrix
//'
//' @details
//' CLR transformation (Faith et al. 2007):
//' \preformatted{
//' Z_row(i,j) = (MI(i,j) - mean(row_i)) / sd(row_i)
//' Z_col(i,j) = (MI(i,j) - mean(col_j)) / sd(col_j)
//' CLR(i,j) = sqrt(max(0, Z_row)^2 + max(0, Z_col)^2)
//' }
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat apply_clr_cpp(const arma::mat& mi_matrix, int n_cores = 1) {
    const uword n = mi_matrix.n_rows;

    // Compute row means and standard deviations
    arma::vec row_means = arma::mean(mi_matrix, 1);
    arma::vec row_sds(n);

    for (uword i = 0; i < n; ++i) {
        row_sds(i) = arma::stddev(mi_matrix.row(i));
        if (row_sds(i) < 1e-10) row_sds(i) = 1.0;  // Avoid division by zero
    }

    arma::mat clr_matrix(n, n, arma::fill::zeros);

#ifdef _OPENMP
    if (n_cores > 1) {
        omp_set_num_threads(n_cores);
    }
#endif

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(n_cores > 1)
#endif
    for (uword i = 0; i < n; ++i) {
        for (uword j = 0; j < n; ++j) {
            // Z-scores
            double z_row = (mi_matrix(i, j) - row_means(i)) / row_sds(i);
            double z_col = (mi_matrix(i, j) - row_means(j)) / row_sds(j);

            // CLR: only positive z-scores contribute
            double z_row_pos = (z_row > 0.0) ? z_row : 0.0;
            double z_col_pos = (z_col > 0.0) ? z_col : 0.0;

            clr_matrix(i, j) = std::sqrt(z_row_pos * z_row_pos + z_col_pos * z_col_pos);
        }
    }

    return clr_matrix;
}


//' Compute MI+CLR similarity matrix (triangular output)
//'
//' Full pipeline: discretize -> MI -> CLR, returning triangular format.
//'
//' @param expr Expression matrix (genes x samples)
//' @param n_bins Number of bins for discretization
//' @param n_cores Number of OpenMP threads
//' @return List with triangular similarity data
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List compute_mi_clr_tri_cpp(const arma::mat& expr, int n_bins, int n_cores = 1) {
    const uword n = expr.n_rows;

    Rcpp::Rcout << "Discretizing expression values..." << std::endl;
    arma::imat expr_disc = discretize_matrix_cpp(expr, n_bins, n_cores);

    Rcpp::Rcout << "Computing mutual information matrix..." << std::endl;
    arma::mat mi_matrix = compute_mi_matrix_cpp(expr_disc, n_bins, n_cores);

    Rcpp::Rcout << "Applying CLR transformation..." << std::endl;
    arma::mat clr_matrix = apply_clr_cpp(mi_matrix, n_cores);

    // Extract upper triangle (excluding diagonal) in column-major order
    uword tri_size = n * (n - 1) / 2;
    arma::vec upper_tri(tri_size);

    uword idx = 0;
    for (uword j = 1; j < n; ++j) {
        for (uword i = 0; i < j; ++i) {
            upper_tri(idx++) = clr_matrix(i, j);
        }
    }

    // Diagonal value (CLR of self = 0 typically, but we'll use the actual value)
    double diag_val = clr_matrix(0, 0);

    return Rcpp::List::create(
        Rcpp::Named("data") = upper_tri,
        Rcpp::Named("n") = n,
        Rcpp::Named("diag_value") = diag_val
    );
}


//' Compute MI+CLR similarity matrix (full matrix output)
//'
//' Full pipeline: discretize -> MI -> CLR, returning full matrix.
//'
//' @param expr Expression matrix (genes x samples)
//' @param n_bins Number of bins for discretization
//' @param n_cores Number of OpenMP threads
//' @return CLR-transformed similarity matrix
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_mi_clr_cpp(const arma::mat& expr, int n_bins, int n_cores = 1) {
    Rcpp::Rcout << "Discretizing expression values..." << std::endl;
    arma::imat expr_disc = discretize_matrix_cpp(expr, n_bins, n_cores);

    Rcpp::Rcout << "Computing mutual information matrix..." << std::endl;
    arma::mat mi_matrix = compute_mi_matrix_cpp(expr_disc, n_bins, n_cores);

    Rcpp::Rcout << "Applying CLR transformation..." << std::endl;
    arma::mat clr_matrix = apply_clr_cpp(mi_matrix, n_cores);

    return clr_matrix;
}
