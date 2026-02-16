// bootstrap_ccs.cpp
// Fast bootstrapped CCS scoring for homeolog selection confidence
//
// For each multi-copy group, this function:
// 1. Subsamples the reference set (80% without replacement)
// 2. Computes column-wise Pearson correlation between paired columns
// 3. Records which candidate wins and the selected candidate's score
//
// Author: Martin Paliocha <martin.paliocha@nmbu.no>

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <random>

using namespace Rcpp;
using namespace arma;

// Compute Pearson correlation between two vectors, handling NAs
static double pearson_cor(const vec& x, const vec& y) {
    uword n = x.n_elem;
    if (n < 3) return NA_REAL;

    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    uword count = 0;

    for (uword i = 0; i < n; ++i) {
        if (std::isfinite(x(i)) && std::isfinite(y(i))) {
            sx += x(i);
            sy += y(i);
            sxx += x(i) * x(i);
            syy += y(i) * y(i);
            sxy += x(i) * y(i);
            ++count;
        }
    }

    if (count < 3) return NA_REAL;

    double mx = sx / static_cast<double>(count);
    double my = sy / static_cast<double>(count);
    double vx = sxx / static_cast<double>(count) - mx * mx;
    double vy = syy / static_cast<double>(count) - my * my;

    if (vx <= 0 || vy <= 0) return NA_REAL;

    return (sxy / static_cast<double>(count) - mx * my) /
           (std::sqrt(vx) * std::sqrt(vy));
}


//' Bootstrap CCS scoring for a single multi-copy group
//'
//' For one group of candidates, runs n_bootstrap iterations. Each iteration
//' subsamples 80% of reference rows, computes CCS for each candidate, and
//' records which candidate wins.
//'
//' @param ref_sp1 Matrix of co-expression vectors for sp1 (n_ref x n_candidates)
//'   Each column is one candidate's co-expression against the reference set.
//' @param ref_sp2 Matrix of co-expression vectors for sp2 (n_ref x n_candidates)
//' @param selected_idx 0-based index of the candidate selected by full scoring
//' @param n_bootstrap Number of bootstrap iterations
//' @param seed Random seed for reproducibility
//' @return List with:
//'   - n_same: number of bootstraps selecting the same candidate
//'   - scores: vector of selected candidate's score across bootstraps
//'
//' @keywords internal
// [[Rcpp::export]]
List bootstrap_group_scores_cpp(const arma::mat& ref_sp1,
                                const arma::mat& ref_sp2,
                                int selected_idx,
                                int n_bootstrap,
                                unsigned int seed) {
    uword n_ref = ref_sp1.n_rows;
    uword n_cand = ref_sp1.n_cols;
    uword subsample_size = std::max(static_cast<uword>(3),
                                    static_cast<uword>(n_ref * 0.8));

    std::mt19937 rng(seed);
    std::vector<uword> indices(n_ref);
    std::iota(indices.begin(), indices.end(), 0);

    int n_same = 0;
    NumericVector scores(n_bootstrap);

    for (int b = 0; b < n_bootstrap; ++b) {
        // Subsample without replacement
        std::shuffle(indices.begin(), indices.end(), rng);

        // Extract subsampled rows
        uvec sub_idx(subsample_size);
        for (uword k = 0; k < subsample_size; ++k) {
            sub_idx(k) = indices[k];
        }

        mat sub_sp1 = ref_sp1.rows(sub_idx);
        mat sub_sp2 = ref_sp2.rows(sub_idx);

        // Compute CCS for each candidate
        double best_score = -2.0;
        int best_cand = -1;

        for (uword c = 0; c < n_cand; ++c) {
            double score = pearson_cor(sub_sp1.col(c), sub_sp2.col(c));
            if (std::isfinite(score) && score > best_score) {
                best_score = score;
                best_cand = static_cast<int>(c);
            }
        }

        if (best_cand == selected_idx) {
            ++n_same;
        }

        // Record score for the selected candidate
        scores[b] = pearson_cor(
            sub_sp1.col(static_cast<uword>(selected_idx)),
            sub_sp2.col(static_cast<uword>(selected_idx))
        );
    }

    return List::create(
        Named("n_same") = n_same,
        Named("scores") = scores
    );
}
