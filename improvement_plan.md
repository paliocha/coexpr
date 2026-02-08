# coexpr Improvement Plan

## Package Overview

The package is well-architected (~6,000 LOC across R, C++, and tests) implementing cross-species co-expression analysis. The core pipeline — similarity matrices (PCC+MR, MI+CLR) → CCS → ORS — is complete and production-ready, with good documentation, C++ performance-critical code, and an efficient `TriSimilarity` S4 class for memory optimization.

---

## Main Areas for Improvement

### 1. Test Coverage Gaps
- **`multicopy.R` is barely tested** — only 42 LOC covering 2 of 5 strategies (`strict` and `all_pairs`). The `best_hit`, `mean`, and `max` strategies have zero test coverage.
- No integration tests exercising the full pipeline end-to-end.
- No tests with real-world-scale data or performance regression tests.

### 2. Missing Multiple Testing Correction
- `test_ors_significance()` computes raw p-values (`1 - ORS`) but offers no FDR/Benjamini-Hochberg adjustment. For genome-scale analyses with thousands of ortholog pairs, this is a significant omission.

### 3. TriSimilarity Subsetting Performance
- The `[` operator and `extractRows()` use R-level loops for element access (`TriSimilarity-class.R` lines ~263-278, ~389-390). For large extractions this will be slow — these should be vectorized or moved to C++.

### 4. Multicopy Aggregation Bug
- `aggregate_by_mean()` in `multicopy.R` appears to always keep the first `gene_sp2` value after summarizing CCS, which is nonsensical for N:M relationships where there are multiple sp2 genes.

### 5. No Alternative Correlation Methods
- Only Pearson is available for the similarity step. Spearman or Kendall rank correlations (useful for non-normal expression data) are not offered, despite being straightforward to add.

### 6. No Caching / Memoization
- Similarity matrices (the most expensive computation) are recalculated each time. There's no built-in caching or checkpointing mechanism for iterative workflows.

### 7. Documentation Gaps
- The vignette uses only the small bundled test data — no guidance on preprocessing real expression matrices (normalization, log-transform, filtering low-expressed genes).
- Limited guidance on MI+CLR parameter tuning (number of bins, sample size requirements).

### 8. Visualization Helpers
- All plotting is left entirely to the user. Given the well-defined output types (ORS distributions, CCS heatmaps, conservation boxplots), lightweight plot functions would improve usability.

### 9. R CMD check Artifacts
- `coexpr.Rcheck/` and `coexpr_0.1.0.tar.gz` are sitting in the parent directory (untracked). The `.Rbuildignore` or `.gitignore` should handle these.

---

## Priority Ranking

| Priority | Area | Impact |
|----------|------|--------|
| **High** | Complete multicopy tests | Correctness risk — untested code paths |
| **High** | Fix `aggregate_by_mean()` bug | Incorrect results for N:M orthologs |
| **High** | Add multiple testing correction | Statistical validity for large-scale use |
| **Medium** | Vectorize TriSimilarity subsetting | Performance for large genomes |
| **Medium** | Add Spearman/Kendall option | Flexibility for non-normal data |
| **Low** | Caching layer | Workflow convenience |
| **Low** | Visualization helpers | Usability polish |
| **Low** | Clean up check artifacts | Repo hygiene |
