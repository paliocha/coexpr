# Architecture: Backend Selection

This document explains why `coexpr` uses different computational backends for different operations.

## Summary

| Operation | Backend | Library | Rationale |
|-----------|---------|---------|-----------|
| Pearson correlation | OpenBLAS/Fortran | `Rfast::cora()` | Maps to BLAS DGEMM |
| Mutual rank transform | Armadillo/C++ | `RcppArmadillo` | Custom algorithm + OpenMP |
| CCS correlation | Base R | `stats::cor()` | Small vectors, overhead not worth it |

## Why a Hybrid Approach?

Neither OpenBLAS/Fortran nor Armadillo/C++ is universally better—they excel at different tasks.

### OpenBLAS/Fortran (via Rfast)

**Best for**: Standard linear algebra operations

- BLAS routines have 30+ years of optimization
- OpenBLAS auto-detects CPU and uses optimal kernels
- Hardware-specific tuning for matrix multiplication, decompositions, eigenvalues

**Limitations**:
- Black box: limited control over algorithm implementation
- Not ideal for custom algorithms requiring non-standard logic

### Armadillo/C++ (via RcppArmadillo)

**Best for**: Custom algorithms with fine-grained control

- Full algorithmic control in readable C++
- Easy OpenMP integration (`#pragma omp parallel for`)
- Expression templates avoid temporary allocations
- Mix linear algebra with custom logic seamlessly

**Limitations**:
- Development and maintenance overhead
- Longer compilation times
- For pure BLAS operations, may be slower than direct OpenBLAS

---

## Operation-Specific Decisions

### 1. Pearson Correlation → OpenBLAS/Fortran

**File**: `R/similarity.R` (uses `Rfast::cora()`)

The correlation matrix calculation is essentially:

```
cor(X) = (X - μ)ᵀ(X - μ) / (n-1)  # with normalization
```

This maps directly to **BLAS Level 3 DGEMM** (matrix-matrix multiplication). Writing this in Armadillo/C++ would be:
- More code to maintain
- Same or slower performance (Armadillo uses BLAS internally anyway)
- No algorithmic benefit

**Decision**: Use Rfast ✓

### 2. Mutual Rank Transform → Armadillo/C++

**File**: `src/mutual_rank.cpp`

The mutual rank algorithm:

```cpp
for each pair (i, j):
    R_ij = rank of S_ij in row i (descending)
    R_ji = rank of S_ji in row j (descending)
    MR_ij = 1 - log(sqrt(R_ij * R_ji)) / log(n)
```

This is **not a standard linear algebra operation**:
- Requires row-wise ranking (custom sorting operation)
- Requires element-wise geometric mean
- Benefits from OpenMP parallelization of the outer loop
- Has two memory strategies (cached vs streaming)

**Decision**: Use RcppArmadillo ✓

### 3. CCS Calculation → Base R

**File**: `R/ccs.R` (uses `stats::cor()`)

CCS computes correlation of two vectors (~1,000–10,000 elements):

```r
cor(coexpr_sp1, coexpr_sp2)
```

For vectors this small:
- BLAS call overhead exceeds computation time
- Base R `cor()` is sufficient
- Parallelization happens at the pair level (via `furrr`), not within correlation

**Decision**: Use Base R ✓

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────┐
│ Expression Matrix (n genes × m samples)             │
└───────────────────────┬─────────────────────────────┘
                        ↓
┌───────────────────────────────────────────────────────┐
│ Rfast::cora() - OpenBLAS/Fortran                      │
│ • DGEMM-based correlation                             │
│ • Hardware-optimized matrix multiplication            │
│ • O(n² × m)                                           │
└───────────────────────┬───────────────────────────────┘
                        ↓
┌───────────────────────────────────────────────────────┐
│ mutual_rank_transform_cached_cpp() - Armadillo/C++    │
│ • Custom ranking algorithm                            │
│ • OpenMP parallel for loop                            │
│ • O(n² log n)                                         │
└───────────────────────┬───────────────────────────────┘
                        ↓
┌───────────────────────────────────────────────────────┐
│ Similarity Matrix (n × n)                             │
└───────────────────────┬───────────────────────────────┘
                        ↓
┌───────────────────────────────────────────────────────┐
│ calculate_ccs() - Base R + furrr                      │
│ • cor() on small vectors                              │
│ • Parallelization at pair level                       │
│ • O(k × r) where k=pairs, r=reference size            │
└───────────────────────────────────────────────────────┘
```

---

## Potential Future Optimizations

These are **not currently needed** but documented for future reference:

### BLAS Backend Selection

R uses whatever BLAS it was compiled with. Users can optimize:

```r
# Check current BLAS
sessionInfo()$BLAS

# Options by platform:
# - Apple Silicon: Accelerate framework (default on macOS)
# - Intel: MKL (via Microsoft R Open or manual configuration)
# - Generic: OpenBLAS (common on Linux)
```

### GPU Acceleration

For very large matrices (n > 50,000 genes):
- Replace `Rfast::cora()` with GPU correlation (e.g., via `torch` or CUDA)
- Only worthwhile when matrix operations dominate runtime

### Sparse Matrices

If many correlations are near zero:
- Sparse representation could reduce memory
- Would require rewriting mutual rank for sparse input
- Not typical for expression data (correlations are dense)

---

## Conclusion

The hybrid architecture uses each backend where it performs best:

| Use Case | Best Choice |
|----------|-------------|
| Standard linear algebra (correlation, SVD) | OpenBLAS/Fortran |
| Custom algorithms with control needs | Armadillo/C++ |
| Small vector operations | Base R |

This is the optimal design for this workload. No architectural changes recommended.
