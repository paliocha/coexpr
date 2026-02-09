# coexpr: Cross-Species Co-Expression Network Analysis

R package implementing the methodology from **Grønvold & Hvidsten (2017)**:
*Cross species comparative transcriptomics using co-expression networks*.

## Motivation

Comparing gene regulation across species is fundamental to understanding how
biological processes evolve. Traditional approaches require matched experimental
conditions (same tissues, developmental stages, treatments) in both species,
which is rarely available and limits the data that can be used.

This package implements an alternative strategy that **compares co-expression
patterns rather than expression levels directly**. Because co-expression
captures regulatory relationships (which genes are controlled together), it can
be compared across species even when the underlying experiments are completely
different. A drought experiment in rice and a developmental time-course in
Arabidopsis both reveal the same regulatory wiring, just from different angles.

## Theoretical Background

### Co-expression as a proxy for regulation

If two genes are consistently co-expressed across many conditions, they likely
share regulatory elements or participate in the same pathway. This
"guilt-by-association" principle is the foundation of co-expression network
analysis (Wolfe et al. 2005). The key insight of Grønvold & Hvidsten is that
**orthologous genes with conserved regulation will have conserved
co-expression partners**, even when the experiments used in each species are
unrelated.

### The reference frame: 1:1 orthologs

Comparing co-expression across species requires a common coordinate system.
The package uses **1:1 orthologous genes** — genes with a single counterpart
in each species — as a shared reference frame. For each gene, its
co-expression pattern is represented as a vector of similarities to this
reference set. Two orthologous genes with conserved regulation will have
correlated reference vectors, even though the underlying expression data come
from different experiments.

### Pipeline overview

The analysis proceeds in four steps:

```
Expression matrices (per species)
        │
        ▼
   ┌─────────────┐
   │  Similarity  │  S_ij = f(E_i, E_j)
   │   matrices   │  within-species co-expression
   └──────┬──────┘
          │
          ▼
   ┌─────────────┐
   │     CCS      │  CCS_ij = cor(S^sp1_r,i, S^sp2_r,j)
   │              │  cross-species co-expression correlation
   └──────┬──────┘
          │
          ▼
   ┌─────────────┐
   │     ORS      │  ORS_ij = rank(CCS_ij) / n
   │              │  statistical ranking of conservation
   └──────┬──────┘
          │
          ▼
   Conserved / diverged orthologs
```

### Step 1: Similarity matrix

For each species, a gene-by-gene similarity matrix captures within-species
co-expression. The package supports several methods:

**PCC+MR** (Pearson correlation + Mutual Rank; recommended):

1. Compute Pearson correlation between all gene pairs:
   *S*<sup>PCC</sup><sub>ij</sub> = cor(*E*<sub>i</sub>, *E*<sub>j</sub>)

2. Transform via mutual rank normalization (Obayashi & Kinoshita 2009):
   *S*<sup>MR</sup><sub>ij</sub> = 1 − log(√(*R*<sub>ij</sub> · *R*<sub>ji</sub>)) / log(*n*)

   where *R*<sub>ij</sub> is the rank of gene *j* in gene *i*'s correlation
   row (descending). The log-mutual-rank emphasises the strongest
   co-expression partners and reduces bias from large co-expression clusters
   ("hub genes").

**Alternative methods**: Spearman correlation (robust to non-normality),
MI+CLR (mutual information with context likelihood ratio normalization;
Faith et al. 2007), or raw Pearson/Spearman without mutual rank. All can be
selected via function parameters. PCC+MR consistently outperforms alternatives
in cross-species benchmarks (Grønvold & Hvidsten, Figure 3).

### Step 2: Co-expression Correlation Score (CCS)

For an ortholog pair (*i* in species 1, *j* in species 2), CCS measures
whether *i* and *j* co-express with the same reference genes:

*CCS*<sub>ij</sub> = cor(*S*<sup>sp1</sup><sub>r,i</sub>, *S*<sup>sp2</sup><sub>r,j</sub>)

where *r* is the set of 1:1 reference orthologs. High CCS means the gene
pair has conserved co-expression partners; CCS near zero means no conserved
regulatory relationship.

### Step 3: Ortholog Rank Score (ORS)

CCS values are not directly comparable across datasets because their magnitude
depends on the number and diversity of samples. ORS solves this by converting
CCS to a rank-based percentile:

*ORS*<sub>ij</sub> = rank(*CCS*<sub>ij</sub>) / *n*

The log-transformed version provides an intuitive interpretation:

*logORS*<sub>ij</sub> = −log<sub>10</sub>(1 + 10<sup>−4</sup> − *ORS*<sub>ij</sub>)

| logORS | Meaning |
|--------|---------|
| > 0 | Above median conservation |
| > 1 | Top 10% — conserved |
| > 2 | Top 1% — highly conserved |
| > 3 | Top 0.1% — extremely conserved |

### Multi-copy orthologs

Gene duplication is pervasive and creates 1:N and N:M ortholog relationships.
The paper shows that duplicated genes exhibit **progressively lower
conservation**: 1:2 orthologs have lower ORS than 1:1, 1:3 lower still, and
so on — consistent with expression divergence after duplication (Ohno 1970).
A notable exception is recent whole-genome duplication (e.g. ~13 Mya in
*Glycine max*), where duplicates retain ancestral expression patterns.

The package provides five strategies for handling multi-copy orthologs:

| Strategy | Description | When to use |
|----------|-------------|-------------|
| `strict` | 1:1 orthologs only | Default; matches the paper |
| `best_hit` | Select highest-CCS copy | Need more gene coverage |
| `max` | Most conserved paralog | Identify retained ancestral function |
| `mean` | Average across paralogs | Assume sub-functionalization |
| `all_pairs` | Keep all N×M combinations | Study duplication-driven divergence |

## Installation

```r
# From GitHub
devtools::install_github("paliocha/coexpr")

# From local source
devtools::install_local("/path/to/coexpr")
```

**Requirements**: R ≥ 4.1, a C++17 compiler (for Rcpp/RcppArmadillo).

## Usage

### Minimal workflow

```r
library(coexpr)

# 1. Similarity matrices (one per species)
sim_sp1 <- calculate_pcc_mr(expr_sp1)
sim_sp2 <- calculate_pcc_mr(expr_sp2)

# 2. Co-expression correlation score
ccs <- calculate_ccs(sim_sp1, sim_sp2, orthologs)

# 3. Ortholog rank score
ors <- calculate_ors(ccs)

# 4. Significance testing (BH-corrected)
ors <- test_ors_significance(ors)
```

### Similarity methods

```r
# PCC+MR (default, recommended)
sim <- calculate_pcc_mr(expr)

# Spearman + MR (robust to non-normal data)
sim <- calculate_pcc_mr(expr, cor_method = "spearman")

# Raw Pearson (no mutual rank)
sim <- calculate_pcc_mr(expr, method = "pcc")

# MI+CLR (mutual information, needs ≥30 samples)
sim <- calculate_mi_clr(expr, n_bins = 10)

# Parallel computation
sim <- calculate_pcc_mr(expr, n_cores = 4)

# Cache to disk (avoids recomputation)
sim <- calculate_pcc_mr(expr, cache_dir = "cache/")
```

### Multi-copy orthologs

```r
# Default: strict 1:1 only
ortho_strict <- handle_multicopy_orthologs(orthologs, strategy = "strict")

# Keep all pairs for divergence analysis
ortho_all <- handle_multicopy_orthologs(orthologs, strategy = "all_pairs")
ccs <- calculate_ccs(sim_sp1, sim_sp2, ortho_all)

# Best-hit requires similarity matrices
ortho_best <- handle_multicopy_orthologs(
  orthologs, strategy = "best_hit",
  similarity_sp1 = sim_sp1, similarity_sp2 = sim_sp2
)
```

### Visualisation

```r
library(ggplot2)

# logORS distribution
ggplot(ors, aes(x = logORS)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = c(1, 2, 3), linetype = "dashed") +
  labs(x = "logORS", y = "Density") +
  theme_minimal()

# Conservation by ortholog type
ggplot(ors, aes(x = type, y = logORS, fill = type)) +
  geom_boxplot() +
  theme_minimal()
```

## Input data

| Input | Format | Description |
|-------|--------|-------------|
| Expression matrix | Numeric matrix (genes × samples) | Log-transformed or VST-normalised. Gene IDs as rownames. |
| Ortholog table | Data frame | Columns `gene_sp1`, `gene_sp2`, and optionally `type` ("1:1", "1:N", "N:1", "N:M"). |

**Sample requirements** (from the paper):
- Diverse samples (different tissues, conditions, studies) outperform biological
  replicates of the same condition.
- Diminishing returns after ~100–200 diverse samples.
- MI+CLR needs at least 3 × `n_bins` samples (default: 30).

## Citation

If you use this package, please cite:

> Grønvold L, Hvidsten TR. Cross species comparative transcriptomics using
> co-expression networks. NMBU Philosophiae Doctor (PhD) Thesis 2017: 101.
> [hdl:11250/2712242](http://hdl.handle.net/11250/2712242)

## References

- Grønvold L, Hvidsten TR (2017). Cross species comparative transcriptomics
  using co-expression networks. *NMBU PhD Thesis*.
- Obayashi T, Kinoshita K (2009). Rank of correlation coefficient as a
  comparable measure for biological significance of gene coexpression.
  *DNA Research* 16(5): 249–260.
  [doi:10.1093/dnares/dsp016](https://doi.org/10.1093/dnares/dsp016)
- Faith JJ et al. (2007). Large-scale mapping and validation of *Escherichia
  coli* transcriptional regulation from a compendium of expression profiles.
  *PLoS Biology* 5: e8.
  [doi:10.1371/journal.pbio.0050008](https://doi.org/10.1371/journal.pbio.0050008)
- Ohno S (1970). *Evolution by Gene Duplication*. Springer.
- Wolfe CJ, Kohane IS, Butte AJ (2005). Systematic survey reveals general
  applicability of "guilt-by-association" within gene coexpression networks.
  *BMC Bioinformatics* 6: 227.
  [doi:10.1186/1471-2105-6-227](https://doi.org/10.1186/1471-2105-6-227)

## License

MIT
