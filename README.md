# Cross-species co-expression network analysis in R

R package implementing the methodology from **Grønvold & Hvidsten**: *Cross species comparative transcriptomics using co-expression networks*.

## Overview

Compare gene regulation across species by analyzing co-expression patterns without requiring directly comparable samples (same tissues, developmental stages, or conditions).

**Key innovation**: Uses 1:1 orthologous genes as a "reference frame" to compare co-expression patterns indirectly.

## Core Concepts

### 1. Similarity Matrix (PCC+MR)

Calculate within-species gene co-expression using **Pearson correlation + Mutual Rank** normalization:

```r
sim_matrix <- calculate_pcc_mr(expression_matrix)
```

- **Input**: Gene expression matrix (genes × samples), log-transformed
- **Output**: Symmetric similarity matrix (genes × genes), values in [0, 1]
- **Method**: Best-performing according to paper (Figure 3)

### 2. Co-expression Correlation Score (CCS)

Measures how well orthologous genes retain the same co-expression partners:

```r
ccs_results <- calculate_ccs(sim_sp1, sim_sp2, orthologs)
```

- **Input**: Two similarity matrices + ortholog table
- **Output**: CCS for each ortholog pair (correlation of co-expression vectors)
- **Reference**: Uses 1:1 orthologs only (unduplicated genes)

### 3. Ortholog Rank Score (ORS)

Statistical significance measure - ranks CCS among all genes:

```r
ors_results <- calculate_ors(ccs_results)
```

- **Output**: logORS values where:
  - `>1` = top 10% (conserved)
  - `>2` = top 1% (highly conserved)
  - `>3` = top 0.1% (extremely conserved)

## Multi-Copy Orthologs: Paper vs. This Package

### What the Paper Does

**Paper's approach**: Excludes all 1:N and N:M orthologs, uses only 1:1 reference

**Key findings** (Figure 6-8):
- Duplicated genes show **progressively lower conservation**: 1:2 < 1:3 < 1:4
- **Exception**: Recent WGD duplicates (~13 Mya in *Glycine max*) show NO divergence
- **Preferential retention**: Retained duplicates tend to be high-ORS genes

### What This Package Adds

Multiple strategies for handling multi-copy orthologs:

```r
orthologs_processed <- handle_multicopy_orthologs(
  orthologs_raw,
  strategy = "strict"  # or "best_hit", "mean", "max", "all_pairs"
)
```

| Strategy | Description | Use Case |
|----------|-------------|----------|
| `strict` (default) | Only 1:1 orthologs | Publication-ready, matches paper |
| `best_hit` | Select most similar copy | Need more gene coverage |
| `max` | Most conserved paralog | Identify conserved duplicates |
| `mean` | Average across paralogs | Assume sub-functionalization |
| `all_pairs` | Keep all combinations | Study expression divergence |

**Recommendation**: Start with `"strict"` for main analysis.

## Installation

```r
# Install dependencies
install.packages(c("dplyr", "tidyr", "purrr", "Rfast", "future", "furrr"))

# Install coexpr from source
devtools::install_local("/path/to/coexpr")
```

## Usage

### Complete Workflow

```r
library(coexpr)

# Step 1: Calculate similarity matrices for each species
expr_at <- load_expression_data("arabidopsis.csv")  # genes × samples
expr_os <- load_expression_data("rice.csv")

sim_at <- calculate_pcc_mr(expr_at)
sim_os <- calculate_pcc_mr(expr_os)

# Step 2: Load and process orthologs
orthologs_raw <- read.csv("orthologs.csv")  # columns: gene_sp1, gene_sp2

# Strategy 1: Strict (only 1:1, matches paper)
orthologs_1to1 <- handle_multicopy_orthologs(
  orthologs_raw,
  strategy = "strict"
)

# Step 3: Calculate CCS
ccs_results <- calculate_ccs(
  sim_sp1 = sim_at,
  sim_sp2 = sim_os,
  orthologs = orthologs_1to1,
  use_only_1to1 = TRUE
)

# Step 4: Calculate ORS
ors_results <- calculate_ors(ccs_results, return_log = TRUE)

# Step 5: Test significance
ors_results <- test_ors_significance(ors_results, alpha = 0.05)

# Step 6: Summarize
summarize_conservation(ors_results, by_type = TRUE)
```

### Visualize Results

```r
library(ggplot2)

# Distribution of logORS
ggplot(ors_results, aes(x = logors)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = c(1, 2, 3), linetype = "dashed") +
  labs(title = "Expression Conservation Across Species",
       x = "logORS", y = "Density") +
  theme_minimal()

# Compare ortholog types (if using multi-copy strategies)
ggplot(ors_results, aes(x = type, y = logors, fill = type)) +
  geom_boxplot() +
  labs(title = "Conservation by Ortholog Type",
       x = "Ortholog Type", y = "logORS") +
  theme_minimal()
```

### Identify Conserved Gene Pairs

```r
# Highly conserved orthologs (top 1%)
conserved <- ors_results %>%
  filter(logors > 2) %>%
  arrange(desc(logors))

# Export for downstream analysis
write.csv(conserved, "conserved_orthologs.csv", row.names = FALSE)
```

## Key Design Decisions

### 1. Native R Pipe (`|>`)

Uses R 4.1+ native pipe instead of `%>%`:

```r
# Modern R style
result <- data |>
  process_step1() |>
  process_step2()
```

### 2. Rfast for Performance

Uses `Rfast` package for matrix operations (10-100x faster than base R):

```r
# Fast correlation calculation
sim <- Rfast::cora(t(expr_matrix))
```

### 3. Vectorized Operations

No explicit loops - everything is vectorized for performance:

```r
# Vectorized rank calculation
ors <- rank(ccs) / length(ccs)
```

## Expected Performance

- **Similarity matrix**: ~1-2 min for 20,000 genes × 100 samples
- **CCS calculation**: ~10 sec for 10,000 ortholog pairs
- **ORS calculation**: ~1 sec

**Bottleneck**: Similarity matrix calculation (do this once, cache results)

## Citation

If you use this package, please cite the original paper:

> Grønvold L, Hvidsten TR. Cross species comparative transcriptomics using co-expression networks. hdl:11250/2712242

## References

**Key papers**:
- Ohno S (1970) Evolution by Gene Duplication
- Obayashi T, Kinoshita K (2009) Rank of correlation coefficient. DNA Research 16:249-260
- Tirosh I, Barkai N (2007) Comparative analysis indicates regulatory neofunctionalization. Genome Biology 8:R50

**Sample requirements** (from paper):
- More diverse samples > replicates from same study
- Diminishing returns after ~100-200 diverse samples
- Paper used 35-2545 samples across 5 plant species

## Development

```bash
# Run tests
Rscript -e "devtools::test()"

# Check package
Rscript -e "devtools::check()"

# Build documentation
Rscript -e "devtools::document()"
```

## License

MIT License
