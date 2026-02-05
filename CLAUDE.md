# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Implementation of cross-species comparative transcriptomics using co-expression networks, based on the manuscript "Cross species comparative transcriptomics using co-expression networks" (Grønvold & Hvidsten).

**Core methodology**: Compare gene regulation across species by calculating co-expression correlation scores (CCS) between orthologous genes, then ranking these using the Ortholog Rank Score (ORS) to identify conserved expression patterns without requiring directly comparable samples.

## Key Algorithms to Implement

### 1. PCC+MR Similarity Matrix (Best performing method from paper)
- Calculate Pearson correlation between gene expression vectors: `S^PCC_ij = cor(E_i, E_j)`
- Transform using log mutual rank normalization: `S^(PCC+MR)_ij = 1 - log(√(R_ij·R_ji)) / log(n)`
- Where `R_ij` is the rank of `S^PCC_ij` in row i (ordered high to low)
- Use `Rfast::Dist()` and `Rfast::colRanks()` for fast matrix operations

### 2. Co-expression Correlation Score (CCS)
- For species pair with 1:1 ortholog reference set r:
- `CCS_ij = cor(S^sp1_r,i, S^sp2_r,j)`
- Only correlate rows/columns corresponding to ortholog pairs
- Use native R pipe `|>` for data flow

### 3. Ortholog Rank Score (ORS)
- `ORS_ij` = proportion of genes in species 2 with CCS ≤ CCS of ortholog pair
- Transform: `logORS_ij = -log_10(1 + 10^-4 - (ORS^sp1_ij + ORS^sp2_ji)/2)`
- ORS > 1 means top 10%, > 2 means top 1%, > 3 means top 0.1%
- Higher ORS indicates more conserved co-expression

## R Development Standards

**Code style**:
- Native pipe operator `|>` (not `%>%`)
- Modern tidyverse syntax (dplyr 1.1+, tidyr 1.3+)
- Use `Rfast` package for matrix calculations (faster than base R)
- Vectorize operations; avoid explicit loops where possible
- Use `future` + `furrr` for parallel processing of species pairs

**Key packages**:
- Data manipulation: `dplyr`, `tidyr`, `purrr`
- Matrix ops: `Rfast` (for fast correlation, ranking, distance)
- Parallel: `future`, `furrr`
- Visualization: `ggplot2`, `ComplexHeatmap`
- Orthology: Interface with Ensembl/biomaRt for ortholog mapping

## Data Structures

**Expression matrices**:
- Genes × samples, with log-transformed expression values
- Store as regular matrices for Rfast compatibility (not tibbles)
- Gene IDs as rownames

**Ortholog tables**:
- Tibbles with columns: `gene_sp1`, `gene_sp2`, `ortholog_type` (1:1, 1:N, N:1, N:M)
- Filter to 1:1 orthologs for reference set

**Similarity matrices**:
- Square symmetric matrices (n_genes × n_genes) per species
- Upper triangle sufficient due to symmetry
- Store as full matrices for correlation calculation

## Mathematical Implementation Notes

**Mutual rank normalization**:
- Requires ranking twice: row-wise and column-wise
- Log transformation emphasizes strong correlations
- Subtract from 1 to normalize to [0,1] range with 1 = most similar

**Handling duplicated genes**:
- Paper shows 1:2, 1:3, 1:4 orthologs have progressively lower ORS
- Recent WGD duplicates (e.g., Glycine max) may show equal conservation
- Consider separate analysis for singleton vs duplicated gene orthologs

**Sample requirements**:
- More diverse samples (different studies) > replicates from same study
- Paper used ~10-1400 samples per species across 5 plant species
- Diminishing returns after ~100-200 diverse samples

## Expected Outputs

1. **Species-pair CCS matrices**: n_orthologs × n_genes for each species pair
2. **ORS values**: Per ortholog pair, measures conservation significance
3. **logORS distributions**: Compare between/within species to validate method
4. **Visualization**: Heatmaps, density plots, boxplots by duplication status

## Validation Checks

- Within-species ORS should be higher than between-species ORS
- 1:1 orthologs should have higher ORS than 1:N orthologs (except recent WGD)
- Median logORS should reflect phylogenetic distance
- Check that PCC+MR outperforms PCC, MI, SCC, CLR methods on your data

## Handling Multi-Copy Orthologs

**Paper's approach**: Uses ONLY 1:1 orthologs as reference (excludes 1:N and N:M)

**Key findings from paper**:
- 1:2, 1:3, 1:4 orthologs show progressively lower ORS (expression divergence)
- Exception: Recent WGD duplicates (~13 Mya in Glycine max) show NO divergence
- Hypothesis: Retained duplicates are preferentially high-ORS genes

**Our implementation** (`handle_multicopy_orthologs()` function):

Provides multiple strategies for incorporating multi-copy orthologs:

1. **`"strict"` (default, matches paper)**:
   - Uses only 1:1 orthologs as reference
   - Excludes all 1:N and N:M relationships
   - Most conservative, best for cross-species comparison

2. **`"best_hit"`**:
   - For 1:N: Select single best copy based on CCS magnitude
   - For N:M: Select best reciprocal pair
   - Good when one copy clearly retains ancestral function

3. **`"mean"`**:
   - For 1:N: Average CCS across all N copies
   - For N:M: Average across all N×M combinations
   - Assumes partial sub-functionalization

4. **`"max"`**:
   - For 1:N: Use maximum CCS (most conserved copy)
   - For N:M: Use maximum across all combinations
   - Identifies the most conserved paralog

5. **`"all_pairs"`**:
   - Keep all ortholog combinations separately
   - Useful for analyzing expression divergence among paralogs
   - Increases dataset size substantially

**Recommendation**: Start with `"strict"` for main analysis, then use `"best_hit"` or `"max"` for extended analysis if you need more gene coverage.

## R Package Structure

```
coexpr/
├── R/                  # R source code
│   ├── similarity.R    # PCC+MR similarity matrix calculation
│   ├── ccs.R          # Co-expression correlation score
│   ├── ors.R          # Ortholog rank score
│   ├── multicopy.R    # Multi-copy ortholog handling
│   └── utils.R        # Helper functions
├── man/               # Documentation (auto-generated)
├── tests/             # Unit tests
│   └── testthat/
├── vignettes/         # Long-form documentation
└── data-raw/          # Scripts to process raw data
```

## Performance Considerations

- Pre-compute similarity matrices once per species (most expensive step)
- Parallelize across ortholog pairs or species combinations
- Use sparse matrix representations if many zero correlations
- Consider chunking large gene sets for memory efficiency
- Cache intermediate results (similarity matrices, ortholog maps)
