# Otitis Media Microbiome — Figure Code

Code repository accompanying the manuscript:

> **[Manuscript Title]**  
> Leichtel, Künstner et al.  
> [Journal], [Year]  
> DOI: [DOI]

---

## Overview

This repository contains the R code used to generate all figures in the manuscript.
The analysis examines microbiome composition in otitis media patients, comparing
bacterial and fungal communities before and after treatment, and identifying microbiome
features associated with treatment response.

---

## Data

Two preprocessed phyloseq objects are provided as starting points for all figure scripts.
Raw sequencing data were processed using nf-core/taxprofiler and nf-core/detaxizer
(Nextflow workflows with Singularity containers), followed by taxonomic profiling with
Bracken/Kraken2. Preprocessing steps included selection of the highest-yield sample per
patient per timepoint, filtering of taxa present in fewer than 10% of samples, and removal
of samples below minimum read depth (≥500 reads for bacteria, ≥200 for fungi). Alpha
diversity was estimated using DivNet and breakaway.

| File | Description |
|------|-------------|
| `data/bacteria_phyloseq_processed.rds` | Filtered bacterial phyloseq object with sample metadata and precomputed alpha diversity estimates |
| `data/fungi_phyloseq_processed.rds` | Filtered fungal phyloseq object with sample metadata and precomputed alpha diversity estimates |
| `data/network_analysis_results.rds` | Precomputed network objects (see note below) |

---

## Network Analysis Note

Co-occurrence networks were constructed using Spearman correlation rather than SpiecEasi.
Although the code in `15_network_analysis.R` attempts SpiecEasi first, it consistently
fell back to Spearman correlation due to computational constraints — consistent with what
is reported in the Methods section of the manuscript. The precomputed
`network_analysis_results.rds` is provided so that `16_Fig4_5.R` (Fig 4 & 5)
can be run directly without repeating the network construction step, which is
computationally intensive.

---

## Code

All scripts should be run from this directory. All scripts source
`_helper.R` which defines shared color palettes and themes.

| Script | Figure | Description |
|--------|--------|-------------|
| `_helper.R` | — | Shared color palettes, ggplot2 theme, and `create_taxa_plot()` function |
| `13_diversity_analysis.R` | Fig 1 | Alpha diversity (DivNet Shannon) and beta diversity (Aitchison PCoA) before vs. after treatment, for bacteria and fungi |
| `14_taxa_DA.R` | Fig 2 & 3 | Taxonomic composition bar plots at Phylum, Genus, and Species level for bacteria (Fig 2) and fungi (Fig 3) |
| `15_network_analysis.R` | Fig 4 & 5 (input) | Co-occurrence network construction via Spearman correlation, separately for bacteria and fungi, before and after treatment |
| `16_advanced_network.R` | Fig 4 & 5 (output) | Genus-level network visualization (Fig 4) and cross-domain bacteria–fungi interaction network (Fig 5) |
| `20_Responder_diversity.R` | Fig 6 | Alpha and beta diversity stratified by treatment responder status and timepoint |
| `21_Responder_taxa.R` | Fig 7 | Differential abundance analysis (ANCOM-BC2) comparing responders vs. non-responders at T1 and T2 |

---

## Dependencies

All analyses were performed in R. The following packages are required:

```r
# CRAN / Bioconductor
tidyverse, phyloseq, patchwork, ggpubr, vegan, microbiome,
ANCOMBC, breakaway, DivNet, SpiecEasi, igraph, ggnetwork,
ggrepel, pals, RColorBrewer, gt, lme4, lmerTest


```
