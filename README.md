# miome <img src="man/figures/logo.png" align="right" height="120" alt="" />

An interactive R Shiny application for end-to-end microbiome analysis:

- **Alpha diversity** — scatter plot visualization + linear mixed effects model (LMM) statistics
- **Beta diversity** — PCoA ordination (2D & 3D interactive) + PERMANOVA and ANOSIM statistics
- **Taxonomy bar plots** — relative abundance across all taxonomic levels
- **Differential abundance analysis** — MaAsLin3

---

## Install

```r
# Install devtools if needed
install.packages("devtools")

# Install miome from GitHub
devtools::install_github("ChihChunChen-Celine/miome")

# Install optional dependency (required for MaAsLin3 tab only)
devtools::install_github("biobakery/maaslin3")

```

## Launch the app
```r
library(miome)
run_app()
```

## Files needed
**Required**  
One taxonomy or microbial gene composition (samples in column, features in row)  
One mapping/metadata file (samples in row, features in column)  
**Optional**  
Precomputed UniFrac distance matrix (.qza or .tsv) - required for UniFrac  


