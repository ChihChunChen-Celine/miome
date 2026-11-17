# miome <img src="man/figures/logo.png" align="right" height="120" alt="" />

An interactive R Shiny application for end-to-end microbiome analysis:

- **Alpha diversity** — scatter plot visualization + linear mixed effects model (LMM) statistics
- **Beta diversity** — PCoA ordination (2D & 3D interactive) + PERMANOVA and ANOSIM statistics
- **Taxonomy bar plots** — relative abundance across all taxonomic levels
- **Differential abundance analysis** — MaAsLin3
- **Targeted measurements** (qPCR, metabolites) — linear (mixed) models + descriptive statistics

---

## Install

```r
# Install devtools if needed
install.packages("devtools")

# Install miome from GitHub
devtools::install_github("ChihChunChen-Celine/miome")

# Install optional dependency (required for MaAsLin3 tab only; needs R >= 4.4)
devtools::install_github("biobakery/maaslin3")
```

## Launch the app
```r
library(miome)
run_app()
```

## Files needed
**Required**

1. **Abundance table** (.tsv/.txt, taxa in rows, samples in columns), raw counts or
   relative abundance:
   - **Shotgun (MetaPhlAn): use the full merged table** (`merged_all.txt`, the output of
     `merge_metaphlan_tables.py`), whose first column holds the whole lineage
     (`k__Bacteria|p__Firmicutes|...|s__Blautia_wexlerae`) and which lists every rank.
     This is what makes phylum-to-genus analysis possible. A species-only export
     (just `Blautia_wexlerae`) still loads, but can only be analysed "as uploaded".
   - **16S (QIIME 2): the collapsed table** (e.g. L7 species level,
     `d__Bacteria;p__...;g__...;s__...`), with or without the `#OTU ID` header. Use a
     rarefied table if you want comparable alpha diversity.
   - HUMAnN pathway / gene tables (e.g. `pathabundance_relab_unstratified.tsv`) load too
     and are analysed as uploaded (no taxonomic levels).
   - Sample names may carry a pipeline suffix (`GI220_kneaddata_combined_Abundance`,
     `GI091_metaphlan`); they are matched to the metadata ID (`GI220`, `GI091`).
2. **Metadata** (.csv, one row per sample). One column must hold the Sample IDs exactly
   as they appear in the abundance table; the app finds it automatically. The metadata
   may list only a subset of the sequenced samples: only samples present in both files
   are analysed, and the app reports which ones were left out.

After uploading, pick your subject/animal ID column and baseline flag column in
**Step 2** of the Data page (leave them at "(none)" if your study has no repeated
measures or baseline samples). Treatment, formulas and covariates are chosen later
on each analysis page. Column names are tidied on loading (`BW (kg)` becomes `BW.kg`).

**Optional**
- Precomputed UniFrac distance matrix (.qza or .tsv) - required for UniFrac.
- **Targeted measurements** (.csv), uploaded on the *Measurements* page. Sequencing
  data is not needed for these: load the metadata alone and every metadata row is
  used, sequenced or not. qPCR results, fatty acids / sterols, bile acids or other per-sample values. One row
  per sample, one column per measurement; the first column holds an ID that also
  appears in the metadata (the Sample ID, or another ID column such as a sequencing
  number). You choose which measurements to test; p-values are adjusted across your
  selection, and descriptive statistics (N, mean, SD, median, min, max, 95% CI) are
  given per group.

## Example data
Small synthetic datasets for trying the app are installed with the package:
```r
system.file("testdata", package = "miome")
```
