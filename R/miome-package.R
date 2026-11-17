#' Miome: Microbiome Explorer
#' @keywords internal
"_PACKAGE"

#' @importFrom utils head tail read.csv read.table write.csv unzip packageVersion
#' @importFrom stats setNames as.dist cmdscale as.formula terms aggregate dist
#' @importFrom stats reformulate median p.adjust complete.cases relevel
#' @importFrom magrittr %>%
NULL

# Column names used inside dplyr / ggplot2 non-standard evaluation.
utils::globalVariables(c(
  "SampleID", "RelAbundance", "Taxon", "sampleid", "max_abund",
  "relative_abundance", "median_abundance", "min_abundance", "max_abundance"
))
