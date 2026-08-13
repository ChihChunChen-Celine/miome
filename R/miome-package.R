#' Miome: Microbiome Explorer
#' @keywords internal
"_PACKAGE"

# ── Import base utils functions used throughout the package ──────────────────
#' @importFrom utils head read.csv read.table write.csv unzip packageVersion
#' @importFrom stats setNames as.dist cmdscale as.formula terms aggregate dist
#' @importFrom stats reformulate
NULL

# ── Silence "no visible binding for global variable" NOTES ───────────────────
# These are column names used inside dplyr/ggplot2 non-standard evaluation.
utils::globalVariables(c(
  "SampleID", "RelAbundance", "Taxon", "clade_name"
))
