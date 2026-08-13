# -----------------------------------------------------------------------------
# config.R
# Central configuration for the 16S Microbiome Explorer
# Covers: PCoA / Beta Diversity / Alpha Diversity / MaAsLin3
# -----------------------------------------------------------------------------

# -- Palette options -----------------------------------------------------------

#' Categorical colour palette options (for grouping variables)
#' @export
CAT_PALETTES <- c("Set1", "Set2", "Set3", "Dark2", "Paired", "Accent")

#' Continuous colour palette options (for numeric variables)
#' @export
CONT_PALETTES <- c("RdBu", "RdYlBu", "Spectral", "YlOrRd", "Blues")

#' Supported beta diversity distance methods
#'
#' Restricted to methods actually implemented in \code{.compute_dist()}
#' (data.R): most are passed straight to \code{vegan::vegdist()}, while
#' \code{"aitchison"} is computed via a CLR transform + Euclidean distance
#' and \code{"unifrac_weighted"} / \code{"unifrac_unweighted"} are computed
#' via \code{phyloseq::UniFrac()} using an uploaded phylogenetic tree.
#' @export
DIST_METHODS <- c(
  "Bray-Curtis" = "bray",
  "Jaccard"     = "jaccard",
  "Euclidean"   = "euclidean",
  "Aitchison"    = "aitchison",
  "UniFrac (weighted)"   = "unifrac_weighted",
  "UniFrac (unweighted)" = "unifrac_unweighted",
  "Manhattan"   = "manhattan",
  "Canberra"    = "canberra",
  "Kulczynski"  = "kulczynski"
)

#' Supported alpha diversity metrics
#' @export
ALPHA_METRICS <- c(
  "Shannon",
  "Simpson",
  "Observed_ASVs",
  "Chao1",
  "Pielou"
)

#' Taxonomic levels available for the taxonomy bar plot
#'
#' Maps a human-readable taxonomic level (shown in the Shiny UI dropdown
#' and used as \code{level} in \code{compute_taxonomy_barplot_data()} /
#' \code{build_taxon_barplot()}, data.R and plots.R) to the taxonomy
#' clade-name prefix used to parse it out of the semicolon-delimited
#' lineage string (e.g. \code{"f__Ruminococcaceae"}). Uses the
#' QIIME2 / SILVA-style prefix convention typical of 16S rRNA taxonomy
#' assignments (\code{d__} for Domain), not MetaPhlAn's \code{k__}.
#' @export
TAXA_PREFIXES <- c(
  "Domain"  = "d",
  "Phylum"  = "p",
  "Class"   = "c",
  "Order"   = "o",
  "Family"  = "f",
  "Genus"   = "g",
  "Species" = "s"
)

#' Taxonomic level choices shown in the Shiny UI
#' @export
TAXA_LEVELS <- names(TAXA_PREFIXES)

#' Default number of top taxa shown individually in the taxonomy bar plot
#' before the remainder are collapsed into \code{"Others"}.
#' @export
TAXA_DEFAULT_TOP_N <- 14L

#' Supported p-value adjustment methods
#' @export
PADJ_METHODS <- c(
  "Benjamini-Hochberg (FDR)" = "BH",
  "Benjamini-Yekutieli (BY)" = "BY",
  "Bonferroni"               = "bonferroni",
  "Holm"                     = "holm",
  "None"                     = "none"
)

# -- Main configuration builder ------------------------------------------------

#' Build a Microbiome Explorer configuration object
#'
#'
#' @section Beta Diversity Parameters:
#' \describe{
#'   \item{dist_method}{Distance metric key passed to \code{vegan::vegdist}.
#'     Choose from \code{DIST_METHODS}.}
#'   \item{default_covariates}{Character vector of metadata columns to
#'     pre-select as covariates in PERMANOVA (Shiny UI default).}
#'   \item{default_strata}{Metadata column to pre-select for restricted
#'     permutations (Shiny UI default). Set \code{NULL} to disable.}
#'   \item{padj_method}{P-value adjustment method for pairwise tests.
#'     Choose from \code{PADJ_METHODS}.}
#' }
#'
#' @section Alpha Diversity Parameters:
#' \describe{
#'   \item{alpha_metrics}{Character vector of metrics to compute.
#'     Choose from \code{ALPHA_METRICS}.}
#'   \item{alpha_export_w_in, alpha_export_h_in}{Alpha diversity plot
#'     export dimensions in inches.}
#' }
#'
#' @section Plot / Export Parameters:
#' \describe{
#'   \item{plot_title}{Title string shown on PCoA plots.}
#'   \item{dist_method_label}{Human-readable label for the distance metric
#'     (used in plot subtitles).}
#'   \item{marker_size, marker_opacity}{Point aesthetics.}
#'   \item{export_dpi, export_w_in, export_h_in}{TIFF export settings.}
#'   \item{export_3d_w_px, export_3d_h_px}{In-browser 3D PNG snapshot size.}
#'   \item{study_design}{Study design label, e.g. \code{"crossover"} or
#'     \code{"parallel"}. Used by app.R to default the "connect paired
#'     samples" checkbox to TRUE for crossover designs.}
#' }
#'
#' @return A list of class \code{pcoa_config}.
#' @export
#' @param dist_method Distance method.
#' @param dist_method_label Human-readable distance label.
#' @param plot_title Plot title.
#' @param study_design Study design string (e.g. "crossover").
#' @param treatment_col Treatment column name.
#' @param timepoint_col Timepoint column name.
#' @param sequence_col Sequence column name.
#' @param subject_id_col Subject/animal ID column name.
#' @param is_baseline_col Baseline flag column name.
#' @param alpha_metrics Character vector of alpha diversity metrics.
#' @param default_covariates Default PERMANOVA covariates.
#' @param default_strata Default strata column.
#' @param padj_method P-value adjustment method.
#' @param marker_size 3D marker size.
#' @param marker_opacity 3D marker opacity.
#' @param export_dpi Export resolution.
#' @param export_w_in Export width (inches).
#' @param export_h_in Export height (inches).
#' @param export_3d_w_px 3D export width (px).
#' @param export_3d_h_px 3D export height (px).
#' @param alpha_export_w_in Alpha plot export width (inches).
#' @param alpha_export_h_in Alpha plot export height (inches).
#' @param maaslin_normalization MaAsLin3 normalization.
#' @param maaslin_transform MaAsLin3 transform.
#' @param maaslin_method MaAsLin3 analysis method.
#' @param maaslin_min_prevalence MaAsLin3 min prevalence.
#' @param maaslin_min_abundance MaAsLin3 min abundance.
#' @param maaslin_max_sig MaAsLin3 max significance.
pcoa_config <- function(
    # -- Column names ---------------------------------------------
  subject_id_col     = "sampleid",
  treatment_col      = "Treatment",
  timepoint_col      = "Timepoint",
  sequence_col       = "Sequence",
  is_baseline_col    = "IsBaseline",
  # -- Study design ----------------------------------------------
  study_design        = "crossover",
  # -- Beta diversity -------------------------------------------
  dist_method        = "bray",
  dist_method_label  = "Bray-Curtis",
  default_covariates = c("Timepoint", "Sequence"),
  default_strata     = "sampleid",
  padj_method        = "BH",
  # -- Alpha diversity ------------------------------------------
  alpha_metrics      = ALPHA_METRICS,
  alpha_export_w_in  = 8,
  alpha_export_h_in  = 6,
  # -- Plot aesthetics ------------------------------------------
  plot_title         = "PCoA",
  marker_size        = 6,
  marker_opacity     = 0.85,
  # -- 2D / TIFF export -----------------------------------------
  export_dpi         = 300,
  export_w_in        = 9,
  export_h_in        = 7,
  # -- 3D PNG snapshot ------------------------------------------
  export_3d_w_px     = 2400,
  export_3d_h_px     = 1800
) {
  # -- Validate column name inputs ------------------------------
  char_args <- list(
    subject_id_col  = subject_id_col,
    treatment_col   = treatment_col,
    timepoint_col   = timepoint_col,
    sequence_col    = sequence_col,
    is_baseline_col = is_baseline_col
  )
  invisible(lapply(names(char_args), function(nm) {
    if (!is.character(char_args[[nm]]) || nchar(char_args[[nm]]) == 0)
      stop(sprintf("'%s' must be a non-empty character string.", nm))
  }))
  # -- Validate alpha metrics -----------------------------------
  bad_metrics <- setdiff(alpha_metrics, ALPHA_METRICS)
  if (length(bad_metrics) > 0)
    warning(sprintf(
      "Unknown alpha metric(s) ignored: %s. Choose from: %s",
      paste(bad_metrics, collapse = ", "),
      paste(ALPHA_METRICS, collapse = ", ")
    ))
  alpha_metrics <- intersect(alpha_metrics, ALPHA_METRICS)
  # -- Validate padj method -------------------------------------
  if (!padj_method %in% c("BH", "bonferroni", "holm", "none",
                          "fdr", "hochberg", "hommel", "BY"))
    stop("'padj_method' must be a valid p.adjust method.")
  # -- Assemble config list -------------------------------------
  structure(
    list(
      subject_id_col     = subject_id_col,
      treatment_col      = treatment_col,
      timepoint_col      = timepoint_col,
      sequence_col       = sequence_col,
      is_baseline_col    = is_baseline_col,
      study_design       = study_design,
      dist_method        = dist_method,
      dist_method_label  = dist_method_label,
      default_covariates = default_covariates,
      default_strata     = default_strata,
      padj_method        = padj_method,
      alpha_metrics      = alpha_metrics,
      alpha_export_w_in  = alpha_export_w_in,
      alpha_export_h_in  = alpha_export_h_in,
      plot_title         = plot_title,
      marker_size        = marker_size,
      marker_opacity     = marker_opacity,
      export_dpi         = export_dpi,
      export_w_in        = export_w_in,
      export_h_in        = export_h_in,
      export_3d_w_px     = export_3d_w_px,
      export_3d_h_px     = export_3d_h_px
    ),
    class = "pcoa_config"
  )
  # -- MaAsLin3 defaults ----------------------------------------------------
  #' @param maaslin_normalization  Default normalisation passed to MaAsLin3
  #'   (\code{"TSS"}, \code{"CLR"}, \code{"CSS"}, \code{"NONE"}).
  #' @param maaslin_transform      Default transform (\code{"LOG"},
  #'   \code{"SQRT"}, \code{"NONE"}).
  #' @param maaslin_method         Analysis method (\code{"LM"},
  #'   \code{"NEGBIN"}, \code{"ZINB"}).
  #' @param maaslin_min_prevalence Minimum prevalence filter (default 0.1).
  #' @param maaslin_min_abundance  Minimum abundance filter (default 0).
  #' @param maaslin_max_sig        FDR threshold for reporting (default 0.25).
}

# -- print method --------------------------------------------------------------

#' @export
print.pcoa_config <- function(x, ...) {
  cat("-- Microbiome Explorer Configuration ----------------------\n")
  cat(sprintf("  Subject ID col  : %s\n",   x$subject_id_col))
  cat(sprintf("  Treatment col   : %s\n",   x$treatment_col))
  cat(sprintf("  Study design    : %s\n",   x$study_design))
  cat(sprintf("  Distance method : %s\n",   x$dist_method_label))
  cat(sprintf("  Covariates      : %s\n",   paste(x$default_covariates,
                                                  collapse = ", ")))
  cat(sprintf("  Strata          : %s\n",   x$default_strata))
  cat(sprintf("  Alpha metrics   : %s\n",   paste(x$alpha_metrics,
                                                  collapse = ", ")))
  cat(sprintf("  p-adj method    : %s\n",   x$padj_method))
  cat("------------------------------------------------------------\n")
  invisible(x)
}
