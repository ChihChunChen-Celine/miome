# -----------------------------------------------------------------------------
# config.R
# Central configuration for the Microbiome Explorer
# Covers: PCoA / Beta Diversity / Alpha Diversity / MaAsLin3
# -----------------------------------------------------------------------------

# -- Palette options -----------------------------------------------------------

#' Categorical colour palette options (for grouping variables)
#'
#' \code{"Miome"} is the package's own palette (see \code{MIOME_PALETTE});
#' the others are RColorBrewer palettes.
#' @export
CAT_PALETTES <- c("Miome", "Set1", "Set2", "Set3", "Dark2", "Paired", "Accent")

#' The miome categorical palette
#'
#' Eight hues in a fixed order, chosen so that neighbouring colours stay
#' distinguishable for readers with colour-vision deficiency.
#' @export
MIOME_PALETTE <- c("#2F6FB0", "#DAA22A", "#B8456B", "#22AE93",
                   "#7A4FA6", "#95BD4C", "#CF6A3C", "#96A3EA")

#' Continuous colour palette options (for numeric variables)
#' @export
CONT_PALETTES <- c("RdBu", "RdYlBu", "Spectral", "YlOrRd", "Blues")

#' Supported beta diversity distance methods
#'
#' Most are passed straight to \code{vegan::vegdist()}. \code{"aitchison"} is
#' computed as a CLR transform + Euclidean distance. The two UniFrac options
#' require a precomputed distance matrix (QIIME 2 \code{.qza} or plain
#' \code{.tsv}) uploaded by the user; no phylogenetic tree is used.
#' @export
DIST_METHODS <- c(
  "Bray-Curtis"          = "bray",
  "Jaccard"              = "jaccard",
  "Euclidean"            = "euclidean",
  "Aitchison"            = "aitchison",
  "UniFrac (weighted)"   = "unifrac_weighted",
  "UniFrac (unweighted)" = "unifrac_unweighted",
  "Manhattan"            = "manhattan",
  "Canberra"             = "canberra",
  "Kulczynski"           = "kulczynski"
)

#' Supported alpha diversity metrics
#' @export
ALPHA_METRICS <- c("Shannon", "Simpson", "Observed_taxa", "Chao1", "Pielou")

#' Taxonomic levels available for taxonomy plots and MaAsLin3
#'
#' Maps a rank name to its lineage prefix letter (e.g. \code{"f"} for
#' \code{"f__Ruminococcaceae"}). Lineages may use \code{";"} (QIIME 2 / SILVA)
#' or \code{"|"} (MetaPhlAn) separators; MetaPhlAn's \code{k__} is accepted as
#' an alias for \code{d__}.
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
#' Column-name arguments may be \code{""} when your metadata has no such
#' column; every function then simply skips the step that needs it.
#'
#' @param subject_id_col Subject / animal ID column (random effect, strata,
#'   paired lines).
#' @param treatment_col Treatment column name.
#' @param timepoint_col Timepoint column name (also orders paired lines).
#' @param sequence_col Sequence column name (crossover designs).
#' @param is_baseline_col Baseline flag column (TRUE/FALSE, Yes/No or 1/0).
#' @param complete_subjects_only Logical; leave out every subject that has
#'   fewer samples than the others (e.g. a drop-out with no final sample).
#' @param study_design Optional label for your records, e.g.
#'   \code{"crossover"}, \code{"longitudinal"} or \code{"cross-sectional"}.
#'   Not used by the app.
#' @param dist_method Distance method key, one of \code{DIST_METHODS}.
#' @param dist_method_label Human-readable distance label for plot titles.
#' @param default_covariates Default PERMANOVA covariates (Shiny UI).
#' @param default_strata Default strata column (Shiny UI); \code{NULL} = none.
#' @param padj_method P-value adjustment method, see \code{PADJ_METHODS}.
#' @param alpha_metrics Alpha diversity metrics, see \code{ALPHA_METRICS}.
#' @param alpha_export_w_in,alpha_export_h_in Alpha plot export size (inches).
#' @param plot_title Title shown on PCoA plots.
#' @param marker_size,marker_opacity Point aesthetics.
#' @param export_dpi,export_w_in,export_h_in 2D TIFF export settings.
#' @param export_3d_w_px,export_3d_h_px 3D PNG/TIFF snapshot size (pixels).
#' @return A list of class \code{pcoa_config}.
#' @export
pcoa_config <- function(
    subject_id_col     = "sampleid",
    treatment_col      = "Treatment",
    timepoint_col      = "Timepoint",
    sequence_col       = "Sequence",
    is_baseline_col    = "IsBaseline",
    complete_subjects_only = FALSE,
    study_design       = "crossover",
    dist_method        = "bray",
    dist_method_label  = "Bray-Curtis",
    default_covariates = c("Timepoint", "Sequence"),
    default_strata     = "sampleid",
    padj_method        = "BH",
    alpha_metrics      = ALPHA_METRICS,
    alpha_export_w_in  = 8,
    alpha_export_h_in  = 6,
    plot_title         = "PCoA",
    marker_size        = 6,
    marker_opacity     = 0.85,
    export_dpi         = 300,
    export_w_in        = 9,
    export_h_in        = 7,
    export_3d_w_px     = 2400,
    export_3d_h_px     = 1800) {

  cols <- list(subject_id_col = subject_id_col, treatment_col = treatment_col,
               timepoint_col = timepoint_col, sequence_col = sequence_col,
               is_baseline_col = is_baseline_col)
  for (nm in names(cols))
    if (!is.character(cols[[nm]]) || length(cols[[nm]]) != 1)
      stop(sprintf("'%s' must be a single character string (\"\" if absent).", nm),
           call. = FALSE)

  bad <- setdiff(alpha_metrics, ALPHA_METRICS)
  if (length(bad))
    warning(sprintf("Unknown alpha metric(s) ignored: %s. Choose from: %s",
                    paste(bad, collapse = ", "),
                    paste(ALPHA_METRICS, collapse = ", ")), call. = FALSE)

  if (!padj_method %in% stats::p.adjust.methods)
    stop("'padj_method' must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), call. = FALSE)

  structure(
    c(cols, list(
      complete_subjects_only = isTRUE(complete_subjects_only),
      study_design       = study_design,
      dist_method        = dist_method,
      dist_method_label  = dist_method_label,
      default_covariates = default_covariates,
      default_strata     = default_strata,
      padj_method        = padj_method,
      alpha_metrics      = intersect(alpha_metrics, ALPHA_METRICS),
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
    )),
    class = "pcoa_config"
  )
}

#' @export
print.pcoa_config <- function(x, ...) {
  show <- function(lbl, v) cat(sprintf("  %-16s: %s\n", lbl,
    if (length(v) && any(nzchar(v))) paste(v, collapse = ", ") else "(none)"))
  cat("-- Microbiome Explorer Configuration ----------------------\n")
  show("Subject ID col",  x$subject_id_col)
  show("Treatment col",   x$treatment_col)
  show("Timepoint col",   x$timepoint_col)
  show("Baseline col",    x$is_baseline_col)
  show("Study design",    x$study_design)
  show("Distance method", x$dist_method_label)
  show("Covariates",      x$default_covariates)
  show("Strata",          x$default_strata)
  show("Alpha metrics",   x$alpha_metrics)
  show("p-adj method",    x$padj_method)
  cat("------------------------------------------------------------\n")
  invisible(x)
}
