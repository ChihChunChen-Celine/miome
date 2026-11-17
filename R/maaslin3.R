# maaslin3.R
# -----------------------------------------------------------------------------
# MaAsLin3 wrappers for the Microbiome Explorer app.
# -----------------------------------------------------------------------------

#' Check a MaAsLin3 fixed-effects formula against the metadata
#'
#' @param formula_str Fixed effects without a leading \code{~}, e.g.
#'   \code{"Treatment + Timepoint"}.
#' @param meta Metadata data frame.
#' @return Invisibly \code{TRUE}; stops with a readable message otherwise.
#' @export
validate_maaslin_formula <- function(formula_str, meta) {
  .parse_fixed_formula(formula_str, meta)
  invisible(TRUE)
}

#' Prepare feature and metadata tables for MaAsLin3
#'
#' Collapses the OTU table to \code{level} (QIIME 2 or MetaPhlAn lineages),
#' converts to relative abundance (\%), drops baseline samples and aligns the
#' metadata. Taxa not resolved down to \code{level} are kept as
#' \code{"<deepest known rank>_unclassified"} features. A table with taxon
#' names only (no lineage) is used as uploaded with \code{level = "Species"}.
#'
#' @param otu_path Path to the OTU/species table (.tsv, taxa as rows).
#' @param meta_enriched Metadata data frame (Sample IDs as row names).
#' @param level One of \code{TAXA_LEVELS}, or \code{NULL} to keep rows as-is.
#' @param config A \code{pcoa_config()} object.
#' @param min_prevalence,min_abundance Kept for API compatibility; filtering
#'   happens inside \code{run_maaslin()}.
#' @return A list: \code{features} (samples x taxa, \%), \code{metadata},
#'   \code{otu_scale_info}.
#' @export
prepare_maaslin_input <- function(otu_path, meta_enriched, level = "Species",
                                  config = pcoa_config(),
                                  min_prevalence = 0.1, min_abundance = 0) {
  mat  <- .harmonise_samples(.read_otu_table(otu_path), rownames(meta_enriched))
  meta <- meta_enriched[!.is_baseline(meta_enriched, config$is_baseline_col), , drop = FALSE]
  ids  <- intersect(colnames(mat), rownames(meta))
  if (!length(ids))
    stop("No overlapping Sample IDs between OTU table and metadata.", call. = FALSE)
  mat <- mat[, ids, drop = FALSE]

  if (!is.null(level)) {
    if (!level %in% names(TAXA_PREFIXES))
      stop(sprintf("Unknown taxonomic level '%s'. Choose from: %s", level,
                   paste(names(TAXA_PREFIXES), collapse = ", ")), call. = FALSE)
    mat <- .collapse_to_level(mat, level)
  }

  # % per sample is computed over ALL taxa. Taxa not resolved down to `level`
  # stay in as "<deepest known rank>_unclassified" features; only rows with no
  # taxonomy at all ("Unclassified") are left out of the feature list.
  features <- .to_relative_abundance(t(mat))
  scale_info <- attr(features, "scale_info")
  features <- features[, colnames(features) != "Unclassified", drop = FALSE]
  features <- features[, colSums(features) > 0, drop = FALSE]   # absent in this subset
  if (!is.null(level) && .has_lineage(colnames(features))) {
    # readable, unique names: "Blautia"; clashes get the parent rank,
    # e.g. "Blautia.uncultured_bacterium"
    toks <- lapply(strsplit(colnames(features), "[;|]"),
                   function(t) sub("^[A-Za-z]__", "", t))
    nm <- vapply(toks, function(t) t[length(t)], character(1))
    for (k in c(2, 99)) {                    # lengthen only the names that clash
      dup <- nm %in% nm[duplicated(nm)]
      nm[dup] <- vapply(toks[dup], function(t) paste(utils::tail(t, k), collapse = "."),
                        character(1))
    }
    colnames(features) <- nm
  }
  if (!ncol(features))
    stop("No taxa left to test at ", level, " level.", call. = FALSE)
  features <- as.data.frame(features, check.names = FALSE)
  meta     <- meta[rownames(features), , drop = FALSE]

  if (anyDuplicated(colnames(meta)))
    stop("Duplicate metadata column names: ",
         paste(unique(colnames(meta)[duplicated(colnames(meta))]), collapse = ", "),
         call. = FALSE)
  chr <- vapply(meta, is.character, logical(1))
  meta[chr] <- lapply(meta[chr], factor)
  single <- vapply(meta, function(x) is.factor(x) && nlevels(droplevels(x)) < 2, logical(1))
  if (any(single)) {
    message("Dropping single-level metadata columns: ",
            paste(names(meta)[single], collapse = ", "))
    meta <- meta[, !single, drop = FALSE]
  }
  message(sprintf("prepare_maaslin_input(): %d samples x %d taxa.",
                  nrow(features), ncol(features)))
  list(features = features, metadata = meta, otu_scale_info = scale_info)
}

# MaAsLin3 before 1.5.5 - including the Bioconductor 1.4.x release - stops with
# "supplied row-names must match list extent" when given more than one random
# effect. Fixed upstream in 1.5.5 ("Fix random effects handling for multiple
# crossed random intercepts", Aug 2026). One random effect works in all versions.
.MAASLIN_MULTI_RANEF <- "1.5.5"
.maaslin_multi_ranef_ok <- function()
  requireNamespace("maaslin3", quietly = TRUE) &&
    utils::packageVersion("maaslin3") >= .MAASLIN_MULTI_RANEF
.maaslin_multi_ranef_msg <- function() sprintf(paste0(
  "Your MaAsLin3 (version %s) can use only one random effect: versions before %s ",
  "stop with an error when given more. Choose one random effect, or update ",
  "MaAsLin3 with devtools::install_github(\"biobakery/maaslin3\") and restart R."),
  utils::packageVersion("maaslin3"), .MAASLIN_MULTI_RANEF)

#' Run MaAsLin3
#'
#' @param features Samples x features data frame from
#'   \code{prepare_maaslin_input()}.
#' @param metadata Metadata data frame aligned to \code{features}.
#' @param formula_str Fixed effects without \code{~}.
#' @param group_effects_var Optional column added as an ANOVA-style
#'   \code{group()} term.
#' @param random_effect Optional column(s); each becomes its own random
#'   intercept, e.g. \code{c("CatID", "Litter")} gives
#'   \code{(1|CatID) + (1|Litter)}. More than one needs MaAsLin3 1.5.5 or
#'   later; older versions get a clear error before anything is fitted.
#' @param use_random Logical; include \code{random_effect}?
#' @param reference Optional named vector/list giving the reference (control)
#'   group of categorical variables, e.g. \code{c(GROUP = "HEALTHY")}. Every
#'   other group is compared with it. Default: the alphabetically first level.
#' @param small_random_effects Passed to MaAsLin3.
#' @param normalization,transform,min_prevalence,min_abundance,max_significance
#'   Passed to \code{maaslin3::maaslin3()}. \code{min_abundance} is in percent.
#' @param padj_method Multiple-testing correction (\code{correction}).
#' @return Results data frame sorted by q-value, with attributes
#'   \code{"formula_used"} and \code{"reference_used"}.
#' @export
run_maaslin <- function(features, metadata, formula_str,
                        group_effects_var    = NULL,
                        random_effect        = NULL,
                        use_random           = FALSE,
                        reference            = NULL,
                        small_random_effects = FALSE,
                        normalization        = "TSS",
                        transform            = "LOG",
                        min_prevalence       = 0.1,
                        min_abundance        = 0,
                        padj_method          = "BH",
                        max_significance     = 0.25) {
  if (!requireNamespace("maaslin3", quietly = TRUE))
    stop("Package 'maaslin3' is required: devtools::install_github('biobakery/maaslin3')",
         call. = FALSE)

  nz <- function(x) { x <- trimws(x[!is.na(x)]); x <- unique(x[nzchar(x)]); if (length(x)) x }
  group_var   <- nz(group_effects_var)[1]
  rand_effect <- if (isTRUE(use_random)) nz(random_effect)
  if (isTRUE(use_random) && is.null(rand_effect))
    stop("'Include random effect' is ticked but no random effect column is chosen.",
         call. = FALSE)

  parsed <- .parse_fixed_formula(formula_str, metadata)
  for (x in c(group_var, rand_effect)) {
    if (!x %in% colnames(metadata))
      stop(sprintf("Column '%s' not found in metadata.", x), call. = FALSE)
    if (x %in% parsed$vars_used)
      stop(sprintf("'%s' is already in the formula; it can't also be a group/random effect.", x),
           call. = FALSE)
  }
  if (!is.null(group_var) && group_var %in% rand_effect)
    stop("Grouping factor and random effect cannot be the same column.", call. = FALSE)
  if (length(rand_effect) > 1 && !.maaslin_multi_ranef_ok())
    stop(.maaslin_multi_ranef_msg(), call. = FALSE)

  metadata <- as.data.frame(metadata[, unique(c(parsed$vars_used, group_var, rand_effect)),
                                     drop = FALSE])
  metadata <- metadata[stats::complete.cases(metadata), , drop = FALSE]
  metadata[] <- lapply(metadata, function(x) {
    if (is.character(x)) x <- factor(x)
    if (is.factor(x)) droplevels(x) else x
  })
  for (v in intersect(names(reference), colnames(metadata))) {
    ref <- as.character(reference[[v]])
    if (is.factor(metadata[[v]]) && length(ref) == 1 && ref %in% levels(metadata[[v]]))
      metadata[[v]] <- stats::relevel(metadata[[v]], ref = ref)
    else if (is.factor(metadata[[v]]))
      warning(sprintf("Reference '%s' is not a level of '%s' (after filtering); using '%s'.",
                      ref, v, levels(metadata[[v]])[1]), call. = FALSE)
  }
  ref_used <- vapply(Filter(is.factor, metadata[setdiff(colnames(metadata), rand_effect)]),
                     function(x) levels(x)[1], character(1))
  bad <- vapply(metadata, function(x) is.factor(x) && nlevels(x) < 2, logical(1))
  if (any(bad))
    stop("After removing missing values these columns have <2 levels: ",
         paste(names(metadata)[bad], collapse = ", "), call. = FALSE)
  for (re in rand_effect) {
    n_lv <- length(unique(metadata[[re]]))
    if (n_lv >= nrow(metadata))
      stop(sprintf(paste0("Random effect '%s' has a separate group for every sample - it ",
                          "needs repeated samples per group."), re), call. = FALSE)
  }

  common <- intersect(rownames(features), rownames(metadata))
  if (length(common) < 3)
    stop("Fewer than 3 samples remain after alignment/NA filtering.", call. = FALSE)
  features <- as.data.frame(features[common, , drop = FALSE], check.names = FALSE)
  metadata <- metadata[common, , drop = FALSE]

  rhs <- trimws(formula_str)
  if (!is.null(group_var))   rhs <- paste0(rhs, " + group(", group_var, ")")
  if (!is.null(rand_effect))
    rhs <- paste0(rhs, " + ", paste0("(1|", rand_effect, ")", collapse = " + "))
  full_formula <- paste("~", rhs)
  message("MaAsLin3: ", full_formula, " | ", nrow(features), " samples x ",
          ncol(features), " features")

  out_dir <- tempfile("maaslin3_")
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  maaslin3::maaslin3(
    input_data = features, input_metadata = metadata, output = out_dir,
    formula = full_formula, small_random_effects = small_random_effects,
    normalization = normalization, transform = transform,
    min_prevalence = min_prevalence, min_abundance = min_abundance,
    max_significance = max_significance, correction = padj_method,
    plot_summary_plot = FALSE, plot_associations = FALSE, save_models = FALSE)

  res_path <- file.path(out_dir, "all_results.tsv")
  if (!file.exists(res_path)) stop("MaAsLin3 produced no results file.", call. = FALSE)
  res <- utils::read.delim(res_path, stringsAsFactors = FALSE, check.names = FALSE,
                           comment.char = "")   # MaAsLin3 quotes its fields
  if (!nrow(res)) stop("MaAsLin3 returned an empty results table.", call. = FALSE)

  res <- .standardise_maaslin_cols(res)
  res <- res[order(res$qval, res$pval, na.last = TRUE), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "formula_used") <- full_formula
  attr(res, "reference_used") <- ref_used
  res
}

# -- Internal: harmonise result column names across MaAsLin versions ----------
.standardise_maaslin_cols <- function(df) {
  rename_map <- list(
    feature  = c("feature", "Feature", "taxon", "Taxon"),
    metadata = c("metadata", "Metadata", "variable", "Variable"),
    value    = c("value", "Value", "level", "Level"),
    coef     = c("coef", "Coef", "coefficient", "Coefficient", "estimate", "Estimate"),
    stderr   = c("stderr", "StdErr", "se", "SE"),
    pval     = c("pval", "pval_individual", "pValue", "p.value", "P.value", "p_value"),
    qval     = c("qval", "qval_individual", "qValue", "q.value", "Q.value",
                 "q_value", "padj", "p_adjusted"))
  for (new in names(rename_map)) {
    if (new %in% colnames(df)) next
    old <- intersect(rename_map[[new]], colnames(df))[1]
    if (!is.na(old)) colnames(df)[colnames(df) == old] <- new else df[[new]] <- NA
  }
  df
}

#' Significant MaAsLin3 results
#'
#' @param results Data frame from \code{run_maaslin()}.
#' @param threshold q-value cut-off.
#' @return Rows with \code{qval < threshold}, or a one-row message.
#' @export
summarise_maaslin_results <- function(results, threshold = 0.25) {
  if (is.null(results) || !nrow(results))
    return(data.frame(message = "No results available."))
  sig <- results[!is.na(results$qval) & results$qval < threshold, , drop = FALSE]
  if (!nrow(sig))
    return(data.frame(message = sprintf("No features significant at q < %.2f.", threshold)))
  num <- intersect(c("coef", "stderr", "pval", "qval"), colnames(sig))
  sig[num] <- lapply(sig[num], round, 6)
  sig
}
