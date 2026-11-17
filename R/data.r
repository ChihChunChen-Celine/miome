# data.R
# -----------------------------------------------------------------------------
# Data loading, alpha diversity, baseline covariates, PCoA and taxonomy
# bar-plot data. File reading, metadata cleaning, counts-vs-percent detection
# and taxonomy parsing all live in utils.R and are shared with every other
# file in the package.
# -----------------------------------------------------------------------------

# -- Internal: read a precomputed distance matrix (.qza or .tsv) ---------------
.read_precomputed_dist <- function(path) {
  if (grepl("\\.qza$", path, ignore.case = TRUE)) {
    tmp_dir <- tempfile("qza_")
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    utils::unzip(path, exdir = tmp_dir)
    path <- list.files(tmp_dir, pattern = "distance-matrix\\.tsv$",
                       recursive = TRUE, full.names = TRUE)[1]
    if (is.na(path)) stop("Could not find distance-matrix.tsv inside the .qza file.")
  }
  mat <- as.matrix(utils::read.table(path, sep = "\t", header = TRUE,
                                     row.names = 1, check.names = FALSE,
                                     quote = "", comment.char = ""))
  stats::as.dist(mat)
}

# -- Internal: distance matrix from a samples x taxa matrix --------------------
.compute_dist <- function(otu_mat, method, precomputed_dist_path = NULL) {

  # UniFrac: values come from the uploaded matrix; otu_mat only supplies IDs.
  if (method %in% c("unifrac_weighted", "unifrac_unweighted")) {
    if (is.null(precomputed_dist_path))
      stop("A precomputed distance matrix (.qza or .tsv) is required for ",
           "UniFrac. Please upload one in the sidebar.", call. = FALSE)
    d_mat  <- as.matrix(.read_precomputed_dist(precomputed_dist_path))
    shared <- intersect(rownames(otu_mat), rownames(d_mat))
    if (!length(shared))
      stop("No overlapping Sample IDs between the OTU table and the distance ",
           "matrix.\n  OTU table IDs (first 5): ",
           paste(head(rownames(otu_mat), 5), collapse = ", "),
           "\n  Distance matrix IDs (first 5): ",
           paste(head(rownames(d_mat), 5), collapse = ", "), call. = FALSE)
    dropped <- setdiff(rownames(otu_mat), shared)
    if (length(dropped))
      warning(sprintf("%d sample(s) absent from the distance matrix were dropped: %s",
                      length(dropped), paste(dropped, collapse = ", ")),
              call. = FALSE)
    d <- stats::as.dist(d_mat[shared, shared])
    attr(d, "scale_info") <- list(note = paste(
      "UniFrac uses a precomputed distance matrix;",
      "OTU table scale is not applicable."))
    return(d)
  }

  # Everything else is computed on relative abundance (%).
  otu_mat    <- .to_relative_abundance(otu_mat)
  scale_info <- attr(otu_mat, "scale_info")

  d <- if (method == "aitchison") {
    # Pseudo-count = half the smallest non-zero relative abundance observed.
    nz     <- otu_mat[otu_mat > 0]
    pseudo <- if (length(nz)) min(nz) / 2 else 0.5
    clr    <- log(replace(otu_mat, otu_mat == 0, pseudo))
    stats::dist(clr - rowMeans(clr))
  } else {
    vegan::vegdist(otu_mat, method = method)
  }
  attr(d, "scale_info") <- scale_info
  d
}

# -- Internal: load and align OTU table + metadata (one read of each file) ----
# Returns otu = samples x features (deepest rank only for all-ranks tables).
.load_otu_meta <- function(otu_path, meta_path, is_baseline_col = "IsBaseline") {
  mat  <- .read_otu_table(otu_path)
  meta <- .read_meta(meta_path, is_baseline_col, sample_ids = colnames(mat))
  mat  <- .harmonise_samples(mat, rownames(meta))
  otu  <- t(.feature_rows(mat))

  shared <- intersect(rownames(otu), rownames(meta))
  if (!length(shared))
    stop("No matching Sample IDs between OTU table and metadata.\n",
         "  OTU table sample names (first 5): ",
         paste(head(rownames(otu), 5), collapse = ", "), "\n",
         "  Metadata row names     (first 5): ",
         paste(head(rownames(meta), 5), collapse = ", "), call. = FALSE)
  # Which samples are analysed: only those present in BOTH files. Analysing a
  # subset of a sequencing run is normal, so this is reported, not an error.
  match_info <- list(
    n_table = nrow(otu), n_meta = nrow(meta), n_matched = length(shared),
    table_only = setdiff(rownames(otu), rownames(meta)),
    meta_only  = setdiff(rownames(meta), rownames(otu)))
  message(sprintf("%d sample(s) analysed: in both the table (%d) and the metadata (%d).",
                  length(shared), nrow(otu), nrow(meta)))
  if (length(match_info$table_only))
    message(sprintf("  %d sequenced sample(s) not in the metadata were left out: %s",
                    length(match_info$table_only),
                    paste(match_info$table_only, collapse = ", ")))
  if (length(match_info$meta_only))
    warning(sprintf("%d metadata sample(s) have no match in the abundance table and were dropped: %s",
                    length(match_info$meta_only),
                    paste(match_info$meta_only, collapse = ", ")), call. = FALSE)

  otu <- otu[shared, , drop = FALSE]
  list(otu = otu[, colSums(otu) > 0, drop = FALSE],     # taxa absent from this subset
       meta = meta[shared, , drop = FALSE],
       match_info = match_info,
       tax_levels = .available_levels(rownames(mat)))
}

# -- Internal: subjects that lack some of their samples ---------------------------
# An animal is "incomplete" when it has fewer analysed samples than the animals
# with the most (e.g. a drop-out without a final sample). They are always
# reported; they are left out of everything only when the configuration asks.
.handle_incomplete_subjects <- function(dat, config) {
  subj <- config$subject_id_col
  mi   <- dat$match_info
  # integer(0), not NULL: with the element absent, mi$incomplete would partially
  # match mi$incomplete_dropped
  mi$incomplete <- integer(0); mi$incomplete_dropped <- FALSE
  mi$subject_unique <- nzchar(subj) && subj %in% colnames(dat$meta) &&
    !anyDuplicated(stats::na.omit(dat$meta[[subj]]))
  # No subject column named? Look for one, so repeated samples are not silently
  # treated as independent (name-based and conservative - see .subject_candidates).
  mi$subject_candidates <- if (!nzchar(subj) || !subj %in% colnames(dat$meta))
    setdiff(.subject_candidates(dat$meta), colnames(dat$meta)[0]) else character(0)
  if (length(mi$subject_candidates))
    message("Column(s) that look like an animal/subject ID with repeated samples: ",
            paste(mi$subject_candidates, collapse = ", "),
            ". Set one in Step 2 if samples come from the same animal more than once.")
  if (nzchar(subj) && subj %in% colnames(dat$meta)) {
    n <- table(dat$meta[[subj]])
    if (length(n) > 1 && max(n) > 1 && any(n < max(n))) {
      inc <- n[n < max(n)]
      mi$incomplete <- stats::setNames(as.integer(inc), names(inc))
      mi$expected_n <- as.integer(max(n))
      txt <- paste(sprintf("%s (%d of %d)", names(inc), as.integer(inc), max(n)), collapse = ", ")
      if (isTRUE(config$complete_subjects_only)) {
        keep     <- !dat$meta[[subj]] %in% names(inc)
        dat$meta <- dat$meta[keep, , drop = FALSE]
        dat$otu  <- dat$otu[keep, , drop = FALSE]
        dat$otu  <- dat$otu[, colSums(dat$otu) > 0, drop = FALSE]
        mi$incomplete_dropped <- TRUE
        mi$n_matched <- nrow(dat$meta)
        message("Left out (missing samples): ", txt)
      } else {
        message("Subject(s) with missing samples are still included: ", txt)
      }
    }
  }
  dat$match_info <- mi
  dat
}

# -- Internal: alpha metrics for a samples x features matrix ------------------
.alpha_table <- function(otu) {
  # Chao1 needs integer counts; vegan errors otherwise -> NA.
  chao1 <- tryCatch(vegan::estimateR(otu)["S.chao1", ],
                    error = function(e) rep(NA_real_, nrow(otu)))
  observed <- rowSums(otu > 0)
  shannon  <- vegan::diversity(otu, index = "shannon")
  data.frame(
    Shannon       = shannon,
    Simpson       = vegan::diversity(otu, index = "simpson"),
    Observed_taxa = observed,
    Chao1         = unname(chao1),
    Pielou        = ifelse(observed > 1, shannon / log(observed), NA_real_),
    row.names     = rownames(otu)
  )[, ALPHA_METRICS, drop = FALSE]
}

# -- Internal: baseline covariates from in-memory data -------------------------
# meta must already contain the alpha columns; d is a dist over (at least) the
# baseline samples. Adds Baseline_<metric> for EVERY alpha metric present plus
# Baseline_PC1 / Baseline_PC2 to post-baseline rows; baseline rows get NA.
.baseline_covariates <- function(meta, d, config) {
  subj <- config$subject_id_col
  base <- config$is_baseline_col
  miss <- setdiff(c(base, subj), colnames(meta))
  if (length(miss))
    stop(sprintf("Missing required column(s) in metadata: %s",
                 paste(miss, collapse = ", ")), call. = FALSE)

  is_base <- .is_baseline(meta, base)
  ids     <- intersect(rownames(meta)[is_base], labels(d))
  if (length(ids) < 2)
    stop("Need at least 2 baseline samples to compute baseline covariates.",
         call. = FALSE)

  metrics <- intersect(ALPHA_METRICS, colnames(meta))
  b <- meta[ids, metrics, drop = FALSE]
  colnames(b) <- paste0("Baseline_", metrics)

  pc <- stats::cmdscale(stats::as.dist(as.matrix(d)[ids, ids]),
                        k = min(length(ids) - 1, 2))
  b$Baseline_PC1 <- pc[, 1]
  b$Baseline_PC2 <- if (ncol(pc) >= 2) pc[, 2] else NA_real_
  b[[subj]] <- meta[ids, subj]

  dup <- unique(b[[subj]][duplicated(b[[subj]])])
  if (length(dup))
    warning(sprintf(paste0("Subject(s) with more than one baseline sample - ",
                           "only the first is used: %s"),
                    paste(dup, collapse = ", ")), call. = FALSE)
  b <- b[!duplicated(b[[subj]]), , drop = FALSE]

  cov_cols <- setdiff(colnames(b), subj)
  meta[intersect(cov_cols, colnames(meta))] <- NULL     # re-upload safe
  idx <- match(meta[[subj]], b[[subj]])
  idx[is_base] <- NA                                     # baseline rows -> NA
  meta[cov_cols] <- b[idx, cov_cols, drop = FALSE]

  message(sprintf("Baseline covariates computed for %d subject(s).", nrow(b)))
  meta
}

# -- Internal: PCoA on a dist object -------------------------------------------
.pcoa_from_dist <- function(d, meta) {
  ids <- labels(d)
  if (length(ids) < 3) stop("Need at least 3 samples for PCoA.", call. = FALSE)
  res <- stats::cmdscale(d, k = min(length(ids) - 1, 10), eig = TRUE)
  eig <- pmax(res$eig, 0)
  pts <- res$points
  pcoa_df <- data.frame(
    SampleID = ids,
    PC1 = pts[, 1], PC2 = pts[, 2],
    PC3 = if (ncol(pts) >= 3) pts[, 3] else 0,
    stringsAsFactors = FALSE
  )
  pcoa_df <- cbind(pcoa_df, meta[ids, setdiff(colnames(meta),
                                              colnames(pcoa_df)), drop = FALSE])
  rownames(pcoa_df) <- NULL
  list(pcoa_df = pcoa_df, var_pct = round(eig / sum(eig) * 100, 1))
}

.new_pcoa_obj <- function(d, meta, config, otu_path, meta_path, dat = NULL) {
  meta <- meta[labels(d), , drop = FALSE]
  p    <- .pcoa_from_dist(d, meta)
  structure(list(
    pcoa_df        = p$pcoa_df,
    var_pct        = p$var_pct,
    meta_cols      = colnames(meta),
    dist           = d,
    dist_mat       = as.matrix(d),
    meta           = meta,
    meta_enriched  = meta,
    config         = config,
    otu_path       = otu_path,
    meta_path      = meta_path,
    otu_scale_info = attr(d, "scale_info"),
    match_info     = dat$match_info,
    tax_levels     = if (is.null(dat$tax_levels)) names(TAXA_PREFIXES) else dat$tax_levels
  ), class = "pcoa_obj")
}

# -- 0. prepare_metadata() ------------------------------------------------------
#' Load the metadata on its own (no sequencing table needed)
#'
#' Used for analyses that do not depend on sequencing, such as targeted
#' measurements: every metadata row is kept, whether or not that sample was
#' sequenced. The baseline flag is parsed and, when
#' \code{config$complete_subjects_only} is set, subjects with fewer rows than
#' the others are left out.
#'
#' @param meta_path Path to the metadata (.csv), or a data frame with Sample
#'   IDs as row names.
#' @param config A \code{pcoa_config()} object.
#' @param sample_ids Optional sample names of an abundance table, used only to
#'   find the Sample ID column.
#' @return A list: \code{meta} (data frame) and \code{info} (row count,
#'   incomplete subjects, subject-ID candidates).
#' @export
prepare_metadata <- function(meta_path, config = pcoa_config(), sample_ids = NULL) {
  meta <- .read_meta(meta_path, config$is_baseline_col, sample_ids = sample_ids)
  dat  <- list(meta = meta, otu = matrix(1, nrow(meta), 1, dimnames = list(rownames(meta), "x")),
               match_info = list(n_meta = nrow(meta), n_matched = nrow(meta)))
  dat  <- suppressMessages(.handle_incomplete_subjects(dat, config))
  list(meta = dat$meta, info = dat$match_info)
}

# -- 1. compute_alpha_diversity() ---------------------------------------------
#' Compute alpha diversity metrics and append them to the metadata
#'
#' @param otu_path  Path to OTU/species table (.tsv, taxa as rows, samples as
#'   columns). QIIME 2, MetaPhlAn and plain tables are accepted.
#' @param meta_path Path to mapping file (.csv, first column = Sample ID), or a
#'   data frame with Sample IDs as row names.
#' @param config    A \code{pcoa_config()} object.
#' @return Metadata data frame with Shannon, Simpson, Observed_taxa, Chao1
#'   and Pielou columns added. Chao1 is \code{NA} for non-count data.
#' @export
compute_alpha_diversity <- function(otu_path, meta_path, config = pcoa_config()) {
  dat <- .load_otu_meta(otu_path, meta_path, config$is_baseline_col)
  cbind(dat$meta[setdiff(colnames(dat$meta), ALPHA_METRICS)],
        .alpha_table(dat$otu)[rownames(dat$meta), , drop = FALSE])
}

# -- 2. prepare_baseline_covariates() -----------------------------------------
#' Compute and merge baseline covariates into post-baseline metadata
#'
#' For every alpha metric adds \code{Baseline_<metric>} (e.g.
#' \code{Baseline_Shannon}, \code{Baseline_Observed_taxa}) and the subject's
#' baseline PCoA coordinates (\code{Baseline_PC1}, \code{Baseline_PC2}) to its
#' post-baseline rows. Baseline rows get \code{NA}.
#'
#' @inheritParams compute_alpha_diversity
#' @param dist_method Distance method for the baseline PCoA.
#' @param precomputed_dist_path Optional precomputed distance matrix (.qza or
#'   .tsv); required for UniFrac.
#' @return Metadata data frame with the covariate columns added.
#' @export
prepare_baseline_covariates <- function(otu_path, meta_path,
                                        dist_method           = "bray",
                                        config                = pcoa_config(),
                                        precomputed_dist_path = NULL) {
  dat  <- .load_otu_meta(otu_path, meta_path, config$is_baseline_col)
  meta <- dat$meta
  if (!any(ALPHA_METRICS %in% colnames(meta)))
    meta <- cbind(meta, .alpha_table(dat$otu)[rownames(meta), , drop = FALSE])
  base <- .is_baseline(meta, config$is_baseline_col)
  d    <- .compute_dist(dat$otu[base, , drop = FALSE], dist_method,
                        precomputed_dist_path)
  .baseline_covariates(meta, d, config)
}

# -- 3. prepare_pcoa() --------------------------------------------------------
#' Load OTU table and metadata and compute a PCoA
#'
#' @inheritParams prepare_baseline_covariates
#' @return A \code{pcoa_obj} list: \code{pcoa_df}, \code{var_pct},
#'   \code{meta_cols}, \code{dist}, \code{dist_mat}, \code{meta},
#'   \code{meta_enriched}, \code{config}, \code{otu_path}, \code{meta_path},
#'   \code{otu_scale_info}, \code{match_info} (which samples were matched)
#'   and \code{tax_levels} (taxonomic levels the table supports).
#' @export
prepare_pcoa <- function(otu_path, meta_path,
                         config                = pcoa_config(),
                         dist_method           = "bray",
                         precomputed_dist_path = NULL) {
  dat <- .load_otu_meta(otu_path, meta_path, config$is_baseline_col)
  d   <- .compute_dist(dat$otu, dist_method, precomputed_dist_path)
  .new_pcoa_obj(d, dat$meta, config, otu_path, meta_path, dat)
}

# -- 4. prepare_all() - master pipeline ---------------------------------------
#' Master data preparation pipeline
#'
#' Reads the files once, then computes alpha diversity, baseline covariates
#' (when the baseline column exists) and the full PCoA on all samples.
#'
#' @inheritParams prepare_baseline_covariates
#' @return A \code{pcoa_obj} whose \code{$meta} / \code{$meta_enriched} hold
#'   the alpha metrics and baseline covariates.
#' @export
prepare_all <- function(otu_path, meta_path,
                        config                = pcoa_config(),
                        dist_method           = "bray",
                        precomputed_dist_path = NULL) {

  dat  <- .load_otu_meta(otu_path, meta_path, config$is_baseline_col)
  dat  <- .handle_incomplete_subjects(dat, config)
  meta <- cbind(dat$meta[setdiff(colnames(dat$meta), ALPHA_METRICS)],
                .alpha_table(dat$otu)[rownames(dat$meta), , drop = FALSE])
  message(sprintf("Alpha diversity computed for %d samples.", nrow(meta)))

  d <- tryCatch(
    .compute_dist(dat$otu, dist_method, precomputed_dist_path),
    error = function(e)
      stop("Distance computation failed: ", conditionMessage(e), call. = FALSE))

  base_col <- config$is_baseline_col
  if (base_col %in% colnames(meta)) {
    meta <- tryCatch(.baseline_covariates(meta, d, config), error = function(e) {
      warning("Baseline covariates skipped: ", conditionMessage(e), call. = FALSE)
      meta
    })
  } else if (nzchar(base_col)) {
    message(sprintf("No '%s' column - baseline covariates skipped.", base_col))
  }

  obj <- .new_pcoa_obj(d, meta, config, otu_path, meta_path, dat)
  message(sprintf(
    "prepare_all() complete: %d samples (%d baseline), %d metadata columns, distance = %s",
    nrow(obj$meta), sum(.is_baseline(obj$meta, base_col)),
    ncol(obj$meta), dist_method))
  obj
}

# -- 5. compute_taxonomy_barplot_data() ---------------------------------------
#' Collapse the OTU table to one taxonomic level for a stacked bar plot
#'
#' Works with QIIME 2 / SILVA (\code{d__...;p__...}) and MetaPhlAn
#' (\code{k__...|p__...}) lineages, including MetaPhlAn tables that list every
#' rank (each level is then read from its own rows, so nothing is double
#' counted). Taxa not resolved down to \code{level} are shown as
#' \code{"<deepest known rank>_unclassified"}. A table with taxon names only
#' (no lineage) can be shown at \code{level = "Species"}, i.e. as uploaded.
#'
#' @inheritParams compute_alpha_diversity
#' @param level     One of \code{TAXA_LEVELS} (e.g. \code{"Family"}).
#' @param group_var Metadata column used for faceting (kept for API
#'   compatibility; faceting happens in \code{build_taxon_barplot()}).
#' @param top_n     Number of taxa shown individually; the rest become
#'   \code{"Others"}.
#' @return Long data frame: \code{SampleID}, \code{Taxon},
#'   \code{RelAbundance} (\%), plus all metadata columns.
#' @export
compute_taxonomy_barplot_data <- function(otu_path, meta_path,
                                          level     = "Family",
                                          group_var = NULL,
                                          top_n     = TAXA_DEFAULT_TOP_N,
                                          config    = pcoa_config()) {
  mat  <- .read_otu_table(otu_path)
  meta <- .read_meta(meta_path, config$is_baseline_col, sample_ids = colnames(mat))
  mat  <- .harmonise_samples(mat, rownames(meta))
  ids  <- intersect(colnames(mat), rownames(meta))
  if (!length(ids))
    stop("No matching Sample IDs between OTU table columns and metadata rows.",
         call. = FALSE)
  # metadata columns that would clash with the plot's own columns
  clash <- colnames(meta) %in% c("SampleID", "Taxon", "RelAbundance")
  colnames(meta)[clash] <- paste0(colnames(meta)[clash], ".meta")

  lvl <- .collapse_to_level(mat[, ids, drop = FALSE], level)
  tot <- colSums(lvl); tot[tot == 0] <- 1
  rel <- sweep(lvl, 2, tot, "/") * 100
  ranked <- names(sort(rowSums(rel), decreasing = TRUE))
  top    <- utils::head(ranked, top_n)
  lab    <- stats::setNames(.short_taxon(rownames(rel)), rownames(rel))
  rel    <- rowsum(rel, ifelse(rownames(rel) %in% top, lab[rownames(rel)], "Others"))

  df <- data.frame(SampleID     = rep(colnames(rel), each = nrow(rel)),
                   Taxon        = rep(rownames(rel), times = ncol(rel)),
                   RelAbundance = as.vector(rel),
                   stringsAsFactors = FALSE)
  meta$SampleID <- rownames(meta)
  df <- merge(df, meta[ids, , drop = FALSE], by = "SampleID", all.x = TRUE)
  attr(df, "top_n_used") <- top_n
  attr(df, "level_used") <- level
  df
}
