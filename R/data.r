# data.R
# -----------------------------------------------------------------------------
# Data loading, alpha diversity, baseline covariate preparation, and PCoA.
# Phylogenetic tree support has been removed; UniFrac is supported via a
# precomputed distance matrix (.qza artifact or plain .tsv).
# -----------------------------------------------------------------------------

# -- Internal helper: read a precomputed distance matrix -----------------------
# Accepts either a QIIME 2 .qza artifact (a renamed ZIP containing
# data/distance-matrix.tsv) or a plain tab-separated matrix file.
.read_precomputed_dist <- function(path) {
  if (grepl("\\.qza$", path, ignore.case = TRUE)) {
    tmp_dir <- tempfile("qza_")
    dir.create(tmp_dir)
    utils::unzip(path, exdir = tmp_dir)
    tsv_path <- list.files(
      tmp_dir,
      pattern   = "distance-matrix\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(tsv_path) == 0)
      stop("Could not find data/distance-matrix.tsv inside: ", path)
    tsv_path <- tsv_path[1]
  } else {
    # Plain .tsv - use as-is
    tsv_path <- path
  }

  mat <- as.matrix(
    read.table(
      tsv_path,
      sep          = "\t",
      header       = TRUE,
      row.names    = 1,
      check.names  = FALSE,
      quote        = "",
      comment.char = ""
    )
  )
  stats::as.dist(mat)
}

# -- Internal helper: detect whether an OTU table is raw counts or percent ----
# Used only to choose an appropriate pseudo-count for the Aitchison/CLR
# transform below - vegan::vegdist() and the alpha diversity metrics don't
# need this, since they're scale-invariant (or, for Chao1, require counts
# and simply return NA if given anything else - see compute_alpha_diversity()).
#
# Detection is per-sample (per row here - otu_mat is samples-as-rows,
# taxa-as-columns throughout data.R), checked against the WHOLE table: if
# most sample rows sum to ~100, the table is treated as percent; otherwise
# it's treated as raw counts. This mirrors the per-column check used in the
# standalone summarize_abundance.R script, just transposed to this table's
# orientation.
.detect_otu_scale <- function(otu_mat, tol = 1) {
  row_sums <- rowSums(otu_mat, na.rm = TRUE)
  if (length(row_sums) == 0 || all(row_sums == 0)) return("counts")
  pct_frac <- mean(abs(row_sums - 100) <= tol)
  if (pct_frac >= 0.5) "percent" else "counts"
}

# -- Internal helper: convert an OTU table to relative abundance (%) ----------
# Applied before every beta-diversity distance calculation (.compute_dist(),
# except the UniFrac branch, which doesn't use otu_mat's values) and before
# MaAsLin3 input prep (maaslin3.R's prepare_maaslin_input()), so both
# consistently work on percent abundance regardless of whether the uploaded
# table was raw counts or already relative abundance.
#
# Detection + conversion is per-SAMPLE (per row here - otu_mat is
# samples-as-rows, taxa-as-columns, matching the orientation used
# throughout data.R and the MaAsLin3 features table): a sample row whose
# sum is already within `tol` of 100 is left untouched (already percent);
# any other row is treated as raw counts and divided by its own row sum,
# then multiplied by 100. Per-sample (not whole-table) detection means a
# table mixing already-normalized and raw-count samples is still handled
# correctly rather than assuming the whole table is one or the other.
#
# NOTE: this is deliberately NOT applied to compute_alpha_diversity().
# Chao1 requires genuine read counts and would silently return NA if fed
# percent data (see compute_alpha_diversity()), so alpha diversity keeps
# working on whatever scale was actually uploaded.
.to_relative_abundance <- function(otu_mat, tol = 1) {
  row_sums <- rowSums(otu_mat, na.rm = TRUE)

  zero_rows <- rownames(otu_mat)[row_sums == 0]
  if (length(zero_rows) > 0)
    warning(sprintf(
      "Sample(s) with a total abundance of 0 - left as 0 after conversion: %s",
      paste(zero_rows, collapse = ", ")
    ))

  is_percent  <- abs(row_sums - 100) <= tol
  n_total     <- length(row_sums)
  n_percent   <- sum(is_percent)
  convert_rows <- rownames(otu_mat)[!is_percent & row_sums > 0]
  n_converted  <- length(convert_rows)

  # message() only reaches the R/Shiny SERVER console or log - never the
  # browser UI, so it's invisible on a deployed app. Attaching the same
  # detection result as an attribute lets callers (prepare_pcoa(),
  # prepare_maaslin_input()) pass it up into the pcoa_obj / prepped-input
  # list, where app.R can render it as an actual on-screen info box - see
  # output$otu_scale_info_ui and the MaAsLin3 info panel in app.R.
  message(sprintf(
    paste0(
      "Relative-abundance check: %d of %d sample(s) already look like ",
      "percent (row sum ~100, tol=%.1f) and were left as-is; %d ",
      "sample(s) looked like raw counts and were converted to percent."
    ),
    n_percent, n_total, tol, n_converted
  ))

  out <- otu_mat
  if (n_converted > 0) {
    out[convert_rows, ] <- sweep(
      otu_mat[convert_rows, , drop = FALSE],
      MARGIN = 1,
      STATS  = row_sums[convert_rows],
      FUN    = "/"
    ) * 100
  }

  attr(out, "scale_info") <- list(
    n_total           = n_total,
    n_already_percent = n_percent,
    n_converted       = n_converted,
    converted_samples = convert_rows,
    tol               = tol
  )
  out
}

# -- Internal helper: compute distance matrix ----------------------------------
# Aitchison (CLR + Euclidean) and standard vegan metrics are computed
# directly from the OTU table.  UniFrac metrics require a precomputed
# distance matrix supplied via precomputed_dist_path; tree support has
# been removed entirely.
.compute_dist <- function(otu_mat,
                          method,
                          precomputed_dist_path = NULL) {

  # -- UniFrac: must come from a precomputed matrix --------------
  # Handled FIRST and separately: the actual distance values come
  # entirely from the externally-supplied matrix, not from otu_mat, which
  # is only used below to match sample IDs - so otu_mat's counts/percent
  # scale is irrelevant here and it's left unconverted.
  if (method %in% c("unifrac_weighted", "unifrac_unweighted")) {
    if (is.null(precomputed_dist_path))
      stop(paste0(
        "A precomputed distance matrix (.qza or .tsv) is required ",
        "for UniFrac metrics. Please upload one in the sidebar."
      ))

    d_mat  <- as.matrix(.read_precomputed_dist(precomputed_dist_path))
    shared <- intersect(rownames(otu_mat), rownames(d_mat))

    if (length(shared) == 0)
      stop(paste0(
        "No overlapping Sample IDs between the OTU table and the ",
        "uploaded distance matrix.\n",
        "  OTU table IDs (first 5) : ",
        paste(head(rownames(otu_mat), 5), collapse = ", "), "\n",
        "  Distance matrix IDs (first 5): ",
        paste(head(rownames(d_mat),   5), collapse = ", ")
      ))

    if (length(shared) < nrow(otu_mat))
      warning(sprintf(
        paste0(
          "%d OTU-table sample(s) are absent from the distance matrix ",
          "and will be dropped: %s"
        ),
        nrow(otu_mat) - length(shared),
        paste(setdiff(rownames(otu_mat), rownames(d_mat)), collapse = ", ")
      ))

    d <- stats::as.dist(d_mat[shared, shared])
    attr(d, "scale_info") <- list(
      n_total = nrow(otu_mat), n_already_percent = NA, n_converted = NA,
      converted_samples = character(0), tol = NA,
      note = "UniFrac uses a precomputed distance matrix; OTU table scale is not applicable."
    )
    return(d)
  }

  # Every metric below is computed directly FROM otu_mat, so make sure
  # it's on a consistent relative-abundance (%) scale first - detects and
  # converts raw counts automatically (.to_relative_abundance()); leaves
  # already-percent data untouched. This means the PCoA ordination and
  # PERMANOVA/ANOSIM/PERMDISP results (all derived from this same distance
  # object - see prepare_pcoa() / prepare_baseline_covariates()) are
  # always computed on percent abundance regardless of what scale was
  # uploaded.
  otu_mat <- .to_relative_abundance(otu_mat)
  scale_info <- attr(otu_mat, "scale_info")

  # -- Aitchison: CLR-transform then Euclidean -------------------
  if (method == "aitchison") {
    otu_clr <- otu_mat

    # A pseudo-count of 0.5 is the standard "half a read" convention for
    # raw counts, but it's meaningless on percentage data - a real
    # low-abundance taxon might genuinely sit at 0.05%, so flatly
    # substituting 0.5% for its zeros would swamp the true signal. Since
    # otu_mat is now always percent (converted above if needed), use the
    # data-driven pseudo-count; .detect_otu_scale() is kept as a safety
    # net in case this function is ever called directly with unconverted
    # data.
    scale_detected <- .detect_otu_scale(otu_mat)

    if (scale_detected == "counts") {
      pseudo <- 0.5
    } else {
      # Percent data: scale the pseudo-count to the data itself - half
      # of the smallest non-zero relative abundance actually observed,
      # so it stays proportionate regardless of how the % was derived.
      nonzero_vals <- otu_mat[otu_mat > 0]
      pseudo <- if (length(nonzero_vals) > 0)
        min(nonzero_vals, na.rm = TRUE) / 2 else 0.5
    }

    message(sprintf(
      paste0(
        "Aitchison distance: OTU table detected as '%s' ",
        "(sample row sums %s ~100); using pseudo-count = %.6g for zeros."
      ),
      scale_detected,
      if (scale_detected == "percent") "close to" else "not close to",
      pseudo
    ))

    otu_clr[otu_clr == 0] <- pseudo
    otu_clr <- t(apply(otu_clr, 1,
                       function(x) log(x) - mean(log(x))))
    d <- stats::dist(otu_clr, method = "euclidean")
    attr(d, "scale_info") <- scale_info
    return(d)
  }

  # -- All other metrics: delegate to vegan ---------------------
  d <- vegan::vegdist(otu_mat, method = method)
  attr(d, "scale_info") <- scale_info
  d
}

# -- Internal helper: load and align OTU + metadata ---------------------------
# Used by all public functions to avoid redundant disk reads.
.load_otu_meta <- function(otu_path,
                           meta_path,
                           is_baseline_col = "IsBaseline") {

  # -- OTU table (taxa as rows, samples as columns) --------------
  otu_raw <- read.table(
    otu_path,
    sep          = "\t",
    header       = TRUE,
    row.names    = 1,
    check.names  = FALSE,
    quote        = "",
    fill         = TRUE
  )

  # Strip MetaPhlAn suffix if present
  colnames(otu_raw) <- gsub("_metaphlan$", "", colnames(otu_raw))
  otu <- t(otu_raw)   # now: samples as rows, taxa as columns

  # -- Metadata --------------------------------------------------
  meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)

  # Trim whitespace from all character columns
  meta <- data.frame(
    lapply(meta, function(x) if (is.character(x)) trimws(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(meta)
  )

  # Coerce is_baseline_col to logical
  if (is_baseline_col %in% colnames(meta)) {
    meta[[is_baseline_col]] <- as.logical(meta[[is_baseline_col]])
  }

  # -- Match samples present in both tables ----------------------
  shared <- intersect(rownames(otu), rownames(meta))

  if (length(shared) == 0)
    stop(paste0(
      "No matching Sample IDs between OTU table and metadata.\n",
      "  OTU table rownames (first 5): ",
      paste(head(rownames(otu),  5), collapse = ", "), "\n",
      "  Metadata rownames  (first 5): ",
      paste(head(rownames(meta), 5), collapse = ", ")
    ))

  if (length(shared) < nrow(meta))
    warning(sprintf(
      "%d metadata sample(s) have no match in OTU table and will be dropped: %s",
      nrow(meta) - length(shared),
      paste(setdiff(rownames(meta), rownames(otu)), collapse = ", ")
    ))

  list(
    otu  = otu [shared, , drop = FALSE],
    meta = meta[shared, , drop = FALSE]
  )
}

# -- 1. compute_alpha_diversity() ---------------------------------------------
#' Compute alpha diversity metrics and append to metadata
#'
#' Computes Shannon, Simpson, Pielou's J, Observed ASVs, and Chao1 for
#' every sample and appends them as new columns to the metadata.
#'
#' @param otu_path  Path to OTU/species table (.tsv, samples as columns).
#' @param meta_path Path to mapping file (.csv, first column = Sample ID).
#' @param config    A \code{pcoa_config()} object. Used for
#'   \code{is_baseline_col} during metadata loading.
#' @return A data frame: original metadata + Shannon, Simpson, Pielou,
#'   Observed_ASVs, Chao1 columns.
#' @importFrom vegan diversity estimateR
#' @export
compute_alpha_diversity <- function(otu_path,
                                    meta_path,
                                    config = pcoa_config()) {

  dat  <- .load_otu_meta(otu_path, meta_path,
                         is_baseline_col = config$is_baseline_col)
  otu  <- dat$otu
  meta <- dat$meta

  # -- Compute metrics -------------------------------------------
  chao1_vals <- tryCatch(
    vegan::estimateR(otu)["S.chao1", ],
    error = function(e) {
      warning("Chao1 computation failed: ", e$message,
              " - filling with NA.")
      stats::setNames(rep(NA_real_, nrow(otu)), rownames(otu))
    }
  )

  observed <- rowSums(otu > 0)
  shannon  <- vegan::diversity(otu, index = "shannon")

  alpha_df <- data.frame(
    Shannon       = shannon,
    Simpson       = vegan::diversity(otu, index = "simpson"),
    Pielou        = ifelse(observed > 1,
                           shannon / log(observed),
                           NA_real_),
    Observed_ASVs = observed,
    Chao1         = chao1_vals,
    row.names     = rownames(otu),
    stringsAsFactors = FALSE
  )

  # Merge into metadata (preserves all meta rows)
  meta_enriched <- cbind(
    meta,
    alpha_df[rownames(meta), ALPHA_METRICS, drop = FALSE]
  )

  meta_enriched
}

# -- 2. prepare_baseline_covariates() -----------------------------------------
#' Compute and merge baseline covariates into post-baseline metadata
#'
#' Calculates per-subject baseline alpha diversity (Shannon, Observed)
#' and baseline PCoA coordinates (PC1, PC2), then merges them as numeric
#' covariates onto post-baseline rows only.
#' Baseline rows are returned with NA in the covariate columns.
#'
#' @param otu_path    Path to OTU table (.tsv).
#' @param meta_path   Path to mapping file (.csv). Should already contain
#'   alpha diversity columns if called after \code{compute_alpha_diversity()}.
#' @param dist_method Distance method for baseline PCoA (default \code{"bray"}).
#' @param config      A \code{pcoa_config()} object. Used for
#'   \code{subject_id_col} and \code{is_baseline_col}.
#' @param precomputed_dist_path Optional path to a precomputed distance
#'   matrix (.qza or .tsv) - required when \code{dist_method} is a UniFrac
#'   variant.
#' @return A data frame with the same rows as the input metadata plus
#'   columns: \code{Baseline_Shannon}, \code{Baseline_Observed},
#'   \code{Baseline_PC1}, \code{Baseline_PC2}.
#' @importFrom vegan vegdist diversity
#' @importFrom stats cmdscale
#' @export
prepare_baseline_covariates <- function(otu_path,
                                        meta_path,
                                        dist_method           = "bray",
                                        config                = pcoa_config(),
                                        precomputed_dist_path = NULL) {

  subject_col <- config$subject_id_col
  base_col    <- config$is_baseline_col

  dat  <- .load_otu_meta(otu_path, meta_path, is_baseline_col = base_col)
  otu  <- dat$otu
  meta <- dat$meta

  # -- Validate required columns ---------------------------------
  missing_cols <- setdiff(c(base_col, subject_col), colnames(meta))
  if (length(missing_cols) > 0)
    stop(sprintf(
      "Missing required column(s) in metadata: %s",
      paste(missing_cols, collapse = ", ")
    ))

  # -- Identify baseline samples ---------------------------------
  baseline_ids <- rownames(meta)[meta[[base_col]] == TRUE]

  if (length(baseline_ids) == 0)
    stop(paste0(
      "No baseline samples found. Check '", base_col, "' column.\n",
      "  Unique values found: ",
      paste(unique(meta[[base_col]]), collapse = ", ")
    ))

  message(sprintf(
    "Found %d baseline sample(s): %s",
    length(baseline_ids),
    paste(baseline_ids, collapse = ", ")
  ))

  otu_base <- otu[baseline_ids, , drop = FALSE]

  # -- STEP 1: Baseline alpha diversity -------------------------
  if ("Shannon" %in% colnames(meta)) {
    # Re-use already computed values - avoids redundant calculation
    baseline_alpha <- data.frame(
      Baseline_Shannon  = meta[baseline_ids, "Shannon"],
      Baseline_Observed = meta[baseline_ids, "Observed_ASVs"],
      stringsAsFactors  = FALSE
    )
  } else {
    baseline_alpha <- data.frame(
      Baseline_Shannon  = vegan::diversity(otu_base, index = "shannon"),
      Baseline_Observed = rowSums(otu_base > 0),
      stringsAsFactors  = FALSE
    )
  }
  baseline_alpha[[subject_col]] <- meta[baseline_ids, subject_col]

  # -- STEP 2: Baseline PCoA coordinates ------------------------
  if (nrow(otu_base) < 2)
    stop("Need at least 2 baseline samples to compute baseline PCoA.")

  dist_base <- .compute_dist(
    otu_base,
    method                = dist_method,
    precomputed_dist_path = precomputed_dist_path
  )

  pcoa_base <- stats::cmdscale(
    dist_base,
    k   = min(nrow(otu_base) - 1, 3),
    eig = TRUE
  )

  baseline_pcoa <- data.frame(
    Baseline_PC1     = pcoa_base$points[, 1],
    Baseline_PC2     = if (ncol(pcoa_base$points) >= 2)
      pcoa_base$points[, 2] else NA_real_,
    stringsAsFactors = FALSE
  )
  baseline_pcoa[[subject_col]] <- meta[baseline_ids, subject_col]

  # -- Merge alpha + PCoA -> one row per subject ------------------
  baseline_covs <- merge(baseline_alpha, baseline_pcoa, by = subject_col)
  baseline_covs <- baseline_covs[
    !duplicated(baseline_covs[[subject_col]]), , drop = FALSE
  ]

  message(sprintf(
    "Baseline covariates computed for %d subject(s): %s",
    nrow(baseline_covs),
    paste(baseline_covs[[subject_col]], collapse = ", ")
  ))

  # -- STEP 3: Merge covariates into post-baseline rows ---------
  covariate_cols <- c("Baseline_Shannon", "Baseline_Observed",
                      "Baseline_PC1",     "Baseline_PC2")

  meta_post            <- meta[meta[[base_col]] == FALSE, , drop = FALSE]
  meta_post$SampleID_tmp <- rownames(meta_post)   # preserve row names

  meta_post <- merge(
    meta_post,
    baseline_covs[, c(subject_col, covariate_cols)],
    by    = subject_col,
    all.x = TRUE   # subjects missing a post-baseline row (e.g. dropouts)
    # are retained via their Phase 1 row
  )

  rownames(meta_post)    <- meta_post$SampleID_tmp
  meta_post$SampleID_tmp <- NULL

  # -- Baseline rows: add NA covariate columns for rbind ---------
  meta_baseline                  <- meta[meta[[base_col]] == TRUE, , drop = FALSE]
  meta_baseline[, covariate_cols] <- NA

  # -- Combine and sort by SampleID -----------------------------
  meta_final <- rbind(meta_post, meta_baseline)
  meta_final <- meta_final[order(rownames(meta_final)), , drop = FALSE]

  meta_final
}

# -- 3. prepare_pcoa() --------------------------------------------------------
#' Load OTU/species table and metadata, compute PCoA
#'
#' @param otu_path  Path to OTU/species table (.tsv, samples as columns).
#' @param meta_path Path to metadata file (.csv, first column = Sample ID).
#' @param config    A \code{pcoa_config} object from \code{pcoa_config()}.
#' @param dist_method Distance method passed to \code{.compute_dist()}.
#' @param precomputed_dist_path Optional path to a precomputed distance
#'   matrix (.qza or .tsv). Required when \code{dist_method} is a UniFrac
#'   variant.
#' @return A list of class \code{pcoa_obj} containing:
#'   \item{pcoa_df}{data frame with PC coordinates + metadata}
#'   \item{var_pct}{variance explained per axis}
#'   \item{meta_cols}{metadata column names}
#'   \item{dist}{dist object}
#'   \item{dist_mat}{distance matrix (for easier subsetting)}
#'   \item{meta}{aligned metadata data frame}
#'   \item{meta_enriched}{metadata with alpha + baseline covariates}
#'   \item{config}{the pcoa_config object}
#' @importFrom vegan vegdist
#' @importFrom stats cmdscale
#' @export
prepare_pcoa <- function(otu_path,
                         meta_path,
                         config                = pcoa_config(),
                         dist_method           = "bray",
                         precomputed_dist_path = NULL) {

  dat  <- .load_otu_meta(otu_path, meta_path,
                         is_baseline_col = config$is_baseline_col)
  otu  <- dat$otu
  meta <- dat$meta

  # -- Compute distance matrix -----------------------------------
  dist_bc <- .compute_dist(
    otu,
    method                = dist_method,
    precomputed_dist_path = precomputed_dist_path
  )

  # Subset OTU/meta to samples actually present in the distance matrix
  # (relevant when a precomputed matrix covers fewer samples than the
  #  OTU table - e.g. upload of a filtered UniFrac .qza)
  dist_samples <- labels(dist_bc)
  shared       <- intersect(rownames(otu), dist_samples)
  otu  <- otu [shared, , drop = FALSE]
  meta <- meta[shared, , drop = FALSE]

  # -- PCoA via classical MDS ------------------------------------
  pcoa_res <- stats::cmdscale(
    dist_bc,
    k   = min(length(shared) - 1, 10),
    eig = TRUE
  )

  # -- Variance explained (set negative eigenvalues to 0) --------
  eig         <- pcoa_res$eig
  eig[eig < 0] <- 0
  var_pct     <- round(eig / sum(eig) * 100, 1)

  # -- Coordinate data frame -------------------------------------
  n_axes  <- ncol(pcoa_res$points)
  pcoa_df <- data.frame(
    SampleID = rownames(pcoa_res$points),
    PC1      = pcoa_res$points[, 1],
    PC2      = pcoa_res$points[, 2],
    PC3      = if (n_axes >= 3) pcoa_res$points[, 3] else 0,
    stringsAsFactors = FALSE
  )

  # Attach metadata columns
  pcoa_df  <- cbind(pcoa_df,
                    meta[pcoa_df$SampleID, , drop = FALSE])
  rownames(pcoa_df) <- NULL

  structure(
    list(
      pcoa_df       = pcoa_df,
      var_pct       = var_pct,
      meta_cols     = colnames(meta),
      dist          = dist_bc,
      dist_mat      = as.matrix(dist_bc),
      meta          = meta,
      meta_enriched = meta,
      config        = config,
      otu_path      = otu_path,
      meta_path     = meta_path,
      # Detection result from .to_relative_abundance() (via .compute_dist())
      # - which samples were already percent vs. converted from raw counts
      # for beta diversity. See app.R's output$otu_scale_info_ui for the
      # on-screen display.
      otu_scale_info = attr(dist_bc, "scale_info")
    ),
    class = "pcoa_obj"
  )
}

# -- 4. prepare_all() - master pipeline ---------------------------------------
#' Master data preparation pipeline
#'
#' Chains together in sequence:
#' \enumerate{
#'   \item \code{compute_alpha_diversity()} - Shannon, Simpson, Pielou,
#'     Observed ASVs, Chao1
#'   \item \code{prepare_baseline_covariates()} - Baseline_Shannon,
#'     Baseline_Observed, Baseline_PC1, Baseline_PC2 (only if the
#'     \code{is_baseline_col} column is present in the metadata)
#'   \item \code{prepare_pcoa()} - full PCoA on all samples
#' }
#'
#' The OTU table is loaded once per sub-function via the shared
#' \code{.load_otu_meta()} helper; a temporary CSV is used to pass the
#' alpha-enriched metadata between steps without re-reading the raw files.
#'
#' @param otu_path    Path to OTU table (.tsv, taxa as rows).
#' @param meta_path   Path to mapping file (.csv, first column = Sample ID).
#' @param config      A \code{pcoa_config()} object.
#' @param dist_method Distance method (default \code{"bray"}).
#'   Use \code{"unifrac_weighted"} or \code{"unifrac_unweighted"} together
#'   with \code{precomputed_dist_path}.
#' @param precomputed_dist_path Optional path to a precomputed distance
#'   matrix (.qza artifact or plain .tsv). Required when \code{dist_method}
#'   is a UniFrac variant.
#' @return A \code{pcoa_obj} with \code{$meta_enriched} containing all
#'   computed covariates, ready for beta diversity statistics and LMM.
#' @export
prepare_all <- function(otu_path,
                        meta_path,
                        config                = pcoa_config(),
                        dist_method           = "bray",
                        precomputed_dist_path = NULL) {

  # -- Step A: alpha diversity -> enriched metadata ---------------
  message("-- Step A: Computing alpha diversity -------------------")

  meta_alpha <- compute_alpha_diversity(
    otu_path  = otu_path,
    meta_path = meta_path,
    config    = config
  )

  # Write alpha-enriched metadata to a temp file so later steps can
  # read it without keeping the full OTU matrix in memory twice.
  # row.names = TRUE preserves SampleIDs as row names on re-read.
  tmp_meta <- tempfile(fileext = ".csv")
  write.csv(meta_alpha, tmp_meta, row.names = TRUE)

  message(sprintf(
    "Alpha diversity computed for %d samples. Columns added: %s",
    nrow(meta_alpha),
    paste(intersect(ALPHA_METRICS, colnames(meta_alpha)), collapse = ", ")
  ))

  # -- Step B: baseline covariates -------------------------------
  # Only executed when the is_baseline_col column exists in the metadata.
  # Adds: Baseline_Shannon, Baseline_Observed, Baseline_PC1, Baseline_PC2
  # to post-baseline rows; baseline rows receive NA for these columns.
  # Subjects missing a post-baseline row (e.g. dropouts / MONA) are
  # handled via all.x = TRUE inside prepare_baseline_covariates().
  message("-- Step B: Computing baseline covariates ---------------")

  base_col <- config$is_baseline_col

  if (base_col %in% colnames(meta_alpha)) {

    meta_enriched <- tryCatch(
      prepare_baseline_covariates(
        otu_path              = otu_path,
        meta_path             = tmp_meta,
        dist_method           = dist_method,
        config                = config,
        precomputed_dist_path = precomputed_dist_path
      ),
      error = function(e) {
        warning(sprintf(
          paste0(
            "Baseline covariate computation failed: %s\n",
            "Continuing without baseline covariates."
          ),
          e$message
        ))
        meta_alpha   # fall back to alpha-only metadata
      }
    )

    # Overwrite temp file with fully enriched metadata for Step C
    write.csv(meta_enriched, tmp_meta, row.names = TRUE)

    message(sprintf(
      "Baseline covariates merged. Final metadata dimensions: %d x %d",
      nrow(meta_enriched), ncol(meta_enriched)
    ))

  } else {
    # No baseline column found - skip and warn
    warning(sprintf(
      paste0(
        "'%s' column not found in metadata. ",
        "Skipping baseline covariate computation. ",
        "Add a %s (TRUE/FALSE) column to your mapping file, or set ",
        "is_baseline_col in pcoa_config() to match your column name."
      ),
      base_col, base_col
    ))
    meta_enriched <- meta_alpha
    write.csv(meta_enriched, tmp_meta, row.names = TRUE)
  }

  # -- Step C: full PCoA on ALL samples -------------------------
  # Baseline samples are included in the ordination; exclusion for
  # statistical tests happens downstream in stats.R via is_baseline_col.
  message("-- Step C: Computing full PCoA -------------------------")

  pcoa_obj <- tryCatch(
    prepare_pcoa(
      otu_path              = otu_path,
      meta_path             = tmp_meta,
      config                = config,
      dist_method           = dist_method,
      precomputed_dist_path = precomputed_dist_path
    ),
    error = function(e) {
      stop(sprintf(
        "PCoA computation failed: %s\nCheck OTU table format.",
        e$message
      ))
    }
  )

  # -- Step D: attach enriched metadata to pcoa_obj --------------
  # $meta_enriched is consumed by:
  #   app.R  -> rv$meta_enriched  (alpha plot, LMM table, metadata tab)
  #   stats.R -> run_alpha_lmm_all_metrics(meta_enriched = ...)
  #   plots.R -> build_alpha_scatter(meta_enriched = ...)
  pcoa_obj$meta_enriched <- meta_enriched

  # -- Step E: clean up temp file --------------------------------
  if (file.exists(tmp_meta)) {
    file.remove(tmp_meta)
    message("Temporary metadata file cleaned up.")
  }

  # -- Step F: summary -------------------------------------------
  n_baseline <- if (base_col %in% colnames(meta_enriched))
    sum(meta_enriched[[base_col]] == TRUE,  na.rm = TRUE) else NA_integer_
  n_post     <- if (base_col %in% colnames(meta_enriched))
    sum(meta_enriched[[base_col]] == FALSE, na.rm = TRUE) else NA_integer_

  message(sprintf(
    paste0(
      "\n-- prepare_all() complete ------------------------------\n",
      "  Samples total         : %d\n",
      "  Baseline samples      : %s\n",
      "  Post-baseline samples : %s\n",
      "  Metadata columns      : %d\n",
      "  Alpha metrics         : %s\n",
      "  Baseline covariates   : %s\n",
      "  Distance method       : %s\n",
      "  Precomputed dist used : %s\n",
      "--------------------------------------------------------"
    ),
    nrow(meta_enriched),
    if (is.na(n_baseline)) "N/A" else n_baseline,
    if (is.na(n_post))     "N/A" else n_post,
    ncol(meta_enriched),
    paste(intersect(ALPHA_METRICS, colnames(meta_enriched)),
          collapse = ", "),
    paste(intersect(
      c("Baseline_Shannon", "Baseline_Observed",
        "Baseline_PC1",     "Baseline_PC2"),
      colnames(meta_enriched)
    ), collapse = ", "),
    dist_method,
    if (!is.null(precomputed_dist_path))
      basename(precomputed_dist_path) else "none"
  ))

  pcoa_obj
}

# -- 5. compute_taxonomy_barplot_data() ---------------------------------------
#' Parse clade names and collapse to a chosen taxonomic level
#'
#' Reads the same OTU/species table already loaded for PCoA (taxa as rows,
#' samples as columns).  Clade names must follow MetaPhlAn semicolon
#' notation, e.g.:
#'   k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;...
#'
#' @param otu_path   Path to OTU table (.tsv, taxa as rows, samples as cols).
#' @param meta_path  Path to mapping file (.csv, first column = Sample ID).
#' @param level      One of \code{TAXA_LEVELS} (e.g. \code{"Family"}).
#' @param group_var  Metadata column to facet by (e.g. \code{"Treatment"}).
#' @param top_n      Number of top taxa to show individually; the rest are
#'   collapsed into \code{"Others"} (default \code{TAXA_DEFAULT_TOP_N}).
#' @param config     A \code{pcoa_config()} object.
#' @return A data frame with columns:
#'   \code{SampleID}, \code{Taxon}, \code{RelAbundance},
#'   plus all metadata columns.
#' @export
compute_taxonomy_barplot_data <- function(otu_path,
                                          meta_path,
                                          level    = "Family",
                                          group_var,
                                          top_n    = TAXA_DEFAULT_TOP_N,
                                          config   = pcoa_config()) {

  # -- Load raw OTU table (taxa as rows) ------------------------
  otu_raw <- read.table(
    otu_path,
    sep          = "\t",
    header       = TRUE,
    row.names    = 1,
    check.names  = FALSE,
    quote        = "",
    fill         = TRUE
  )
  colnames(otu_raw) <- gsub("_metaphlan$", "", colnames(otu_raw))

  # -- Load and clean metadata -----------------------------------
  meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
  meta <- data.frame(
    lapply(meta, function(x) if (is.character(x)) trimws(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(meta)
  )
  if (config$is_baseline_col %in% colnames(meta)) {
    meta[[config$is_baseline_col]] <- as.logical(meta[[config$is_baseline_col]])
  }

  # -- Align samples ---------------------------------------------
  shared <- intersect(colnames(otu_raw), rownames(meta))
  if (length(shared) == 0)
    stop(paste0(
      "No matching Sample IDs between OTU table columns and metadata rows.\n",
      "  OTU table columns (first 5): ",
      paste(head(colnames(otu_raw), 5), collapse = ", "), "\n",
      "  Metadata rows     (first 5): ",
      paste(head(rownames(meta),    5), collapse = ", ")
    ))
  otu_raw <- otu_raw[, shared, drop = FALSE]
  meta    <- meta[shared,     , drop = FALSE]

  # -- Internal: extract one taxonomic level from a clade string -
  .extract_level <- function(clade, prefix) {
    parts <- unlist(strsplit(as.character(clade), ";"))

    if (prefix == "s") {
      # Species: return "Genus species" or "Genus unclassified"
      genus_part   <- parts[grep("^g__", parts)]
      species_part <- parts[grep("^s__", parts)]
      genus   <- if (length(genus_part)   > 0)
        sub("g__", "", genus_part[1])   else NA_character_
      species <- if (length(species_part) > 0)
        sub("s__", "", species_part[1]) else NA_character_
      if (is.na(species) || species == "" || species == "_") {
        return(if (!is.na(genus))
          paste0(genus, " unclassified") else "Unclassified")
      }
      return(species)
    }

    target <- parts[grep(paste0("^", prefix, "__"), parts)]
    if (length(target) == 0) return(NA_character_)
    sub(paste0(prefix, "__"), "", target[1])
  }

  prefix <- TAXA_PREFIXES[[level]]

  # -- Wide -> long -----------------------------------------------
  feat_long <- data.frame(
    clade_name = rownames(otu_raw),
    otu_raw,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    tidyr::pivot_longer(
      cols      = -clade_name,
      names_to  = "SampleID",
      values_to = "Abundance"
    )

  # -- Assign taxonomic label ------------------------------------
  feat_long[[level]] <- vapply(
    feat_long$clade_name,
    .extract_level,
    FUN.VALUE = character(1),
    prefix    = prefix
  )

  # Clades that don't resolve to this rank (e.g. reads classified only
  # to a higher level, with no "p__..." tag etc.) are kept as their own
  # "Unclassified" bucket rather than dropped - dropping them would
  # silently shrink that sample's total abundance below 100%.
  feat_long[[level]][is.na(feat_long[[level]])] <- "Unclassified"

  # -- Attach metadata -------------------------------------------
  meta_df         <- meta
  meta_df$SampleID <- rownames(meta)
  feat_long <- merge(feat_long, meta_df, by = "SampleID", all.x = TRUE)

  # -- Identify top N taxa across all samples --------------------
  taxon_totals <- tapply(feat_long$Abundance, feat_long[[level]], sum)
  taxon_totals <- sort(taxon_totals, decreasing = TRUE)
  top_taxa     <- names(taxon_totals)[seq_len(min(top_n, length(taxon_totals)))]

  feat_long$Taxon <- ifelse(
    feat_long[[level]] %in% top_taxa,
    feat_long[[level]],
    "Others"
  )

  # -- Summarise to sample  taxon -------------------------------
  # Sum abundance per SampleID  Taxon (multiple clade rows may map
  # to the same taxon after collapsing "Others"/"Unclassified").
  # IMPORTANT: aggregate only on SampleID + Taxon, then re-attach
  # metadata afterward. Using "Abundance ~ ." here would implicitly
  # group by every metadata column too, and aggregate()'s default
  # na.action = na.omit silently drops any row with an NA in ANY of
  # those columns - e.g. wiping out baseline samples regardless of
  # the exclude_baseline setting, if any unrelated metadata column
  # happens to be NA for them.
  agg <- stats::aggregate(
    Abundance ~ SampleID + Taxon,
    data = feat_long[, c("SampleID", "Taxon", "Abundance")],
    FUN  = sum
  )

  # Re-attach metadata (one row per sample) now that aggregation is done
  agg <- merge(agg, meta_df, by = "SampleID", all.x = TRUE)

  # -- Relative abundance (per sample) --------------------------
  sample_totals        <- tapply(agg$Abundance, agg$SampleID, sum)
  agg$RelAbundance     <- agg$Abundance /
    sample_totals[agg$SampleID] * 100

  # Tag the parameters actually used to build this table so callers
  # (app.R) can detect a stale plot - e.g. the user changed the
  # "Top N taxa" input after clicking "Build taxonomy plot" but
  # hasn't re-clicked it, so the displayed data no longer matches
  # the current UI selection.
  attr(agg, "top_n_used") <- top_n
  attr(agg, "level_used") <- level

  agg
}
