# maaslin.R
# ─────────────────────────────────────────────────────────────────────────────
# MaAsLin3 wrapper functions for the Microbiome Explorer Shiny app.
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. validate_maaslin_formula() ────────────────────────────────────────────
#' @export
validate_maaslin_formula <- function(formula_str, meta) {
  if (is.null(formula_str) || !nzchar(trimws(formula_str)))
    stop("Formula is empty. Please enter at least one fixed effect.",
         call. = FALSE)
  fml <- tryCatch(
    stats::as.formula(paste("~", formula_str)),
    error = function(e)
      stop("Cannot parse formula '", formula_str, "':\n  ", e$message,
           call. = FALSE)
  )
  term_labels <- attr(stats::terms(fml), "term.labels")
  simple_terms <- term_labels[!grepl("[:*^]", term_labels)]
  simple_terms <- gsub("^I\\((.+)\\)$", "\\1", simple_terms)
  missing_cols <- setdiff(simple_terms, colnames(meta))
  if (length(missing_cols) > 0)
    stop(sprintf(
      "Formula term(s) not found in metadata: %s\nAvailable columns: %s",
      paste(missing_cols, collapse = ", "),
      paste(colnames(meta), collapse = ", ")
    ), call. = FALSE)
  invisible(TRUE)
}

# ── 2. prepare_maaslin_input() ───────────────────────────────────────────────
#' @export
prepare_maaslin_input <- function(otu_path,
                                  meta_enriched,
                                  level           = "Species",
                                  config          = pcoa_config(),
                                  min_prevalence  = 0.1,
                                  min_abundance   = 0) {
  base_col <- config$is_baseline_col

  otu_raw <- read.table(
    otu_path,
    sep = "\t", header = TRUE, row.names = 1,
    check.names = FALSE, quote = "", fill = TRUE
  )
  colnames(otu_raw) <- gsub("_metaphlan$", "", colnames(otu_raw))

  meta_use <- meta_enriched
  if (base_col %in% colnames(meta_use)) {
    meta_use <- meta_use[
      !is.na(meta_use[[base_col]]) & meta_use[[base_col]] == FALSE,
      , drop = FALSE
    ]
  }

  shared <- intersect(colnames(otu_raw), rownames(meta_use))
  if (length(shared) == 0)
    stop("No overlapping Sample IDs between OTU table and metadata.",
         call. = FALSE)

  otu_use  <- otu_raw[, shared, drop = FALSE]
  meta_use <- meta_use[shared, , drop = FALSE]

  # ── Collapse taxonomy (same as before) ────────────────────────
  if (!is.null(level) && level != "Species") {
    prefix <- TAXA_PREFIXES[[level]]

    # Always build the feature name from the FULL resolved lineage up to
    # and including the target rank (e.g. "Firmicutes_Clostridia_
    # Lachnospiraceae_uncultured"), rather than keying off the single
    # tag at that rank. This keeps every feature distinct - two
    # "uncultured" genera from different families no longer collapse
    # into one shared feature when aggregate() sums by taxon below, and
    # a clade with no tag at all is no longer silently dropped as long
    # as an ancestor rank is resolved.
    .extract_level_local <- function(clade, pfx) {
      parts <- trimws(unlist(strsplit(as.character(clade), ";")))

      .rank_val <- function(p) {
        tgt <- parts[grep(paste0("^", p, "__"), parts)]
        if (length(tgt) == 0) return(NA_character_)
        v <- sub(paste0(p, "__"), "", tgt[1])
        if (v == "" || v == "_") NA_character_ else v
      }

      rank_order <- unname(TAXA_PREFIXES)
      target_i   <- match(pfx, rank_order)
      ranks_upto <- rank_order[seq_len(target_i)]

      vals <- vapply(ranks_upto, .rank_val, character(1))
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) return(NA_character_)
      paste(vals, collapse = "_")
    }
    taxa_labels <- vapply(rownames(otu_use), .extract_level_local,
                          FUN.VALUE = character(1), pfx = prefix)
    keep        <- !is.na(taxa_labels)
    otu_use     <- otu_use[keep, , drop = FALSE]
    taxa_labels <- taxa_labels[keep]
    otu_df <- data.frame(taxon = taxa_labels, otu_use, check.names = FALSE)
    otu_agg <- stats::aggregate(. ~ taxon, data = otu_df, FUN = sum)
    rownames(otu_agg) <- otu_agg$taxon
    otu_agg$taxon     <- NULL
    features <- as.data.frame(t(otu_agg), check.names = FALSE)
  } else {
    clean_names <- vapply(rownames(otu_use), function(x) {
      parts <- unlist(strsplit(x, ";"))
      trimws(parts[length(parts)])
    }, FUN.VALUE = character(1))
    rownames(otu_use) <- make.unique(clean_names)
    features <- as.data.frame(t(otu_use), check.names = FALSE)
  }

  # ── Convert to relative abundance (%) ──────────────────────────
  # Ensures MaAsLin3 always receives percent abundance regardless of
  # whether the uploaded OTU table was raw counts or already percent -
  # detects and converts raw counts automatically (.to_relative_abundance()
  # in data.R, samples-as-rows / taxa-as-columns, matching `features`
  # here); already-percent samples are left untouched. Applied BEFORE
  # MaAsLin3's own `normalization` option runs, so the behavior is
  # consistent regardless of what normalization the user selects in the
  # UI (e.g. still correct even if they pick "NONE").
  features_mat <- .to_relative_abundance(as.matrix(features))
  # as.data.frame() below drops custom attributes, so capture the
  # detection result now and return it separately for app.R to display.
  otu_scale_info <- attr(features_mat, "scale_info")
  features <- as.data.frame(features_mat, check.names = FALSE)

  common   <- intersect(rownames(features), rownames(meta_use))
  features <- features[common, , drop = FALSE]
  meta_use <- meta_use[common, , drop = FALSE]

  # ── Sanity checks (INSIDE the function now) ───────────────────
  dup_cols <- colnames(meta_use)[duplicated(colnames(meta_use))]
  if (length(dup_cols) > 0)
    stop("Duplicate metadata column names detected: ",
         paste(unique(dup_cols), collapse = ", "), call. = FALSE)

  bad_cols <- names(meta_use)[
    !vapply(meta_use,
            function(x) is.atomic(x) || is.factor(x),
            logical(1))
  ]
  if (length(bad_cols) > 0)
    stop("Metadata columns are not plain vectors: ",
         paste(bad_cols, collapse = ", "), call. = FALSE)

  char_cols <- vapply(meta_use, is.character, logical(1))
  meta_use[char_cols] <- lapply(meta_use[char_cols], factor)

  single_lvl <- vapply(meta_use, function(x) {
    is.factor(x) && nlevels(droplevels(x)) < 2
  }, logical(1))
  if (any(single_lvl)) {
    message("Dropping single-level metadata columns: ",
            paste(names(meta_use)[single_lvl], collapse = ", "))
    meta_use <- meta_use[, !single_lvl, drop = FALSE]
  }

  # ── Diagnostic dump — very useful for debugging ───────────────
  message("prepare_maaslin_input() finished: ",
          nrow(features), " samples × ", ncol(features), " taxa; ",
          "metadata cols: ",
          paste(colnames(meta_use), collapse = ", "))

  list(features = features, metadata = meta_use, otu_scale_info = otu_scale_info)
}

# ── 3. run_maaslin() ─────────────────────────────────────────────────────────
#' @export
run_maaslin <- function(features,
                        metadata,
                        formula_str,
                        group_effects_var    = NULL,
                        random_effect        = NULL,
                        use_random           = FALSE,
                        small_random_effects = FALSE,
                        normalization        = "TSS",
                        transform            = "LOG",
                        min_prevalence       = 0.1,
                        min_abundance        = 0,
                        padj_method          = "BH",
                        max_significance     = 0.25) {

  if (!requireNamespace("maaslin3", quietly = TRUE))
    stop("Package 'maaslin3' is required.", call. = FALSE)

  if (is.null(formula_str) || !nzchar(trimws(formula_str)))
    stop("Formula is empty.", call. = FALSE)

  group_var <- if (!is.null(group_effects_var) &&
                   nzchar(trimws(group_effects_var)))
    trimws(group_effects_var) else NULL
  rand_effect <- if (isTRUE(use_random) &&
                     !is.null(random_effect) &&
                     nzchar(trimws(random_effect)))
    trimws(random_effect) else NULL

  if (!is.null(group_var) && !group_var %in% colnames(metadata))
    stop(sprintf("Grouping factor '%s' not found in metadata.", group_var),
         call. = FALSE)
  if (!is.null(rand_effect) && !rand_effect %in% colnames(metadata))
    stop(sprintf("Random effect column '%s' not found in metadata.",
                 rand_effect), call. = FALSE)
  if (!is.null(group_var) && !is.null(rand_effect) &&
      identical(group_var, rand_effect))
    stop("Grouping factor and random effect cannot be the same column.",
         call. = FALSE)

  fixed_term_labels <- tryCatch(
    attr(stats::terms(stats::as.formula(paste("~", formula_str))),
         "term.labels"),
    error = function(e) character(0)
  )

  if (!is.null(group_var) && group_var %in% fixed_term_labels)
    stop(sprintf("'%s' used as grouping factor is also in formula.",
                 group_var), call. = FALSE)
  if (!is.null(rand_effect) && rand_effect %in% fixed_term_labels)
    stop(sprintf("'%s' used as random effect is also in formula.",
                 rand_effect), call. = FALSE)

  # ── Trim metadata to only columns we need (INSIDE now) ────────
  needed <- unique(c(
    fixed_term_labels,
    if (!is.null(group_var))   group_var   else character(),
    if (!is.null(rand_effect)) rand_effect else character()
  ))
  missing_needed <- setdiff(needed, colnames(metadata))
  if (length(missing_needed))
    stop("Metadata missing columns: ",
         paste(missing_needed, collapse = ", "), call. = FALSE)
  metadata <- metadata[, needed, drop = FALSE]

  complete_rows <- stats::complete.cases(metadata)
  if (any(!complete_rows)) {
    message(sprintf("Dropping %d samples with NA in modeling columns.",
                    sum(!complete_rows)))
    metadata <- metadata[complete_rows, , drop = FALSE]
    features <- features[rownames(metadata), , drop = FALSE]
  }

  # Force character → factor once more, and drop unused levels.
  # (Trimming above may have re-introduced character columns if the
  # source was a tibble, and subsetting rows may leave unused levels.)
  for (nm in colnames(metadata)) {
    x <- metadata[[nm]]
    if (is.character(x)) x <- factor(x)
    if (is.factor(x))    x <- droplevels(x)
    metadata[[nm]] <- x
  }

  # Re-check: any factor with < 2 levels after row drops is fatal
  bad_lvl <- vapply(metadata, function(x) {
    is.factor(x) && nlevels(x) < 2
  }, logical(1))
  if (any(bad_lvl))
    stop("After NA filtering these factor columns have <2 levels: ",
         paste(names(metadata)[bad_lvl], collapse = ", "),
         ". Adjust your formula or filters.", call. = FALSE)

  # Coerce to plain base data.frame (not tibble / data.table).
  # This is the single most common fix for the generic
  # "argument 1 is not a vector" error inside MaAsLin3, because
  # MaAsLin3 uses base-R subsetting internally and tibble's
  # [.tbl_df returns a 1-column tibble instead of a vector.
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  features <- as.data.frame(features, stringsAsFactors = FALSE,
                            check.names = FALSE)
  rownames(metadata) <- rownames(metadata)   # keep row names intact
  rownames(features) <- rownames(features)

  # Align order once more
  common <- intersect(rownames(features), rownames(metadata))
  if (length(common) < 3)
    stop("Fewer than 3 samples remain after alignment/NA filtering.",
         call. = FALSE)
  features <- features[common, , drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]

  # ── Build unified lme4-style formula ──────────────────────────
  rhs <- trimws(formula_str)
  if (!is.null(group_var))
    rhs <- paste0(rhs, " + group(", group_var, ")")
  if (!is.null(rand_effect))
    rhs <- paste0(rhs, " + (1|", rand_effect, ")")
  full_formula <- paste("~", rhs)

  # ── Diagnostic dump right before calling MaAsLin3 ─────────────
  message("── MaAsLin3 pre-flight ──")
  message("  formula:    ", full_formula)
  message("  n_samples:  ", nrow(features))
  message("  n_features: ", ncol(features))
  message("  metadata columns and types:")
  for (nm in colnames(metadata)) {
    x <- metadata[[nm]]
    message(sprintf("    %-25s %s (%s)",
                    nm,
                    class(x)[1],
                    if (is.factor(x))
                      paste(levels(x), collapse = "/")
                    else
                      paste(range(x, na.rm = TRUE), collapse = "..")))
  }

  # ── Temp output dir ───────────────────────────────────────────
  out_dir <- tempfile("maaslin3_")
  dir.create(out_dir, recursive = TRUE)
  on.exit(
    tryCatch(unlink(out_dir, recursive = TRUE),
             error = function(e) NULL),
    add = TRUE
  )

  # ── Run MaAsLin3 ──────────────────────────────────────────────
  fit <- maaslin3::maaslin3(
    input_data           = features,
    input_metadata       = metadata,
    output               = out_dir,
    formula              = full_formula,
    small_random_effects = small_random_effects,
    normalization        = normalization,
    transform            = transform,
    min_prevalence       = min_prevalence,
    min_abundance        = min_abundance,
    max_significance     = max_significance,
    correction           = padj_method,
    plot_summary_plot    = FALSE,
    plot_associations    = FALSE,
    save_models          = FALSE
  )

  # ── Read results file ─────────────────────────────────────────
  res_path <- file.path(out_dir, "all_results.tsv")
  if (!file.exists(res_path))
    stop("MaAsLin3 produced no results file.", call. = FALSE)
  res <- read.table(res_path, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
  if (nrow(res) == 0)
    stop("MaAsLin3 returned an empty results table.", call. = FALSE)

  res <- .standardise_maaslin_cols(res)
  res <- res[order(res$qval, res$pval, na.last = TRUE), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "formula_used") <- full_formula
  res
}

# ── Internal: standardise column names ───────────────────────────────────────
.standardise_maaslin_cols <- function(df) {
  rename_map <- list(
    feature  = c("feature", "Feature", "taxon", "Taxon"),
    metadata = c("metadata", "Metadata", "variable", "Variable"),
    value    = c("value", "Value", "level", "Level"),
    coef     = c("coef", "Coef", "coefficient", "Coefficient",
                 "estimate", "Estimate"),
    stderr   = c("stderr", "StdErr", "se", "SE"),
    pval     = c("pval", "pval_individual", "pValue", "p.value",
                 "P.value", "p_value"),
    qval     = c("qval", "qval_individual", "qValue", "q.value",
                 "Q.value", "q_value", "padj", "p_adjusted")
  )
  current <- colnames(df)
  for (new_name in names(rename_map)) {
    if (new_name %in% current) next
    for (old_name in rename_map[[new_name]]) {
      if (old_name %in% current) {
        colnames(df)[colnames(df) == old_name] <- new_name
        break
      }
    }
  }
  for (col in names(rename_map)) {
    if (!col %in% colnames(df)) df[[col]] <- NA
  }
  df
}

# ── 4. summarise_maaslin_results() ───────────────────────────────────────────
#' @export
summarise_maaslin_results <- function(results, threshold = 0.25) {
  if (is.null(results) || nrow(results) == 0)
    return(data.frame(message = "No results available."))
  sig <- results[!is.na(results$qval) & results$qval < threshold, ,
                 drop = FALSE]
  if (nrow(sig) == 0)
    return(data.frame(message = sprintf(
      "No features significant at FDR < %.2f.", threshold
    )))
  num_cols <- intersect(c("coef", "stderr", "pval", "qval"), colnames(sig))
  sig[num_cols] <- lapply(sig[num_cols], function(x) round(x, 6))
  sig[order(sig$qval, sig$pval, na.last = TRUE), , drop = FALSE]
}
