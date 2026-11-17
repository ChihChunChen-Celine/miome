# descriptive_summary.R
# -----------------------------------------------------------------------------
# Descriptive relative-abundance summary (median / min / max / range / mean / SD / 95% CI) of every
# taxon, per user-chosen group, for every taxonomic level - listed from the
# highest rank (Phylum) to the lowest (Species) on ONE Excel sheet.
#
#   1. summarize_taxa_by_group()           - computes the numbers (long data frame)
#   2. export_taxa_descriptive_workbook()  - computes + writes the formatted .xlsx
#
# Lineage parsing, layout detection, file reading and the counts-vs-percent
# rule are shared helpers in utils.R.
# -----------------------------------------------------------------------------

# -- Internal: OTU table from a path or a data frame / matrix -----------------
.read_taxa_matrix <- function(otu) {
  if (is.character(otu) && length(otu) == 1) return(.read_otu_table(otu))
  if (!is.data.frame(otu) && !is.matrix(otu))
    stop("'otu' must be a file path or a data frame / matrix (taxa as rows, ",
         "samples as columns).", call. = FALSE)
  mat <- as.matrix(otu)
  colnames(mat) <- sub("_metaphlan$", "", colnames(mat))
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (anyNA(mat)) {
    warning("Non-numeric or missing values in the OTU table were set to 0.", call. = FALSE)
    mat[is.na(mat)] <- 0
  }
  mat
}

.lineage_keys <- function(lin, sep)
  apply(lin, 1, function(r) paste(r[!is.na(r)], collapse = sep))

# Relative abundance (%) for a taxa x samples matrix (shared rule in utils.R)
.relativize_taxa_by_sample <- function(mat) {
  out <- t(.to_relative_abundance(t(mat)))
  attr(out, "scale_info") <- NULL
  out
}

# Format a number the way Excel's TEXT(x,"0.00") does
.fmt2 <- function(x) formatC(round(x + 1e-9, 2), format = "f", digits = 2)

# -- 1. summarize_taxa_by_group() ---------------------------------------------
#' Descriptive relative-abundance summary per group, for every taxonomic level
#'
#' For each requested level (Phylum -> Species) computes the median, minimum,
#' maximum, range, mean, standard deviation and 95\% CI of the mean of each taxon's relative abundance (\%) within every
#' group defined by one or more metadata columns.
#'
#' Two table layouts are supported (see \code{table_type}):
#' \itemize{
#'   \item \strong{terminal}: rows are the deepest classification only (e.g.
#'     QIIME 2 L7). Higher levels are built by summing the rows beneath.
#'   \item \strong{all_ranks}: a row for every rank (un-filtered MetaPhlAn).
#'     Each level is taken from the rows at exactly that rank.
#' }
#'
#' @param otu OTU/taxa table: file path (.tsv, taxa as rows) or data frame /
#'   matrix.
#' @param meta Metadata: file path (.csv, first column = Sample ID) or data
#'   frame with Sample IDs as row names.
#' @param group_columns One or more metadata columns defining the groups;
#'   several are pasted together with \code{group_sep}.
#' @param sort_columns Columns that set the ORDER of the groups (default:
#'   \code{group_columns}).
#' @param group_levels Optional exact order of the final group labels.
#' @param group_sep Separator for pasted group labels.
#' @param tax_levels Levels to summarise, any of \code{names(TAXA_PREFIXES)}.
#' @param exclude_baseline Logical; drop baseline samples.
#' @param is_baseline_col Baseline flag column.
#' @param table_type \code{"auto"}, \code{"terminal"} or \code{"all_ranks"}.
#' @param drop_all_zero Logical; drop taxa that are zero in every sample.
#' @param digits Decimal places for median / min / max.
#' @param lineage_sep Separator used in the output taxon labels.
#' @return A list: \code{summary} (long data frame), \code{group_levels},
#'   \code{tax_levels}, \code{table_type}.
#' @export
summarize_taxa_by_group <- function(otu,
                                    meta,
                                    group_columns,
                                    sort_columns     = group_columns,
                                    group_levels     = NULL,
                                    group_sep        = "_",
                                    tax_levels       = c("Phylum", "Class", "Order",
                                                         "Family", "Genus", "Species"),
                                    exclude_baseline = FALSE,
                                    is_baseline_col  = "IsBaseline",
                                    table_type       = c("auto", "terminal", "all_ranks"),
                                    drop_all_zero    = FALSE,
                                    digits           = 3,
                                    lineage_sep      = "|") {
  table_type <- match.arg(table_type)
  all_levels <- names(TAXA_PREFIXES)
  bad_lv <- setdiff(tax_levels, all_levels)
  if (length(bad_lv))
    stop(sprintf("Unknown taxonomic level(s): %s. Choose from: %s",
                 paste(bad_lv, collapse = ", "), paste(all_levels, collapse = ", ")),
         call. = FALSE)
  tax_levels <- all_levels[all_levels %in% tax_levels]
  if (missing(group_columns) || !length(group_columns))
    stop("'group_columns' must name at least one metadata column.", call. = FALSE)

  mat  <- .read_taxa_matrix(otu)
  meta <- .read_meta(meta, is_baseline_col, sample_ids = colnames(mat))
  mat  <- .harmonise_samples(mat, rownames(meta))
  need <- unique(c(group_columns, sort_columns))
  miss <- setdiff(need, colnames(meta))
  if (length(miss))
    stop(sprintf("Column(s) not found in metadata: %s\nAvailable columns: %s",
                 paste(miss, collapse = ", "), paste(colnames(meta), collapse = ", ")),
         call. = FALSE)

  shared <- intersect(rownames(meta), colnames(mat))
  if (!length(shared))
    stop("No matching Sample IDs between OTU table and metadata.\n",
         "  OTU table columns (first 5): ", paste(head(colnames(mat), 5), collapse = ", "),
         "\n  Metadata rows     (first 5): ", paste(head(rownames(meta), 5), collapse = ", "),
         call. = FALSE)
  meta <- meta[shared, , drop = FALSE]

  if (isTRUE(exclude_baseline)) {
    meta <- meta[!.is_baseline(meta, is_baseline_col), , drop = FALSE]
    message(sprintf("Baseline samples excluded. N remaining = %d", nrow(meta)))
  }
  keep <- stats::complete.cases(meta[, need, drop = FALSE])
  if (any(!keep))
    message(sprintf("%d sample(s) dropped due to NA in: %s",
                    sum(!keep), paste(need, collapse = ", ")))
  meta <- meta[keep, , drop = FALSE]
  if (!nrow(meta)) stop("No samples remaining after filtering.", call. = FALSE)

  # -- Group labels and their order (radix = locale-independent) -------
  meta  <- meta[do.call(order, c(unname(as.list(meta[sort_columns])),
                                 list(method = "radix"))), , drop = FALSE]
  grp   <- do.call(paste, c(unname(lapply(meta[group_columns], as.character)),
                            list(sep = group_sep)))
  found <- unique(grp)
  if (is.null(group_levels)) {
    group_levels <- found
  } else {
    group_levels <- as.character(group_levels)
    not_listed   <- setdiff(found, group_levels)
    if (length(not_listed))
      stop(sprintf("'group_levels' is missing group(s) present in the data: %s",
                   paste(not_listed, collapse = ", ")), call. = FALSE)
    extra <- setdiff(group_levels, found)
    if (length(extra))
      warning(sprintf("'group_levels' entries not found in the data were ignored: %s",
                      paste(extra, collapse = ", ")), call. = FALSE)
    group_levels <- group_levels[group_levels %in% found]
  }

  samples <- rownames(meta)
  mat     <- mat[, samples, drop = FALSE]
  grp_idx <- split(seq_along(samples), factor(grp, levels = group_levels))
  message(sprintf("Groups (%d): %s", length(group_levels),
                  paste(sprintf("%s [n=%d]", group_levels, lengths(grp_idx)),
                        collapse = "; ")))

  # -- Parse lineages (strain/SGB rows repeat their species -> dropped) -
  is_strain <- grepl("(^|[;|])t__", rownames(mat))
  if (any(is_strain)) {
    message(sprintf("Ignoring %d strain-level (t__) row(s) to avoid double counting.",
                    sum(is_strain)))
    mat <- mat[!is_strain, , drop = FALSE]
  }
  lin   <- .lineage_matrix(rownames(mat))
  if (all(is.na(lin))) {
    # Names only (e.g. a species-level MetaPhlAn export without lineage):
    # list the taxa as uploaded, as one "Species" block.
    message("The table has taxon names only (no lineage): taxa are listed as ",
            "uploaded under 'Species'; higher levels are not available.")
    lin[, "Species"] <- rownames(mat)
    tax_levels <- "Species"
    table_type <- "terminal"
  }
  depth <- .lineage_depth(lin)
  if (table_type == "auto") {
    table_type <- .detect_layout(lin, lineage_sep)
    message(sprintf("Table layout detected: %s.", table_type))
  }

  # Terminal tables: % over every row once, so unresolved reads stay in the
  # denominator instead of inflating the classified taxa.
  rel_terminal <- if (table_type == "terminal") .relativize_taxa_by_sample(mat)

  out_list  <- list()
  done_lvls <- character(0)
  for (lv in tax_levels) {
    ti  <- match(lv, all_levels)
    sel <- !is.na(lin[, ti]) & (table_type == "terminal" | depth == ti)
    if (!any(sel)) { message(sprintf("No rows resolve to %s level - skipped.", lv)); next }
    keys <- .lineage_keys(lin[sel, seq_len(ti), drop = FALSE], lineage_sep)

    if (table_type == "terminal") {
      lvl_mat <- rowsum(rel_terminal[sel, , drop = FALSE], group = keys, reorder = FALSE)
      lost <- 1 - sum(lvl_mat) / sum(rel_terminal)
      if (is.finite(lost) && lost > 0.001)
        message(sprintf(paste0("%s: %.1f%% of total relative abundance is not ",
                               "classified down to this level and is not listed."),
                        lv, lost * 100))
    } else {
      lvl_mat <- .relativize_taxa_by_sample(
        rowsum(mat[sel, , drop = FALSE], group = keys, reorder = FALSE))
    }
    lvl_mat <- lvl_mat[sort(rownames(lvl_mat), method = "radix"), , drop = FALSE]
    if (isTRUE(drop_all_zero)) {
      lvl_mat <- lvl_mat[apply(lvl_mat, 1, max) > 0, , drop = FALSE]
      if (!nrow(lvl_mat)) next
    }

    for (g in group_levels) {
      sub <- lvl_mat[, grp_idx[[g]], drop = FALSE]
      # mean, SD and t-based 95% CI of the mean (NA with a single sample)
      n_ <- ncol(sub); m_ <- rowMeans(sub)
      s_ <- if (n_ > 1) apply(sub, 1, stats::sd) else rep(NA_real_, nrow(sub))
      h_ <- if (n_ > 1) stats::qt(0.975, n_ - 1) * s_ / sqrt(n_) else rep(NA_real_, nrow(sub))
      out_list[[length(out_list) + 1]] <- data.frame(
        Level            = lv,
        Taxon            = rownames(lvl_mat),
        Group            = g,
        N_samples        = length(grp_idx[[g]]),
        median_abundance = round(apply(sub, 1, stats::median), digits),
        min_abundance    = round(apply(sub, 1, min), digits),
        max_abundance    = round(apply(sub, 1, max), digits),
        mean_abundance   = round(m_, digits),
        sd_abundance     = round(s_, digits),
        ci95_low         = round(m_ - h_, digits),
        ci95_high        = round(m_ + h_, digits),
        stringsAsFactors = FALSE)
    }
    done_lvls <- c(done_lvls, lv)
  }

  if (!length(out_list))
    stop("Nothing to summarise: no taxa resolved at the requested level(s). ",
         "Check that clade names use rank tags like 'p__Firmicutes'.", call. = FALSE)

  summary_df <- do.call(rbind, out_list)
  summary_df$range <- paste0(.fmt2(summary_df$min_abundance), "-",
                             .fmt2(summary_df$max_abundance))
  summary_df <- summary_df[, c("Level", "Taxon", "Group", "N_samples", "median_abundance",
                               "min_abundance", "max_abundance", "range", "mean_abundance",
                               "sd_abundance", "ci95_low", "ci95_high")]
  rownames(summary_df) <- NULL
  list(summary = summary_df, group_levels = group_levels,
       tax_levels = done_lvls, table_type = table_type)
}

# -- 2. export_taxa_descriptive_workbook() ------------------------------------
#' Export the per-group, all-level taxa summary to a formatted Excel sheet
#'
#' One sheet; each taxonomic level forms a block (highest rank first) with one
#' merged header per group spanning median / min / max / range / mean / SD / 95\% CI, and one row
#' per taxon labelled with its full lineage.
#'
#' @inheritParams summarize_taxa_by_group
#' @param outfile Path of the .xlsx file to write.
#' @param sheet_name Worksheet name (max 31 characters).
#' @param range_as_formula Logical; write \code{range} as a live Excel formula
#'   instead of text. Off by default: a formula shows as an empty cell until
#'   Excel recalculates it, which it does not do in Protected View or previews.
#' @param font_name,font_size Base font for the sheet.
#' @param zoom Sheet zoom in percent.
#' @param taxon_col_width Width of the Taxon column.
#' @param ... Passed to \code{summarize_taxa_by_group()}.
#' @return Invisibly, the summary list plus \code{outfile}.
#' @export
export_taxa_descriptive_workbook <- function(otu,
                                             meta,
                                             group_columns,
                                             outfile          = "taxa_summary_Descriptive.xlsx",
                                             sheet_name       = "ALL_LEVELS_abundance_by_group",
                                             range_as_formula = FALSE,
                                             font_name        = "Aptos Narrow",
                                             font_size        = 11,
                                             zoom             = 70,
                                             taxon_col_width  = 26.14,
                                             ...) {
  res <- summarize_taxa_by_group(otu = otu, meta = meta, group_columns = group_columns, ...)
  sm      <- res$summary
  groups  <- res$group_levels
  n_g     <- length(groups)
  stat_hd <- c("median_abundance", "min_abundance", "max_abundance", "range",
               "mean_abundance", "sd_abundance", "ci95_low", "ci95_high")
  k       <- length(stat_hd)
  n_cols  <- 1 + k * n_g

  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontSize = font_size, fontName = font_name)
  openxlsx::addWorksheet(wb, sheet_name, zoom = zoom)
  st_center <- openxlsx::createStyle(halign = "center")
  st_num    <- openxlsx::createStyle(halign = "center", numFmt = "0.00")
  st_top    <- openxlsx::createStyle(border = "top",    borderStyle = "thin")
  st_bottom <- openxlsx::createStyle(border = "bottom", borderStyle = "thin")

  # openxlsx adds 0.71 padding to widths, so subtract it to match the reference
  openxlsx::setColWidths(wb, sheet_name, cols = 1, widths = taxon_col_width - 0.71)
  openxlsx::setColWidths(wb, sheet_name, cols = 2:n_cols,
                         widths = rep(c(20.71, 17.14, 17.71, 13.43, 18.71, 15.71, 12.43, 12.43),
                                      times = n_g) - 0.71)

  row <- 1
  for (lv in res$tax_levels) {
    d     <- sm[sm$Level == lv, , drop = FALSE]
    taxa  <- unique(d$Taxon)
    hdr_row <- row + 1; col_row <- row + 2
    first_row <- row + 3; last_row <- row + 2 + length(taxa)

    openxlsx::writeData(wb, sheet_name, lv, startCol = 1, startRow = row, colNames = FALSE)
    openxlsx::writeData(wb, sheet_name, t(c("Taxon", rep(stat_hd, times = n_g))),
                        startCol = 1, startRow = col_row, colNames = FALSE)
    openxlsx::writeData(wb, sheet_name, taxa, startCol = 1, startRow = first_row,
                        colNames = FALSE)

    for (gi in seq_len(n_g)) {
      c0 <- 2 + (gi - 1) * k
      openxlsx::writeData(wb, sheet_name, groups[gi], startCol = c0, startRow = hdr_row,
                          colNames = FALSE)
      openxlsx::mergeCells(wb, sheet_name, cols = c0:(c0 + k - 1), rows = hdr_row)
      dg <- d[d$Group == groups[gi], , drop = FALSE]
      dg <- dg[match(taxa, dg$Taxon), , drop = FALSE]
      openxlsx::writeData(wb, sheet_name,
                          as.matrix(dg[, c("median_abundance", "min_abundance", "max_abundance")]),
                          startCol = c0, startRow = first_row, colNames = FALSE, rowNames = FALSE)
      if (isTRUE(range_as_formula)) {
        r <- first_row:last_row
        openxlsx::writeFormula(wb, sheet_name, startCol = c0 + 3, startRow = first_row,
          x = sprintf('TEXT(%s%d,"0.00") & "-" & TEXT(%s%d,"0.00")',
                      openxlsx::int2col(c0 + 1), r, openxlsx::int2col(c0 + 2), r))
      } else {
        openxlsx::writeData(wb, sheet_name, dg$range, startCol = c0 + 3,
                            startRow = first_row, colNames = FALSE)
      }
      openxlsx::writeData(wb, sheet_name,
                          as.matrix(dg[, c("mean_abundance", "sd_abundance", "ci95_low", "ci95_high")]),
                          startCol = c0 + 4, startRow = first_row, colNames = FALSE, rowNames = FALSE)
    }

    openxlsx::addStyle(wb, sheet_name, st_center, rows = c(hdr_row, col_row),
                       cols = 2:n_cols, gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(wb, sheet_name, st_num, rows = first_row:last_row,
                       cols = 2:n_cols, gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(wb, sheet_name, st_top, rows = hdr_row, cols = 1:n_cols,
                       gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(wb, sheet_name, st_bottom, rows = c(col_row, last_row),
                       cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
    row <- last_row + 2
  }

  openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)
  message(sprintf("Workbook written to: %s", outfile))
  res$outfile <- outfile
  invisible(res)
}
