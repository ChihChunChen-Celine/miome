# descriptive_summary.R
# -----------------------------------------------------------------------------
# Descriptive relative-abundance summary (median / min / max / range) of every
# taxon, per user-chosen group, for every taxonomic level - listed from the
# highest rank (Phylum) to the lowest (Species) on ONE Excel sheet.
#
#   1. summarize_taxa_by_group()           - computes the numbers (long data frame)
#   2. export_taxa_descriptive_workbook()  - computes + writes the formatted .xlsx
#
# Like every file under R/, this file only DEFINES functions - nothing is read
# or executed when the package loads.
#
# Compatible with:
#   config.R -> TAXA_PREFIXES (rank names + prefixes), pcoa_config()$is_baseline_col
#   data.R   -> .to_relative_abundance() (same counts -> percent rule used for
#               beta diversity and MaAsLin3)
#   app.R    -> rv$obj$otu_path, rv$meta_enriched can be passed straight in
# -----------------------------------------------------------------------------

# -- Internal: read an OTU/taxa table -> numeric matrix (taxa rows, samples cols)
# Accepts a file path (read exactly as data.R does) or a ready-made data frame /
# matrix with taxa as rows and samples as columns.
.read_taxa_matrix <- function(otu) {
  if (is.character(otu) && length(otu) == 1) {
    tab <- read.table(
      otu,
      sep         = "\t",
      header      = TRUE,
      row.names   = 1,
      check.names = FALSE,
      quote       = "",
      fill        = TRUE
    )
  } else if (is.data.frame(otu) || is.matrix(otu)) {
    tab <- as.data.frame(otu, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop("'otu' must be a file path or a data frame / matrix (taxa as rows, ",
         "samples as columns).", call. = FALSE)
  }

  colnames(tab) <- gsub("_metaphlan$", "", colnames(tab))

  mat <- as.matrix(tab)
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (anyNA(mat)) {
    warning("Non-numeric or missing values in the OTU table were set to 0.",
            call. = FALSE)
    mat[is.na(mat)] <- 0
  }
  mat
}

# -- Internal: read metadata (path or data frame) -> data frame, IDs = rownames
.read_meta_frame <- function(meta) {
  if (is.character(meta) && length(meta) == 1) {
    meta <- read.csv(meta, row.names = 1, stringsAsFactors = FALSE)
  } else if (!is.data.frame(meta)) {
    stop("'meta' must be a file path (.csv, first column = Sample ID) or a ",
         "data frame with Sample IDs as row names.", call. = FALSE)
  }
  data.frame(
    lapply(meta, function(x) if (is.character(x)) trimws(x) else x),
    stringsAsFactors = FALSE,
    row.names        = rownames(meta),
    check.names      = FALSE
  )
}

# -- Internal: split every clade string into one token per taxonomic rank -----
# Returns a character matrix (one row per clade, one column per rank in
# TAXA_PREFIXES order: Domain ... Species). Each cell holds the ORIGINAL tag
# (e.g. "p__Firmicutes", "k__Bacteria") or NA if that rank is absent/empty.
# Accepts ";" or "|" as the separator, and "k__" as an alias for "d__".
.lineage_matrix <- function(clades) {
  rank_letters <- unname(TAXA_PREFIXES)
  out <- matrix(
    NA_character_, nrow = length(clades), ncol = length(rank_letters),
    dimnames = list(NULL, names(TAXA_PREFIXES))
  )
  for (i in seq_along(clades)) {
    parts <- trimws(strsplit(as.character(clades[i]), "[;|]")[[1]])
    for (p in parts) {
      m <- regmatches(p, regexec("^([A-Za-z])__(.*)$", p))[[1]]
      if (length(m) != 3) next
      rk <- tolower(m[2])
      if (rk == "k") rk <- "d"
      j <- match(rk, rank_letters)
      if (is.na(j)) next
      val <- trimws(m[3])
      if (!nzchar(val) || val == "_") next
      out[i, j] <- paste0(m[2], "__", val)
    }
  }
  out
}

# -- Internal: join the resolved tokens of each row into one lineage string ---
.lineage_keys <- function(lin, sep) {
  apply(lin, 1, function(r) paste(r[!is.na(r)], collapse = sep))
}

# -- Internal: relative-abundance (%) for a taxa-by-samples matrix ------------
# Thin wrapper around data.R's .to_relative_abundance() (samples-as-rows), so
# the counts-vs-percent rule is identical to the one used everywhere else in
# the package: samples whose total is already ~100 are left as-is, all others
# are treated as raw counts and rescaled to 100.
.relativize_taxa_by_sample <- function(mat, quiet = FALSE) {
  conv <- if (quiet) suppressMessages(.to_relative_abundance(t(mat)))
  else .to_relative_abundance(t(mat))
  out <- t(conv)
  attr(out, "scale_info") <- NULL
  out
}

# -- Internal: format a number the way Excel's TEXT(x,"0.00") does -----------
.fmt2 <- function(x) formatC(round(x + 1e-9, 2), format = "f", digits = 2)

# -- 1. summarize_taxa_by_group() ---------------------------------------------
#' Descriptive relative-abundance summary per group, for every taxonomic level
#'
#' For each requested taxonomic level (Phylum -> Species), computes the
#' median, minimum, maximum and range of each taxon's relative abundance (%)
#' within every group defined by one or more metadata columns.
#'
#' Relative abundance is obtained with the same rule used elsewhere in the
#' package (\code{.to_relative_abundance()} in data.R): a sample whose total is
#' already ~100 is treated as percent and left alone; any other sample is
#' treated as raw counts and rescaled to 100.
#'
#' Two table layouts are supported (see \code{table_type}):
#' \itemize{
#'   \item \strong{terminal}: rows are the deepest classification only (e.g. a
#'     species-only table, or a QIIME2 L7 export). Higher levels are built by
#'     summing the rows beneath each taxon.
#'   \item \strong{all_ranks}: the table already contains a row for every rank
#'     (typical of an un-filtered MetaPhlAn merged table). Each level is then
#'     taken from the rows at exactly that rank, so nothing is double counted.
#' }
#'
#' @param otu  OTU/taxa table: a file path (.tsv, taxa as rows, samples as
#'   columns - same file the rest of the app uses) or a data frame / matrix.
#' @param meta Metadata: a file path (.csv, first column = Sample ID) or a
#'   data frame with Sample IDs as row names (e.g. \code{rv$meta_enriched}).
#' @param group_columns Character vector of one or more metadata columns that
#'   define the groups. Multiple columns are pasted together with
#'   \code{group_sep} to form the group label (e.g.
#'   \code{c("Treatment", "Timepoint")} -> \code{"Treatment1_D042(Initial)"}).
#' @param sort_columns Metadata columns that control the ORDER of the groups
#'   (default: same as \code{group_columns}). Use e.g.
#'   \code{c("Timepoint", "Treatment")} to list every treatment within the
#'   first timepoint before moving on to the next timepoint.
#' @param group_levels Optional character vector giving the exact order of the
#'   final group labels; overrides \code{sort_columns}.
#' @param group_sep Separator used when pasting several group columns
#'   (default \code{"_"}).
#' @param tax_levels Taxonomic levels to summarise, any of
#'   \code{names(TAXA_PREFIXES)}. Always emitted from highest to lowest rank.
#'   Default: Phylum, Class, Order, Family, Genus, Species.
#' @param exclude_baseline Logical; drop samples flagged \code{TRUE} in
#'   \code{is_baseline_col} (default \code{FALSE}: the groups you choose decide
#'   what is included).
#' @param is_baseline_col Baseline flag column (default \code{"IsBaseline"}).
#' @param table_type \code{"auto"} (default; detected from the table),
#'   \code{"terminal"} or \code{"all_ranks"}.
#' @param drop_all_zero Logical; remove taxa that are zero in every selected
#'   sample (default \code{FALSE}).
#' @param digits Decimal places kept for median / min / max (default 3).
#' @param lineage_sep Separator used in the output taxon labels (default
#'   \code{"|"}, e.g. \code{k__Bacteria|p__Firmicutes}).
#' @return A list with:
#'   \item{summary}{Long data frame: Level, Taxon, Group, N_samples,
#'     median_abundance, min_abundance, max_abundance, range.}
#'   \item{group_levels}{Ordered group labels.}
#'   \item{tax_levels}{Levels actually produced, highest -> lowest.}
#'   \item{table_type}{Layout used (\code{"terminal"} or \code{"all_ranks"}).}
#' @importFrom stats median
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

  # -- Validate levels ---------------------------------------------
  all_levels <- names(TAXA_PREFIXES)
  bad_lv <- setdiff(tax_levels, all_levels)
  if (length(bad_lv) > 0)
    stop(sprintf("Unknown taxonomic level(s): %s. Choose from: %s",
                 paste(bad_lv, collapse = ", "),
                 paste(all_levels, collapse = ", ")), call. = FALSE)
  tax_levels <- all_levels[all_levels %in% tax_levels]   # high -> low order

  if (missing(group_columns) || length(group_columns) == 0)
    stop("'group_columns' must name at least one metadata column.",
         call. = FALSE)

  # -- Load data ---------------------------------------------------
  mat  <- .read_taxa_matrix(otu)
  meta <- .read_meta_frame(meta)

  need <- unique(c(group_columns, sort_columns))
  miss <- setdiff(need, colnames(meta))
  if (length(miss) > 0)
    stop(sprintf("Column(s) not found in metadata: %s\nAvailable columns: %s",
                 paste(miss, collapse = ", "),
                 paste(colnames(meta), collapse = ", ")), call. = FALSE)

  # -- Align samples -----------------------------------------------
  shared <- intersect(rownames(meta), colnames(mat))
  if (length(shared) == 0)
    stop(paste0(
      "No matching Sample IDs between OTU table and metadata.\n",
      "  OTU table columns (first 5): ",
      paste(head(colnames(mat), 5), collapse = ", "), "\n",
      "  Metadata rows     (first 5): ",
      paste(head(rownames(meta), 5), collapse = ", ")), call. = FALSE)
  meta <- meta[shared, , drop = FALSE]

  # -- Optional baseline exclusion + drop NA group values ----------
  if (isTRUE(exclude_baseline) && is_baseline_col %in% colnames(meta)) {
    b    <- as.logical(meta[[is_baseline_col]])
    meta <- meta[!(b %in% TRUE), , drop = FALSE]
    message(sprintf("Baseline samples excluded. N remaining = %d", nrow(meta)))
  }
  keep <- stats::complete.cases(meta[, need, drop = FALSE])
  if (any(!keep))
    message(sprintf("%d sample(s) dropped due to NA in: %s",
                    sum(!keep), paste(need, collapse = ", ")))
  meta <- meta[keep, , drop = FALSE]
  if (nrow(meta) == 0)
    stop("No samples remaining after filtering.", call. = FALSE)

  # -- Group labels and their order --------------------------------
  # method = "radix" -> locale-independent (C-locale) ordering of text.
  ord  <- do.call(order, c(unname(as.list(meta[sort_columns])),
                           list(method = "radix")))
  meta <- meta[ord, , drop = FALSE]
  grp  <- do.call(paste, c(unname(lapply(meta[group_columns], as.character)),
                           list(sep = group_sep)))
  found <- unique(grp)

  if (is.null(group_levels)) {
    group_levels <- found
  } else {
    group_levels <- as.character(group_levels)
    not_listed   <- setdiff(found, group_levels)
    if (length(not_listed) > 0)
      stop(sprintf("'group_levels' is missing group(s) present in the data: %s",
                   paste(not_listed, collapse = ", ")), call. = FALSE)
    extra <- setdiff(group_levels, found)
    if (length(extra) > 0)
      warning(sprintf("'group_levels' entries not found in the data were ignored: %s",
                      paste(extra, collapse = ", ")), call. = FALSE)
    group_levels <- group_levels[group_levels %in% found]
  }

  samples  <- rownames(meta)
  mat      <- mat[, samples, drop = FALSE]
  grp_idx  <- split(seq_along(samples), factor(grp, levels = group_levels))
  message(sprintf("Groups (%d): %s", length(group_levels),
                  paste(sprintf("%s [n=%d]", group_levels,
                                lengths(grp_idx)), collapse = "; ")))

  # -- Parse lineages ----------------------------------------------
  clades <- rownames(mat)

  # MetaPhlAn strain/SGB rows (t__...) repeat their parent species' abundance;
  # keeping them would double count.
  is_strain <- grepl("(^|[;|])t__", clades)
  if (any(is_strain)) {
    message(sprintf("Ignoring %d strain-level (t__) row(s) to avoid double counting.",
                    sum(is_strain)))
    mat    <- mat[!is_strain, , drop = FALSE]
    clades <- rownames(mat)
  }

  lin   <- .lineage_matrix(clades)
  depth <- apply(lin, 1, function(r) {
    w <- which(!is.na(r)); if (length(w)) max(w) else 0L
  })

  # -- Detect table layout -----------------------------------------
  if (table_type == "auto") {
    tok_list    <- lapply(seq_len(nrow(lin)), function(i) lin[i, !is.na(lin[i, ])])
    keys_full   <- vapply(tok_list, paste, character(1), collapse = lineage_sep)
    parent_keys <- vapply(tok_list, function(t)
      if (length(t) >= 2) paste(t[-length(t)], collapse = lineage_sep)
      else NA_character_, character(1))
    has_parent  <- parent_keys[!is.na(parent_keys)] %in% keys_full
    table_type  <- if (length(has_parent) > 0 && mean(has_parent) >= 0.9)
      "all_ranks" else "terminal"
    message(sprintf("Table layout detected: %s.", table_type))
  }

  # Terminal tables: convert to % ONCE over every row, so the denominator is
  # each sample's whole community (reads not resolved down to a level stay in
  # the denominator instead of inflating the classified taxa).
  rel_terminal <- if (table_type == "terminal") .relativize_taxa_by_sample(mat)
  else NULL

  # -- Summarise level by level (highest -> lowest) ----------------
  out_list   <- list()
  done_lvls  <- character(0)

  for (lv in tax_levels) {
    ti     <- match(lv, all_levels)
    tok_ok <- !is.na(lin[, ti])
    sel    <- if (table_type == "all_ranks") tok_ok & depth == ti else tok_ok

    if (!any(sel)) {
      message(sprintf("No rows resolve to %s level - skipped.", lv))
      next
    }

    keys <- .lineage_keys(lin[sel, seq_len(ti), drop = FALSE], lineage_sep)

    if (table_type == "terminal") {
      lvl_mat <- rowsum(rel_terminal[sel, , drop = FALSE], group = keys,
                        reorder = FALSE)
      lost <- 1 - sum(lvl_mat) / sum(rel_terminal)
      if (is.finite(lost) && lost > 0.001)
        message(sprintf(
          "%s: %.1f%% of total relative abundance is not classified down to this level and is not listed.",
          lv, lost * 100))
    } else {
      lvl_mat <- rowsum(mat[sel, , drop = FALSE], group = keys, reorder = FALSE)
      # report the counts-vs-percent check once, not once per level
      lvl_mat <- .relativize_taxa_by_sample(lvl_mat,
                                            quiet = length(done_lvls) > 0)
    }

    lvl_mat <- lvl_mat[sort(rownames(lvl_mat), method = "radix"), , drop = FALSE]

    if (isTRUE(drop_all_zero)) {
      lvl_mat <- lvl_mat[apply(lvl_mat, 1, max) > 0, , drop = FALSE]
      if (nrow(lvl_mat) == 0) next
    }

    for (g in group_levels) {
      idx <- grp_idx[[g]]
      sub <- lvl_mat[, idx, drop = FALSE]
      out_list[[length(out_list) + 1]] <- data.frame(
        Level            = lv,
        Taxon            = rownames(lvl_mat),
        Group            = g,
        N_samples        = length(idx),
        median_abundance = round(apply(sub, 1, stats::median), digits),
        min_abundance    = round(apply(sub, 1, min), digits),
        max_abundance    = round(apply(sub, 1, max), digits),
        stringsAsFactors = FALSE
      )
    }
    done_lvls <- c(done_lvls, lv)
  }

  if (length(out_list) == 0)
    stop("Nothing to summarise: no taxa resolved at the requested level(s). ",
         "Check that clade names use rank tags like 'p__Firmicutes'.",
         call. = FALSE)

  summary_df <- do.call(rbind, out_list)
  summary_df$range <- paste0(.fmt2(summary_df$min_abundance), "-",
                             .fmt2(summary_df$max_abundance))
  rownames(summary_df) <- NULL

  list(
    summary      = summary_df,
    group_levels = group_levels,
    tax_levels   = done_lvls,
    table_type   = table_type
  )
}

# -- 2. export_taxa_descriptive_workbook() ------------------------------------
#' Export the per-group, all-level taxa summary to a formatted Excel sheet
#'
#' Writes ONE sheet in which each taxonomic level forms a block, listed from
#' the highest rank to the lowest. Every block has the level name, one merged
#' header per group spanning four columns (median_abundance, min_abundance,
#' max_abundance, range), and one row per taxon labelled with its full lineage,
#' sorted alphabetically so children stay grouped under their parents.
#' The \code{range} column is a live Excel formula
#' (\code{=TEXT(min,"0.00") & "-" & TEXT(max,"0.00")}).
#'
#' @inheritParams summarize_taxa_by_group
#' @param outfile Path of the .xlsx file to write.
#' @param sheet_name Worksheet name (max 31 characters).
#' @param range_as_formula Logical; write \code{range} as an Excel formula
#'   (default \code{TRUE}, as in the reference layout) or as plain text.
#' @param font_name,font_size Base font for the whole sheet.
#' @param zoom Sheet zoom in percent (default 70).
#' @param taxon_col_width Width of the Taxon column (column A).
#' @param ... Further arguments passed to \code{summarize_taxa_by_group()}
#'   (\code{sort_columns}, \code{group_levels}, \code{group_sep},
#'   \code{tax_levels}, \code{exclude_baseline}, \code{is_baseline_col},
#'   \code{table_type}, \code{drop_all_zero}, \code{digits},
#'   \code{lineage_sep}).
#' @return Invisibly, the list returned by \code{summarize_taxa_by_group()}
#'   plus \code{outfile}.
#' @importFrom openxlsx createWorkbook addWorksheet writeData writeFormula
#' @importFrom openxlsx mergeCells createStyle addStyle setColWidths
#' @importFrom openxlsx saveWorkbook modifyBaseFont int2col
#' @export
#' @examples
#' \dontrun{
#' # From the Shiny session objects:
#' export_taxa_descriptive_workbook(
#'   otu           = rv$obj$otu_path,
#'   meta          = rv$meta_enriched,
#'   group_columns = c("Treatment", "Timepoint"),
#'   sort_columns  = c("Timepoint", "Treatment"),
#'   outfile       = "taxa_summary_Descriptive.xlsx"
#' )
#' }
export_taxa_descriptive_workbook <- function(otu,
                                             meta,
                                             group_columns,
                                             outfile          = "taxa_summary_Descriptive.xlsx",
                                             sheet_name       = "ALL_LEVELS_abundance_by_group",
                                             range_as_formula = TRUE,
                                             font_name        = "Aptos Narrow",
                                             font_size        = 11,
                                             zoom             = 70,
                                             taxon_col_width  = 26.14,
                                             ...) {

  res <- summarize_taxa_by_group(
    otu           = otu,
    meta          = meta,
    group_columns = group_columns,
    ...
  )

  sm      <- res$summary
  groups  <- res$group_levels
  n_g     <- length(groups)
  n_cols  <- 1 + 4 * n_g
  stat_hd <- c("median_abundance", "min_abundance", "max_abundance", "range")

  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontSize = font_size, fontName = font_name)
  openxlsx::addWorksheet(wb, sheet_name, zoom = zoom)

  st_center <- openxlsx::createStyle(halign = "center")
  st_num    <- openxlsx::createStyle(halign = "center", numFmt = "0.00")
  st_top    <- openxlsx::createStyle(border = "top",    borderStyle = "thin")
  st_bottom <- openxlsx::createStyle(border = "bottom", borderStyle = "thin")

  # Column widths: Taxon column + the repeating 4-column group pattern.
  # Values are the raw widths stored in the reference file; openxlsx adds a
  # 0.71 padding to whatever it is given, so subtract it to land on the same
  # width.
  openxlsx::setColWidths(wb, sheet_name, cols = 1,
                         widths = taxon_col_width - 0.71)
  openxlsx::setColWidths(
    wb, sheet_name, cols = 2:n_cols,
    widths = rep(c(20.71, 17.14, 17.71, 13.43), times = n_g) - 0.71
  )

  row <- 1
  for (lv in res$tax_levels) {

    d     <- sm[sm$Level == lv, , drop = FALSE]
    taxa  <- unique(d$Taxon)
    n_tax <- length(taxa)

    hdr_row   <- row + 1          # merged group labels
    col_row   <- row + 2          # Taxon / median / min / max / range
    first_row <- row + 3
    last_row  <- row + 2 + n_tax

    # Level title (column A)
    openxlsx::writeData(wb, sheet_name, lv, startCol = 1, startRow = row,
                        colNames = FALSE)

    # Column-header row
    openxlsx::writeData(
      wb, sheet_name,
      t(c("Taxon", rep(stat_hd, times = n_g))),
      startCol = 1, startRow = col_row, colNames = FALSE
    )

    # Taxon labels
    openxlsx::writeData(wb, sheet_name, taxa, startCol = 1,
                        startRow = first_row, colNames = FALSE)

    for (gi in seq_len(n_g)) {
      c0 <- 2 + (gi - 1) * 4

      # Merged group label across the 4 stat columns
      openxlsx::writeData(wb, sheet_name, groups[gi], startCol = c0,
                          startRow = hdr_row, colNames = FALSE)
      openxlsx::mergeCells(wb, sheet_name, cols = c0:(c0 + 3), rows = hdr_row)

      dg <- d[d$Group == groups[gi], , drop = FALSE]
      dg <- dg[match(taxa, dg$Taxon), , drop = FALSE]

      # median / min / max
      openxlsx::writeData(
        wb, sheet_name,
        as.matrix(dg[, c("median_abundance", "min_abundance", "max_abundance")]),
        startCol = c0, startRow = first_row, colNames = FALSE, rowNames = FALSE
      )

      # range
      if (isTRUE(range_as_formula)) {
        rows_i <- first_row:last_row
        f <- sprintf('TEXT(%s%d,"0.00") & "-" & TEXT(%s%d,"0.00")',
                     openxlsx::int2col(c0 + 1), rows_i,
                     openxlsx::int2col(c0 + 2), rows_i)
        openxlsx::writeFormula(wb, sheet_name, x = f,
                               startCol = c0 + 3, startRow = first_row)
      } else {
        openxlsx::writeData(wb, sheet_name, dg$range, startCol = c0 + 3,
                            startRow = first_row, colNames = FALSE)
      }
    }

    # -- Formatting for this block ---------------------------------
    all_cols <- 1:n_cols
    val_cols <- 2:n_cols
    # centred headers + numbers
    openxlsx::addStyle(wb, sheet_name, st_center, rows = c(hdr_row, col_row),
                       cols = val_cols, gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(wb, sheet_name, st_num, rows = first_row:last_row,
                       cols = val_cols, gridExpand = TRUE, stack = TRUE)
    # thin rules: above the group header, under the column header,
    # and under the last taxon
    openxlsx::addStyle(wb, sheet_name, st_top, rows = hdr_row,
                       cols = all_cols, gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(wb, sheet_name, st_bottom, rows = c(col_row, last_row),
                       cols = all_cols, gridExpand = TRUE, stack = TRUE)

    row <- last_row + 2           # one blank row between blocks
  }

  openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)
  message(sprintf("Workbook written to: %s", outfile))

  res$outfile <- outfile
  invisible(res)
}
