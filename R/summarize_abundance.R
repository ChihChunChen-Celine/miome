# summarize_abundance.R
# -----------------------------------------------------------------------------
# Summarize taxon relative abundance by group, and export counts + relative
# abundance + metadata to a single Excel workbook.
#
# NOTE ON PACKAGE STRUCTURE: files under R/ are sourced whenever the package
# is loaded (devtools::load_all(), devtools::document(), library(miome)).
# They may therefore only DEFINE functions - never read a file, call
# library(), or run analysis code at the top level. (That's what broke
# devtools::document() previously: the old version of this file read
# "C26-29_mapping_R.txt" and a hardcoded counts path immediately on source,
# which don't exist outside the original author's machine.) Every function
# below takes file paths / data frames as arguments and does nothing until
# you call it - use export_abundance_workbook() as the entry point.
#
# This works from an OTU/taxonomy table shaped like a QIIME2 "L7-species"
# export (taxa as rows, read with skip = 1, row.names = 1) rather than the
# transposed samples-as-rows matrix used in data.R, so it has its own small
# loading step instead of reusing .load_otu_meta().
# -----------------------------------------------------------------------------

# -- Internal helper: detect percent vs. counts, convert only if needed ------
# Per-column (per-sample) detection: a column whose sum is within `tol` of
# 100 is left as-is (already percent); any other column is treated as raw
# counts and converted to percent by dividing by its own column sum. This
# means a table mixing already-normalized and raw-count columns is still
# handled correctly rather than assuming the whole table is one or the
# other.
#' @keywords internal
.detect_and_relativize <- function(counts_numeric, tol = 1) {
  col_sums <- colSums(counts_numeric, na.rm = TRUE)

  zero_cols <- names(col_sums)[col_sums == 0]
  if (length(zero_cols) > 0)
    warning(sprintf(
      paste0(
        "Sample column(s) with a total of 0 (empty/all-NA) - will remain ",
        "0 after conversion: %s"
      ),
      paste(zero_cols, collapse = ", ")
    ))

  is_percent <- abs(col_sums - 100) <= tol

  message(sprintf(
    paste0(
      "detect_and_relativize(): %d of %d sample column(s) already look ",
      "like percent (colSum ~100, tol=%.1f) and were left as-is; %d ",
      "column(s) looked like raw counts and were converted to percent."
    ),
    sum(is_percent), length(col_sums), tol, sum(!is_percent & col_sums > 0)
  ))

  out <- counts_numeric
  convert_cols <- names(col_sums)[!is_percent & col_sums > 0]
  if (length(convert_cols) > 0) {
    out[convert_cols] <- sweep(
      counts_numeric[convert_cols],
      MARGIN = 2,
      STATS  = col_sums[convert_cols],
      FUN    = "/"
    ) * 100
  }
  out
}

# -- 1. summarize_abundance() -------------------------------------------------
#' Summarize taxon relative abundance by group
#'
#' Converts an OTU/taxonomy count table to relative abundance (percent) -
#' auto-detecting and skipping sample columns that are already percent (see
#' \code{.detect_and_relativize()}) - filters out very rare taxa, and
#' computes median/min/max relative abundance per taxon within each level
#' of \code{group_column}.
#'
#' @param counts     Data frame of taxon abundances, taxa as rows (row
#'   names = taxon/clade string), samples as columns. Typically read with
#'   \code{read.delim(..., skip = 1, row.names = 1, check.names = FALSE)}
#'   from a QIIME2 taxa-by-level export.
#' @param metadata   Data frame of sample metadata. Must contain a
#'   \code{sampleid} column matching \code{counts}' column names, plus
#'   \code{group_column}.
#' @param group_column Metadata column to group by (default
#'   \code{"Treatment"}).
#' @param remove_first_row Logical; drop the first row of \code{counts}
#'   before processing (e.g. a stray header/comment row). Default
#'   \code{TRUE}.
#' @param pct_tol    Tolerance (percentage points) around 100 for a sample
#'   column to be treated as "already percent" rather than raw counts.
#'   Default 1.
#' @param outfile    Optional path to write the summary as a CSV.
#'   \code{NULL} (default) skips writing.
#' @return A data frame with columns: the grouping column, \code{Taxon},
#'   \code{median_abundance}, \code{min_abundance}, \code{max_abundance},
#'   \code{range}.
#' @importFrom dplyr mutate across everything select filter rowwise ungroup group_by summarise left_join c_across where
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom utils write.csv
#' @importFrom stats median
#' @export
summarize_abundance <- function(counts,
                                metadata,
                                group_column      = "Treatment",
                                remove_first_row  = TRUE,
                                pct_tol           = 1,
                                outfile           = NULL) {

  missing_cols <- setdiff(c("sampleid", group_column), colnames(metadata))
  if (length(missing_cols) > 0)
    stop(sprintf(
      "metadata is missing column(s): %s\nAvailable metadata columns: %s",
      paste(missing_cols, collapse = ", "),
      paste(colnames(metadata), collapse = ", ")
    ))

  if (remove_first_row) counts <- counts[-1, ]

  counts_numeric <- counts %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(.))) %>%
    as.data.frame()
  rownames(counts_numeric) <- rownames(counts)

  relative_abundances <- .detect_and_relativize(counts_numeric, tol = pct_tol)

  relative_abundances <- relative_abundances %>%
    tibble::rownames_to_column("Taxon")

  relative_abundances <- relative_abundances %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      max_abund = max(dplyr::c_across(dplyr::where(is.numeric)), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(max_abund >= 0.01) %>%
    dplyr::select(-max_abund)

  long_data <- relative_abundances %>%
    tidyr::pivot_longer(
      cols      = -Taxon,
      names_to  = "sampleid",
      values_to = "relative_abundance"
    )

  long_data <- long_data %>%
    dplyr::left_join(
      metadata[, c("sampleid", group_column)],
      by = "sampleid"
    )

  summary_stats <- long_data %>%
    dplyr::group_by(.data[[group_column]], Taxon) %>%
    dplyr::summarise(
      median_abundance = stats::median(relative_abundance, na.rm = TRUE),
      min_abundance    = min(relative_abundance, na.rm = TRUE),
      max_abundance    = max(relative_abundance, na.rm = TRUE),
      .groups          = "drop"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(median_abundance, min_abundance, max_abundance), round, 3
      ),
      range = paste0(min_abundance, "-", max_abundance)
    )

  if (!is.null(outfile)) utils::write.csv(summary_stats, outfile, row.names = FALSE)

  summary_stats
}

# -- 2. get_relative_abundances() ---------------------------------------------
#' Convert an OTU/taxonomy count table to relative abundance only
#'
#' Same conversion logic as \code{summarize_abundance()} (auto-detects and
#' skips sample columns already in percent) but skips the group-summary
#' step and returns the full per-sample relative-abundance table. Used for
#' the "Relative Abundances" sheet in \code{export_abundance_workbook()}.
#'
#' @inheritParams summarize_abundance
#' @return A data frame: \code{Taxon} column plus one relative-abundance
#'   column per sample.
#' @importFrom dplyr mutate across everything
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @export
get_relative_abundances <- function(counts,
                                    remove_first_row = TRUE,
                                    pct_tol           = 1) {

  if (remove_first_row) counts <- counts[-1, ]

  counts_numeric <- counts %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(.))) %>%
    as.data.frame()
  rownames(counts_numeric) <- rownames(counts)

  relative_abundances <- .detect_and_relativize(counts_numeric, tol = pct_tol)

  relative_abundances %>%
    tibble::rownames_to_column("Taxon")
}

# -- 3. export_abundance_workbook() -------------------------------------------
#' Read an OTU/taxonomy table + metadata, summarize abundance, and export
#' an Excel workbook (Counts / Relative Abundances / Metadata sheets)
#'
#' This is the interactive entry point for the abundance-summary workflow.
#' It performs the file reads that used to run automatically when this
#' file was sourced - which is what broke \code{devtools::load_all()} /
#' \code{devtools::document()}, since files under \code{R/} must define
#' functions only and must never read files or execute analysis code at
#' the top level. Call this function explicitly with your own file paths
#' after loading the package.
#'
#' @param counts_path   Path to the taxa-by-level counts \code{.tsv}
#'   (e.g. a QIIME2 \code{L7-species-table.tsv} export). Read with
#'   \code{skip = 1, row.names = 1}.
#' @param metadata_path Path to the tab-delimited mapping file. Must
#'   contain a \code{sampleid} column matching the counts table's sample
#'   columns, plus \code{group_column} and \code{id_column}.
#' @param group_column  Metadata column to group by in the summary
#'   (default \code{"Treatment"}).
#' @param id_column     Metadata column used as a human-readable sample
#'   label in the Excel sheets (default \code{"Internal_ID"}).
#' @param remove_first_row Logical; drop the first row of the counts table
#'   before processing (default \code{TRUE}).
#' @param pct_tol       Tolerance (percentage points) around 100 for a
#'   sample column to be treated as already-percent rather than raw counts
#'   (default 1).
#' @param summary_csv   Optional path to also write the group summary as a
#'   CSV. \code{NULL} (default) skips writing.
#' @param workbook_path Path for the output \code{.xlsx} workbook
#'   (default \code{"species_counts_relabun.xlsx"}).
#' @return Invisibly, a list with \code{summary}, \code{relative_abundance},
#'   and \code{metadata} data frames (the same content written to
#'   \code{workbook_path}).
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom utils read.delim
#' @export
export_abundance_workbook <- function(counts_path,
                                      metadata_path,
                                      group_column      = "Treatment",
                                      id_column         = "Internal_ID",
                                      remove_first_row  = TRUE,
                                      pct_tol           = 1,
                                      summary_csv       = NULL,
                                      workbook_path     = "species_counts_relabun.xlsx") {

  metadata <- utils::read.delim(
    metadata_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE
  )

  counts <- utils::read.delim(
    counts_path, sep = "\t", header = TRUE, skip = 1,
    row.names = 1, check.names = FALSE
  )

  missing_id_cols <- setdiff(c("sampleid", id_column), colnames(metadata))
  if (length(missing_id_cols) > 0)
    stop(sprintf(
      "metadata is missing column(s): %s\nAvailable metadata columns: %s",
      paste(missing_id_cols, collapse = ", "),
      paste(colnames(metadata), collapse = ", ")
    ))

  summary_stats <- summarize_abundance(
    counts            = counts,
    metadata          = metadata,
    group_column      = group_column,
    remove_first_row  = remove_first_row,
    pct_tol           = pct_tol,
    outfile           = summary_csv
  )

  rel_abund <- get_relative_abundances(
    counts, remove_first_row = FALSE, pct_tol = pct_tol
  )

  counts_with_taxa <- counts %>% tibble::rownames_to_column("Taxon")

  id_mapping <- metadata %>%
    dplyr::select(sampleid, dplyr::all_of(id_column)) %>%
    tibble::column_to_rownames("sampleid") %>%
    t() %>%
    as.data.frame()

  wb <- openxlsx::createWorkbook()

  # ---------- COUNTS SHEET ----------
  openxlsx::addWorksheet(wb, "Counts")
  openxlsx::writeData(
    wb, "Counts", t(c("Taxon", colnames(counts))),
    startRow = 1, colNames = FALSE
  )
  openxlsx::writeData(
    wb, "Counts",
    cbind(data.frame(label = id_column), id_mapping),
    startRow = 2, colNames = FALSE
  )
  openxlsx::writeData(
    wb, "Counts", counts_with_taxa, startRow = 3, colNames = FALSE
  )

  # ---------- RELATIVE ABUNDANCES SHEET ----------
  openxlsx::addWorksheet(wb, "Relative Abundances")
  openxlsx::writeData(
    wb, "Relative Abundances", t(c("Taxon", colnames(counts))),
    startRow = 1, colNames = FALSE
  )
  openxlsx::writeData(
    wb, "Relative Abundances",
    cbind(data.frame(label = id_column), id_mapping),
    startRow = 2, colNames = FALSE
  )
  openxlsx::writeData(
    wb, "Relative Abundances", rel_abund, startRow = 3, colNames = FALSE
  )

  # ---------- METADATA SHEET ----------
  openxlsx::addWorksheet(wb, "Metadata")
  openxlsx::writeData(wb, "Metadata", metadata)

  openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)

  message(sprintf("Workbook written to: %s", workbook_path))

  invisible(list(
    summary            = summary_stats,
    relative_abundance = rel_abund,
    metadata           = metadata
  ))
}
