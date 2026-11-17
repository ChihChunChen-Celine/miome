# summarize_abundance.R
# -----------------------------------------------------------------------------
# Taxon relative-abundance summary by group, and a 3-sheet Excel workbook
# (Counts / Relative Abundances / Metadata).
#
# For new work prefer summarize_taxa_by_group() / export_taxa_descriptive_
# workbook() (descriptive_summary.R), which cover every taxonomic level.
# These functions are kept so existing scripts keep working.
# -----------------------------------------------------------------------------

# -- Internal: numeric copy of a counts table, optional first-row removal ----
# remove_first_row = "auto" drops the first row only when it is not numeric
# (e.g. a stray "#q2:types" row); TRUE / FALSE force the choice.
.clean_counts <- function(counts, remove_first_row = "auto") {
  first_bad <- any(is.na(suppressWarnings(as.numeric(unlist(counts[1, ])))))
  drop <- if (identical(remove_first_row, "auto")) first_bad else isTRUE(remove_first_row)
  if (drop) {
    message("Dropping first row of the counts table: ", rownames(counts)[1])
    counts <- counts[-1, , drop = FALSE]
  }
  out <- as.data.frame(lapply(counts, function(x) as.numeric(as.character(x))),
                       check.names = FALSE, row.names = rownames(counts))
  out[is.na(out)] <- 0
  out
}

# Relative abundance with the package-wide rule (utils.R), taxa as rows.
.relab_taxa_rows <- function(counts_numeric) {
  rel <- t(.to_relative_abundance(t(as.matrix(counts_numeric))))
  attr(rel, "scale_info") <- NULL
  as.data.frame(rel, check.names = FALSE)
}

#' Summarize taxon relative abundance by group
#'
#' Converts to relative abundance (\%) using the same counts-vs-percent rule
#' as the rest of the package, drops taxa whose maximum is below 0.01\%, and
#' computes median/min/max per taxon within each level of \code{group_column}.
#'
#' @param counts Data frame of taxon abundances, taxa as rows (row names =
#'   taxon), samples as columns.
#' @param metadata Data frame with a \code{sampleid} column matching the
#'   sample columns, plus \code{group_column}.
#' @param group_column Metadata column to group by.
#' @param remove_first_row \code{"auto"} (drop the first row only if it is
#'   non-numeric), \code{TRUE} or \code{FALSE}.
#' @param pct_tol Kept for compatibility (unused).
#' @param outfile Optional CSV path for the summary.
#' @return Data frame: group column, \code{Taxon}, \code{median_abundance},
#'   \code{min_abundance}, \code{max_abundance}, \code{range}.
#' @export
summarize_abundance <- function(counts, metadata, group_column = "Treatment",
                                remove_first_row = "auto", pct_tol = 1,
                                outfile = NULL) {
  miss <- setdiff(c("sampleid", group_column), colnames(metadata))
  if (length(miss))
    stop(sprintf("metadata is missing column(s): %s\nAvailable columns: %s",
                 paste(miss, collapse = ", "), paste(colnames(metadata), collapse = ", ")),
         call. = FALSE)

  rel <- .relab_taxa_rows(.clean_counts(counts, remove_first_row))
  rel <- rel[apply(rel, 1, max) >= 0.01, , drop = FALSE]

  long <- data.frame(Taxon = rep(rownames(rel), times = ncol(rel)),
                     sampleid = rep(colnames(rel), each = nrow(rel)),
                     relative_abundance = unlist(rel, use.names = FALSE),
                     stringsAsFactors = FALSE)
  long[[group_column]] <- metadata[[group_column]][match(long$sampleid, metadata$sampleid)]

  out <- long %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_column, "Taxon")))) %>%
    dplyr::summarise(median_abundance = round(stats::median(relative_abundance), 3),
                     min_abundance    = round(min(relative_abundance), 3),
                     max_abundance    = round(max(relative_abundance), 3),
                     .groups = "drop") %>%
    dplyr::mutate(range = paste0(min_abundance, "-", max_abundance)) %>%
    as.data.frame()

  if (!is.null(outfile)) utils::write.csv(out, outfile, row.names = FALSE)
  out
}

#' Convert a taxonomy count table to relative abundance (\%)
#'
#' @inheritParams summarize_abundance
#' @return Data frame: \code{Taxon} column plus one column per sample.
#' @export
get_relative_abundances <- function(counts, remove_first_row = "auto", pct_tol = 1) {
  rel <- .relab_taxa_rows(.clean_counts(counts, remove_first_row))
  cbind(Taxon = rownames(rel), rel, stringsAsFactors = FALSE)
}

#' Summarize abundance and export a Counts / Relative Abundances / Metadata
#' workbook
#'
#' @param counts_path Taxa-by-level counts \code{.tsv} (e.g. QIIME 2
#'   \code{L7} export); comment lines and a \code{#OTU ID} header are handled.
#' @param metadata_path Tab-delimited mapping file with a \code{sampleid}
#'   column, \code{group_column} and \code{id_column}.
#' @param group_column Metadata column to group by in the summary.
#' @param id_column Metadata column used as a readable sample label.
#' @param remove_first_row See \code{summarize_abundance()}.
#' @param pct_tol Kept for compatibility (unused).
#' @param summary_csv Optional CSV path for the group summary.
#' @param workbook_path Output \code{.xlsx} path.
#' @return Invisibly, a list with \code{summary}, \code{relative_abundance}
#'   and \code{metadata}.
#' @export
export_abundance_workbook <- function(counts_path, metadata_path,
                                      group_column = "Treatment",
                                      id_column = "Internal_ID",
                                      remove_first_row = "auto", pct_tol = 1,
                                      summary_csv = NULL,
                                      workbook_path = "species_counts_relabun.xlsx") {
  metadata <- utils::read.delim(metadata_path, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE)
  miss <- setdiff(c("sampleid", id_column), colnames(metadata))
  if (length(miss))
    stop(sprintf("metadata is missing column(s): %s\nAvailable columns: %s",
                 paste(miss, collapse = ", "), paste(colnames(metadata), collapse = ", ")),
         call. = FALSE)

  # Clean ONCE so every sheet and the summary use the same rows.
  counts <- .clean_counts(as.data.frame(.read_otu_table(counts_path), check.names = FALSE),
                          remove_first_row)
  summary_stats <- summarize_abundance(counts, metadata, group_column,
                                       remove_first_row = FALSE, outfile = summary_csv)
  rel_abund <- get_relative_abundances(counts, remove_first_row = FALSE)

  # Label row aligned to the counts table's sample columns (not metadata order).
  ids <- metadata[[id_column]][match(colnames(counts), metadata$sampleid)]
  header <- t(c("Taxon", colnames(counts)))
  label  <- t(c(id_column, ifelse(is.na(ids), "", as.character(ids))))

  wb <- openxlsx::createWorkbook()
  for (sh in list(list("Counts", cbind(Taxon = rownames(counts), counts)),
                  list("Relative Abundances", rel_abund))) {
    openxlsx::addWorksheet(wb, sh[[1]])
    openxlsx::writeData(wb, sh[[1]], header, startRow = 1, colNames = FALSE)
    openxlsx::writeData(wb, sh[[1]], label,  startRow = 2, colNames = FALSE)
    openxlsx::writeData(wb, sh[[1]], sh[[2]], startRow = 3, colNames = FALSE)
  }
  openxlsx::addWorksheet(wb, "Metadata")
  openxlsx::writeData(wb, "Metadata", metadata)
  openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
  message(sprintf("Workbook written to: %s", workbook_path))

  invisible(list(summary = summary_stats, relative_abundance = rel_abund,
                 metadata = metadata))
}
