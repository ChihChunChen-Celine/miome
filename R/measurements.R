# measurements.R
# -----------------------------------------------------------------------------
# Targeted measurements (qPCR, fatty acids / sterols, bile acids, ...):
#   1. read_measurements()        - read a samples x measurements table
#   2. link_measurements()        - attach it to the analysed metadata
#   3. run_measurement_models()   - LMM / LM per selected measurement
#   4. summarize_measurements()   - descriptive statistics per group
#
# The models use the same engine as alpha diversity (.fit_models_all in
# stats.R), so formula, random effects, reference group, baseline handling and
# the two cross-feature p-value adjustments behave identically.
# -----------------------------------------------------------------------------

#' Read a table of targeted measurements
#'
#' One row per sample, one column per measurement (e.g. qPCR targets,
#' fatty acids, bile acids). The first column holds an ID that also appears in
#' the metadata (the Sample ID or any other ID column such as a sequencing
#' number). Column names are tidied exactly like metadata names
#' (\code{"myristic acid (14:0) (ug/mg)"} becomes
#' \code{myristic.acid.14.0.ug.mg}); the names as typed are kept as labels.
#'
#' @param path Path to a .csv (or tab-separated) file.
#' @return Data frame of the numeric columns, IDs as row names, with attribute
#'   \code{"labels"} (original column names).
#' @export
read_measurements <- function(path) {
  x   <- .read_meta(path)
  lab <- attr(x, "labels")
  num <- vapply(x, is.numeric, logical(1))
  # text columns that are numbers apart from a few entries such as "<LOD", "n.d."
  for (v in names(x)[!num]) {
    val <- suppressWarnings(as.numeric(x[[v]]))
    if (sum(!is.na(val)) >= max(3, 0.5 * sum(!is.na(x[[v]])))) {
      message(sprintf("'%s': %d non-numeric entr%s (e.g. '%s') treated as missing.", v,
                      sum(is.na(val) & !is.na(x[[v]])),
                      if (sum(is.na(val) & !is.na(x[[v]])) == 1) "y" else "ies",
                      x[[v]][is.na(val) & !is.na(x[[v]])][1]))
      x[[v]] <- val; num[v] <- TRUE
    }
  }
  if (any(!num))
    message("Ignoring non-numeric column(s): ", paste(names(x)[!num], collapse = ", "))
  x <- x[, num, drop = FALSE]
  x <- x[, colSums(!is.na(x)) > 0, drop = FALSE]
  if (!ncol(x)) stop("No numeric measurement columns found in the file.", call. = FALSE)
  attr(x, "labels") <- lab[colnames(x)]
  x
}

#' Attach measurements to the metadata
#'
#' Finds how the measurement IDs relate to the metadata: they may be the
#' Sample IDs themselves or the values of any metadata column that identifies
#' one sample per row (e.g. a sequencing number). Only samples present in
#' \code{meta} are kept, so the measurements are analysed on the same samples
#' as everything else.
#'
#' @param meas Data frame from \code{read_measurements()}.
#' @param meta Metadata data frame (Sample IDs as row names), e.g.
#'   \code{pcoa_obj$meta}.
#' @return A list: \code{data} (metadata + measurement columns), \code{features},
#'   \code{labels}, and \code{link} (how the files were matched, and which rows
#'   or samples had no partner).
#' @export
link_measurements <- function(meas, meta) {
  ids  <- rownames(meas)
  cand <- c(list(`Sample ID` = rownames(meta)), lapply(meta, function(v) trimws(as.character(v))))
  score <- vapply(cand, function(v) {
    hit <- .match_ids(ids, v)
    # every matched value must point to exactly one metadata row
    if (anyDuplicated(v[v %in% hit[!is.na(hit)]])) 0 else sum(!is.na(hit))
  }, numeric(1))
  if (max(score) == 0)
    stop("The IDs in the measurement file (first column: ", paste(utils::head(ids, 3), collapse = ", "),
         ", ...) do not match the Sample IDs or any ID column of the metadata.", call. = FALSE)
  by   <- names(score)[which.max(score)]
  hit  <- .match_ids(ids, cand[[by]])
  rows <- match(hit[!is.na(hit)], cand[[by]])
  m    <- meas[!is.na(hit), , drop = FALSE]

  out <- meta[rows, , drop = FALSE]
  clash <- intersect(colnames(m), colnames(out))
  if (length(clash))
    message(sprintf("%d measurement column(s) also exist in the metadata; the values from the measurement file are used: %s",
                    length(clash), paste(utils::head(clash, 6), collapse = ", ")))
  out[clash] <- NULL
  out <- cbind(out, m)
  rownames(out) <- rownames(meta)[rows]

  link <- list(by = by, n_rows = length(ids), n_linked = nrow(out),
               unlinked_rows = ids[is.na(hit)],
               samples_without = setdiff(rownames(meta), rownames(out)), clash = clash)
  message(sprintf("Measurements linked through '%s': %d of %d rows matched an analysed sample.",
                  by, link$n_linked, link$n_rows))
  lab <- attr(meas, "labels")
  if (is.null(lab)) lab <- stats::setNames(colnames(meas), colnames(meas))
  list(data = out, features = colnames(m), labels = lab[colnames(m)], link = link)
}

# Baseline value of each feature per subject -> "Baseline_<feature>" on the
# subject's post-baseline rows (same idea as the alpha-diversity covariates).
.add_feature_baselines <- function(data, features, subject_id_col, is_baseline_col) {
  if (!nzchar(subject_id_col) || !subject_id_col %in% colnames(data)) return(data)
  base <- .is_baseline(data, is_baseline_col)
  if (!any(base) || all(base)) return(data)
  b    <- data[base, , drop = FALSE]
  b    <- b[!duplicated(b[[subject_id_col]]), , drop = FALSE]
  idx  <- match(data[[subject_id_col]], b[[subject_id_col]])
  idx[base] <- NA
  for (f in features) {
    v <- b[[f]][idx]
    if (any(!is.na(v))) data[[paste0("Baseline_", f)]] <- v
  }
  data
}

#' Linear (mixed) models for selected measurements
#'
#' Fits one model per selected measurement, exactly like the alpha-diversity
#' models: fixed-effects formula, optional random intercept(s) (none = ordinary
#' linear model), treatment reference group, baseline samples excluded and -
#' when a subject and baseline column are given - each animal's baseline value
#' of that measurement added as a covariate.
#'
#' P-values are adjusted across \strong{the measurements you selected}, so
#' select the ones you planned to test before looking at results.
#'
#' @param data Data frame from \code{link_measurements()$data}.
#' @param features Measurement columns to analyse.
#' @param fixed_formula Fixed effects without \code{~}, e.g. \code{"Group + Sex"}.
#' @param random_effect Column(s) for random intercepts, or \code{NULL}.
#' @param treatment_ref,treatment_col Reference group and its column.
#' @param is_baseline_col,subject_id_col Baseline flag and subject columns
#'   (\code{""} if the study has none).
#' @param baseline_adjust Logical; add \code{Baseline_<measurement>} when
#'   baseline samples exist.
#' @param transform \code{"none"} or \code{"log10"}. With \code{"log10"},
#'   zeros are replaced by half the smallest positive value of that
#'   measurement; estimates are then differences on the log10 scale.
#' @param padj_method p-value adjustment method.
#' @param fdr_target \code{"auto"}, \code{"main"} or \code{"interaction"}; see
#'   \code{run_alpha_lmm_all_metrics()}.
#' @param labels Optional named vector of display names for the measurements.
#' @return Data frame, one row per model term per measurement, with
#'   \code{p_adjusted} (within one model), \code{p_adj_across_features} (one
#'   comparison adjusted over the selected measurements) and
#'   \code{p_adj_across_features_and_groups} (all treatment comparisons of all
#'   selected measurements as one family).
#' @export
run_measurement_models <- function(data, features, fixed_formula,
                                   random_effect   = NULL,
                                   treatment_ref   = NULL,
                                   treatment_col   = "Treatment",
                                   is_baseline_col = "",
                                   subject_id_col  = "",
                                   baseline_adjust = TRUE,
                                   transform       = c("none", "log10"),
                                   padj_method     = "BH",
                                   fdr_target      = c("auto", "main", "interaction"),
                                   labels          = NULL) {
  transform  <- match.arg(transform)
  fdr_target <- match.arg(fdr_target)
  features   <- unique(features)
  miss <- setdiff(features, colnames(data))
  if (length(miss)) stop("Measurement(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!length(features)) stop("Select at least one measurement.", call. = FALSE)
  notnum <- features[!vapply(data[features], is.numeric, logical(1))]
  if (length(notnum)) stop("Not numeric: ", paste(notnum, collapse = ", "), call. = FALSE)

  skipped <- character(0)
  if (transform == "log10") {
    for (f in features) {
      x <- data[[f]]
      if (any(x < 0, na.rm = TRUE)) { skipped <- c(skipped, f); next }
      pos <- x[!is.na(x) & x > 0]
      if (!length(pos)) { skipped <- c(skipped, f); next }
      x[!is.na(x) & x == 0] <- min(pos) / 2
      data[[f]] <- log10(x)
    }
    if (length(skipped))
      warning("log10 not possible (negative or no positive values); left out: ",
              paste(skipped, collapse = ", "), call. = FALSE)
  }
  use <- setdiff(features, skipped)
  if (!length(use)) stop("No measurement left to model after the transformation.", call. = FALSE)
  if (isTRUE(baseline_adjust))
    data <- .add_feature_baselines(data, use, subject_id_col, is_baseline_col)

  out <- .fit_models_all(data, use, fixed_formula, random_effect, treatment_ref,
                         treatment_col, is_baseline_col, padj_method,
                         add_baseline = isTRUE(baseline_adjust),
                         fdr_target = fdr_target)
  att <- attributes(out)[c("fdr_target", "fdr_terms")]
  names(out)[names(out) == "Metric"] <- "Feature"
  names(out)[names(out) == "p_adj_across_metrics"] <- "p_adj_across_features"
  names(out)[names(out) == "p_adj_across_metrics_and_groups"] <- "p_adj_across_features_and_groups"
  out$Label     <- if (is.null(labels)) out$Feature else
    ifelse(is.na(labels[out$Feature]), out$Feature, labels[out$Feature])
  out$Transform <- transform
  out <- out[, c("Feature", "Label", setdiff(names(out), c("Feature", "Label")))]
  attr(out, "fdr_target") <- att$fdr_target
  attr(out, "fdr_terms")  <- att$fdr_terms
  attr(out, "n_features") <- length(use)
  attr(out, "skipped")    <- skipped
  out
}

#' Descriptive statistics of measurements per group
#'
#' For every selected measurement and every group: number of samples, mean,
#' standard deviation, median, minimum, maximum and the 95\% confidence
#' interval of the mean (t-distribution). Values are summarised as uploaded
#' (no transformation).
#'
#' @param data Data frame from \code{link_measurements()$data}.
#' @param features Measurement columns to summarise.
#' @param group_columns One or more metadata columns defining the groups;
#'   several are joined with \code{group_sep}.
#' @param group_sep Separator for joined group labels.
#' @param exclude_baseline Logical; drop baseline samples first.
#' @param is_baseline_col Baseline flag column.
#' @param labels Optional named vector of display names.
#' @param digits Significant digits kept.
#' @return Long data frame: Feature, Label, Group, N, N_missing, mean, sd,
#'   median, min, max, ci95_low, ci95_high.
#' @export
summarize_measurements <- function(data, features, group_columns, group_sep = "_",
                                   exclude_baseline = FALSE, is_baseline_col = "",
                                   labels = NULL, digits = 4) {
  miss <- setdiff(c(features, group_columns), colnames(data))
  if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!length(features) || !length(group_columns))
    stop("Select at least one measurement and one grouping column.", call. = FALSE)
  if (isTRUE(exclude_baseline)) data <- data[!.is_baseline(data, is_baseline_col), , drop = FALSE]
  data <- data[stats::complete.cases(data[, group_columns, drop = FALSE]), , drop = FALSE]
  if (!nrow(data)) stop("No samples left to summarise.", call. = FALSE)
  grp <- do.call(paste, c(unname(lapply(data[group_columns], as.character)), list(sep = group_sep)))
  lv  <- unique(grp[do.call(order, c(unname(as.list(data[group_columns])), list(method = "radix")))])

  one <- function(x) {
    n <- sum(!is.na(x)); x <- x[!is.na(x)]
    m <- if (n) mean(x) else NA_real_
    s <- if (n > 1) stats::sd(x) else NA_real_
    h <- if (n > 1) stats::qt(0.975, n - 1) * s / sqrt(n) else NA_real_
    c(N = n, mean = m, sd = s, median = if (n) stats::median(x) else NA_real_,
      min = if (n) min(x) else NA_real_, max = if (n) max(x) else NA_real_,
      ci95_low = m - h, ci95_high = m + h)
  }
  out <- do.call(rbind, lapply(features, function(f) do.call(rbind, lapply(lv, function(g) {
    x <- data[[f]][grp == g]
    st <- one(x)
    data.frame(Feature = f,
               Label = if (is.null(labels) || is.na(labels[f])) f else unname(labels[f]),
               Group = g, N = st[["N"]], N_missing = sum(is.na(x)),
               t(signif(st[-1], digits)), stringsAsFactors = FALSE)
  }))))
  rownames(out) <- NULL
  out
}

# -- export_measurements_descriptive_workbook() ----------------------------------
#' Export the measurement descriptive statistics to a formatted Excel sheet
#'
#' Same layout as the taxa summary: one row per measurement and the groups side
#' by side, each under a merged header, with N, median, min, max, range, mean,
#' standard deviation and the 95\% CI of the mean.
#'
#' @inheritParams summarize_measurements
#' @param outfile Path of the .xlsx file to write.
#' @param sheet_name Worksheet name (max 31 characters).
#' @param font_name,font_size Base font for the sheet.
#' @param zoom Sheet zoom in percent.
#' @param ... Passed to \code{summarize_measurements()}.
#' @return Invisibly, the long-format summary (as \code{summarize_measurements()}).
#' @export
export_measurements_descriptive_workbook <- function(data, features, group_columns,
                                                     outfile    = "measurement_descriptives.xlsx",
                                                     sheet_name = "Descriptives_by_group",
                                                     font_name  = "Aptos Narrow",
                                                     font_size  = 11,
                                                     zoom       = 70,
                                                     ...) {
  sm <- summarize_measurements(data, features, group_columns, ...)
  groups <- unique(sm$Group)
  feats  <- unique(sm$Feature)
  num    <- function(x) ifelse(is.na(x), "", format(x, trim = TRUE, scientific = FALSE, drop0trailing = TRUE))
  # " to " rather than "-": measurements can be negative
  sm$range <- ifelse(is.na(sm$min), "", paste(vapply(sm$min, num, ""), "to", vapply(sm$max, num, "")))
  stat_hd <- c("N", "median", "min", "max", "range", "mean", "sd", "ci95_low", "ci95_high")
  k <- length(stat_hd); n_g <- length(groups); n_cols <- 1 + k * n_g

  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontSize = font_size, fontName = font_name)
  openxlsx::addWorksheet(wb, sheet_name, zoom = zoom)
  st_center <- openxlsx::createStyle(halign = "center")
  st_top    <- openxlsx::createStyle(border = "top",    borderStyle = "thin")
  st_bottom <- openxlsx::createStyle(border = "bottom", borderStyle = "thin")
  openxlsx::setColWidths(wb, sheet_name, cols = 1, widths = 42)
  openxlsx::setColWidths(wb, sheet_name, cols = 2:n_cols,
                         widths = rep(c(7, 12, 12, 12, 20, 12, 12, 12, 12), times = n_g))

  first_row <- 4; last_row <- 3 + length(feats)
  openxlsx::writeData(wb, sheet_name, paste("Grouped by:", paste(group_columns, collapse = " + ")),
                      startCol = 1, startRow = 1, colNames = FALSE)
  openxlsx::writeData(wb, sheet_name, t(c("Measurement", rep(stat_hd, times = n_g))),
                      startCol = 1, startRow = 3, colNames = FALSE)
  openxlsx::writeData(wb, sheet_name, sm$Label[match(feats, sm$Feature)], startCol = 1,
                      startRow = first_row, colNames = FALSE)
  for (gi in seq_len(n_g)) {
    c0 <- 2 + (gi - 1) * k
    openxlsx::writeData(wb, sheet_name, groups[gi], startCol = c0, startRow = 2, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet_name, cols = c0:(c0 + k - 1), rows = 2)
    dg <- sm[sm$Group == groups[gi], , drop = FALSE]
    dg <- dg[match(feats, dg$Feature), , drop = FALSE]
    openxlsx::writeData(wb, sheet_name, dg[, c("N", "median", "min", "max")], startCol = c0,
                        startRow = first_row, colNames = FALSE)
    openxlsx::writeData(wb, sheet_name, dg$range, startCol = c0 + 4, startRow = first_row,
                        colNames = FALSE)
    openxlsx::writeData(wb, sheet_name, dg[, c("mean", "sd", "ci95_low", "ci95_high")],
                        startCol = c0 + 5, startRow = first_row, colNames = FALSE)
  }
  openxlsx::addStyle(wb, sheet_name, st_center, rows = 2:last_row, cols = 2:n_cols,
                     gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, sheet_name, st_top, rows = 2, cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, sheet_name, st_bottom, rows = c(3, last_row), cols = 1:n_cols,
                     gridExpand = TRUE, stack = TRUE)
  openxlsx::freezePane(wb, sheet_name, firstActiveRow = first_row, firstActiveCol = 2)
  openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)
  message(sprintf("Workbook written to: %s", outfile))
  invisible(sm[, setdiff(names(sm), "range")])
}
