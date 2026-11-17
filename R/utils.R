# utils.R
# -----------------------------------------------------------------------------
# Shared internal helpers. Every other file calls these instead of keeping its
# own copy, so a fix in one place applies everywhere.
# -----------------------------------------------------------------------------

# -- Baseline flag ------------------------------------------------------------
# Parses TRUE/FALSE, Yes/No, Y/N, 1/0, T/F (any case) into a logical vector.
# Anything unrecognised becomes NA.
.as_flag <- function(x) {
  if (is.logical(x)) return(x)
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(s))
  out[s %in% c("true", "t", "yes", "y", "1")]  <- TRUE
  out[s %in% c("false", "f", "no", "n", "0")]  <- FALSE
  out
}

# NA-safe "is this row a baseline sample?" -> never returns NA.
# Use  meta[.is_baseline(meta, col), ]  and  meta[!.is_baseline(meta, col), ]
# everywhere instead of  meta[meta[[col]] == TRUE, ]  /  == FALSE.
.is_baseline <- function(meta, col) {
  if (is.null(col) || !col %in% colnames(meta))
    return(rep(FALSE, nrow(meta)))
  .as_flag(meta[[col]]) %in% TRUE
}

# -- Columns that look like a subject / animal ID with repeated samples ----------
# Only columns whose NAME looks like an ID are considered, and only when their
# values actually repeat. A value-pattern search alone is useless here: in real
# metadata, measurement columns repeat by chance, two animals can share an
# owner, and several different animals can share a pet name. This stays silent
# unless the file really carries an ID-style column with repeats; it can never
# replace naming the column in Step 2.
.ID_NAME_PATTERN <- paste0("(^|[^a-z])(id|ids)([^a-z]|$)|subject|animal|patient|",
                           "catid|dogid|microchip|chip|pairedid|cageid|litter")

.subject_candidates <- function(meta) {
  nm  <- colnames(meta)
  nml <- tolower(gsub("[^A-Za-z]", "", nm))
  looks_id <- grepl(.ID_NAME_PATTERN, tolower(nm)) | grepl("^(cat|dog|pet|animal)id$", nml)
  out <- vapply(nm, function(v) {
    x <- as.character(meta[[v]]); x <- x[!is.na(x) & nzchar(x)]
    length(x) > 1 && anyDuplicated(x) > 0 && length(unique(x)) > 1
  }, logical(1))
  nm[looks_id & out]
}

# -- Pick a default for a selectInput -----------------------------------------
# Replaces the repeated  if (x %in% choices) x else choices[1]  pattern.
.pick <- function(preferred, choices) {
  hit <- intersect(preferred, choices)
  if (length(hit)) hit[1] else choices[1]
}

# -- Robust OTU / taxa table reader --------------------------------------------
# Returns a numeric matrix, taxa as rows, samples as columns. Handles:
#   * QIIME2 / biom exports ("# Constructed from biom file", "#OTU ID" header)
#   * MetaPhlAn merged tables ("#mpa_v..." line, "_metaphlan" suffixes)
#   * non-numeric annotation columns (e.g. "taxonomy", "NCBI_tax_id")
#   * blank / NA cells (set to 0 with a warning)
.read_otu_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  # drop leading comment lines, but keep a "#OTU ID"-style header
  while (length(lines) > 1 && startsWith(lines[1], "#") &&
         startsWith(lines[2], "#")) lines <- lines[-1]
  if (startsWith(lines[1], "#")) {
    if (grepl("\t", lines[1])) lines[1] <- sub("^#", "", lines[1])
    else lines <- lines[-1]
  }
  tab <- utils::read.table(text = lines, sep = "\t", header = TRUE,
                           row.names = 1, check.names = FALSE, quote = "",
                           comment.char = "", fill = TRUE)
  colnames(tab) <- sub("_metaphlan$", "", colnames(tab))

  num_ok <- vapply(tab, function(v)
    is.numeric(v) || all(is.na(suppressWarnings(as.numeric(v))) == is.na(v)),
    logical(1))
  meta_like <- tolower(colnames(tab)) %in% c("ncbi_tax_id", "taxonomy",
                                             "clade_taxid")
  drop <- !num_ok | meta_like
  if (any(drop))
    message("Ignoring non-sample column(s): ",
            paste(colnames(tab)[drop], collapse = ", "))
  mat <- as.matrix(tab[, !drop, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (anyNA(mat)) {
    warning("Empty or non-numeric cells in the OTU table were set to 0.",
            call. = FALSE)
    mat[is.na(mat)] <- 0
  }
  mat
}

# -- Match table sample names to metadata Sample IDs -------------------------------
# Pipelines often add a suffix to sample names ("GI220_kneaddata_combined_Abundance",
# "GI091_metaphlan"). A table name is matched to a metadata ID when it is that
# ID exactly, or starts with it followed by a separator (_ . - or space); the
# longest such ID wins. Returns, for each table name, the ID or NA.
.match_ids <- function(table_names, ids) {
  ids <- unique(as.character(ids[!is.na(ids)]))
  out   <- ids[match(table_names, ids)]
  exact <- !is.na(out)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  for (i in which(!exact)) {
    # byte-wise, so accented text elsewhere in the metadata cannot break it
    cand <- ids[startsWith(table_names[i], ids)]
    if (length(cand)) {
      raw  <- charToRaw(table_names[i])
      nxt  <- vapply(nchar(cand, type = "bytes"), function(k)
        if (k < length(raw)) rawToChar(raw[k + 1]) else "", character(1))
      cand <- cand[nxt %in% c("_", ".", "-", " ")]
    }
    if (length(cand)) out[i] <- cand[which.max(nchar(cand, type = "bytes"))]
  }
  # a prefix match is dropped when that ID is claimed more than once
  # (e.g. "S1_R1" and "S1_R2", or an exact "S1" is also present)
  out[!exact & out %in% out[duplicated(out) & !is.na(out)]] <- NA
  out
}

# Rename the table's sample columns to the metadata IDs they match.
.harmonise_samples <- function(mat, meta_ids) {
  hit <- .match_ids(colnames(mat), meta_ids)
  ren <- !is.na(hit) & hit != colnames(mat)
  if (any(ren)) {
    message(sprintf("%d sample name(s) shortened to match the metadata, e.g. '%s' -> '%s'.",
                    sum(ren), colnames(mat)[which(ren)[1]], hit[which(ren)[1]]))
    colnames(mat)[ren] <- hit[ren]
  }
  mat
}

# -- Locale-proof handling of non-ASCII text ---------------------------------------
# R's table parser rewrites characters it cannot show in the current locale
# (e.g. the micro sign becomes the literal text "<ce><bc>" on a "C" locale), so
# results would differ between computers. Non-ASCII characters are therefore
# swapped for ASCII tokens ("@@U03BC@@") before parsing and restored afterwards.
.ascii_encode <- function(x) {
  x   <- as.character(x)
  bad <- !is.na(x) & grepl("[^\001-\177]", x, useBytes = TRUE)
  x[bad] <- vapply(x[bad], function(s) {
    if (!validUTF8(s)) s <- iconv(s, "latin1", "UTF-8", sub = "?")
    cp <- utf8ToInt(s)
    paste(ifelse(cp > 127, sprintf("@@U%04X@@", cp), intToUtf8(cp, multiple = TRUE)),
          collapse = "")
  }, character(1), USE.NAMES = FALSE)
  x
}
.ascii_decode <- function(x) {
  hit <- !is.na(x) & grepl("@@U[0-9A-F]{4,6}@@", x)
  x[hit] <- vapply(x[hit], function(s) {
    m  <- gregexpr("@@U[0-9A-F]{4,6}@@", s)
    tk <- regmatches(s, m)[[1]]
    regmatches(s, m) <- list(vapply(tk, function(t)
      intToUtf8(strtoi(gsub("@@U|@@", "", t), 16L)), character(1)))
    s
  }, character(1), USE.NAMES = FALSE)
  Encoding(x[hit]) <- "UTF-8"
  x
}

# -- Tidy column names ------------------------------------------------------------
# "BW (kg)" -> "BW.kg", "Age at presentation (years)" -> "Age.at.presentation.years",
# a second "Sex" column -> "Sex.1". Greek letters are spelled out (micro -> "u",
# alpha -> "alpha", beta -> "beta"), other non-ASCII characters are dropped, so
# the names are valid in formulas and identical on every computer.
.clean_names <- function(x) {
  x <- .ascii_encode(gsub("[\r\n]+", " ", x))
  greek <- c("03B1" = "alpha", "03B2" = "beta", "03B3" = "gamma", "03B4" = "delta",
             "03BC" = "u", "00B5" = "u", "0394" = "Delta", "03C9" = "omega")
  for (g in names(greek)) x <- gsub(sprintf("@@U%s@@", g), greek[[g]], x, fixed = TRUE)
  x <- gsub("@@U[0-9A-F]{4,6}@@", "", x)
  x <- make.names(trimws(x), unique = FALSE)
  x <- sub("\\.+$", "", gsub("\\.{2,}", ".", x))
  x[!nzchar(x)] <- "X"
  make.unique(x, sep = ".")
}

# -- Metadata reader (path or data frame) -> data frame, Sample IDs as rownames
# * sample_ids (the abundance table's sample names): the metadata column that
#   matches them best is used as the Sample ID column (normally the first one).
# * column names are tidied (.clean_names); duplicated names are reported.
# * "28.6%"-style text columns become numbers.
# * the baseline flag accepts TRUE/FALSE, Yes/No, Y/N, 1/0.
.read_meta <- function(meta, is_baseline_col = NULL, sample_ids = NULL) {
  if (is.character(meta) && length(meta) == 1) {
    # Read as text first: copes with Excel's UTF-8 "BOM" mark, non-UTF-8
    # characters (e.g. a micro sign) and tab- as well as comma-separated files.
    lines <- readLines(meta, warn = FALSE)
    lines[1] <- sub("^\xef\xbb\xbf", "", lines[1], useBytes = TRUE)
    lines <- .ascii_encode(lines)               # shield non-ASCII text (see above)
    n_tab <- lengths(regmatches(lines[1], gregexpr("\t", lines[1])))
    n_com <- lengths(regmatches(lines[1], gregexpr(",", lines[1], fixed = TRUE)))
    raw <- utils::read.table(text = lines, sep = if (n_tab > n_com) "\t" else ",",
                             header = TRUE, quote = "\"", fill = TRUE,
                             comment.char = "", check.names = FALSE,
                             stringsAsFactors = FALSE, na.strings = c("", "NA"),
                             strip.white = TRUE)
    orig <- trimws(gsub("[\r\n]+", " ", colnames(raw)))
    dup  <- unique(orig[duplicated(orig) & nzchar(orig)])
    colnames(raw) <- .clean_names(colnames(raw))
    labels <- stats::setNames(.ascii_decode(gsub(" +", " ", orig)), colnames(raw))  # as typed
    dup    <- .ascii_decode(dup)
    chr    <- vapply(raw, is.character, logical(1))
    raw[chr] <- lapply(raw[chr], .ascii_decode)       # restore accents in the cells
    raw <- raw[rowSums(!is.na(raw)) > 0, colSums(!is.na(raw)) > 0 | seq_along(raw) == 1,
               drop = FALSE]
    if (length(dup))
      message("Metadata has repeated column name(s): ", paste(dup, collapse = ", "),
              ". The 2nd copy is renamed with '.1' (e.g. ",
              .clean_names(c(dup[1], dup[1]))[2], ") - check you pick the right one.")

    id_col <- 1L
    if (!is.null(sample_ids)) {
      hits <- vapply(raw, function(x)
        sum(!is.na(.match_ids(sample_ids, trimws(as.character(x))))), numeric(1))
      if (max(hits) > 0 && hits[1] < max(hits)) {
        id_col <- which.max(hits)
        message(sprintf("Sample IDs taken from metadata column '%s' (best match to the table).",
                        colnames(raw)[id_col]))
      }
    }
    ids <- trimws(as.character(raw[[id_col]]))
    if (anyNA(ids) || any(!nzchar(ids)))
      stop("The Sample ID column of the metadata has empty cells.", call. = FALSE)
    if (anyDuplicated(ids))
      stop("Sample IDs in the metadata are not unique: ",
           paste(unique(ids[duplicated(ids)]), collapse = ", "),
           ". Each row must be one sample.", call. = FALSE)
    meta <- raw[, -id_col, drop = FALSE]
    rownames(meta) <- ids
    attr(meta, "labels") <- labels[colnames(meta)]
  }
  lab <- attr(meta, "labels")
  meta[] <- lapply(meta, function(x) {
    if (!is.character(x)) return(x)
    x <- trimws(x)
    v <- x[!is.na(x) & nzchar(x)]
    if (length(v) && all(grepl("^-?[0-9.]+ ?%$", v)))          # "28.60%" -> 28.6
      return(suppressWarnings(as.numeric(sub(" ?%$", "", x))))
    x
  })
  if (!is.null(is_baseline_col) && is_baseline_col %in% colnames(meta)) {
    flag <- .as_flag(meta[[is_baseline_col]])
    raw_b <- trimws(as.character(meta[[is_baseline_col]]))
    bad  <- is.na(flag) & !is.na(raw_b) & nzchar(raw_b)
    if (any(bad))
      warning(sprintf("'%s' has unrecognised value(s) treated as NA: %s",
                      is_baseline_col,
                      paste(unique(meta[[is_baseline_col]][bad]),
                            collapse = ", ")), call. = FALSE)
    meta[[is_baseline_col]] <- flag
  }
  attr(meta, "labels") <- lab
  meta
}

# Sample names of an abundance table without reading the whole file
.peek_sample_ids <- function(path) {
  lines <- readLines(path, n = 5, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  hdr <- lines[grepl("\t", lines)][1]
  if (is.na(hdr)) return(character(0))
  sub("_metaphlan$", "", trimws(strsplit(hdr, "\t", fixed = TRUE)[[1]][-1]))
}

# -- Counts vs. relative abundance ---------------------------------------------
# Everything downstream works in percent. Rule per sample (x is samples x taxa):
#   * all values whole numbers AND sum not ~100 -> raw counts  -> divide by sum x 100
#   * sum ~1 (proportions, e.g. HUMAnN "relab")  -> multiply by 100
#   * anything else (already percent; MetaPhlAn with an unclassified share,
#     so < 100)                                 -> left as uploaded
.to_relative_abundance <- function(x, tol = 1) {
  rs        <- rowSums(x)
  whole     <- apply(x, 1, function(r) all(abs(r - round(r)) < 1e-8))
  is_counts <- whole & abs(rs - 100) > tol & rs > 0
  is_prop   <- !is_counts & abs(rs - 1) <= tol / 100
  if (any(rs == 0))
    warning("Sample(s) with total abundance 0 left as 0: ",
            paste(rownames(x)[rs == 0], collapse = ", "), call. = FALSE)
  out <- x
  out[is_counts, ] <- x[is_counts, , drop = FALSE] / rs[is_counts] * 100
  out[is_prop, ]   <- x[is_prop, , drop = FALSE] * 100
  attr(out, "scale_info") <- list(
    n_total           = nrow(x),
    n_already_percent = sum(!is_counts & !is_prop),
    n_converted       = sum(is_counts),
    n_proportion      = sum(is_prop),
    converted_samples = rownames(x)[is_counts]
  )
  out
}

# -- Lineage parsing (moved here from descriptive_summary.R) -------------------
# One row per clade, one column per rank (Domain..Species); cells hold the
# original tag ("p__Firmicutes") or NA. Accepts ";" or "|" and k__ for d__.
.lineage_matrix <- function(clades) {
  rank_letters <- unname(TAXA_PREFIXES)
  parts <- strsplit(as.character(clades), "[;|]")
  out <- t(vapply(parts, function(p) {
    p   <- trimws(p)
    rk  <- tolower(substr(p, 1, 1)); rk[rk == "k"] <- "d"
    val <- sub("^[A-Za-z]__", "", p)
    ok  <- grepl("^[A-Za-z]__", p) & nzchar(val) & val != "_" &
           rk %in% rank_letters
    row <- rep(NA_character_, length(rank_letters))
    row[match(rk[ok], rank_letters)] <- p[ok]
    row
  }, character(length(rank_letters))))
  dimnames(out) <- list(NULL, names(TAXA_PREFIXES))
  out
}

.lineage_depth <- function(lin)
  apply(lin, 1, function(r) { w <- which(!is.na(r)); if (length(w)) max(w) else 0L })

# "terminal" = one row per feature (QIIME2 L7, ASV table);
# "all_ranks" = a row for every rank (un-filtered MetaPhlAn) -> summing would
# double count, so each level must be read from rows at exactly that depth.
.detect_layout <- function(lin, sep = "|") {
  keys   <- apply(lin, 1, function(r) paste(r[!is.na(r)], collapse = sep))
  parent <- apply(lin, 1, function(r) {
    r <- r[!is.na(r)]; if (length(r) >= 2) paste(r[-length(r)], collapse = sep)
    else NA_character_ })
  hp <- parent[!is.na(parent)] %in% keys
  if (length(hp) && mean(hp) >= 0.9) "all_ranks" else "terminal"
}

# Does the table carry rank tags (p__, g__, ...) at all?  A species-only table
# such as a MetaPhlAn export reduced to "Collinsella_intestinalis" does not.
.has_lineage <- function(clades) any(!is.na(.lineage_matrix(clades)))

# Taxonomic levels a table can be summarised at.
.available_levels <- function(clades) {
  lin <- .lineage_matrix(clades[!grepl("(^|[;|])t__", clades)])
  if (all(is.na(lin))) return("Species")
  names(TAXA_PREFIXES)[seq_len(max(.lineage_depth(lin)))]
}

# -- THE single taxonomy collapser ----------------------------------------------
# mat: taxa x samples numeric matrix with clade strings as rownames.
# Returns a taxa x samples matrix at `level`; rownames = full lineage up to
# that rank (so equal names in different families stay separate).
#  * Rows not resolved down to `level` keep their deepest known rank and are
#    labelled "<that rank>_unclassified" (e.g. g__Streptococcus;__ at Species
#    level -> "Streptococcus_unclassified"), as MetaPhlAn does.
#  * Rows with no rank tag at all are pooled as "Unclassified".
#  * Tables without any lineage (names only) are used as they are, and only
#    at the "Species" (as-uploaded) level.
# Used by: taxonomy bar plot, MaAsLin3 input.
.collapse_to_level <- function(mat, level, sep = "|") {
  strain <- grepl("(^|[;|])t__", rownames(mat))     # MetaPhlAn SGB rows
  mat    <- mat[!strain, , drop = FALSE]
  lin    <- .lineage_matrix(rownames(mat))
  ti     <- match(level, names(TAXA_PREFIXES))

  if (all(is.na(lin))) {
    if (level != "Species")
      stop("This table has taxon names only (no lineage such as ",
           "'k__Bacteria|p__...|g__...'), so it can only be shown as uploaded ",
           "(choose 'Species'). Upload the full MetaPhlAn / QIIME 2 table to ",
           "summarise at ", level, " level.", call. = FALSE)
    attr(mat, "layout") <- "no_lineage"
    return(mat)
  }

  layout <- .detect_layout(lin, sep)
  depth  <- .lineage_depth(lin)
  sel    <- if (layout == "all_ranks") depth == ti else rep(TRUE, nrow(mat))
  letter <- unname(TAXA_PREFIXES)[ti]
  keys   <- apply(lin[sel, seq_len(ti), drop = FALSE], 1, function(r) {
    known <- r[!is.na(r)]
    if (!length(known)) return("Unclassified")
    if (!is.na(r[ti])) return(paste(known, collapse = sep))
    last <- sub("^[A-Za-z]__", "", known[length(known)])
    paste(c(known, paste0(letter, "__", last, "_unclassified")), collapse = sep)
  })
  out <- rowsum(mat[sel, , drop = FALSE], group = keys)
  attr(out, "layout") <- layout
  out
}

# Rows to use for alpha / beta diversity: on an all-ranks table only the
# deepest rows, otherwise higher-rank rows inflate richness and distances.
.feature_rows <- function(mat) {
  strain <- grepl("(^|[;|])t__", rownames(mat))
  lin    <- .lineage_matrix(rownames(mat)[!strain])
  if (all(is.na(lin)) || .detect_layout(lin) == "terminal") return(mat)
  d    <- .lineage_depth(lin)
  keep <- !strain
  keep[!strain] <- d == max(d)
  message(sprintf("All-ranks table detected: using the %d deepest-rank rows ",
                  sum(keep)), "for alpha/beta diversity.")
  mat[keep, , drop = FALSE]
}

# Short plot label for a lineage key: the last rank's name, e.g. "Blautia".
# If two different taxa end with the same name (e.g. "uncultured" in two
# families), the parent rank is prefixed so they stay distinguishable.
.short_taxon <- function(keys) {
  toks <- lapply(strsplit(keys, "[;|]"), function(t) sub("^[A-Za-z]__", "", t))
  last <- vapply(toks, function(t) t[length(t)], character(1))
  dup  <- last %in% last[duplicated(last)]
  last[dup] <- vapply(toks[dup], function(t)
    paste(utils::tail(t, 2), collapse = " / "), character(1))
  last
}

# -- TIFF writer used by every export_*_tiff() ----------------------------------
.save_tiff <- function(p, outfile, width, height, dpi) {
  ggplot2::ggsave(outfile, p, device = "tiff", width = width, height = height,
                  units = "in", dpi = dpi, compression = "lzw")
  invisible(outfile)
}

# -- Safe download file-name part ------------------------------------------------
.fname <- function(...) paste(gsub("[^A-Za-z0-9]", "_", c(...)), collapse = "_")
