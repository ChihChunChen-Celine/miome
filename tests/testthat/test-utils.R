# Regression tests for the shared helpers (bugs A3, A4, A5 in the review)

test_that("baseline flags parse beginner-style values and are never NA", {
  expect_identical(.as_flag(c("Yes", "no", "1", "0", "TRUE", "f", NA, "x")),
                   c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, NA, NA))
  m <- data.frame(b = c(TRUE, NA, FALSE))
  expect_identical(.is_baseline(m, "b"), c(TRUE, FALSE, FALSE))
  expect_identical(.is_baseline(m, "absent"), c(FALSE, FALSE, FALSE))
  expect_identical(.is_baseline(m, ""), c(FALSE, FALSE, FALSE))
})

test_that("reader handles QIIME 2 '#OTU ID' header and annotation columns", {
  f <- tempfile()
  writeLines(c("# Constructed from biom file", "#OTU ID\tS1\tS2\ttaxonomy",
               "t1\t5\t3\td__B", "t2\t1\t\td__B"), f)
  x <- suppressWarnings(suppressMessages(.read_otu_table(f)))
  expect_identical(colnames(x), c("S1", "S2"))
  expect_identical(rownames(x), c("t1", "t2"))
  expect_equal(x["t2", "S2"], 0)
})

test_that("reader handles MetaPhlAn comment line, NCBI column and suffix", {
  f <- tempfile()
  writeLines(c("#mpa_vJan21", "clade_name\tNCBI_tax_id\tA_metaphlan\tB_metaphlan",
               "k__Bacteria\t2\t100\t100", "k__Bacteria|p__F\t2|1\t60\t70"), f)
  y <- suppressMessages(.read_otu_table(f))
  expect_identical(colnames(y), c("A", "B"))
  expect_equal(nrow(y), 2)
})

test_that("low-depth count tables are converted; relative tables are not", {
  cnt <- matrix(c(10, 20, 5, 3, 0, 40), 2, byrow = TRUE, dimnames = list(c("a", "b"), NULL))
  expect_equal(unname(rowSums(.to_relative_abundance(cnt))), c(100, 100))
  pct <- matrix(c(50.5, 49.5, 0, 10.1, 80, 9.9), 2, byrow = TRUE,
                dimnames = list(c("a", "b"), NULL))
  expect_equal(attr(.to_relative_abundance(pct), "scale_info")$n_converted, 0)
  expect_equal(attr(.to_relative_abundance(pct / 100), "scale_info")$n_converted, 0)
  expect_equal(attr(.to_relative_abundance(pct / 100), "scale_info")$n_proportion, 2)
  stacked <- pct * 3
  expect_equal(attr(.to_relative_abundance(stacked), "scale_info")$n_converted, 0)
})

test_that("collapse works for terminal QIIME tables and keeps totals", {
  q <- matrix(c(5, 5, 3, 1, 2, 4, 1, 0, 7, 7), 5, byrow = TRUE, dimnames = list(
    c("d__B;p__F;c__C;o__O;f__L;g__X;s__x1", "d__B;p__F;c__C;o__O;f__L;g__Y;s__uncultured",
      "d__B;p__F;c__C;o__O;f__M;g__Z;s__uncultured", "d__B;p__F;c__C;o__O;__;__;__",
      "Unassigned;__;__;__;__;__;__"), c("S1", "S2")))
  fam <- .collapse_to_level(q, "Family")
  expect_equal(attr(fam, "layout"), "terminal")
  # resolved only to order -> kept under its deepest known rank, not dropped
  expect_equal(fam["d__B|p__F|c__C|o__O|f__O_unclassified", "S1"], 1)
  expect_equal(fam["Unclassified", "S1"], 7)                 # no taxonomy at all
  expect_equal(sum(fam[, "S1"]), sum(q[, "S1"]))
  sp <- .collapse_to_level(q, "Species")
  expect_equal(sum(grepl("uncultured", rownames(sp))), 2)
  expect_equal(colSums(sp), colSums(q))
})

test_that("tables with taxon names only (no lineage) are used as uploaded", {
  m <- matrix(1:6, 3, dimnames = list(c("Blautia_wexlerae", "GGB3602_SGB66069", "Clostridium_sp_AT4"),
                                      c("A", "B")))
  expect_false(.has_lineage(rownames(m)))
  expect_identical(.available_levels(rownames(m)), "Species")
  expect_identical(unclass(.collapse_to_level(m, "Species"))[, ], m)
  expect_error(.collapse_to_level(m, "Family"), "names only")
  expect_identical(.available_levels("d__B;p__F;c__C;o__O;f__L;g__X;s__x"), names(TAXA_PREFIXES))
  expect_identical(suppressMessages(.feature_rows(m)), m)
})

test_that("column names are tidied and stay unique", {
  expect_identical(.clean_names(c("BW (kg)", "Age at presentation\n(years)", "Sex", "Sex",
                                  "Dry Matter %", "E. Coli Log DNA", "")),
                   c("BW.kg", "Age.at.presentation.years", "Sex", "Sex.1", "Dry.Matter",
                     "E.Coli.Log.DNA", "X"))
})

test_that("metadata reader copes with Excel-style files", {
  f   <- system.file("testdata", "metadata_messy.csv", package = "miome")
  ids <- .peek_sample_ids(system.file("testdata", "metaphlan_species_only.txt", package = "miome"))
  expect_identical(ids, sprintf("D%02d", 1:12))
  expect_message(m <- .read_meta(f, sample_ids = ids), "repeated column name")
  expect_identical(rownames(m), sprintf("D%02d", 1:10))    # ID column found (2nd column, BOM, blank row)
  expect_identical(colnames(m), c("Lab.NO", "Group", "Sex", "Sex.1", "BW.kg",
                                  "Age.at.presentation.years", "Dry.Matter", "Note"))
  expect_type(m$Dry.Matter, "double")                      # "21.09%" -> 21.09
  expect_identical(unique(m$Sex.1), c("FS", "MC"))         # trailing blank trimmed
  expect_identical(m$Note[1], "ok, fine")
  # tab-separated metadata and duplicated IDs
  tsv <- tempfile(fileext = ".txt")
  writeLines(c("SampleID\tGroup", "a\tx", "b\ty"), tsv)
  expect_identical(.read_meta(tsv)$Group, c("x", "y"))
  writeLines(c("SampleID,Group", "a,x", "a,y"), tsv)
  expect_error(.read_meta(tsv), "not unique")
})

test_that("all-ranks MetaPhlAn tables are not double counted (bug A1)", {
  a <- matrix(c(100, 100, 100, 100, 60, 30, 40, 70, 60, 30, 40, 70), 6, byrow = TRUE,
              dimnames = list(c("k__B", "k__B|p__F", "k__B|p__F|c__C1", "k__B|p__F|c__C2",
                                "k__B|p__F|c__C1|o__O1", "k__B|p__F|c__C2|o__O2"), c("A", "B")))
  cl <- .collapse_to_level(a, "Class")
  expect_equal(attr(cl, "layout"), "all_ranks")
  expect_equal(unname(colSums(cl)), c(100, 100))
  ph <- .collapse_to_level(a, "Phylum")
  expect_equal(nrow(ph), 1)
  expect_equal(ph[1, 1], 100)
  fr <- suppressMessages(.feature_rows(a))
  expect_true(all(grepl("o__", rownames(fr))))
})

test_that("short taxon labels stay distinct for duplicated names", {
  expect_identical(.short_taxon(c("k__B|p__F|c__C1", "d__B;p__G")), c("C1", "G"))
  lab <- .short_taxon(c("f__A;g__uncultured", "f__B;g__uncultured"))
  expect_equal(length(unique(lab)), 2)
})

test_that("table sample names with pipeline suffixes are matched to metadata IDs", {
  ids <- c("GI22", "GI220", "GI221", "S1")
  expect_identical(.match_ids(c("GI220_kneaddata_combined_Abundance", "GI221", "GI22.profile",
                                "GI2200_x", "other"), ids),
                   c("GI220", "GI221", "GI22", NA, NA))
  # two columns claiming one ID are ambiguous; an exact match always wins
  expect_identical(.match_ids(c("S1_R1", "S1_R2"), ids), c(NA_character_, NA_character_))
  expect_identical(.match_ids(c("S1", "S1_rep"), ids), c("S1", NA))
  m <- matrix(1:4, 2, dimnames = list(c("a", "b"), c("GI220_Abundance", "zzz")))
  expect_identical(colnames(suppressMessages(.harmonise_samples(m, ids))), c("GI220", "zzz"))
})

test_that("proportions (rows summing to 1) are shown as percent", {
  p <- matrix(c(0.2, 0.3, 0.5, 0.1, 0.1, 0.8), 2, byrow = TRUE, dimnames = list(c("a", "b"), NULL))
  r <- .to_relative_abundance(p)
  expect_equal(unname(rowSums(r)), c(100, 100))
  expect_equal(attr(r, "scale_info")$n_proportion, 2)
  expect_equal(attr(r, "scale_info")$n_converted, 0)
  pct93 <- matrix(c(50, 43, 0), 1, dimnames = list("a", NULL)) + 0.001     # MetaPhlAn with unclassified
  expect_equal(sum(.to_relative_abundance(pct93)), sum(pct93))             # left as uploaded
})

test_that("subject-ID columns are detected by name plus repeats, without false alarms", {
  me <- data.frame(CatID = c("A", "A", "B", "B"), Microchip = c("1", "1", "2", "2"),
                   PetName = c("Mia", "Mia", "Mia", "Rex"),         # 3 different animals, same name
                   Owner = c("Smith", "Smith", "Jones", "Lee"),     # 2 dogs, 1 owner
                   Shannon = c(3.1, 3.1, 2.0, 2.5),                 # repeats by chance
                   Treatment = c("CON", "TRT", "CON", "TRT"),
                   stringsAsFactors = FALSE)
  expect_identical(.subject_candidates(me), c("CatID", "Microchip"))
  # one row per animal -> nothing flagged even though the ID column exists
  expect_identical(.subject_candidates(data.frame(SubjectID = c("a", "b", "c"))), character(0))
  expect_identical(.subject_candidates(data.frame(Treatment = c("x", "x", "y"))), character(0))
})
