# End-to-end regression tests on the bundled synthetic datasets

td   <- function(f) system.file("testdata", f, package = "miome")
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("prepare_all builds baseline covariates for EVERY metric (bug A2)", {
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv")))
  expect_true(all(paste0("Baseline_", ALPHA_METRICS) %in% colnames(obj$meta)))
  expect_true(all(is.na(obj$meta$Baseline_Shannon[.is_baseline(obj$meta, "IsBaseline")])))
  expect_true(file.exists(obj$meta_path))                # bug B4
  expect_equal(nrow(obj$meta), 36)
})

test_that("QIIME and MetaPhlAn inputs give the same family composition (bug A1)", {
  q <- quiet(compute_taxonomy_barplot_data(td("qiime_counts.tsv"), td("metadata.csv"),
                                           "Family", top_n = 5))
  m <- quiet(compute_taxonomy_barplot_data(td("metaphlan_allranks.tsv"), td("metadata.csv"),
                                           "Family", top_n = 5))
  expect_setequal(unique(q$Taxon), unique(m$Taxon))
  expect_true(all(abs(tapply(m$RelAbundance, m$SampleID, sum) - 100) < 1e-6))
  expect_false("Unclassified" %in% m$Taxon)
})

test_that("works when the metadata has no baseline column at all", {
  me <- utils::read.csv(td("metadata.csv"))
  me$IsBaseline <- NULL
  f <- tempfile(fileext = ".csv"); utils::write.csv(me, f, row.names = FALSE)
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), f))
  expect_false(any(grepl("^Baseline_", colnames(obj$meta))))
  res <- quiet(run_global_test(obj, "Treatment", permutations = 99))
  expect_equal(res$N_samples, 36)
})

test_that("LMM for all metrics uses each metric's own baseline", {
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv")))
  res <- quiet(run_alpha_lmm_all_metrics(obj$meta, "Timepoint + Treatment", "sampleid", "CON"))
  used <- mapply(function(m, f) grepl(paste0("`Baseline_", m, "`"), f, fixed = TRUE),
                 res$Metric, res$Formula)
  expect_true(all(tapply(used, res$Metric, any)[setdiff(ALPHA_METRICS, "Chao1")]))
  expect_equal(sum(!is.na(res$p_adj_across_metrics)),
               length(unique(res$Metric[res$Term != "ERROR"])))
  inter <- quiet(run_alpha_lmm_interaction_all(obj$meta, "Treatment * Timepoint", "sampleid", "CON"))
  expect_false(any(grepl("Baseline_", inter$Formula)))
  expect_true(all(grepl(":", inter$Term[!is.na(inter$p_adj_across_metrics)])))
})

test_that("PERMANOVA global, typed formula and pairwise run", {
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv")))
  g <- quiet(run_global_test(obj, "Treatment", covariates = c("Timepoint", "Sequence"),
                             strata = "sampleid", permutations = 99))
  expect_equal(g$Variable, c("Timepoint", "Sequence", "Treatment"))   # covariates shown
  expect_equal(g$Role, c("covariate", "covariate", "primary variable"))
  expect_equal(unique(g$Strata), "sampleid")
  expect_message(run_global_test(obj, "Treatment", covariates = "Sequence", strata = "sampleid",
                                 permutations = 99),
                 "PERMANOVA: distance ~ Sequence \\+ Treatment \\| strata: sampleid")
  me <- obj; me$meta$uid <- rownames(me$meta)       # every block a single sample
  u <- quiet(run_global_test(me, "Treatment", strata = "uid", permutations = 99))
  expect_match(u$Strata, "NOT applied")
  f <- quiet(run_global_test(obj, "Treatment", formula_rhs = "Treatment * Timepoint",
                             permutations = 99))
  expect_true(any(grepl(":", f$Variable)))
  p <- quiet(run_pairwise_test(obj, "Treatment", "anosim", permutations = 99))
  expect_equal(nrow(p), 1)
})

test_that("plots build, paired lines follow time order (bug B1) and jitter is shared (B2)", {
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv")))
  p2 <- build_2d(obj, "Treatment", connect_by = "sampleid")
  path <- ggplot2::layer_data(p2, 1)
  first <- obj$pcoa_df[obj$pcoa_df$sampleid == "C01", ]
  first <- first[order(first$Timepoint), ]
  expect_equal(path$x[1:3], first$PC1)

  pa <- build_alpha_scatter(obj$meta, "Observed_taxa", "Treatment", "Timepoint",
                            connect_by = "sampleid", exclude_baseline = TRUE)
  lines <- ggplot2::layer_data(pa, 1); pts <- ggplot2::layer_data(pa, 3)
  expect_equal(lines$x, pts$x)                         # lines end on the points
  expect_equal(pts$y, lines$y)
  expect_true(all(pts$y == round(pts$y)))              # no vertical jitter

  tb <- quiet(compute_taxonomy_barplot_data(td("qiime_counts.tsv"), td("metadata.csv"),
                                            "Genus", top_n = 20))
  bp <- build_taxon_barplot(tb, "Genus", "Treatment", exclude_baseline = TRUE)
  fills <- unique(ggplot2::layer_data(bp, 1)$fill)
  expect_equal(length(fills), length(unique(tb$Taxon)))   # bug B3: no repeats
  expect_s3_class(build_3d(obj, "Treatment", connect_by = "sampleid"), "plotly")
})

test_that("TIFF exports write files", {
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv")))
  f <- tempfile(fileext = ".tiff")
  export_2d_tiff(obj, "Treatment", f)
  expect_gt(file.size(f), 1000)
  f3 <- tempfile(fileext = ".tiff")
  quiet(export_3d_tiff(obj, "Treatment", f3))
  expect_identical(readBin(f3, "raw", 2), as.raw(c(0x49, 0x49)))   # real TIFF, bug B9
})

test_that("descriptive workbook and MaAsLin input work on MetaPhlAn tables", {
  out <- tempfile(fileext = ".xlsx")
  res <- quiet(export_taxa_descriptive_workbook(td("metaphlan_allranks.tsv"), td("metadata.csv"),
                                                "Treatment", outfile = out,
                                                exclude_baseline = TRUE))
  expect_true(file.exists(out))
  expect_equal(res$table_type, "all_ranks")
  obj <- quiet(prepare_all(td("metaphlan_allranks.tsv"), td("metadata.csv")))
  mi <- quiet(prepare_maaslin_input(td("metaphlan_allranks.tsv"), obj$meta, "Genus"))
  expect_equal(nrow(mi$features), 24)                      # 36 - 12 baseline
  expect_equal(ncol(mi$features), 30)
})

test_that("export_abundance_workbook aligns ID labels with sample columns (bug A6)", {
  cnt <- tempfile(fileext = ".tsv")
  writeLines(c("# Constructed from biom file", "#OTU ID\tS2\tS1", "tA\t5\t1", "tB\t5\t3"), cnt)
  md <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(sampleid = c("S1", "S2"), Internal_ID = c("cat1", "cat2"),
                                Treatment = c("A", "B")), md, sep = "\t", row.names = FALSE,
                     quote = FALSE)
  wb <- tempfile(fileext = ".xlsx")
  res <- quiet(export_abundance_workbook(cnt, md, workbook_path = wb))
  sheet <- openxlsx::read.xlsx(wb, "Counts", colNames = FALSE)
  expect_identical(unlist(sheet[1, 2:3], use.names = FALSE), c("S2", "S1"))
  expect_identical(unlist(sheet[2, 2:3], use.names = FALSE), c("cat2", "cat1"))
  expect_equal(nrow(res$relative_abundance), 2)            # first taxon kept
})

test_that("LMM accepts several random effects", {
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv")))
  r1 <- quiet(run_alpha_lmm(obj$meta, "Shannon", "Timepoint + Treatment", "sampleid", "CON"))
  r2 <- quiet(run_alpha_lmm(obj$meta, "Shannon", "Timepoint + Treatment",
                            c("sampleid", "Room"), "CON"))
  expect_match(r2$Formula[1], "(1|`sampleid`) + (1|`Room`)", fixed = TRUE)
  expect_identical(r1$Term, r2$Term)                         # same fixed effects reported
  all <- quiet(run_alpha_lmm_all_metrics(obj$meta, "Timepoint + Treatment",
                                         c("sampleid", "Room"), "CON"))
  expect_false(any(all$Term == "ERROR"))
  expect_error(run_alpha_lmm(obj$meta, "Shannon", "Room + Treatment", c("sampleid", "Room")),
               "both the fixed effects and the random effects")
  expect_error(run_alpha_lmm(obj$meta, "Shannon", "Treatment", character(0)),
               "at least one random effect")
  me <- obj$meta; me$uid <- rownames(me)                      # one group per sample
  expect_error(quiet(run_alpha_lmm(me, "Shannon", "Treatment", "uid")),
               "fewer groups than samples")
})

test_that("MaAsLin3 runs with one or several random effects", {
  skip_if_not_installed("maaslin3")
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv")))
  mi  <- quiet(prepare_maaslin_input(td("qiime_counts.tsv"), obj$meta, "Family"))
  r1 <- quiet(run_maaslin(mi$features, mi$metadata, "Treatment + Timepoint",
                          random_effect = "sampleid", use_random = TRUE))
  expect_identical(attr(r1, "formula_used"), "~ Treatment + Timepoint + (1|sampleid)")
  # Lactobacillaceae was simulated as enriched under TRT
  top <- r1[r1$metadata == "Treatment" & grepl("Lactobacillaceae", r1$feature), ]
  expect_true(any(top$coef > 0 & top$qval < 0.05, na.rm = TRUE))

  # Two random effects: needs MaAsLin3 >= 1.5.5. Older versions (e.g. the
  # Bioconductor 1.4.0 release) get a clear error from miome before fitting.
  two <- function() quiet(run_maaslin(mi$features, mi$metadata, "Treatment + Timepoint",
                                      random_effect = c("sampleid", "Room"), use_random = TRUE))
  if (.maaslin_multi_ranef_ok()) {
    r2 <- two()
    expect_identical(attr(r2, "formula_used"),
                     "~ Treatment + Timepoint + (1|sampleid) + (1|Room)")
    expect_gt(nrow(r2), 0)
    expect_true(any(r2$metadata == "Treatment"))
  } else {
    expect_error(two(), "can use only one random effect")
  }
  # the message itself, whatever version is installed here
  local_mocked_bindings(.maaslin_multi_ranef_ok = function() FALSE)
  expect_error(two(), "can use only one random effect.*install_github")
  expect_s3_class(quiet(run_maaslin(mi$features, mi$metadata, "Treatment + Timepoint",
                                    random_effect = "sampleid", use_random = TRUE)), "data.frame")
  expect_error(run_maaslin(mi$features, mi$metadata, "Treatment", random_effect = character(0),
                           use_random = TRUE), "no random effect column")
  expect_error(run_maaslin(mi$features, mi$metadata, "Treatment + Room",
                           random_effect = c("sampleid", "Room"), use_random = TRUE),
               "already in the formula")
})


# ---- Real-world formats: names-only shotgun table, subset metadata, dropouts ----

test_that("species-only MetaPhlAn table + subset metadata run end to end", {
  otu <- td("metaphlan_species_only.txt"); md <- td("metadata_messy.csv")
  obj <- quiet(prepare_all(otu, md, pcoa_config(subject_id_col = "", is_baseline_col = "")))
  expect_equal(obj$match_info$n_matched, 10)
  expect_identical(obj$match_info$table_only, c("D11", "D12"))
  expect_identical(obj$tax_levels, "Species")
  expect_true(all(is.na(obj$meta$Chao1)))                       # relative abundance
  tb <- quiet(compute_taxonomy_barplot_data(otu, md, "Species", top_n = 5))
  expect_false("Unclassified" %in% tb$Taxon)
  expect_true("Blautia_wexlerae" %in% tb$Taxon)
  expect_error(quiet(compute_taxonomy_barplot_data(otu, md, "Family")), "names only")
  mi <- quiet(prepare_maaslin_input(otu, obj$meta, "Species", obj$config))
  expect_equal(ncol(mi$features), 19)        # 20 taxa, one present only in a left-out sample
  ds <- quiet(summarize_taxa_by_group(otu, obj$meta, "Group"))
  expect_identical(ds$tax_levels, "Species")
  expect_true("Blautia_wexlerae" %in% ds$summary$Taxon)
  g <- quiet(run_global_test(obj, "Group", covariates = c("Sex.1", "BW.kg"), permutations = 99))
  expect_identical(g$Variable, c("Sex.1", "BW.kg", "Group"))
})

test_that("QIIME 2 L7 table keeps taxa that are unresolved at the chosen level", {
  otu <- td("qiime_L7_rarefied.tsv"); md <- td("metadata_crossover_dropout.csv")
  cfg <- pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline")
  obj <- quiet(prepare_all(otu, md, cfg))
  expect_false(any(is.na(obj$meta$Chao1)))                      # counts
  expect_equal(sum(!is.na(obj$meta$Baseline_Shannon)), 15)      # all post-baseline rows
  for (lv in c("Genus", "Species")) {
    mi <- quiet(prepare_maaslin_input(otu, obj$meta, lv, cfg))
    # everything except the 'Unassigned' row is kept
    raw <- .read_otu_table(otu)[, rownames(mi$features)]
    expect_equal(unname(rowSums(mi$features)),
                 unname(100 - raw["Unassigned;__;__;__;__;__;__", ] / colSums(raw) * 100),
                 info = lv)
    expect_false(anyDuplicated(colnames(mi$features)) > 0, info = lv)
  }
  mi <- quiet(prepare_maaslin_input(otu, obj$meta, "Species", cfg))
  expect_true("Streptococcus_unclassified" %in% colnames(mi$features))
  tb <- quiet(compute_taxonomy_barplot_data(otu, md, "Genus", top_n = 20, config = cfg))
  expect_true("Lachnospiraceae_unclassified" %in% tb$Taxon)
  expect_true(all(abs(tapply(tb$RelAbundance, tb$SampleID, sum) - 100) < 1e-6))
})

test_that("restricted permutations survive a drop-out and flag untestable terms", {
  cfg <- pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline")
  obj <- quiet(prepare_all(td("qiime_L7_rarefied.tsv"), td("metadata_crossover_dropout.csv"), cfg))
  g <- quiet(run_global_test(obj, "Treatment", covariates = c("Timepoint", "Sequence"),
                             strata = "CatID", permutations = 999))
  expect_match(g$Strata[1], "1 single-sample block kept fixed")
  expect_true(is.na(g$p_value[g$Variable == "Sequence"]))       # constant within cat
  expect_match(g$Note[g$Variable == "Sequence"], "does not vary within CatID")
  expect_false(is.na(g$p_value[g$Variable == "Treatment"]))
  expect_equal(unique(g$Permutations), 2^7 - 1)                 # 7 cats with 2 samples
  expect_match(g$Note[g$Variable == "Treatment"], "exact test: all 128 possible permutations")
  g2 <- quiet(run_global_test(obj, "Treatment", covariates = c("Timepoint", "Sequence"),
                              strata = "CatID", permutations = 99))
  expect_identical(g$p_value, g2$p_value)                       # exact -> same every run
  a <- quiet(run_global_test(obj, "Sex", "anosim", strata = "CatID", permutations = 99))
  expect_true(is.na(a$p_value))                                 # Sex never changes within a cat
  p <- quiet(run_pairwise_test(obj, "Treatment", strata = "CatID", permutations = 99))
  expect_equal(nrow(p), 1)
  expect_match(p$Strata, "kept fixed")
  expect_equal(p$Permutations, 2^7 - 1)
})

test_that("cross-metric adjustment matches treatment terms exactly", {
  cfg <- pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline")
  obj <- quiet(prepare_all(td("qiime_L7_rarefied.tsv"), td("metadata_crossover_dropout.csv"), cfg))
  me <- obj$meta
  me$Treatment2 <- ifelse(me$Sex == "F", "x", "y")              # name starts like 'Treatment'
  r <- quiet(run_alpha_lmm_all_metrics(me, "Treatment2 + Timepoint + Treatment", "CatID", "CON"))
  hit <- unique(r$Term[!is.na(r$p_adj_across_metrics)])
  expect_identical(hit, "TreatmentTEST")
  ri <- quiet(run_alpha_lmm_all_metrics(me, "Treatment * Sex", "CatID", "CON", fdr_target = "auto"))
  expect_identical(attr(ri, "fdr_target"), "interaction")
  expect_true(all(grepl(":", ri$Term[!is.na(ri$p_adj_across_metrics)])))
  rm <- quiet(run_alpha_lmm_all_metrics(me, "Treatment + Sex", "CatID", "CON", fdr_target = "auto"))
  expect_identical(attr(rm, "fdr_target"), "main")
})

# ---- Pathway tables, incomplete subjects, cross-sectional models, reference groups ----

test_that("HUMAnN pathway table: suffixed sample names, proportions, no taxonomy", {
  otu <- td("humann_pathabundance_relab.tsv"); md <- td("metadata_messy.csv")
  obj <- quiet(prepare_all(otu, md, pcoa_config(subject_id_col = "", is_baseline_col = "")))
  expect_equal(obj$match_info$n_matched, 10)
  expect_identical(obj$match_info$table_only, "D11_kneaddata_combined_Abundance")
  expect_identical(obj$tax_levels, "Species")
  expect_equal(obj$otu_scale_info$n_proportion, 10)
  tb <- quiet(compute_taxonomy_barplot_data(otu, obj$meta, "Species", top_n = 5))
  expect_true(all(abs(tapply(tb$RelAbundance, tb$SampleID, sum) - 100) < 1e-6))
  expect_true("PWY-0001: glycolysis I" %in% tb$Taxon)
  mi <- quiet(prepare_maaslin_input(otu, obj$meta, "Species", obj$config))
  expect_equal(dim(mi$features), c(10L, 16L))
  expect_equal(unname(rowSums(mi$features)), rep(100, 10), tolerance = 1e-5)
  ds <- quiet(summarize_taxa_by_group(otu, obj$meta, "Group"))
  expect_equal(nrow(ds$summary), 32)                       # 16 pathways x 2 groups
  expect_true(max(ds$summary$max_abundance) > 1)           # percent, not 0-1
})

test_that("animals with missing samples are reported and can be left out", {
  otu <- td("qiime_L7_rarefied.tsv"); md <- td("metadata_crossover_dropout.csv")
  keep <- quiet(prepare_all(otu, md, pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline")))
  expect_equal(nrow(keep$meta), 23)
  expect_identical(keep$match_info$incomplete, c(K08 = 2L))
  expect_false(keep$match_info$incomplete_dropped)
  drop <- quiet(prepare_all(otu, md, pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline",
                                                  complete_subjects_only = TRUE)))
  expect_equal(nrow(drop$meta), 21)
  expect_false("K08" %in% drop$meta$CatID)
  expect_true(drop$match_info$incomplete_dropped)
  g <- quiet(run_global_test(drop, "Treatment", strata = "CatID", permutations = 99))
  expect_identical(g$Strata, "CatID")                      # no single-sample block left
  none <- quiet(prepare_all(otu, md, pcoa_config(subject_id_col = "", is_baseline_col = "IsBaseline",
                                                  complete_subjects_only = TRUE)))
  expect_equal(nrow(none$meta), 23)                        # no subject column -> nothing to do
})

test_that("alpha models fall back to a linear model without a random effect", {
  obj <- quiet(prepare_all(td("metaphlan_species_only.txt"), td("metadata_messy.csv"),
                           pcoa_config(subject_id_col = "", is_baseline_col = "")))
  r <- quiet(run_alpha_lmm_all_metrics(obj$meta, "Group + BW.kg", random_effect = NULL,
                                       treatment_ref = "Control", treatment_col = "Group",
                                       is_baseline_col = ""))
  expect_false(any(r$Term == "ERROR"))
  expect_false(any(grepl("(1|", r$Formula, fixed = TRUE)))
  expect_true("GroupCase" %in% r$Term)                     # Case vs reference Control
  expect_true(all(is.na(r$df)))
  expect_equal(sum(!is.na(r$p_adj_across_metrics)), length(unique(r$Metric)))
  one <- stats::lm(Shannon ~ Group + BW.kg, transform(obj$meta, Group = relevel(factor(Group), "Control")))
  expect_equal(r$Estimate[r$Metric == "Shannon"], unname(stats::coef(one)))
})

test_that("MaAsLin3 reference group can be chosen", {
  skip_if_not_installed("maaslin3")
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv")))
  mi  <- quiet(prepare_maaslin_input(td("qiime_counts.tsv"), obj$meta, "Family"))
  a <- quiet(run_maaslin(mi$features, mi$metadata, "Treatment + Timepoint"))
  b <- quiet(run_maaslin(mi$features, mi$metadata, "Treatment + Timepoint",
                         reference = c(Treatment = "TRT")))
  expect_identical(unname(attr(a, "reference_used")["Treatment"]), "CON")
  expect_identical(unique(a$value[a$metadata == "Treatment"]), "TRT")
  expect_identical(unname(attr(b, "reference_used")["Treatment"]), "TRT")
  expect_identical(unique(b$value[b$metadata == "Treatment"]), "CON")
  fa <- a[a$metadata == "Treatment", ]; fb <- b[b$metadata == "Treatment", ]
  k <- intersect(paste(fa$feature, fa$model), paste(fb$feature, fb$model))
  expect_equal(fa$coef[match(k, paste(fa$feature, fa$model))],
               -fb$coef[match(k, paste(fb$feature, fb$model))], tolerance = 1e-6)   # mirror image
  expect_warning(run_maaslin(mi$features, mi$metadata, "Treatment", reference = c(Treatment = "nope")),
                 "not a level")
})

test_that("both cross-metric adjustments are reported", {
  obj <- quiet(prepare_all(td("metaphlan_species_only.txt"), td("metadata_messy.csv"),
                           pcoa_config(subject_id_col = "", is_baseline_col = "")))
  me <- obj$meta
  me$G3 <- rep(c("Control", "MUO", "NIBD"), length.out = nrow(me))     # 3 groups
  r3 <- quiet(run_alpha_lmm_all_metrics(me, "G3", random_effect = NULL,
                                        treatment_ref = "Control", treatment_col = "G3",
                                        is_baseline_col = ""))
  expect_setequal(attr(r3, "fdr_terms"), c("G3MUO", "G3NIBD"))
  hit <- !is.na(r3$p_adj_across_metrics)
  # B: each comparison adjusted over the metrics only
  for (trm in c("G3MUO", "G3NIBD")) {
    k <- r3$Term == trm & hit
    expect_equal(r3$p_adj_across_metrics[k], stats::p.adjust(r3$p_value[k], "BH"), info = trm)
  }
  # A: one family over both comparisons and all metrics -> never smaller than B
  expect_equal(r3$p_adj_across_metrics_and_groups[hit],
               stats::p.adjust(r3$p_value[hit], "BH"))
  # usually larger, but BH is not monotone in family size, so no >= assertion here
  expect_false(identical(r3$p_adj_across_metrics_and_groups[hit], r3$p_adj_across_metrics[hit]))
  # two groups -> the two columns agree
  r2 <- quiet(run_alpha_lmm_all_metrics(obj$meta, "Group", random_effect = NULL,
                                        treatment_ref = "Control", treatment_col = "Group",
                                        is_baseline_col = ""))
  k <- !is.na(r2$p_adj_across_metrics)
  expect_equal(r2$p_adj_across_metrics[k], r2$p_adj_across_metrics_and_groups[k])
})

test_that("Observed_taxa replaces Observed_ASVs", {
  expect_true("Observed_taxa" %in% ALPHA_METRICS)
  expect_false("Observed_ASVs" %in% ALPHA_METRICS)
  obj <- quiet(prepare_all(td("qiime_L7_rarefied.tsv"), td("metadata_crossover_dropout.csv"),
                           pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline")))
  expect_true(all(c("Observed_taxa", "Baseline_Observed_taxa") %in% colnames(obj$meta)))
})

test_that("subject candidates are reported when Step 2 names no subject", {
  obj <- quiet(prepare_all(td("qiime_L7_rarefied.tsv"), td("metadata_crossover_dropout.csv"),
                           pcoa_config(subject_id_col = "", is_baseline_col = "IsBaseline")))
  expect_identical(obj$match_info$subject_candidates, "CatID")
  named <- quiet(prepare_all(td("qiime_L7_rarefied.tsv"), td("metadata_crossover_dropout.csv"),
                             pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline")))
  expect_identical(named$match_info$subject_candidates, character(0))
})

# ---- Targeted measurements (qPCR, metabolites) ----

test_that("measurement files are read, tidied and linked through any metadata ID column", {
  obj <- quiet(prepare_all(td("metaphlan_species_only.txt"), td("metadata_messy.csv"),
                           pcoa_config(subject_id_col = "", is_baseline_col = "")))
  expect_message(m <- read_measurements(td("measurements_targeted.csv")), "treated as missing")
  expect_identical(colnames(m), c("myristic.acid.14.0.ug.mg", "Cholic.acid.CA", "Blautia.Log.DNA",
                                  "Dysbiosis.index", "X7alphaHSDH.activity", "coprostanol.ug.mg"))
  expect_identical(unname(attr(m, "labels")[2]), "Cholic acid (CA)")
  expect_true(all(vapply(m, is.numeric, logical(1))))
  expect_true(is.na(m$X7alphaHSDH.activity[1]))                       # "<LOD"
  L <- quiet(link_measurements(m, obj$meta))
  expect_identical(L$link$by, "Lab.NO")                               # not the Sample ID
  expect_equal(nrow(L$data), 9)
  expect_identical(L$link$unlinked_rows, "SEQ_012")
  expect_identical(L$link$samples_without, "D10")
  expect_identical(rownames(L$data), sprintf("D%02d", 1:9))
  expect_equal(L$data["D03", "Blautia.Log.DNA"], m["SEQ_003", "Blautia.Log.DNA"])
  expect_error(quiet(link_measurements(m, data.frame(x = 1:3, row.names = c("q", "r", "s")))),
               "do not match")
})

test_that("measurement models adjust p-values over the selected measurements only", {
  obj <- quiet(prepare_all(td("metaphlan_species_only.txt"), td("metadata_messy.csv"),
                           pcoa_config(subject_id_col = "", is_baseline_col = "")))
  L <- quiet(link_measurements(quiet(read_measurements(td("measurements_targeted.csv"))), obj$meta))
  f4 <- L$features[1:4]
  r4 <- quiet(run_measurement_models(L$data, f4, "Group", NULL, "Control", "Group", labels = L$labels))
  k <- r4$Term == "GroupCase"
  expect_equal(sum(k), 4)
  expect_equal(r4$p_adj_across_features[k], stats::p.adjust(r4$p_value[k], "BH"))
  expect_equal(r4$p_adj_across_features[k], r4$p_adj_across_features_and_groups[k])   # 2 groups
  one <- stats::lm(Dysbiosis.index ~ Group, transform(L$data, Group = relevel(factor(Group), "Control")))
  expect_equal(r4$Estimate[r4$Feature == "Dysbiosis.index"], unname(stats::coef(one)))
  expect_equal(unique(r4$N), 9)
  expect_identical(r4$Label[r4$Feature == "Cholic.acid.CA"][1], "Cholic acid (CA)")
  r2 <- quiet(run_measurement_models(L$data, f4[1:2], "Group", NULL, "Control", "Group"))
  k2 <- r2$Term == "GroupCase"
  expect_equal(r2$p_value[k2], r4$p_value[k][1:2])                    # same models
  expect_equal(r2$p_adj_across_features[k2], stats::p.adjust(r2$p_value[k2], "BH"))
  # log10: zeros -> half the smallest positive value; negatives cannot be logged
  lg <- suppressWarnings(quiet(run_measurement_models(
    L$data, c("coprostanol.ug.mg", "Dysbiosis.index"), "Group", NULL, "Control", "Group",
    transform = "log10")))
  expect_identical(attr(lg, "skipped"), "Dysbiosis.index")
  x <- L$data$coprostanol.ug.mg; x[x == 0] <- min(x[x > 0]) / 2
  byhand <- stats::lm(log10(x) ~ relevel(factor(L$data$Group), "Control"))
  expect_equal(lg$Estimate[lg$Feature == "coprostanol.ug.mg"], unname(stats::coef(byhand)))
  expect_error(run_measurement_models(L$data, character(0), "Group", NULL, "Control", "Group"),
               "at least one")
})

test_that("measurement models use each animal's baseline value in repeated-measures studies", {
  cfg <- pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline")
  obj <- quiet(prepare_all(td("qiime_L7_rarefied.tsv"), td("metadata_crossover_dropout.csv"), cfg))
  set.seed(2); me <- obj$meta
  meas <- data.frame(acid = stats::rnorm(nrow(me), 10) + (me$Treatment == "TEST") * 2,
                     sterol = stats::rexp(nrow(me)), row.names = rownames(me))
  L <- quiet(link_measurements(meas, me))
  expect_identical(L$link$by, "Sample ID")
  r <- quiet(run_measurement_models(L$data, c("acid", "sterol"), "Timepoint + Treatment", "CatID",
                                    "CON", "Treatment", is_baseline_col = "IsBaseline",
                                    subject_id_col = "CatID"))
  expect_match(r$Formula[r$Feature == "acid"][1], "`acid` ~ `Baseline_acid` + ", fixed = TRUE)
  expect_match(r$Formula[r$Feature == "sterol"][1], "`Baseline_sterol`", fixed = TRUE)
  expect_equal(unique(r$N), 15)                                        # post-baseline samples
  expect_lt(r$p_value[r$Feature == "acid" & r$Term == "TreatmentTEST"], 0.05)
  off <- quiet(run_measurement_models(L$data, "acid", "Timepoint + Treatment", "CatID", "CON",
                                      "Treatment", is_baseline_col = "IsBaseline",
                                      subject_id_col = "CatID", baseline_adjust = FALSE))
  expect_false(grepl("Baseline_", off$Formula[1]))
})

test_that("descriptive statistics per group are correct", {
  obj <- quiet(prepare_all(td("metaphlan_species_only.txt"), td("metadata_messy.csv"),
                           pcoa_config(subject_id_col = "", is_baseline_col = "")))
  L <- quiet(link_measurements(quiet(read_measurements(td("measurements_targeted.csv"))), obj$meta))
  d <- summarize_measurements(L$data, c("Blautia.Log.DNA", "X7alphaHSDH.activity"), "Group",
                              labels = L$labels, digits = 8)
  expect_identical(names(d), c("Feature", "Label", "Group", "N", "N_missing", "mean", "sd", "median",
                               "min", "max", "ci95_low", "ci95_high"))
  expect_equal(nrow(d), 4)
  x <- L$data$Blautia.Log.DNA[L$data$Group == "Control"]; e <- d[d$Feature == "Blautia.Log.DNA" & d$Group == "Control", ]
  h <- stats::qt(0.975, length(x) - 1) * stats::sd(x) / sqrt(length(x))
  expect_equal(c(e$N, e$mean, e$sd, e$median, e$min, e$max, e$ci95_low, e$ci95_high),
               c(length(x), mean(x), stats::sd(x), stats::median(x), min(x), max(x), mean(x) - h, mean(x) + h))
  m <- d[d$Feature == "X7alphaHSDH.activity" & d$Group == "Control", ]
  expect_equal(c(m$N, m$N_missing), c(4, 1))                           # "<LOD" counted as missing
  two <- summarize_measurements(L$data, "Blautia.Log.DNA", c("Group", "Sex.1"))
  expect_true(all(grepl("_", two$Group)))
  expect_error(summarize_measurements(L$data, "Blautia.Log.DNA", character(0)), "at least one")
})

# ---- Measurements without sequencing data ----

test_that("metadata loads on its own and measurements use every metadata row", {
  cfg <- pcoa_config(subject_id_col = "", is_baseline_col = "")
  M   <- quiet(prepare_metadata(td("metadata_messy.csv"), cfg))
  raw <- utils::read.csv(td("metadata_messy.csv"), check.names = FALSE, fileEncoding = "UTF-8-BOM")
  expect_equal(nrow(M$meta), sum(nzchar(trimws(raw[[1]]))))           # every row, no table needed
  L <- quiet(link_measurements(quiet(read_measurements(td("measurements_targeted.csv"))), M$meta))
  # with sequencing, only the sequenced samples were linked; now none is lost
  obj <- quiet(prepare_all(td("metaphlan_species_only.txt"), td("metadata_messy.csv"), cfg))
  L2  <- quiet(link_measurements(quiet(read_measurements(td("measurements_targeted.csv"))), obj$meta))
  expect_gte(L$link$n_linked, L2$link$n_linked)
  both <- intersect(rownames(L$data), rownames(L2$data))
  expect_equal(L$data[both, "Blautia.Log.DNA"], L2$data[both, "Blautia.Log.DNA"])
  res <- quiet(run_measurement_models(L$data, "Blautia.Log.DNA", "Group", NULL, "Control", "Group"))
  expect_false(any(res$Term == "ERROR"))
  expect_equal(unique(res$N), sum(!is.na(L$data$Blautia.Log.DNA) & !is.na(L$data$Group)))
})

test_that("prepare_metadata reports and optionally leaves out incomplete subjects", {
  cfg <- pcoa_config(subject_id_col = "CatID", is_baseline_col = "IsBaseline")
  M   <- quiet(prepare_metadata(td("metadata_crossover_dropout.csv"), cfg))
  expect_true(length(M$info[["incomplete"]]) >= 1)
  cfg$complete_subjects_only <- TRUE
  M2  <- quiet(prepare_metadata(td("metadata_crossover_dropout.csv"), cfg))
  expect_false(any(names(M$info[["incomplete"]]) %in% M2$meta$CatID))
  expect_lt(nrow(M2$meta), nrow(M$meta))
})

test_that("the Miome palette gives distinct colours for any number of groups", {
  expect_identical(get_cat_pal("Miome", 3), MIOME_PALETTE[1:3])
  expect_identical(get_cat_pal("Miome", 3, reverse = TRUE), rev(MIOME_PALETTE[1:3]))
  expect_equal(anyDuplicated(get_cat_pal("Miome", 15)), 0)
  expect_length(get_cat_pal("Set2", 5), 5)
})

test_that("taxa summary has mean, SD and 95% CI, and the workbook shows the range as text", {
  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv"),
                           pcoa_config(subject_id_col = "sampleid", is_baseline_col = "IsBaseline")))
  ds <- quiet(summarize_taxa_by_group(td("qiime_counts.tsv"), obj$meta, "Treatment",
                                      tax_levels = "Phylum", digits = 8))$summary
  expect_true(all(c("mean_abundance", "sd_abundance", "ci95_low", "ci95_high") %in% names(ds)))
  r <- ds[1, ]; n <- r$N_samples
  expect_equal(r$ci95_high - r$mean_abundance, stats::qt(0.975, n - 1) * r$sd_abundance / sqrt(n),
               tolerance = 1e-6)
  expect_true(all(ds$min_abundance <= ds$mean_abundance & ds$mean_abundance <= ds$max_abundance))
  # check one taxon by hand from the percent table
  raw <- utils::read.delim(td("qiime_counts.tsv"), check.names = FALSE, row.names = 1, skip = 1)
  raw <- raw[, vapply(raw, is.numeric, logical(1)), drop = FALSE]
  pct <- sweep(as.matrix(raw), 2, colSums(raw), "/") * 100
  g   <- rownames(obj$meta)[obj$meta$Treatment == r$Group]
  x   <- colSums(pct[grepl(paste0("p__", sub(".*p__", "", r$Taxon), "(;|$)"), rownames(pct)), g, drop = FALSE])
  expect_equal(c(r$mean_abundance, r$sd_abundance), c(mean(x), stats::sd(x)), tolerance = 1e-6)

  f <- tempfile(fileext = ".xlsx")
  quiet(export_taxa_descriptive_workbook(td("qiime_counts.tsv"), obj$meta, "Treatment", outfile = f,
                                         tax_levels = "Phylum"))
  w <- openxlsx::read.xlsx(f, colNames = FALSE, skipEmptyRows = FALSE)
  expect_identical(unname(unlist(w[3, 2:9])), c("median_abundance", "min_abundance", "max_abundance",
                   "range", "mean_abundance", "sd_abundance", "ci95_low", "ci95_high"))
  expect_match(w[4, 5], "^[0-9.]+-[0-9.]+$")                           # range is text, not a formula
})

test_that("no incomplete subjects means an empty list entry, not a partial match", {
  M <- quiet(prepare_metadata(td("metadata.csv"), pcoa_config(subject_id_col = "sampleid",
                                                              is_baseline_col = "IsBaseline")))
  expect_length(M$info$incomplete, 0)
  expect_false(M$info$subject_unique)
  messy <- quiet(prepare_metadata(td("metadata_messy.csv"), pcoa_config(subject_id_col = "", is_baseline_col = "")))
  one_each <- names(messy$meta)[vapply(messy$meta, function(x) !anyNA(x) && !anyDuplicated(x), logical(1))][1]
  M2 <- quiet(prepare_metadata(td("metadata_messy.csv"), pcoa_config(subject_id_col = one_each,
                                                                     is_baseline_col = "")))
  expect_true(M2$info$subject_unique)                                  # an ID that never repeats
})

test_that("measurement descriptives export side by side and match the long table", {
  M <- quiet(prepare_metadata(td("metadata_messy.csv"), pcoa_config(subject_id_col = "", is_baseline_col = "")))
  L <- quiet(link_measurements(quiet(read_measurements(td("measurements_targeted.csv"))), M$meta))
  f <- tempfile(fileext = ".xlsx"); feats <- c("Blautia.Log.DNA", "Dysbiosis.index", "X7alphaHSDH.activity")
  long <- quiet(export_measurements_descriptive_workbook(L$data, feats, "Group", outfile = f, labels = L$labels))
  w <- openxlsx::read.xlsx(f, colNames = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
  g <- unique(long$Group)
  expect_equal(dim(w), c(3 + length(feats), 1 + 9 * length(g)))
  expect_identical(unname(unlist(w[2, c(2, 11)])), g[1:2])              # groups side by side
  expect_identical(unname(unlist(w[3, 1:10])), c("Measurement", "N", "median", "min", "max", "range",
                                                 "mean", "sd", "ci95_low", "ci95_high"))
  e <- long[long$Feature == "Dysbiosis.index" & long$Group == g[2], ]
  expect_equal(as.numeric(unlist(w[5, c(11:14, 16:19)])),
               c(e$N, e$median, e$min, e$max, e$mean, e$sd, e$ci95_low, e$ci95_high))
  expect_identical(w[5, 15], paste(format(e$min, drop0trailing = TRUE), "to", format(e$max, drop0trailing = TRUE)))
  expect_identical(w[5, 1], "Dysbiosis index")                         # label as typed
})

test_that("custom colours: fixed colour per group and custom gradients", {
  lv <- c("CON", "TRT", "BASE")
  expect_identical(.level_colours("Set2", lv), stats::setNames(get_cat_pal("Set2", 3), lv))
  x <- .level_colours(c(TRT = "#FF0000"), lv)
  expect_identical(x[["TRT"]], "#FF0000")
  expect_identical(unname(x[c("CON", "BASE")]), MIOME_PALETTE[c(1, 3)])      # others keep their palette colour
  y <- .level_colours(structure(c(CON = "#000000", nonsense = "red", GONE = "#123456"), base = "Set1"), lv)
  expect_identical(y[["CON"]], "#000000")
  expect_identical(unname(y[c("TRT", "BASE")]), get_cat_pal("Set1", 3)[2:3])
  expect_identical(.cont_ramp(c("#000000", "#FFFFFF"), n = 3), c("#000000", "#7F7F7F", "#FFFFFF"))
  expect_identical(.cont_ramp(c("#000000", "#FFFFFF"), reverse = TRUE, n = 3)[1], "#FFFFFF")
  expect_length(.cont_ramp("RdBu"), 11)

  obj <- quiet(prepare_all(td("qiime_counts.tsv"), td("metadata.csv"),
                           pcoa_config(subject_id_col = "sampleid", is_baseline_col = "IsBaseline")))
  fills <- function(p) { b <- ggplot2::ggplot_build(p); unique(b$data[[length(b$data)]][c("fill", "group")]) }
  p <- build_2d(obj, "Treatment", palette = c(TRT = "#FF0000"))
  sc <- p$scales$get_scales("fill")$palette(3)
  expect_true("#FF0000" %in% sc)
  expect_s3_class(build_2d(obj, "Shannon", palette = c("#2F6FB0", "#F2F2F2", "#B8456B"), is_cont = TRUE), "ggplot")
  expect_s3_class(build_3d(obj, "Shannon", palette = c("#2F6FB0", "#B8456B"), is_cont = TRUE), "plotly")
  expect_s3_class(build_3d(obj, "Treatment", palette = c(TRT = "#FF0000")), "plotly")
  a <- build_alpha_scatter(obj$meta_enriched, "Shannon", "Treatment", "Treatment", palette = c(CON = "#00FF00"))
  expect_true("#00FF00" %in% ggplot2::ggplot_build(a)$data[[length(ggplot2::ggplot_build(a)$data)]]$colour)
})
