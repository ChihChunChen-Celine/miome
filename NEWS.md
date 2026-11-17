# miome 0.8.1

* **MaAsLin3 with two or more random effects** no longer ends in MaAsLin3's
  own error "supplied row-names must match list extent". That error comes from
  MaAsLin3 itself: versions before 1.5.5, including the Bioconductor release
  1.4.0, cannot handle more than one random effect (fixed upstream in 1.5.5).
  With such a version miome now stops before fitting and says so, and the app
  shows the same note as soon as a second random effect is chosen. One random
  effect works with every version. To use more, update MaAsLin3:
  `devtools::install_github("biobakery/maaslin3")`.
* The MaAsLin3 test no longer fails `devtools::check()` on computers with an
  older MaAsLin3; it checks the clear error there instead.

# miome 0.8.0

## Your own plot colours

* **A colour for each group.** Under every plot (alpha diversity, PCoA,
  measurements) "Choose a colour for each group" opens one colour picker per
  group of the "Colour by" column. A group keeps its colour in every plot and
  in the TIFF downloads; groups you do not touch keep their palette colour.
  "Back to the palette colours" undoes all choices. The taxonomy bar plot has
  the same per taxon.
* **Your own gradient.** PCoA, numeric colouring: choose "My own colours" as
  gradient and pick the low and high colour, optionally a middle colour.
* In R, `palette` accepts the same: a named vector such as
  `c(Control = "#1B9E77", MUO = "#D95F02")`, or for gradients an unnamed
  vector of colours from low to high.

# miome 0.7.2

* **Measurements: descriptive statistics as an Excel sheet** in the same layout
  as the taxa summary: one row per measurement, the groups side by side under
  merged headers, with N, median, min, max, range, mean, SD and 95% CI of the
  mean. New function `export_measurements_descriptive_workbook()`. The .csv
  download (one row per measurement and group) is kept.

# miome 0.7.1

## Bug fixes

* **Settings / result text boxes stayed empty** (0.7.0 only): "Settings used"
  after MaAsLin3, the test settings and PERMDISP on the Beta statistics page.
  The numbers were computed but never displayed.
* **False warning "Animals with missing samples are still included: ."** It
  appeared after every load when no animal was incomplete, ticked or not.
* **Taxa summary (Excel): the range column was empty** when the file was opened
  in Protected View or a preview, because it was written as a formula. It is
  now written as text (`export_taxa_descriptive_workbook(range_as_formula =)`
  defaults to `FALSE`).

## New

* Taxa summary (Excel) gains **mean, standard deviation and the 95% CI of the
  mean** per group, after median / min / max / range.
* The app tells you when the column chosen as subject / animal ID never repeats
  (e.g. a sequencing number): it is then not offered as random effect, strata
  or for paired lines.
* Alpha diversity and Measurements stop with a clear message when the formula
  does not contain the column chosen under "Groups to compare".

## Wording

* "Leave out animals with missing samples" is now "Leave out animals that miss
  a timepoint", shown only when a subject ID is set, with an explanation.
* "Treatment column" is now "Groups to compare (column)", with an explanation;
  the random-effect box says "Leave empty if each animal has one sample".

# miome 0.7.0

## Measurements no longer need sequencing data

* **Load the metadata on its own.** The abundance table is now optional: with
  only a metadata file loaded, the Measurements page works in full (models,
  descriptive statistics, plot). The sequencing pages say what is missing.
* **Measurements use every metadata row**, whether or not that sample was
  sequenced. Before, a sample with a qPCR or metabolite value but no
  sequencing data was left out. When sequencing data is loaded as well, alpha
  diversity and baseline columns are still available as covariates (empty for
  rows that were not sequenced).
* New function `prepare_metadata()`.
* Unchanged, as agreed: entries such as `<LOD` stay missing (no substitution),
  and each measurement file is analysed and p-adjusted on its own.

## New layout

* A **navigation rail** on the left replaces the long sidebar: Data, Alpha
  diversity, PCoA, Beta statistics, Taxonomy, MaAsLin3, Measurements.
* Each page shows its **settings next to its results**. The two scroll
  separately, so results stay in view while you change settings.
* Explanations are folded behind short questions ("How do I write this?");
  rarely changed MaAsLin3 settings sit under "More settings".
* New colours and type; each analysis has its own colour.
* **Plot colours** are chosen once, at the bottom of the navigation, for every
  plot. New default palette `"Miome"` (`MIOME_PALETTE`), checked for
  colour-vision deficiency; the RColorBrewer palettes remain available.

## Changes

* Alpha diversity has its own "P-value adjustment" setting (it used to borrow
  the one from the beta-diversity settings) and its own "Join samples with a
  line" tick box; so does the measurement plot.
* PCoA: ticking "This column is a number" shows a separate gradient choice.
* The metadata preview moved to the Data page; decimals are rounded for
  display (the download keeps full precision).

# miome 0.6.0

## New: Measurements tab (targeted data)

* Upload a table of **targeted measurements** - qPCR, fatty acids / sterols,
  bile acids or any other per-sample values - after loading your study.
  * The first column is an ID; the app works out whether it is the Sample ID
    or another metadata column (e.g. a sequencing number) and reports which
    rows and samples had no partner.
  * **Choose the measurements to analyse.** p-values are adjusted across
    exactly the ones you select (`p_adj_across_features`, and
    `p_adj_across_features_and_groups` for 3+ groups).
  * Same model options as alpha diversity: fixed-effects formula, optional
    random effect(s), reference group, baseline value of each measurement as
    covariate when baseline samples exist. Optional log10 transform.
  * **Descriptive statistics** per group (one or more grouping columns): N,
    mean, SD, median, min, max and the 95% CI of the mean.
  * Box/point plot of any measurement.
* New functions: `read_measurements()`, `link_measurements()`,
  `run_measurement_models()`, `summarize_measurements()`.

## Changes

* Model tables (alpha and measurements) gain an `N` column: the number of
  samples each model actually used.
* Column names with Greek letters are now the same on every computer: the
  micro sign becomes `u` (so "ug/mg"), alpha becomes `alpha`, beta becomes
  `beta` (a column typed as "7[alpha]HSDH activity" is `X7alphaHSDH.activity`).
  Other accented characters are dropped from column names. The names as typed
  are kept as labels in tables and plots, and accented text inside the cells
  is preserved.

# miome 0.5.0

## Results that change compared with 0.4.0

* `Observed_ASVs` is renamed **`Observed_taxa`** (it counts the rows of the
  uploaded table, which are species or pathways, not ASVs). The baseline
  covariate becomes `Baseline_Observed_taxa`.
* The alpha table's `p_adjusted_cross_metric` column is replaced by **two**
  columns, so both ways of correcting are available side by side. They are
  identical when the treatment has two groups:
  * `p_adj_across_metrics` - one comparison (e.g. MUO vs Control) adjusted
    over the alpha metrics only.
  * `p_adj_across_metrics_and_groups` - every treatment term of every metric
    adjusted as a single family (the old behaviour). Usually the larger of the
    two, though Benjamini-Hochberg is not monotone in family size.

## New

* The app warns when a column that looks like an animal/subject ID has
  repeated values but no subject is set in Step 2, so repeated samples are
  not silently treated as independent. Detection is by column name
  (CatID, SubjectID, Microchip, PairedID, ...) plus actual repeats, and stays
  silent otherwise - it cannot replace naming the column in Step 2.
* README states that shotgun studies need the full MetaPhlAn merged table.

# miome 0.4.0

Checked against real project files (MetaPhlAn species table + cross-sectional
metadata; QIIME 2 L7 table + crossover metadata with a drop-out).

## Results that change compared with 0.2.0

* **Restricted permutations now survive a drop-out.** One animal with a single
  sample used to switch the restriction off for the whole PERMANOVA/ANOSIM
  (silently unrestricted). That sample now simply stays in place and all other
  animals are shuffled within themselves. For crossover studies with a
  drop-out the Treatment p-value can change a lot.
* **Exact test when feasible.** If at most 10,000 within-animal shuffles
  exist, all of them are used: the p-value is exact and identical every run.
* **Terms that never change within an animal** (Sequence, Sex, baseline
  covariates) get no p-value under restricted permutations, with a note,
  because shuffling within animals cannot test them.
* **Taxa not resolved to the chosen level are kept**, labelled
  `<deepest known rank>_unclassified` (e.g. `Streptococcus_unclassified`), in
  the taxonomy plot and MaAsLin3. Previously they were pooled as
  "Unclassified" and dropped from MaAsLin3 (at Species level a 16S table lost
  13-62% of each sample).
* Taxa that are absent from every analysed sample (metadata is a subset of
  the run) are removed before distances and MaAsLin3. Aitchison distances
  change slightly.
* **Metadata column names are tidier**: `BW (kg)` is now `BW.kg` (was
  `BW..kg.`), `Age at presentation (years)` is `Age.at.presentation.years`.
* The cross-metric adjusted p-value matches treatment terms exactly: a column
  whose name merely starts like the treatment column (`Group2` next to
  `Group`) is no longer pulled in.

## New

* **Pathway / gene tables (HUMAnN)** are supported: sample names with a
  pipeline suffix (`GI220_kneaddata_combined_Abundance`) are matched to the
  metadata ID (`GI220`), proportions (each sample sums to 1) are shown as
  percent, and the rows are analysed as uploaded.
* **Animals with missing samples** (e.g. a drop-out) are reported after
  loading. Tick "Leave out animals with missing samples" in Step 2 to exclude
  them from every analysis (`pcoa_config(complete_subjects_only = TRUE)`).
* **Reference (control) group for MaAsLin3**: a dropdown appears for every
  categorical variable in the formula; the choice is shown in the info panel
  (`run_maaslin(reference = c(GROUP = "HEALTHY"))`).
* **Alpha diversity for one-sample-per-animal studies**: leave the random
  effect empty and an ordinary linear model is fitted.
* The taxonomy plot now uses exactly the samples of the loaded analysis.
* **Tables with taxon names only** (a MetaPhlAn export reduced to
  `Collinsella_intestinalis`, no `k__...|s__` lineage) now work in the
  taxonomy plot, descriptive summary and MaAsLin3 as "Taxa as uploaded".
  Before, the plot showed one "Unclassified" bar and MaAsLin3 got 0 features.
* **Sample-matching report** after loading: how many samples are analysed,
  which sequenced samples are not in the metadata, and which metadata rows
  have no sequencing data.
* The Sample ID column is found automatically (the metadata column that
  matches the table's sample names), so it need not be the first column.
* Metadata reader handles Excel's UTF-8 mark, tab-separated files, blank
  rows, repeated column names (reported, 2nd copy gets `.1`), duplicated
  Sample IDs (clear error) and percent text such as `28.60%` (read as number).
* **Random seed** for the beta-diversity tests, shown in the info panel.
* Pairwise tests report the strata and permutations actually used.

# miome 0.2.0

## Results that change compared with 0.1.2

* **MetaPhlAn input now works.** Lineages separated by `|` and the `k__`
  prefix are recognised everywhere. Previously the taxonomy plot showed a
  single "Unclassified" bar and MaAsLin3 received no features.
* **No more double counting on all-ranks tables** (un-filtered MetaPhlAn):
  the taxonomy plot, MaAsLin3, and alpha/beta diversity now use only the rows
  at the right rank. Alpha richness and distances from such tables will differ.
* **Baseline covariates for every alpha metric.** `Baseline_Simpson`,
  `Baseline_Observed_taxa`, `Baseline_Chao1` and `Baseline_Pielou` are now
  created and used by the LMM; before, only Shannon was baseline-adjusted.
  `Baseline_Observed` is renamed `Baseline_Observed_taxa`.
* **Low-depth count tables are now normalised.** A count table whose largest
  value is <= 100 was previously treated as already-percent. One rule is now
  used everywhere: whole-number samples that don't already sum to 100 are
  counts.

## Bug fixes

* Blank, `Yes/No`, `Y/N` or `1/0` values in the baseline column no longer
  crash the baseline step (which then silently skipped all covariates).
* QIIME 2 `biom convert` tables with a `#OTU ID` header and a `taxonomy`
  column are read correctly.
* `export_abundance_workbook()`: sample labels now line up with the sample
  columns, and all sheets contain the same taxa. `remove_first_row` now
  defaults to `"auto"` (only drops a non-numeric first row).
* Paired lines in PCoA plots follow time order instead of PC1 order.
* Alpha diversity plot: paired lines end on the points, and integer metrics
  (Observed, Chao1) are no longer jittered vertically.
* Taxonomy bar plot: no repeated colours above 12 taxa; two different taxa
  with the same name (e.g. "uncultured") keep separate labels.
* `prepare_all()` no longer returns a `meta_path` pointing at a deleted
  temporary file.
* 3D TIFF export: fixed a call to a non-existent `magick` function (the true
  3D path never worked) and no longer writes PNG data into a `.tiff` file.
* App: "Run alpha diversity LMM" no longer silently does nothing when there
  is no baseline column; the beta-diversity info panel and file name now
  record a typed formula; PERMDISP uses the same samples as the test; the
  MaAsLin3 help text no longer claims a baseline covariate is added.
* Global PERMANOVA now lists **every fitted term** (covariates first, with a
  `Role` column), prints the model to the R console, and the `Strata` column
  says when the requested strata could not be applied (single-sample blocks).
* "(none)" can be chosen again in every dropdown after a column was picked
  (e.g. "Restrict permutations within" for unrestricted permutations).
* MaAsLin3 results are read correctly (quoted column names in
  `all_results.tsv` are handled).
* The LMM skips metrics with no values (e.g. Chao1 on percent data) instead
  of returning an "ERROR" row.

## New

* **Step 2 - Match your columns** in the app: choose the subject/animal ID
  and baseline flag columns, which are needed while the data loads (the app
  guesses from common names; either can be left as "(none)"). The treatment
  column is now chosen in the Alpha tab next to its reference level, and the
  formula/covariate boxes start empty for you to fill in. A new "Order lines
  by" dropdown controls the order of paired lines on PCoA plots.
* The alpha diversity LMM accepts **several random effects**: each chosen
  column becomes its own random intercept, e.g. `(1|Subject) + (1|Litter)`.
  `random_effect` in `run_alpha_lmm()` and friends now takes a vector.
  The same applies to the MaAsLin3 tab and `run_maaslin()`.
* The cross-metric adjusted p-value targets the treatment interaction when
  the LMM formula contains `*` or `:`.
* `run_alpha_lmm_all_metrics()` gains `add_baseline` and `fdr_term_pattern`;
  `run_alpha_lmm_interaction_all()` is now a shortcut for it.
* Test suite in `tests/testthat` with small synthetic datasets in
  `inst/testdata`.

## Internal

* Shared helpers in `R/utils.R` replace five copies of the table reader,
  three metadata cleaners, three taxonomy parsers and three
  counts-vs-percent rules. The pipeline reads each file once (no temporary
  CSV round trip).
* Roxygen documentation fixed; `devtools::check()` reports no package
  warnings. `rlang`, `tibble`, `tidyr` and `ape` dropped from dependencies.
