# app.R
# -----------------------------------------------------------------------------
# Shiny app. Analysis code lives in data.R / stats.R / plots.R / export.R /
# maaslin3.R / descriptive_summary.R / measurements.R, the look and layout in
# ui.R - this file only wires inputs to them.
# -----------------------------------------------------------------------------

# -- Small UI helpers -----------------------------------------------------------
.info <- function(...) shiny::tags$div(class = "info-box", ...)
.warn <- function(...) shiny::tags$div(class = "warn-box", ...)
# "(none)" needs a real value: Shiny's dropdown treats "" as "nothing chosen",
# so an empty-string option can never be re-selected once a column is picked.
.NONE <- "__none__"
.none <- function(cols, label = "(none)") c(stats::setNames(.NONE, label), cols)
.nz   <- function(x) if (is.null(x) || identical(x, .NONE)) "" else x

# Guess a metadata column from common names (case-insensitive); "" if none.
.guess_col <- function(cols, candidates) {
  hit <- cols[match(tolower(candidates), tolower(cols))]
  hit <- hit[!is.na(hit)]
  if (length(hit)) hit[1] else .NONE
}

.TREAT_NAMES <- c("Treatment", "Group", "Diet", "Arm", "Condition")
.TIME_NAMES  <- c("Timepoint", "Time", "Period", "Day", "Visit", "Phase")

# DT table with rounding and traffic-light p-value colouring.
.dt_table <- function(tbl, round_cols = NULL, colour = list(), page = 10,
                      filter = "none", dom = "tip") {
  dt <- DT::datatable(tbl, rownames = FALSE, filter = filter,
                      options = list(scrollX = TRUE, pageLength = page, dom = dom))
  rc <- intersect(round_cols, colnames(tbl))
  if (length(rc)) dt <- DT::formatRound(dt, rc, digits = 4)
  for (col in intersect(names(colour), colnames(tbl)))
    dt <- DT::formatStyle(dt, col, backgroundColor = DT::styleInterval(
      colour[[col]], c("#DCEBD5", "#F8ECCB", "white")))
  dt
}

# Turn the scale-detection list from .to_relative_abundance() into a sentence.
.format_scale_info <- function(info) {
  if (is.null(info)) return(NULL)
  if (!is.null(info$note)) return(info$note)
  if (isTRUE(info$n_proportion > 0) && info$n_converted == 0)
    sprintf("%d sample(s) were proportions (each sums to 1) - shown as percent.",
            info$n_proportion)
  else if (info$n_converted == 0)
    sprintf("All %d sample(s) were already relative abundance - used as uploaded.",
            info$n_total)
  else if (info$n_already_percent == 0)
    sprintf("All %d sample(s) looked like raw counts - converted to percent.",
            info$n_total)
  else
    sprintf("%d of %d sample(s) were already percent; %d raw-count sample(s) were converted: %s",
            info$n_already_percent, info$n_total, info$n_converted,
            paste(info$converted_samples, collapse = ", "))
}

.today <- function() format(Sys.Date(), "%Y%m%d")

#' Launch the Microbiome Explorer Shiny app
#'
#' @return A Shiny app object (opens in your browser).
#' @export
run_app <- function() {

  # ============================== UI ==========================================
  pcs <- stats::setNames(1:3, paste0("PC", 1:3))
  ui <- shiny::fluidPage(
    title = "miome",
    shiny::tags$head(shiny::tags$style(shiny::HTML(.miome_css())),
                     shiny::tags$script(shiny::HTML(.miome_js()))),
    shiny::tags$div(class = "mi-app",
      .nav_rail(),
      shiny::tags$main(class = "mi-main",

        .page("data", active = TRUE,
          side = shiny::tagList(
            .step(1, "Upload"),
            shiny::fileInput("meta_file", "Metadata (.csv, one row per sample)", accept = ".csv"),
            shiny::fileInput("otu_file", "Abundance table (.tsv, samples as columns)",
                             accept = c(".tsv", ".txt")),
            .help("Which abundance table?",
                  "16S: the QIIME 2 species (L7) table. Shotgun: the full MetaPhlAn merged table ",
                  "(merged_all), or a HUMAnN pathway table. Leave it empty if you only want to ",
                  "analyse measurements such as qPCR or metabolites."),
            shiny::uiOutput("col_map_ui"),
            .step(3, "Distance"),
            shiny::selectInput("dist_method", "Distance between samples:", choices = DIST_METHODS,
                               selected = "bray"),
            shiny::conditionalPanel(
              "input.dist_method == 'unifrac_weighted' || input.dist_method == 'unifrac_unweighted'",
              shiny::fileInput("unifrac_dist_file", "UniFrac distance matrix (.qza or .tsv)",
                               accept = c(".qza", ".tsv", ".txt")),
              shiny::uiOutput("unifrac_warning_ui")),
            .run_btn("load_btn", "Load data")),
          .card("What was loaded", shiny::uiOutput("match_info_ui"),
                shiny::uiOutput("otu_scale_info_sidebar_ui")),
          .card("Metadata", tools = shiny::downloadButton("dl_meta", "Download (.csv)"),
                hint = "Your metadata as the app reads it. With sequencing data, alpha diversity and baseline columns are added.",
                DT::DTOutput("meta_table"))),

        .page("alpha", side = shiny::uiOutput("alpha_side_ui"),
          .card("Plot", shiny::uiOutput("alpha_controls_ui"),
                shiny::plotOutput("alpha_plot", height = "460px")),
          .card("Model results", tools = shiny::downloadButton("dl_alpha_lmm", "Download (.csv)"),
                hint = "One model per alpha metric, fitted when you press Run alpha diversity model.",
                shiny::uiOutput("alpha_lmm_info_ui"), shiny::uiOutput("alpha_lmm_padj_info_ui"),
                DT::DTOutput("alpha_lmm_table"))),

        .page("pcoa", side = shiny::uiOutput("pcoa_side_ui"),
          .card(NULL, shiny::tabsetPanel(type = "pills",
            shiny::tabPanel("2D", shiny::plotOutput("plot2d", height = "620px")),
            shiny::tabPanel("3D, interactive", plotly::plotlyOutput("plot3d", height = "620px"))))),

        .page("beta", side = shiny::uiOutput("beta_side_ui"),
          .card("Test result", tools = shiny::downloadButton("dl_stat", "Download (.csv)"),
                hint = "Appears when you press Run beta diversity test.",
                shiny::uiOutput("otu_scale_info_ui"),
                shiny::verbatimTextOutput("stat_info"),
                DT::DTOutput("stat_table")),
          .card("PERMDISP: is the spread the same in each group?",
                hint = paste("If PERMDISP is significant, read PERMANOVA with caution: the groups may",
                             "differ in how spread out they are rather than in where they sit."),
                shiny::verbatimTextOutput("permdisp_out"))),

        .page("taxonomy", side = shiny::uiOutput("taxon_side_ui"),
          .card("Composition per sample",
                tools = shiny::downloadButton("dl_taxon", "Download plot (TIFF)"),
                shiny::plotOutput("taxon_plot", height = "640px"))),

        .page("maaslin", side = shiny::uiOutput("maaslin_side_ui"),
          .card("Settings used", hint = "Appears when MaAsLin3 has finished.",
                shiny::verbatimTextOutput("maaslin_info")),
          .card("Significant hits",
                hint = "Results with a q-value below the maximum significance set on the left.",
                DT::DTOutput("maaslin_sig_table")),
          .card("All results", tools = shiny::downloadButton("dl_maaslin", "Download (.csv)"),
                DT::DTOutput("maaslin_table"))),

        .page("meas",
          side = shiny::tagList(
            .section("Measurement file"),
            shiny::fileInput("meas_file", "Table (.csv: one row per sample, one column per measurement)",
                             accept = c(".csv", ".tsv", ".txt")),
            .help("How is it matched to my metadata?",
                  "The first column must hold an ID that is also in your metadata: the Sample ID ",
                  "or another ID column such as a sequencing number. The app finds the match. ",
                  "Entries such as <LOD are kept as missing."),
            shiny::uiOutput("meas_controls_ui")),
          .card("Matched samples", shiny::uiOutput("meas_link_ui")),
          .card("Plot", shiny::uiOutput("meas_plot_controls_ui"),
                shiny::plotOutput("meas_plot", height = "400px")),
          .card("Model results", tools = shiny::downloadButton("dl_meas_model", "Download (.csv)"),
                hint = "One model per selected measurement, fitted when you press Run measurement models.",
                shiny::uiOutput("meas_model_info_ui"), DT::DTOutput("meas_model_table")),
          .card("Descriptive statistics per group",
                tools = shiny::tagList(
                  shiny::downloadButton("dl_meas_desc_xlsx", "Download Excel (groups side by side)"),
                  shiny::downloadButton("dl_meas_desc", "Download (.csv)")),
                hint = paste("N, median, minimum, maximum, mean, standard deviation and the 95% confidence",
                             "interval of the mean, for the selected measurements as uploaded (not",
                             "transformed). Choose the grouping columns on the left."),
                DT::DTOutput("meas_desc_table")))
      ))
  )

  # ============================== SERVER ======================================
  server <- function(input, output, session) {

    # obj / meta_enriched exist only with sequencing data; meta_all (every
    # metadata row, sequenced or not) and cfg exist as soon as metadata is loaded.
    rv            <- shiny::reactiveValues(obj = NULL, meta_cols = NULL, meta_enriched = NULL,
                                           cfg = NULL, meta_all = NULL, meta_info = NULL)
    stat_store    <- shiny::reactiveValues(result = NULL, params = NULL)
    maaslin_store <- shiny::reactiveValues(result = NULL, params = NULL)
    dist_labels   <- stats::setNames(names(DIST_METHODS), DIST_METHODS)

    # -- Step 2: match metadata columns ----------------------------------------
    # Same reader as the analysis, so the names offered here are the names used.
    meta_header <- shiny::reactive({
      if (is.null(input$meta_file)) return(NULL)
      ids <- if (!is.null(input$otu_file))
        tryCatch(.peek_sample_ids(input$otu_file$datapath), error = function(e) NULL)
      tryCatch(colnames(suppressMessages(.read_meta(input$meta_file$datapath,
                                                    sample_ids = ids))),
               error = function(e) {
                 shiny::showNotification(paste("Metadata:", conditionMessage(e)),
                                         type = "error", duration = NULL)
                 character(0)
               })
    })

    output$col_map_ui <- shiny::renderUI({
      cols <- meta_header()
      if (is.null(cols))
        return(shiny::tagList(.step(2, "Match your columns"),
                              shiny::tags$p(class = "mi-hint", "Appears once the metadata is uploaded.")))
      if (!length(cols))
        return(.warn("Could not read the metadata file. Check it is a .csv with Sample IDs in the first column."))
      pick <- function(id, label, cands)
        shiny::selectInput(id, label, choices = .none(cols), selected = .guess_col(cols, cands))
      shiny::tagList(
        .step(2, "Match your columns"),
        pick("map_subject", "Subject / animal ID:",
             c("sampleid", "CatID", "SubjectID", "Subject", "AnimalID", "Animal", "PatientID", "ID")),
        pick("map_baseline", "Baseline flag (TRUE/FALSE, Yes/No):",
             c("IsBaseline", "Baseline", "is_baseline")),
        shiny::conditionalPanel(sprintf("input.map_subject != '%s'", .NONE),
          shiny::checkboxInput("complete_only",
                               "Leave out animals that miss a timepoint (e.g. a drop-out)", FALSE),
          .help("What does this do?",
                "Only for studies where each animal is sampled several times. An animal with ",
                "fewer rows in the metadata than the others (say 2 of 3 timepoints) is ",
                "'incomplete'. Unticked: its available samples are still used. Ticked: all of ",
                "its samples are left out of every analysis. This is not about samples without ",
                "sequencing data or without a measurement; those are skipped automatically ",
                "wherever the value is needed.")),
        .help("What are these for?",
              "Baseline samples are left out of the statistics, and each subject's baseline ",
              "value becomes a Baseline_<metric> covariate. Leave both at (none) if your study ",
              "has no repeated measures or no baseline samples."))
    })

    output$unifrac_warning_ui <- shiny::renderUI({
      if (is.null(input$unifrac_dist_file))
        .warn("UniFrac needs a distance matrix computed in QIIME 2. Upload the .qza or .tsv file above.")
    })

    # -- Load & compute ---------------------------------------------------------
    shiny::observeEvent(input$load_btn, {
      if (is.null(input$meta_file)) {
        shiny::showNotification("Please upload the metadata file first.", type = "error")
        return()
      }
      has_seq    <- !is.null(input$otu_file)
      is_unifrac <- has_seq && input$dist_method %in% c("unifrac_weighted", "unifrac_unweighted")
      if (is_unifrac && is.null(input$unifrac_dist_file)) {
        shiny::showNotification("Please upload a precomputed UniFrac distance matrix first.",
                                type = "error")
        return()
      }

      subj <- .nz(input$map_subject)
      cfg <- pcoa_config(
        subject_id_col     = subj,
        treatment_col      = "",
        timepoint_col      = "",
        sequence_col       = "",
        is_baseline_col    = .nz(input$map_baseline),
        complete_subjects_only = isTRUE(input$complete_only),
        study_design       = "",
        default_covariates = character(0),
        default_strata     = if (nzchar(subj)) subj else NULL,
        dist_method        = input$dist_method,
        dist_method_label  = dist_labels[[input$dist_method]])

      tryCatch({
        obj <- if (has_seq)
          shiny::withProgress(message = "Computing diversity and PCoA...", value = 0.5,
            prepare_all(input$otu_file$datapath, input$meta_file$datapath, config = cfg,
                        dist_method = input$dist_method,
                        precomputed_dist_path = if (is_unifrac) input$unifrac_dist_file$datapath))
        # Every metadata row, whether sequenced or not (used by Measurements).
        ma <- prepare_metadata(input$meta_file$datapath, cfg,
                               sample_ids = if (has_seq) .peek_sample_ids(input$otu_file$datapath))
        meta_all <- ma$meta
        if (!is.null(obj)) {   # alpha diversity etc. for the sequenced rows, NA for the others
          extra <- setdiff(colnames(obj$meta_enriched), colnames(meta_all))
          meta_all[extra] <- obj$meta_enriched[match(rownames(meta_all), rownames(obj$meta_enriched)),
                                               extra, drop = FALSE]
        }
        rv$obj <- obj; rv$cfg <- cfg; rv$meta_all <- meta_all; rv$meta_info <- ma$info
        rv$meta_cols     <- if (is.null(obj)) colnames(meta_all) else obj$meta_cols
        rv$meta_enriched <- if (!is.null(obj)) obj$meta_enriched
        stat_store$result <- stat_store$params <- NULL
        maaslin_store$result <- maaslin_store$params <- NULL
        shiny::showNotification(
          if (has_seq) "Data loaded. Choose an analysis on the left."
          else "Metadata loaded. Open Measurements to analyse measured values.", type = "message")
      }, error = function(e) {
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = NULL)
      })
    })

    output$otu_scale_info_sidebar_ui <- shiny::renderUI({
      msg <- .format_scale_info(rv$obj$otu_scale_info); shiny::req(msg)
      .info(shiny::tags$b("Beta diversity input scale: "), msg)
    })

    # Which samples are being analysed (metadata may be a subset of the run)
    output$match_info_ui <- shiny::renderUI({
      if (is.null(rv$meta_all))
        return(.empty(shiny::tags$b("Start here"),
                      "Upload your metadata on the left, add the abundance table if you have ",
                      "sequencing data, then press Load data."))
      seq <- !is.null(rv$obj)
      mi  <- if (seq) rv$obj$match_info else rv$meta_info
      lst <- function(x) paste(c(utils::head(x, 15), if (length(x) > 15) "..."), collapse = ", ")
      shiny::tagList(
        if (!seq)
          .info(shiny::tags$b(sprintf("%d metadata rows loaded", nrow(rv$meta_all))),
                " - no abundance table, so only Measurements is available. Add an abundance ",
                "table and press Load data again for the sequencing analyses.")
        else
          .info(shiny::tags$b(sprintf("%d samples analysed", mi$n_matched)),
                sprintf(" - found in both the abundance table (%d samples) and the metadata (%d rows).",
                        mi$n_table, mi$n_meta),
                if (length(mi$table_only)) shiny::tagList(shiny::tags$br(),
                  sprintf("%d sequenced sample(s) not in your metadata are left out: %s",
                          length(mi$table_only), lst(mi$table_only)))),
        if (length(mi$meta_only))
          .warn(sprintf("%d metadata row(s) have no matching sample in the abundance table and are left out of the sequencing analyses: %s. Check the Sample IDs are spelled the same in both files. Measurements still use these rows.",
                        length(mi$meta_only), lst(mi$meta_only))),
        if (length(mi$subject_candidates))
          .warn("No subject / animal ID is set in Step 2, but ",
                shiny::tags$b(paste(mi$subject_candidates, collapse = ", ")),
                if (length(mi$subject_candidates) > 1) " repeat values." else " repeats values.",
                " If one of them is your animal ID (the same animal sampled more than once), ",
                "set it in Step 2 and press Load data again - otherwise repeated samples ",
                "are treated as independent in every test. Ignore this if each row is a ",
                "different animal."),
        if (isTRUE(mi$subject_unique))
          .warn(shiny::tags$b(rv$cfg$subject_id_col), " is set as the subject / animal ID in Step 2, ",
                "but every row has a different value, so no animal is sampled more than once. ",
                "If each row is a different animal, set it to (none) and press Load data again. ",
                "If animals are sampled repeatedly, choose the column that holds the animal's ",
                "name or number (the same value on each of its rows)."),
        if (length(mi[["incomplete"]])) {
          inc <- mi[["incomplete"]]
          who <- paste(sprintf("%s (%d of %d samples)", names(inc), inc,
                               mi$expected_n), collapse = ", ")
          if (isTRUE(mi$incomplete_dropped))
            .info(shiny::tags$b("Left out - missing samples: "), who)
          else
            .warn("Animals with missing samples are still included: ", who,
                  ". To leave them out of every analysis, tick \"Leave out animals with ",
                  "that miss a timepoint\" in Step 2 and press Load data again.")
        },
        if (seq && length(rv$obj$tax_levels) == 1)
          .info("The table's rows are names without a taxonomic lineage (species-only or ",
                "pathway/gene table), so plots, the descriptive summary and MaAsLin3 use the ",
                "rows as uploaded. For phylum-genus summaries upload the full table with ",
                "lineages (e.g. k__Bacteria|p__...|s__...)."))
    })
    output$otu_scale_info_ui <- shiny::renderUI({
      msg <- .format_scale_info(rv$obj$otu_scale_info); shiny::req(msg)
      .info(shiny::tags$b("Input scale detected: "), msg, shiny::tags$br(),
            "The PCoA plot, PERMANOVA/ANOSIM and PERMDISP all use this same distance object.")
    })

    # Non-baseline levels of the treatment column chosen in the Alpha tab
    treat_levels <- shiny::reactive({
      shiny::req(rv$obj)
      tcol <- .nz(input$lmm_treat_col); me <- rv$meta_enriched
      if (!nzchar(tcol) || !tcol %in% colnames(me)) return(character(0))
      x <- me[!.is_baseline(me, rv$obj$config$is_baseline_col), tcol]
      sort(unique(as.character(x[!is.na(x)])))
    })

    output$treat_ref_ui <- shiny::renderUI({
      tr <- treat_levels()
      if (length(tr))
        shiny::selectInput("treat_ref", "Reference (control) group:", choices = tr)
      else
        .warn("Choose the column with the groups to compare above.")
    })

    # -- Settings column of each sequencing page ---------------------------------
    need_seq <- function()
      if (is.null(rv$meta_all))
        .empty(shiny::tags$b("No data yet"), "Load your files on the Data page first.")
      else
        .empty(shiny::tags$b("Needs sequencing data"),
               "Only metadata is loaded. Add an abundance table on the Data page and press ",
               "Load data again to use this analysis. Measurements work without it.")

    # Things every settings column needs.
    ctx <- shiny::reactive({
      shiny::req(rv$obj)
      cols <- rv$meta_cols; cfg <- rv$obj$config
      lv <- intersect(TAXA_LEVELS, rv$obj$tax_levels)
      if (!length(lv)) lv <- TAXA_LEVELS
      list(cols = cols, cfg = cfg, me_cols = colnames(rv$meta_enriched), lv = lv,
           # the subject column, but only when animals really are sampled more than once
           subj = if (nzchar(cfg$subject_id_col) && cfg$subject_id_col %in% cols &&
                      anyDuplicated(stats::na.omit(rv$meta_enriched[[cfg$subject_id_col]])))
             cfg$subject_id_col,
           # "Species" only for name-only tables
           lv_choices = if (length(lv) == 1) stats::setNames(lv, "Rows as uploaded (no taxonomy)") else lv)
    })
    # Is column `col` one of the variables of a typed formula?
    in_formula <- function(col, fml)
      col %in% tryCatch(all.vars(stats::as.formula(paste("~", fml))), error = function(e) col)
    multi <- function(placeholder) list(placeholder = placeholder, plugins = list("remove_button"),
                                        maxOptions = 500)

    output$alpha_side_ui <- shiny::renderUI({
      if (is.null(rv$obj)) return(need_seq())
      x <- ctx()
      shiny::tagList(
        .section("Model"),
        shiny::textInput("alpha_fixed_effects", "Fixed effects (formula, no leading '~'):",
                         value = "", placeholder = "e.g. Treatment * Timepoint"),
        .help("How do I write this?",
              "List the columns joined by +, e.g. Treatment + Timepoint. Use * for an ",
              "interaction, e.g. Treatment * Timepoint. Column names with spaces need ",
              "backticks. Do not include Baseline_<metric>: each metric's own baseline value ",
              "is added automatically when available."),
        shiny::selectizeInput("lmm_random_effect", "Random effect(s) (animal / subject ID):",
                              choices = x$me_cols, selected = x$subj, multiple = TRUE,
                              options = multi("Leave empty if each animal has one sample")),
        .help("When do I need one?",
              "A random effect tells the model which samples come from the same animal. ",
              "Several samples per animal (repeated timepoints, crossover): add the column ",
              "with the animal ID. Each animal has only one sample: leave the box empty and ",
              "an ordinary linear model is fitted. More than one column gives one random ",
              "intercept each, e.g. (1|Subject) + (1|Litter)."),
        shiny::selectInput("lmm_treat_col", "Groups to compare (column):", choices = .none(x$me_cols),
                           selected = .guess_col(x$me_cols, .TREAT_NAMES)),
        .help("What is this for?",
              "The column whose groups you want to compare, e.g. Treatment or Group. It must ",
              "also be in the formula above. Every other group is compared with the ",
              "reference group chosen below, and the p-values adjusted across metrics are ",
              "those of this column. Other columns in the formula are adjusted for, not tested."),
        shiny::uiOutput("treat_ref_ui"),
        shiny::selectInput("alpha_padj", "P-value adjustment:", choices = PADJ_METHODS,
                           selected = x$cfg$padj_method),
        .run_btn("run_alpha_lmm", "Run alpha diversity model"))
    })

    output$pcoa_side_ui <- shiny::renderUI({
      if (is.null(rv$obj)) return(need_seq())
      x <- ctx(); cols <- x$cols
      shiny::tagList(
        .section("Colours"),
        shiny::selectInput("color_var", "Colour by:", choices = cols,
                           selected = .pick(.guess_col(cols, .TREAT_NAMES), cols)),
        shiny::checkboxInput("is_cont", "This column is a number (colour as a gradient)", FALSE),
        shiny::conditionalPanel("input.is_cont",
          shiny::selectInput("cont_palette", "Gradient:",
                             choices = c(CONT_PALETTES, "My own colours" = "custom")),
          shiny::conditionalPanel("input.cont_palette == 'custom'",
            shiny::tags$div(class = "mi-cols", style = "margin-bottom:6px",
              .grad_input("grad_low", "Low", "#2F6FB0"), .grad_input("grad_mid", "Middle", "#F2F2F2"),
              .grad_input("grad_high", "High", "#B8456B")),
            shiny::checkboxInput("grad_use_mid", "Use the middle colour", TRUE))),
        shiny::uiOutput("pcoa_cols_ui"),
        if (!is.null(x$subj)) shiny::tagList(
          .section("Paired samples"),
          shiny::checkboxInput("connect_pairs",
                               sprintf("Join samples of the same %s with a line", x$subj), FALSE),
          shiny::conditionalPanel("input.connect_pairs",
            shiny::selectInput("pair_order", "Order along the line (e.g. timepoint):",
                               choices = .none(cols, "Sample ID"), selected = .NONE))),
        .section("2D axes"),
        shiny::fluidRow(
          shiny::column(6, shiny::selectInput("pc_x", "X axis:", choices = pcs, selected = 1)),
          shiny::column(6, shiny::selectInput("pc_y", "Y axis:", choices = pcs, selected = 2))),
        .section("Download"),
        shiny::downloadButton("dl_2d", "2D (TIFF)"), shiny::downloadButton("dl_3d", "3D (TIFF)"))
    })

    output$beta_side_ui <- shiny::renderUI({
      if (is.null(rv$obj)) return(need_seq())
      x <- ctx(); cols <- x$cols
      shiny::tagList(
        .section("Test"),
        shiny::selectInput("stat_var", "Primary grouping variable:", choices = cols,
                           selected = .pick(.guess_col(cols, .TREAT_NAMES), cols)),
        shiny::selectInput("stat_method", "Test:",
                           choices = c("PERMANOVA" = "permanova", "ANOSIM" = "anosim")),
        shiny::radioButtons("stat_scope", "Scope:", inline = TRUE,
                            choices = c("Global" = "global", "Pairwise" = "pairwise")),
        shiny::conditionalPanel("input.stat_method == 'permanova'",
          shiny::selectizeInput("covariates", "Covariates (entered before the primary variable):",
                                choices = cols, selected = NULL, multiple = TRUE,
                                options = multi("Click to add covariates...")),
          .help("Which covariates?",
                "For crossover designs include the period (Timepoint) and Sequence columns, ",
                "and optionally a Baseline_* covariate."),
          shiny::textInput("beta_formula", "Or type a full model formula:", value = "",
                           placeholder = "e.g. Treatment * Timepoint"),
          .help("When to type a formula?",
                "Leave blank to use the dropdowns. A typed formula overrides them and allows ",
                "interactions; every model term is then reported. Global PERMANOVA only.")),
        shiny::conditionalPanel("input.stat_method == 'anosim'",
          .warn("ANOSIM does not support covariates. Use PERMANOVA for a covariate-adjusted analysis.")),
        .section("Permutations"),
        shiny::selectInput("strata_var", "Shuffle only within (strata):", choices = .none(cols),
                           selected = if (is.null(x$subj)) .NONE else x$subj),
        .help("What does this do?",
              "For repeated measures (several samples per animal) choose the subject ID, so ",
              "samples are only shuffled within each animal. Leave at (none) when every animal ",
              "has one sample. Terms that never change within an animal (e.g. Sequence, Sex) ",
              "cannot be tested this way and show no p-value."),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput("stat_perm", "Permutations:", value = 999,
                                               min = 99, max = 99999, step = 100)),
          shiny::column(6, shiny::numericInput("stat_seed", "Random seed:", value = 123,
                                               min = 1, step = 1))),
        .help("Why a seed?",
              "p-values come from random shuffles, so they vary slightly between runs unless ",
              "the seed is fixed (same seed = same p-values). Use 9999 permutations for ",
              "results near 0.05. When few shuffles are possible (small paired studies) the ",
              "exact test is used automatically."),
        shiny::selectInput("padj_method", "P-value adjustment:", choices = PADJ_METHODS,
                           selected = x$cfg$padj_method),
        .run_btn("run_stat", "Run beta diversity test"))
    })

    output$taxon_side_ui <- shiny::renderUI({
      if (is.null(rv$obj)) return(need_seq())
      x <- ctx(); cols <- x$cols; lv <- x$lv
      shiny::tagList(
        .section("Bar plot"),
        shiny::selectInput("taxon_level", "Taxonomic level:", choices = x$lv_choices,
                           selected = .pick("Family", lv)),
        shiny::selectInput("taxon_group_var", "Panel per:", choices = cols,
                           selected = .pick(.guess_col(cols, .TREAT_NAMES), cols)),
        shiny::numericInput("taxon_top_n", "Number of taxa shown (the rest become Others):",
                            value = TAXA_DEFAULT_TOP_N, min = 3, max = 30, step = 1),
        shiny::checkboxInput("taxon_exclude_baseline", "Exclude baseline samples", TRUE),
        shiny::selectInput("taxon_x_label_var", "Label bars with:",
                           choices = c("Sample ID" = "SampleID", cols), selected = "SampleID"),
        shiny::selectInput("taxon_x_order_var", "Order bars by:",
                           choices = c("Sample ID (alphabetical)" = "SampleID", cols),
                           selected = "SampleID"),
        .run_btn("run_taxon", "Build taxonomy plot"),
        shiny::tags$p(class = "mi-hint", style = "margin-top:6px",
                      "Press again after changing the level or the number of taxa. Labels and order update on their own."),
        shiny::uiOutput("taxon_cols_ui"),
        .section("Summary per group (Excel)"),
        shiny::selectizeInput("desc_group_cols", "Group by (one or more columns):", choices = cols,
                              selected = NULL, multiple = TRUE, options = multi("Click to add...")),
        shiny::selectizeInput("desc_sort_cols", "Order groups by (optional):", choices = cols,
                              selected = NULL, multiple = TRUE,
                              options = multi("Blank = same order as grouping columns")),
        shiny::checkboxGroupInput("desc_levels", "Taxonomic levels:", choices = x$lv_choices,
                                  selected = if (length(lv) == 1) lv else setdiff(lv, "Domain"),
                                  inline = TRUE),
        shiny::checkboxInput("desc_exclude_baseline", "Exclude baseline samples", FALSE),
        shiny::tags$p(class = "mi-hint", style = "margin-top:6px",
                      "One sheet: each level with median, min, max, range, mean, standard deviation and 95% CI of the mean of relative abundance (%) per group."),
        shiny::downloadButton("dl_desc", "Download summary (.xlsx)"))
    })

    output$maaslin_side_ui <- shiny::renderUI({
      if (is.null(rv$obj)) return(need_seq())
      x <- ctx(); me_cols <- x$me_cols
      shiny::tagList(
        .section("Model"),
        shiny::selectInput("maaslin_level", "Taxonomic level:", choices = x$lv_choices,
                           selected = .pick("Genus", x$lv)),
        shiny::textInput("maaslin_formula", "Fixed effects (formula, no leading '~'):",
                         value = "", placeholder = "e.g. Treatment + Timepoint"),
        .help("How do I write this?",
              "List the columns joined by +. Use * for an interaction, e.g. Treatment * ",
              "Timepoint. Column names with spaces need backticks. No baseline covariate is ",
              "added automatically: type e.g. Baseline_Shannon into the formula if you want one."),
        shiny::uiOutput("maaslin_ref_ui"),
        shiny::checkboxInput("maaslin_use_random", "Include random effect",
                             value = !is.null(x$subj)),
        shiny::conditionalPanel("input.maaslin_use_random",
          shiny::selectizeInput("maaslin_random_effect", "Random effect(s) (animal / subject ID):",
                                choices = me_cols, selected = x$subj, multiple = TRUE,
                                options = multi("Click to add one or more...")),
          # only shown when the installed MaAsLin3 cannot take a second one
          if (requireNamespace("maaslin3", quietly = TRUE) && !.maaslin_multi_ranef_ok())
            shiny::conditionalPanel(
              "input.maaslin_random_effect && input.maaslin_random_effect.length > 1",
              .warn(.maaslin_multi_ranef_msg()))),
        .more("More settings",
          shiny::selectInput("maaslin_group_var", "Grouping factor (optional, ANOVA-style group() term):",
                             choices = .none(me_cols), selected = .NONE),
          shiny::checkboxInput("maaslin_small_random", "Random effect group sizes are small", FALSE),
          shiny::tags$p(class = "mi-hint",
                        "When ticked, MaAsLin3 replaces the random effect with a fixed effect in the prevalence model for groups with few observations."),
          shiny::selectInput("maaslin_normalization", "Normalization:", choices = c("TSS", "CLR", "NONE")),
          shiny::selectInput("maaslin_transform", "Transform:", choices = c("LOG", "PLOG", "NONE")),
          shiny::numericInput("maaslin_min_prevalence", "Min prevalence (fraction of samples):",
                              value = 0.1, min = 0, max = 1, step = 0.05),
          shiny::numericInput("maaslin_min_abundance", "Min abundance (% relative abundance):",
                              value = 0, min = 0, step = 0.01),
          shiny::selectInput("maaslin_padj_method", "P-value adjustment:", choices = PADJ_METHODS,
                             selected = x$cfg$padj_method),
          shiny::numericInput("maaslin_max_sig", "Max significance (q-value) reported:",
                              value = 0.25, min = 0, max = 1, step = 0.05)),
        .run_btn("run_maaslin_btn", "Run MaAsLin3"))
    })

    # Colours: one palette for every plot; PCoA can switch to a gradient.
    # Colours the user fixed for single groups (group -> colour) lie on top.
    cc <- shiny::reactiveValues(map = character(0), ver = 0)
    shiny::observeEvent(input$set_colour, {
      x <- input$set_colour
      if (is.list(x) && length(x$level) == 1 && isTRUE(.is_colour(x$value))) {
        m <- cc$map; m[[as.character(x$level)]] <- x$value; cc$map <- m
      }
    })
    shiny::observeEvent(input$reset_colours, { cc$map <- character(0); cc$ver <- cc$ver + 1 })
    base_pal <- shiny::reactive(if (is.null(input$palette)) "Miome" else input$palette)
    cat_pal  <- shiny::reactive(if (length(cc$map)) structure(cc$map, base = base_pal()) else base_pal())
    grad     <- function(id, default) if (isTRUE(.is_colour(input[[id]]))) input[[id]] else default
    pcoa_pal <- shiny::reactive(
      if (!isTRUE(input$is_cont)) cat_pal()
      else if (is.null(input$cont_palette)) CONT_PALETTES[1]
      else if (identical(input$cont_palette, "custom"))
        c(grad("grad_low", "#2F6FB0"), if (isTRUE(input$grad_use_mid)) grad("grad_mid", "#F2F2F2"),
          grad("grad_high", "#B8456B"))
      else input$cont_palette)
    # Pickers for the groups of one column. Redrawn when the column, the palette
    # or "back to palette colours" changes - not on every pick, which would fold
    # the list shut while it is being used.
    editor <- function(values, what = "group") {
      cc$ver; base_pal()
      lv <- sort(unique(as.character(values[!is.na(values)])))
      .colour_editor(.level_colours(shiny::isolate(cat_pal()), lv, isTRUE(input$rev_pal)), what)
    }
    output$pcoa_cols_ui <- shiny::renderUI({
      shiny::req(rv$obj, input$color_var %in% colnames(rv$obj$meta))
      if (!isTRUE(input$is_cont)) editor(rv$obj$meta[[input$color_var]])
    })
    output$alpha_cols_ui <- shiny::renderUI({
      shiny::req(rv$meta_enriched, input$alpha_color_var %in% colnames(rv$meta_enriched))
      editor(rv$meta_enriched[[input$alpha_color_var]])
    })
    output$meas_cols_ui <- shiny::renderUI({
      m <- meas(); shiny::req(m, input$meas_plot_color %in% colnames(m$data))
      editor(m$data[[input$meas_plot_color]])
    })
    output$taxon_cols_ui <- shiny::renderUI({
      td <- taxon_data_r(); shiny::req(td)
      editor(setdiff(td$Taxon, "Others"), "taxon")
    })

    output$nav_status <- shiny::renderUI({
      shiny::req(rv$meta_all)
      shiny::tags$div(class = "mi-status",
        shiny::tags$b(if (is.null(rv$obj)) sprintf("%d metadata rows", nrow(rv$meta_all))
                      else sprintf("%d samples", nrow(rv$meta_enriched))),
        if (is.null(rv$obj)) "no sequencing data" else rv$obj$config$dist_method_label)
    })

    # ======================= Targeted measurements =================================
    meas <- shiny::reactive({
      shiny::req(input$meas_file, rv$meta_all)
      tryCatch(link_measurements(read_measurements(input$meas_file$datapath), rv$meta_all),
               error = function(e) {
                 shiny::showNotification(paste("Measurements:", conditionMessage(e)),
                                         type = "error", duration = NULL)
                 NULL
               })
    })
    meas_choices <- shiny::reactive({
      m <- meas(); shiny::req(m)
      stats::setNames(m$features, ifelse(is.na(m$labels), m$features, m$labels))
    })

    output$meas_link_ui <- shiny::renderUI({
      if (is.null(rv$meta_all))
        return(.empty(shiny::tags$b("Metadata first"),
                      "Load your metadata on the Data page. An abundance table is not needed."))
      if (is.null(input$meas_file))
        return(.empty(shiny::tags$b("No measurements yet"),
                      "Upload a measurement table on the left to start."))
      m <- meas(); shiny::req(m); l <- m$link
      lst <- function(x) paste(c(utils::head(x, 12), if (length(x) > 12) "..."), collapse = ", ")
      shiny::tagList(
        .info(shiny::tags$b(sprintf("%d samples with measurements", l$n_linked)),
              sprintf(" - %d measurement columns, matched to the metadata through '%s'.",
                      length(m$features), l$by),
              if (length(l$unlinked_rows)) shiny::tagList(shiny::tags$br(),
                sprintf("%d row(s) of the file are not in your metadata and are left out: %s",
                        length(l$unlinked_rows), lst(l$unlinked_rows))),
              if (length(l$samples_without)) shiny::tagList(shiny::tags$br(),
                sprintf("%d metadata row(s) have no measurement: %s",
                        length(l$samples_without), lst(l$samples_without)))),
        if (length(l$clash))
          .info(sprintf("%d column(s) are in both files (%s); the values from the measurement file are used here.",
                        length(l$clash), lst(l$clash))))
    })

    output$meas_controls_ui <- shiny::renderUI({
      m <- meas(); shiny::req(m)
      cfg <- rv$cfg; cols <- setdiff(colnames(m$data), m$features)
      has_base <- any(.is_baseline(m$data, cfg$is_baseline_col)) && nzchar(cfg$subject_id_col)
      shiny::tagList(
        .section("Measurements to analyse"),
        shiny::selectizeInput("meas_features", NULL,
                              choices = meas_choices(), selected = NULL, multiple = TRUE,
                              options = list(placeholder = "Click to add, or use Select all",
                                             plugins = list("remove_button"))),
        shiny::tags$div(style = "margin:-4px 0 12px",
          shiny::actionButton("meas_all", "Select all", class = "btn-xs"),
          shiny::actionButton("meas_none", "Clear", class = "btn-xs")),
        .help("Which should I select?",
              "p-values are adjusted across the measurements you select here, so choose the ",
              "ones you planned to test before looking at results. Leave out totals or ratios ",
              "calculated from other columns unless you intend to report them."),
        .section("Model"),
        shiny::selectInput("meas_transform", "Transform before modelling:",
                           choices = c("None (use values as uploaded)" = "none",
                                       "log10 (concentrations)" = "log10")),
        .help("When to use log10?",
              "Use log10 for concentrations (fatty acids, sterols, bile acids), which are ",
              "usually right-skewed. Keep 'None' for values already on a log scale (qPCR Log DNA) ",
              "or indices that can be negative."),
        shiny::textInput("meas_formula", "Fixed effects (formula, no leading '~'):",
                         value = "", placeholder = "e.g. Group + Sex"),
        shiny::selectizeInput("meas_random", "Random effect(s) (animal ID):", choices = cols,
                              selected = if (nzchar(cfg$subject_id_col) && cfg$subject_id_col %in% cols &&
                                             anyDuplicated(stats::na.omit(m$data[[cfg$subject_id_col]])))
                                cfg$subject_id_col,
                              multiple = TRUE,
                              options = list(placeholder = "Leave empty if each animal has one sample",
                                             plugins = list("remove_button"))),
        shiny::selectInput("meas_treat_col", "Groups to compare (column):", choices = .none(cols),
                           selected = .guess_col(cols, .TREAT_NAMES)),
        .help("What is this for?",
              "The column whose groups you want to compare, e.g. Group or Treatment. It must ",
              "also be in the formula above. Every other group is compared with the ",
              "reference group chosen below, and the p-values adjusted across measurements ",
              "are those of this column. Other columns in the formula (e.g. Sex) are adjusted ",
              "for, not tested."),
        shiny::uiOutput("meas_ref_ui"),
        if (has_base) shiny::checkboxInput("meas_baseline_adjust",
                                           "Adjust for each animal's baseline value", TRUE),
        shiny::selectInput("meas_padj", "P-value adjustment:", choices = PADJ_METHODS,
                           selected = cfg$padj_method),
        .run_btn("run_meas", "Run measurement models"),
        .section("Descriptive statistics"),
        shiny::selectizeInput("meas_group_cols", "Group by (one or more columns):", choices = cols,
                              selected = NULL, multiple = TRUE,
                              options = list(placeholder = "e.g. Group, or Treatment + Timepoint",
                                             plugins = list("remove_button"))),
        shiny::checkboxInput("meas_desc_excl_base", "Exclude baseline samples", FALSE))
    })
    shiny::observeEvent(input$meas_all, shiny::updateSelectizeInput(
      session, "meas_features", choices = meas_choices(), selected = unname(meas_choices())))
    shiny::observeEvent(input$meas_none, shiny::updateSelectizeInput(
      session, "meas_features", choices = meas_choices(), selected = character(0)))

    meas_levels <- shiny::reactive({
      m <- meas(); shiny::req(m); tcol <- .nz(input$meas_treat_col)
      if (!nzchar(tcol) || !tcol %in% colnames(m$data)) return(character(0))
      x <- m$data[!.is_baseline(m$data, rv$cfg$is_baseline_col), tcol]
      sort(unique(as.character(x[!is.na(x)])))
    })
    output$meas_ref_ui <- shiny::renderUI({
      lv <- meas_levels()
      if (length(lv)) shiny::selectInput("meas_treat_ref", "Reference (control) group:", choices = lv)
      else .warn("Choose the column with the groups to compare above.")
    })

    meas_model <- shiny::eventReactive(input$run_meas, {
      m <- meas(); shiny::req(m); cfg <- rv$cfg
      tcol <- .nz(input$meas_treat_col)
      stop_msg <- if (!length(input$meas_features)) "Select at least one measurement first."
        else if (!nzchar(tcol) || is.null(input$meas_treat_ref)) "Choose the column with the groups to compare and its reference group first."
        else if (!nzchar(trimws(.nz(input$meas_formula)))) "Type the fixed effects first, e.g. Group + Sex."
        else if (!in_formula(tcol, input$meas_formula)) sprintf("The formula does not contain '%s', the column chosen under 'Groups to compare'. Add it to the formula or choose another column.", tcol)
      if (!is.null(stop_msg)) { shiny::showNotification(stop_msg, type = "error"); return(NULL) }
      md <- m$data[!.is_baseline(m$data, cfg$is_baseline_col), , drop = FALSE]
      rep_cols <- if (nzchar(cfg$subject_id_col) && cfg$subject_id_col %in% colnames(md) &&
                      anyDuplicated(md[[cfg$subject_id_col]])) cfg$subject_id_col
                  else .subject_candidates(md)
      if (!length(input$meas_random) && length(rep_cols))
        shiny::showNotification(sprintf(paste0(
          "'%s' has repeated values, so some animals appear more than once. Add it as a ",
          "random effect, otherwise repeated samples are treated as independent."),
          rep_cols[1]), type = "warning", duration = 15)
      tryCatch(
        shiny::withProgress(message = "Fitting models...", value = 0.4,
          suppressWarnings(run_measurement_models(
            m$data, input$meas_features, input$meas_formula,
            random_effect = as.character(input$meas_random),
            treatment_ref = input$meas_treat_ref, treatment_col = tcol,
            is_baseline_col = cfg$is_baseline_col, subject_id_col = cfg$subject_id_col,
            baseline_adjust = !identical(input$meas_baseline_adjust, FALSE),
            transform = input$meas_transform, padj_method = input$meas_padj,
            fdr_target = "auto", labels = m$labels))),
        error = function(e) {
          shiny::showNotification(paste("Measurement models failed:", conditionMessage(e)),
                                  type = "error", duration = NULL)
          NULL
        })
    })

    output$meas_model_info_ui <- shiny::renderUI({
      res <- meas_model(); shiny::req(res)
      terms <- attr(res, "fdr_terms"); n_f <- attr(res, "n_features"); mth <- res$Transform[1]
      n_err <- sum(res$Term == "ERROR"); skipped <- attr(res, "skipped")
      shiny::tagList(
        .info(shiny::HTML(sprintf(paste0(
          "Model per measurement: %s. Transform: %s.<br>",
          "<b>p_value</b>: raw. <b>p_adjusted</b>: across the terms of one measurement's model.<br>",
          "<b>p_adj_across_features</b>: one comparison (e.g. %s) adjusted across the <b>%d</b> ",
          "selected measurements.<br>",
          "<b>p_adj_across_features_and_groups</b>: all treatment comparisons of all %d ",
          "measurements as one set.<br>%s"),
          res$Formula[res$Term != "ERROR"][1],
          if (mth == "log10") "log10 (estimates are differences on the log10 scale; zeros replaced by half the smallest positive value)" else "none",
          if (length(terms)) terms[1] else "treatment", n_f, n_f,
          if (length(terms) > 1) sprintf("Your treatment has %d comparisons (%s), so the two cross-measurement columns differ. Report whichever you decided on beforehand.",
                                         length(terms), paste(terms, collapse = ", "))
          else "With two groups there is one comparison, so the two cross-measurement columns are equal."))),
        if (n_err || length(skipped))
          .warn(sprintf("%d measurement(s) could not be modelled%s. See rows with Term = ERROR.",
                        n_err + length(skipped),
                        if (length(skipped)) paste0(" (log10 impossible for: ", paste(skipped, collapse = ", "), ")") else "")))
    })
    output$meas_model_table <- DT::renderDT({
      tbl <- meas_model(); shiny::req(tbl)
      .dt_table(tbl[, setdiff(names(tbl), "Feature")],
                c("Estimate", "SE", "df", "t_value", "p_value", "p_adjusted",
                  "p_adj_across_features", "p_adj_across_features_and_groups"),
                list(p_value = c(0.05, 0.1), p_adj_across_features = c(0.05, 0.25),
                     p_adj_across_features_and_groups = c(0.05, 0.25)),
                page = 20, filter = "top")
    })
    output$dl_meas_model <- shiny::downloadHandler(
      filename = function() sprintf("measurement_models_%s.csv", .today()),
      content  = function(file) utils::write.csv(meas_model(), file, row.names = FALSE))

    meas_desc <- shiny::reactive({
      m <- meas(); shiny::req(m, length(input$meas_features) > 0, length(input$meas_group_cols) > 0)
      summarize_measurements(m$data, input$meas_features, input$meas_group_cols,
                             exclude_baseline = isTRUE(input$meas_desc_excl_base),
                             is_baseline_col = rv$cfg$is_baseline_col, labels = m$labels)
    })
    output$meas_desc_table <- DT::renderDT({
      d <- meas_desc()
      .dt_table(d[, setdiff(names(d), "Feature")], page = 15, filter = "top")
    })
    output$dl_meas_desc <- shiny::downloadHandler(
      filename = function() sprintf("measurement_descriptives_by_%s_%s.csv",
                                    .fname(input$meas_group_cols), .today()),
      content  = function(file) utils::write.csv(meas_desc(), file, row.names = FALSE))
    output$dl_meas_desc_xlsx <- shiny::downloadHandler(
      filename = function() sprintf("measurement_descriptives_by_%s_%s.xlsx",
                                    .fname(input$meas_group_cols), .today()),
      content  = function(file) {
        m <- meas(); shiny::req(m, length(input$meas_features) > 0, length(input$meas_group_cols) > 0)
        export_measurements_descriptive_workbook(
          m$data, input$meas_features, input$meas_group_cols, outfile = file,
          exclude_baseline = isTRUE(input$meas_desc_excl_base),
          is_baseline_col = rv$cfg$is_baseline_col, labels = m$labels)
      })

    output$meas_plot_controls_ui <- shiny::renderUI({
      m <- meas(); shiny::req(m); cols <- setdiff(colnames(m$data), m$features)
      subj <- rv$cfg$subject_id_col
      shiny::tagList(
        shiny::fluidRow(
          shiny::column(4, shiny::selectInput("meas_plot_feature", "Measurement:", choices = meas_choices())),
          shiny::column(4, shiny::selectInput("meas_plot_x", "X axis grouping:", choices = cols,
                                              selected = .pick(.guess_col(cols, .TREAT_NAMES), cols))),
          shiny::column(4, shiny::selectInput("meas_plot_color", "Colour by:", choices = cols,
                                              selected = .pick(.guess_col(cols, .TREAT_NAMES), cols)))),
        if (nzchar(subj) && subj %in% cols)
          shiny::checkboxInput("meas_connect", sprintf("Join samples of the same %s with a line", subj), FALSE),
        shiny::uiOutput("meas_cols_ui"))
    })
    output$meas_plot <- shiny::renderPlot({
      m <- meas(); shiny::req(m, input$meas_plot_feature, input$meas_plot_x, input$meas_plot_color)
      f   <- input$meas_plot_feature
      lab <- if (is.na(m$labels[f])) f else m$labels[[f]]
      build_alpha_scatter(m$data, f, input$meas_plot_x, input$meas_plot_color,
                          palette = cat_pal(), rev_pal = isTRUE(input$rev_pal),
                          connect_by = if (isTRUE(input$meas_connect)) rv$cfg$subject_id_col,
                          exclude_baseline = FALSE, is_baseline_col = rv$cfg$is_baseline_col,
                          title = lab, y_label = lab)
    })

    shiny::observeEvent(input$pc_x, {
      shiny::req(input$pc_y)
      if (input$pc_x == input$pc_y)
        shiny::updateSelectInput(session, "pc_y", selected = setdiff(c("1", "2", "3"), input$pc_x)[1])
    })

    # pcoa_obj whose timepoint column = the "Order lines by" choice
    plot_obj <- shiny::reactive({
      shiny::req(rv$obj)
      o <- rv$obj
      o$config$timepoint_col <- if (isTRUE(input$connect_pairs)) .nz(input$pair_order) else ""
      o
    })

    connect_col <- shiny::reactive({
      shiny::req(rv$obj)
      s <- rv$obj$config$subject_id_col
      if (isTRUE(input$connect_pairs) && nzchar(s)) s else NULL
    })

    # -- PCoA plots ---------------------------------------------------------------
    output$plot3d <- plotly::renderPlotly({
      shiny::req(rv$obj, input$color_var)
      build_3d(plot_obj(), input$color_var, pcoa_pal(), isTRUE(input$is_cont),
               isTRUE(input$rev_pal), connect_col())
    })

    output$plot2d <- shiny::renderPlot({
      shiny::req(rv$obj, input$color_var, input$pc_x, input$pc_y)
      shiny::validate(shiny::need(input$pc_x != input$pc_y, "Please choose two different PC axes."))
      build_2d(plot_obj(), input$color_var, as.integer(input$pc_x), as.integer(input$pc_y),
               pcoa_pal(), isTRUE(input$is_cont), isTRUE(input$rev_pal), connect_col())
    })

    # -- Beta diversity statistics (params frozen at click) ----------------------
    shiny::observeEvent(input$run_stat, {
      shiny::req(rv$obj, input$stat_var)
      obj <- rv$obj
      obj$config$padj_method <- input$padj_method
      perm <- input$covariates
      covs <- if (input$stat_method == "permanova" && length(perm)) perm else NULL
      strata <- if (nzchar(.nz(input$strata_var))) input$strata_var else NULL
      fml <- trimws(.nz(input$beta_formula))
      fml <- if (input$stat_method == "permanova" && input$stat_scope == "global" && nzchar(fml)) fml

      seed <- if (is.null(input$stat_seed) || is.na(input$stat_seed)) 123 else input$stat_seed
      tryCatch({
        set.seed(seed)
        res <- shiny::withProgress(message = "Running permutation test...", value = 0.5, {
          if (input$stat_scope == "global")
            run_global_test(obj, input$stat_var, input$stat_method, covs, strata,
                            input$stat_perm, formula_rhs = fml)
          else
            run_pairwise_test(obj, input$stat_var, input$stat_method, covs, strata,
                              input$stat_perm)
        })
        stat_store$result <- res
        stat_store$params <- list(method = input$stat_method, scope = input$stat_scope,
                                  var = input$stat_var,
                                  covariates = if (is.null(fml)) covs,
                                  formula = fml, strata = strata,
                                  padj_method = input$padj_method, perms = input$stat_perm,
                                  seed = seed)
      }, error = function(e)
        shiny::showNotification(paste("Stat test failed:", conditionMessage(e)), type = "error"))
    })

    output$stat_info <- shiny::renderText({
      p <- stat_store$params; shiny::req(p)
      model <- if (!is.null(p$formula)) paste("formula:", p$formula)
               else paste(c(p$covariates, p$var), collapse = " + ")
      sprintf("Test      : %s\nScope     : %s\nModel     : %s\nStrata    : %s\np-adjust  : %s\nPerms     : %s (requested)\nSeed      : %s",
              toupper(p$method), p$scope, model, if (is.null(p$strata)) "none" else p$strata,
              p$padj_method, p$perms, p$seed)
    })

    output$stat_table <- DT::renderDT({
      shiny::req(stat_store$result)
      .dt_table(stat_store$result,
                c("Statistic_R2", "Statistic_R", "Statistic", "F_value", "p_value", "p_adjusted"),
                list(p_value = c(0.05, 0.1)))
    })

    output$permdisp_out <- shiny::renderPrint({
      p <- stat_store$params
      shiny::req(stat_store$result, p, rv$obj, p$method == "permanova")
      tryCatch({
        extra <- if (!is.null(p$formula)) .parse_fixed_formula(p$formula, rv$obj$meta)$vars_used
        md  <- .filter_meta(rv$obj$meta, unique(c(p$var, p$covariates, p$strata, extra)),
                            rv$obj$config$is_baseline_col, TRUE)
        grp <- factor(md[[p$var]])
        if (nlevels(grp) < 2) return(cat("Not enough groups for PERMDISP.\n"))
        set.seed(p$seed)
        bd <- vegan::betadisper(.sub_dist(rv$obj, rownames(md)), grp)
        cat("-- PERMDISP: Homogeneity of Multivariate Dispersion --\n")
        cat(sprintf("Grouping variable : %s\nSamples           : %d\nPermutations      : %s\n\n",
                    p$var, nrow(md), p$perms))
        print(vegan::permutest(bd, permutations = p$perms))
        cat("\nGroup dispersions (average distance to centroid):\n")
        print(round(tapply(bd$distances, bd$group, mean), 4))
      }, error = function(e) cat("PERMDISP could not be computed:\n", conditionMessage(e), "\n"))
    })

    # -- Alpha diversity -----------------------------------------------------------
    output$alpha_controls_ui <- shiny::renderUI({
      shiny::req(rv$meta_enriched, rv$obj)
      cfg <- rv$obj$config; cols <- rv$meta_cols
      metrics <- intersect(ALPHA_METRICS, colnames(rv$meta_enriched))
      shiny::tagList(
        shiny::fluidRow(
          shiny::column(4, shiny::selectInput("alpha_metric", "Alpha diversity metric:",
                                              choices = metrics, selected = .pick("Shannon", metrics))),
          shiny::column(4, shiny::selectInput("alpha_x_var", "X axis grouping:", choices = cols,
                                              selected = .pick(.guess_col(cols, .TREAT_NAMES), cols))),
          shiny::column(4, shiny::selectInput("alpha_color_var", "Colour by:", choices = cols,
                                              selected = .pick(.guess_col(cols, .TIME_NAMES), cols)))),
        shiny::checkboxInput("alpha_exclude_baseline", "Exclude baseline samples from plot", TRUE),
        if (nzchar(cfg$subject_id_col) && cfg$subject_id_col %in% cols)
          shiny::checkboxInput("alpha_connect",
                               sprintf("Join samples of the same %s with a line", cfg$subject_id_col),
                               isTRUE(shiny::isolate(input$connect_pairs))),
        shiny::uiOutput("alpha_cols_ui"),
        shiny::downloadButton("dl_alpha", "Download plot (TIFF)"))
    })

    alpha_plot_args <- shiny::reactive({
      shiny::req(rv$meta_enriched, rv$obj, input$alpha_metric, input$alpha_x_var,
                 input$alpha_color_var)
      list(meta_enriched = rv$meta_enriched, metric = input$alpha_metric,
           x_var = input$alpha_x_var, color_var = input$alpha_color_var,
           palette = cat_pal(), rev_pal = isTRUE(input$rev_pal),
           connect_by = if (isTRUE(input$alpha_connect) && nzchar(rv$obj$config$subject_id_col))
             rv$obj$config$subject_id_col,
           exclude_baseline = isTRUE(input$alpha_exclude_baseline),
           is_baseline_col = rv$obj$config$is_baseline_col)
    })

    output$alpha_plot <- shiny::renderPlot(do.call(build_alpha_scatter, alpha_plot_args()))


    alpha_lmm_result <- shiny::eventReactive(input$run_alpha_lmm, {
      shiny::req(rv$meta_enriched, rv$obj)
      cfg <- rv$obj$config
      tcol <- .nz(input$lmm_treat_col)
      if (!nzchar(tcol) || is.null(input$treat_ref)) {
        shiny::showNotification("Choose the treatment column and reference level first.",
                                type = "error")
        return(NULL)
      }
      # Warn when samples look paired but the model treats them as independent.
      subj <- cfg$subject_id_col
      md   <- rv$meta_enriched[!.is_baseline(rv$meta_enriched, cfg$is_baseline_col), , drop = FALSE]
      rep_cols <- if (nzchar(subj) && subj %in% colnames(md) && anyDuplicated(md[[subj]])) subj
                  else .subject_candidates(md)
      if (!length(input$lmm_random_effect) && length(rep_cols))
        shiny::showNotification(sprintf(paste0(
          "'%s' has repeated values, so some animals appear more than once. Add it as a ",
          "random effect, otherwise repeated samples are treated as independent."),
          rep_cols[1]), type = "warning", duration = 15)
      if (!nzchar(trimws(.nz(input$alpha_fixed_effects)))) {
        shiny::showNotification("Type the fixed effects first, e.g. Treatment * Timepoint.",
                                type = "error")
        return(NULL)
      }
      if (!in_formula(tcol, input$alpha_fixed_effects)) {
        shiny::showNotification(sprintf("The formula does not contain '%s', the column chosen under 'Groups to compare'. Add it to the formula or choose another column.", tcol),
                                type = "error", duration = 12)
        return(NULL)
      }
      tryCatch(
        shiny::withProgress(message = "Fitting models...", value = 0.3,
          run_alpha_lmm_all_metrics(
            rv$meta_enriched, input$alpha_fixed_effects, as.character(input$lmm_random_effect),
            input$treat_ref, tcol, cfg$is_baseline_col,
            padj_method = if (nzchar(.nz(input$alpha_padj))) input$alpha_padj else "BH",
            fdr_target = "auto")),
        error = function(e) {
          shiny::showNotification(paste("Alpha model failed:", conditionMessage(e)), type = "error")
          NULL
        })
    })

    output$alpha_lmm_info_ui <- shiny::renderUI({
      shiny::req(rv$obj, input$alpha_fixed_effects)
      re <- input$lmm_random_effect
      has_base <- any(grepl("^Baseline_", colnames(rv$meta_enriched)))
      .info(shiny::HTML(sprintf(paste0(
        "Model: Alpha ~ %s%s%s. Reference: %s.<br>%s"),
        if (has_base) "Baseline_&lt;metric&gt; + " else "",
        input$alpha_fixed_effects,
        if (length(re)) paste0(" + ", paste(sprintf("(1|%s)", re), collapse = " + "))
        else " &nbsp;(no random effect: ordinary linear model)",
        if (is.null(input$treat_ref)) "(not set)" else input$treat_ref,
        if (has_base) paste0(
          "Baseline <b>samples</b> are excluded from the model rows; each subject's ",
          "baseline <b>value</b> for that metric is added as the first covariate.") else "")))
    })

    output$alpha_lmm_padj_info_ui <- shiny::renderUI({
      shiny::req(alpha_lmm_result())
      m <- if (nzchar(.nz(input$alpha_padj))) input$alpha_padj else "BH"
      res   <- alpha_lmm_result()
      inter <- identical(attr(res, "fdr_target"), "interaction")
      terms <- attr(res, "fdr_terms")
      what  <- if (inter) "Treatment <i>interaction</i>" else "Treatment main-effect"
      n_cmp <- length(terms)
      .info(shiny::HTML(sprintf(paste0(
        "<b>p_value</b>: raw p-value from that metric's own model.<br>",
        "<b>p_adjusted</b>: %s-adjusted across all terms of one metric's model.<br>",
        "<b>p_adj_across_metrics</b>: %s-adjusted across the metrics for that one ",
        "comparison (the %s rows of every metric form one set).<br>",
        "<b>p_adj_across_metrics_and_groups</b>: every %s row of every metric adjusted ",
        "as one single set (usually, but not always, the larger p-value).<br>%s"),
        m, m, if (n_cmp) terms[1] else "treatment", what,
        if (n_cmp > 1)
          sprintf(paste0("Your treatment has %d comparisons (%s), so the two columns differ: ",
                         "the first treats each comparison as its own question, the second ",
                         "also corrects for the other comparisons. Report whichever you ",
                         "decided on beforehand, not the smaller one."),
                  n_cmp, paste(terms, collapse = ", "))
        else "With two treatment groups there is one comparison, so both columns are equal.")))
    })

    output$alpha_lmm_table <- DT::renderDT({
      tbl <- alpha_lmm_result(); shiny::req(tbl)
      .dt_table(tbl, c("Estimate", "SE", "df", "t_value", "p_value", "p_adjusted",
                       "p_adj_across_metrics", "p_adj_across_metrics_and_groups"),
                list(p_value = c(0.05, 0.1), p_adj_across_metrics = c(0.05, 0.25),
                     p_adj_across_metrics_and_groups = c(0.05, 0.25)),
                page = 20)
    })

    output$meta_table <- DT::renderDT({
      shiny::req(rv$meta_all)
      md  <- if (is.null(rv$obj)) rv$meta_all else rv$meta_enriched
      dec <- names(md)[vapply(md, function(x) is.numeric(x) && any(x %% 1 != 0, na.rm = TRUE), logical(1))]
      .dt_table(md, round_cols = dec, page = 10, filter = "top", dom = "lfrtip")
    })

    # -- Downloads: alpha / stats / metadata -------------------------------------
    output$dl_alpha <- shiny::downloadHandler(
      filename = function() sprintf("alpha_%s_%s.tiff",
        .fname(input$alpha_metric, input$alpha_x_var, "by", input$alpha_color_var), .today()),
      content = function(file) {
        a <- alpha_plot_args()
        export_alpha_tiff(a$meta_enriched, file, rv$obj$config, a$metric, a$x_var,
                          a$color_var, a$palette, a$rev_pal, a$connect_by,
                          a$exclude_baseline, a$is_baseline_col)
      })

    output$dl_alpha_lmm <- shiny::downloadHandler(
      filename = function() sprintf("alpha_lmm_results_%s.csv", .today()),
      content  = function(file) utils::write.csv(alpha_lmm_result(), file, row.names = FALSE))

    output$dl_stat <- shiny::downloadHandler(
      filename = function() {
        p <- stat_store$params
        if (is.null(p)) return("beta_diversity_stats.csv")
        sprintf("beta_%s_%s.csv", .fname(p$method, p$scope,
                                         if (!is.null(p$formula)) "formula" else p$var), .today())
      },
      content = function(file) utils::write.csv(stat_store$result, file, row.names = FALSE))

    output$dl_meta <- shiny::downloadHandler(
      filename = function() sprintf("enriched_metadata_%s.csv", .today()),
      content  = function(file) utils::write.csv(if (is.null(rv$obj)) rv$meta_all else rv$meta_enriched,
                                                 file, row.names = TRUE))

    # -- Taxonomy --------------------------------------------------------------------
    taxon_data_r <- shiny::eventReactive(input$run_taxon, {
      shiny::req(rv$obj, input$taxon_level)
      shiny::withProgress(message = sprintf("Parsing %s abundances...", input$taxon_level),
                          value = 0.4,
        tryCatch(compute_taxonomy_barplot_data(rv$obj$otu_path, rv$obj$meta,
                                               input$taxon_level, input$taxon_group_var,
                                               input$taxon_top_n, rv$obj$config),
                 error = function(e) {
                   shiny::showNotification(paste("Taxonomy plot failed:", conditionMessage(e)),
                                           type = "error")
                   NULL
                 }))
    })

    taxon_args <- shiny::reactive({
      td <- taxon_data_r()
      shiny::req(td, rv$obj, input$taxon_group_var)
      shiny::validate(shiny::need(
        identical(attr(td, "level_used"), input$taxon_level) &&
          identical(as.numeric(attr(td, "top_n_used")), as.numeric(input$taxon_top_n)),
        "Taxonomic level or Top N has changed - click \"Build taxonomy plot\" to apply."))
      list(taxon_data = td,
           level = if (length(rv$obj$tax_levels) == 1) "Feature" else input$taxon_level,
           group_var = input$taxon_group_var,
           palette = cat_pal(), rev_pal = isTRUE(input$rev_pal), top_n = input$taxon_top_n,
           exclude_baseline = isTRUE(input$taxon_exclude_baseline),
           is_baseline_col = rv$obj$config$is_baseline_col,
           x_label_var = input$taxon_x_label_var, x_order_var = input$taxon_x_order_var)
    })

    output$taxon_plot <- shiny::renderPlot(do.call(build_taxon_barplot, taxon_args()))

    output$dl_taxon <- shiny::downloadHandler(
      filename = function() sprintf("taxonomy_%s_%s.tiff",
                                    .fname(input$taxon_level, "by", input$taxon_group_var), .today()),
      content = function(file) do.call(export_taxon_barplot_tiff,
                                       c(taxon_args(), list(outfile = file))))

    output$dl_desc <- shiny::downloadHandler(
      filename = function() sprintf("taxa_summary_by_%s_%s.xlsx",
        if (length(input$desc_group_cols)) .fname(input$desc_group_cols) else "groups", .today()),
      content = function(file) {
        shiny::req(rv$obj)
        if (!length(input$desc_group_cols) || !length(input$desc_levels)) {
          shiny::showNotification("Choose at least one 'Group by' column and one taxonomic level.",
                                  type = "error")
          shiny::req(FALSE)
        }
        tryCatch(shiny::withProgress(message = "Building descriptive summary...", value = 0.5,
          export_taxa_descriptive_workbook(
            rv$obj$otu_path, rv$meta_enriched, input$desc_group_cols, outfile = file,
            sort_columns = if (length(input$desc_sort_cols)) input$desc_sort_cols
                           else input$desc_group_cols,
            tax_levels = input$desc_levels,
            exclude_baseline = isTRUE(input$desc_exclude_baseline),
            is_baseline_col = rv$obj$config$is_baseline_col)),
          error = function(e) {
            shiny::showNotification(paste("Descriptive summary failed:", conditionMessage(e)),
                                    type = "error", duration = NULL)
            shiny::req(FALSE)
          })
      })

    # -- MaAsLin3: reference (control) group for each categorical fixed effect -------
    maaslin_cat_vars <- shiny::reactive({
      shiny::req(rv$meta_enriched, rv$obj)
      fml <- trimws(.nz(input$maaslin_formula))
      if (!nzchar(fml)) return(list())
      vars <- tryCatch(all.vars(stats::as.formula(paste("~", fml))), error = function(e) character(0))
      me   <- rv$meta_enriched[!.is_baseline(rv$meta_enriched, rv$obj$config$is_baseline_col), , drop = FALSE]
      vars <- intersect(vars, colnames(me))
      lv   <- lapply(stats::setNames(vars, vars), function(v) {
        x <- me[[v]]
        if (is.character(x) || is.factor(x) || is.logical(x)) sort(unique(as.character(x[!is.na(x)])))
      })
      Filter(function(l) length(l) >= 2 && length(l) <= 50, lv)
    })
    ref_id <- function(v) paste0("maaslin_ref_", gsub("[^A-Za-z0-9]", "_", v))

    output$maaslin_ref_ui <- shiny::renderUI({
      cats <- maaslin_cat_vars()
      if (!length(cats)) return(NULL)
      shiny::tagList(
        lapply(names(cats), function(v) {
          cur <- shiny::isolate(input[[ref_id(v)]])
          shiny::selectInput(ref_id(v), sprintf("Reference (control) group for %s:", v),
                             choices = cats[[v]],
                             selected = if (!is.null(cur) && cur %in% cats[[v]]) cur else cats[[v]][1])
        }),
        .info("Every other group is compared with the reference group."))
    })

    maaslin_refs <- shiny::reactive({
      cats <- maaslin_cat_vars()
      refs <- lapply(names(cats), function(v) input[[ref_id(v)]])
      names(refs) <- names(cats)
      Filter(function(r) !is.null(r) && nzchar(r), refs)
    })

    # -- MaAsLin3 -----------------------------------------------------------------------
    shiny::observeEvent(input$run_maaslin_btn, {
      shiny::req(rv$obj, rv$meta_enriched, input$maaslin_formula)
      tryCatch({
        validate_maaslin_formula(input$maaslin_formula, rv$meta_enriched)
        prepped <- NULL
        res <- shiny::withProgress(message = "Running MaAsLin3 (this can take a while)...",
                                   value = 0.2, {
          prepped <- prepare_maaslin_input(rv$obj$otu_path, rv$meta_enriched,
                                           input$maaslin_level, rv$obj$config)
          shiny::incProgress(0.3, detail = "Fitting models...")
          run_maaslin(prepped$features, prepped$metadata, input$maaslin_formula,
                      group_effects_var = .nz(input$maaslin_group_var),
                      random_effect = input$maaslin_random_effect,
                      use_random = isTRUE(input$maaslin_use_random),
                      reference = maaslin_refs(),
                      small_random_effects = isTRUE(input$maaslin_small_random),
                      normalization = input$maaslin_normalization,
                      transform = input$maaslin_transform,
                      min_prevalence = input$maaslin_min_prevalence,
                      min_abundance = input$maaslin_min_abundance,
                      padj_method = input$maaslin_padj_method,
                      max_significance = input$maaslin_max_sig)
        })
        maaslin_store$result <- res
        maaslin_store$params <- list(
          level = input$maaslin_level, formula_used = attr(res, "formula_used"),
          reference = attr(res, "reference_used"),
          small_random = isTRUE(input$maaslin_small_random),
          normalization = input$maaslin_normalization, transform = input$maaslin_transform,
          min_prevalence = input$maaslin_min_prevalence, min_abundance = input$maaslin_min_abundance,
          padj_method = input$maaslin_padj_method, max_significance = input$maaslin_max_sig,
          n_features = ncol(prepped$features), n_samples = nrow(prepped$features),
          otu_scale_info = prepped$otu_scale_info)
        shiny::showNotification("MaAsLin3 finished!", type = "message")
      }, error = function(e)
        shiny::showNotification(paste("MaAsLin3 failed:", conditionMessage(e)),
                                type = "error", duration = NULL))
    })

    output$maaslin_info <- shiny::renderText({
      p <- maaslin_store$params; shiny::req(p)
      scale_msg <- .format_scale_info(p$otu_scale_info)
      refs <- if (length(p$reference)) paste(names(p$reference), "=", p$reference, collapse = "; ")
              else "(no categorical variable)"
      sprintf(paste0("Formula            : %s\nReference groups   : %s\nTaxonomic level    : %s\nSmall random fx    : %s\n",
                     "Normalization      : %s\nTransform          : %s\nMin prevalence     : %s\n",
                     "Min abundance (%%)  : %s\np-adjust           : %s\nMax significance   : %s\n",
                     "Features x Samples : %d x %d\nInput scale        : %s"),
              p$formula_used, refs, p$level, p$small_random, p$normalization, p$transform,
              p$min_prevalence, p$min_abundance, p$padj_method, p$max_significance,
              p$n_features, p$n_samples, if (is.null(scale_msg)) "(not available)" else scale_msg)
    })

    output$maaslin_table <- DT::renderDT({
      shiny::req(maaslin_store$result)
      .dt_table(maaslin_store$result, c("coef", "stderr", "pval", "qval"),
                list(qval = c(0.05, 0.25)), page = 15, filter = "top")
    })

    output$maaslin_sig_table <- DT::renderDT({
      shiny::req(maaslin_store$result, maaslin_store$params)
      .dt_table(summarise_maaslin_results(maaslin_store$result,
                                          maaslin_store$params$max_significance), page = 15)
    })

    output$dl_maaslin <- shiny::downloadHandler(
      filename = function() sprintf("maaslin3_results_%s.csv", .today()),
      content  = function(file) utils::write.csv(maaslin_store$result, file, row.names = FALSE))

    # -- PCoA downloads -------------------------------------------------------------------
    output$dl_2d <- shiny::downloadHandler(
      filename = function() sprintf("pcoa_2D_%s_PC%s-PC%s.tiff", .fname(input$color_var),
                                    input$pc_x, input$pc_y),
      content = function(file)
        export_2d_tiff(plot_obj(), input$color_var, file, as.integer(input$pc_x),
                       as.integer(input$pc_y), pcoa_pal(), isTRUE(input$is_cont),
                       isTRUE(input$rev_pal), connect_col()))

    output$dl_3d <- shiny::downloadHandler(
      filename = function() sprintf("pcoa_3D_%s.tiff", .fname(input$color_var)),
      content = function(file)
        export_3d_tiff(plot_obj(), input$color_var, file, pcoa_pal(), isTRUE(input$is_cont),
                       isTRUE(input$rev_pal), connect_col()))

    # Build every settings column right after loading, so its inputs hold their
    # defaults even before the page is opened.
    for (id in c("alpha_side_ui", "pcoa_side_ui", "beta_side_ui", "taxon_side_ui",
                 "maaslin_side_ui", "meas_controls_ui"))
      shiny::outputOptions(output, id, suspendWhenHidden = FALSE)
  }

  shiny::shinyApp(ui = ui, server = server)
}
