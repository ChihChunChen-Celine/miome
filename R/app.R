# app.R
#' Launch the PCoA / Microbiome Explorer Shiny App
#'
#' @return Launches a Shiny app in your browser.
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarLayout
#' @importFrom shiny sidebarPanel mainPanel tabsetPanel tabPanel
#' @importFrom shiny fileInput actionButton uiOutput renderUI
#' @importFrom shiny observeEvent reactive reactiveValues req
#' @importFrom shiny validate need showNotification updateSelectInput
#' @importFrom shiny downloadButton downloadHandler renderPlot
#' @importFrom shiny plotOutput tagList tags br hr selectInput
#' @importFrom shiny checkboxInput checkboxGroupInput conditionalPanel
#' @importFrom shiny numericInput radioButtons verbatimTextOutput
#' @importFrom shiny withProgress setProgress
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom DT DTOutput renderDT datatable
#' @export
run_app <- function() {

  # -- UI --------------------------------------------------------------------
  ui <- shiny::fluidPage(
    shiny::titlePanel("Microbiome Explorer - PCoA + Alpha Diversity"),

    # -- Inline CSS for warning / info boxes ------------------------------
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .warn-box {
          background-color: #fff3cd;
          border: 1px solid #ffc107;
          border-radius: 4px;
          padding: 8px 12px;
          margin-top: 6px;
          font-size: 12px;
        }
        .info-box {
          background-color: #d1ecf1;
          border: 1px solid #bee5eb;
          border-radius: 4px;
          padding: 8px 12px;
          margin-top: 6px;
          font-size: 12px;
        }
      "))
    ),

    shiny::sidebarLayout(

      # -- Sidebar -----------------------------------------------------------
      shiny::sidebarPanel(
        width = 3,

        # -- Step 1: Upload ----------------------------------------------
        shiny::tags$b("Step 1 - Upload data"),
        shiny::fileInput(
          "otu_file",
          "Species/OTU table (.tsv, samples as columns)",
          accept = c(".tsv", ".txt")
        ),
        shiny::fileInput(
          "meta_file",
          "Metadata (.csv, first column = Sample ID)",
          accept = c(".csv")
        ),

        # -- Precomputed UniFrac distance matrix (optional) --------------
        shiny::fileInput(
          "unifrac_dist_file",
          "Precomputed UniFrac distance matrix (.qza or .tsv) - required for UniFrac",
          accept = c(".qza", ".tsv", ".txt")
        ),
        # Warn if UniFrac selected but no matrix uploaded
        shiny::conditionalPanel(
          condition = "input.dist_method == 'unifrac_weighted' ||
                       input.dist_method == 'unifrac_unweighted'",
          shiny::uiOutput("unifrac_warning_ui")
        ),

        # -- Distance metric ---------------------------------------------
        shiny::selectInput(
          "dist_method", "Distance metric:",
          choices  = DIST_METHODS,
          selected = "bray"
        ),
        shiny::actionButton(
          "load_btn", "Load & Compute",
          class = "btn-primary btn-block"
        ),
        shiny::tags$hr(),

        # -- Dynamic post-load UI -----------------------------------------
        shiny::uiOutput("post_load_ui")
      ),

      # -- Main panel --------------------------------------------------------
      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          id = "main_tabs",

          # Tab 1: Alpha Diversity
          shiny::tabPanel(
            "Alpha Diversity",
            shiny::br(),
            shiny::uiOutput("alpha_controls_ui"),
            shiny::br(),
            shiny::plotOutput("alpha_plot", height = "500px"),
            shiny::br(),
            shiny::tags$b("Linear Mixed Model - Treatment effect"),
            shiny::uiOutput("alpha_lmm_info_ui"),
            shiny::br(),
            shiny::uiOutput("alpha_lmm_padj_info_ui"),
            shiny::br(),
            shiny::downloadButton(
              "dl_alpha_lmm", "Download LMM results (.csv)"
            ),
            shiny::br(), shiny::br(),
            DT::DTOutput("alpha_lmm_table")
          ),

          # Tab 2: 3D interactive PCoA
          shiny::tabPanel(
            "3D interactive",
            plotly::plotlyOutput("plot3d", height = "700px")
          ),

          # Tab 3: 2D PCoA
          shiny::tabPanel(
            "2D preview",
            shiny::plotOutput("plot2d", height = "700px")
          ),

          # Tab 4: Beta Diversity Statistics
          shiny::tabPanel(
            "Beta Diversity Stats",
            shiny::br(),
            shiny::uiOutput("otu_scale_info_ui"),
            shiny::br(),
            shiny::verbatimTextOutput("stat_info"),
            shiny::br(),
            shiny::tags$b("Test result"),
            shiny::br(), shiny::br(),
            shiny::downloadButton(
              "dl_stat", "Download stat results (.csv)"
            ),
            shiny::br(), shiny::br(),
            DT::DTOutput("stat_table"),
            shiny::br(),
            shiny::tags$b("PERMDISP - Homogeneity of dispersion"),
            shiny::tags$div(
              class = "info-box",
              "If PERMDISP is significant, interpret PERMANOVA with caution
               (difference may reflect dispersion, not centroid location)."
            ),
            shiny::br(),
            shiny::verbatimTextOutput("permdisp_out")
          ),

          # Tab 5: Enriched metadata preview
          shiny::tabPanel(
            "Metadata Preview",
            shiny::br(),
            shiny::tags$b(
              "Enriched metadata (alpha diversity + baseline covariates merged)"
            ),
            shiny::br(), shiny::br(),
            shiny::downloadButton(
              "dl_meta", "Download enriched metadata (.csv)"
            ),
            shiny::br(), shiny::br(),
            DT::DTOutput("meta_table")
          ),

          # Tab 6: Taxonomy bar plots
          shiny::tabPanel(
            "Taxonomy",
            shiny::br(),
            shiny::uiOutput("taxon_controls_ui"),
            shiny::br(),
            shiny::plotOutput("taxon_plot", height = "650px"),
            shiny::br(),
            shiny::downloadButton("dl_taxon", "Download taxonomy plot (TIFF)")
          ),

          # Tab 7: MaAsLin3
          shiny::tabPanel(
            "MaAsLin3",
            shiny::br(),
            shiny::verbatimTextOutput("maaslin_info"),
            shiny::br(),
            shiny::tags$b("All results"),
            shiny::br(), shiny::br(),
            shiny::downloadButton(
              "dl_maaslin", "Download full results (.csv)"
            ),
            shiny::br(), shiny::br(),
            DT::DTOutput("maaslin_table"),
            shiny::br(),
            shiny::tags$b("Significant hits"),
            shiny::tags$div(
              class = "info-box",
              "Filtered to q-value below the Max significance threshold
                 set in the sidebar."
            ),
            shiny::br(),
            DT::DTOutput("maaslin_sig_table")
          )
        )
      )
    )
  )

  # -- SERVER --------------------------------------------------------------
  server <- function(input, output, session) {

    rv <- shiny::reactiveValues(
      obj           = NULL,
      meta_cols     = NULL,
      meta_enriched = NULL
    )

    stat_store <- shiny::reactiveValues(
      result = NULL,
      params = NULL
    )

    maaslin_store <- shiny::reactiveValues(
      result = NULL,
      params = NULL
    )

    # -- Distance label lookup ------------------------------------------
    dist_labels <- stats::setNames(names(DIST_METHODS), DIST_METHODS)

    # -- Format the otu_scale_info list (data.R / .to_relative_abundance())
    # into a user-facing message. Shared by the sidebar box and the Beta
    # Diversity Stats tab so both stay in sync automatically.
    .format_scale_info <- function(info) {
      if (is.null(info)) return(NULL)
      if (!is.null(info$note)) return(info$note)   # UniFrac case

      if (info$n_converted == 0) {
        sprintf(
          paste0(
            "All %d sample(s) were already in percent (relative abundance) ",
            "- used as uploaded, no conversion applied."
          ),
          info$n_total
        )
      } else if (info$n_already_percent == 0) {
        sprintf(
          paste0(
            "All %d sample(s) looked like raw counts - converted to ",
            "percent before beta diversity / MaAsLin3."
          ),
          info$n_total
        )
      } else {
        sprintf(
          paste0(
            "%d of %d sample(s) were already percent (used as-is); %d ",
            "looked like raw counts and were converted to percent: %s"
          ),
          info$n_already_percent, info$n_total, info$n_converted,
          paste(info$converted_samples, collapse = ", ")
        )
      }
    }

    # -- OTU scale info: sidebar (visible immediately after load) -------
    output$otu_scale_info_sidebar_ui <- shiny::renderUI({
      shiny::req(rv$obj)
      msg <- .format_scale_info(rv$obj$otu_scale_info)
      shiny::req(msg)
      shiny::tags$div(
        class = "info-box",
        shiny::tags$b("Beta diversity input scale: "), msg
      )
    })

    # -- OTU scale info: Beta Diversity Stats tab ------------------------
    output$otu_scale_info_ui <- shiny::renderUI({
      shiny::req(rv$obj)
      msg <- .format_scale_info(rv$obj$otu_scale_info)
      shiny::req(msg)
      shiny::tags$div(
        class = "info-box",
        shiny::tags$b("Input scale detected: "), msg,
        shiny::tags$br(),
        "PCoA plot, PERMANOVA/ANOSIM, and PERMDISP all use this same
         distance object, so this applies to all of them."
      )
    })

    # -- UniFrac upload warning -----------------------------------------
    # Shown when a UniFrac metric is selected but no matrix has been
    # uploaded yet; disappears as soon as the user supplies a file.
    output$unifrac_warning_ui <- shiny::renderUI({
      if (is.null(input$unifrac_dist_file)) {
        shiny::tags$div(
          class = "warn-box",
          "\u26a0 No UniFrac distance matrix uploaded.
           Please upload a precomputed .qza or .tsv file."
        )
      }
    })

    # -- Load & compute on button click --------------------------------
    shiny::observeEvent(input$load_btn, {
      shiny::req(input$otu_file, input$meta_file)

      tryCatch({
        cfg <- pcoa_config(
          dist_method       = input$dist_method,
          plot_title        = "PCoA",
          dist_method_label = dist_labels[[input$dist_method]]
        )

        # Guard: UniFrac requires a precomputed distance matrix
        if (input$dist_method %in%
            c("unifrac_weighted", "unifrac_unweighted") &&
            is.null(input$unifrac_dist_file)) {
          shiny::showNotification(
            paste(
              "Please upload a precomputed UniFrac distance matrix",
              "(.qza or .tsv) before computing."
            ),
            type = "error"
          )
          return()
        }

        # Resolve the precomputed-distance path (NULL for all other metrics)
        precomputed_path <- if (
          input$dist_method %in%
          c("unifrac_weighted", "unifrac_unweighted") &&
          !is.null(input$unifrac_dist_file)
        ) {
          input$unifrac_dist_file$datapath
        } else {
          NULL
        }

        # Run full pipeline: alpha -> baseline covariates -> PCoA
        obj <- prepare_all(
          otu_path              = input$otu_file$datapath,
          meta_path             = input$meta_file$datapath,
          config                = cfg,
          dist_method           = input$dist_method,
          precomputed_dist_path = precomputed_path
        )

        rv$obj           <- obj
        rv$obj$otu_path  <- input$otu_file$datapath
        rv$obj$meta_path <- input$meta_file$datapath
        rv$meta_cols     <- obj$meta_cols
        rv$meta_enriched <- obj$meta_enriched

        # Clear stats from a previous dataset
        stat_store$result <- NULL
        stat_store$params <- NULL
        maaslin_store$result <- NULL
        maaslin_store$params <- NULL

        shiny::showNotification(
          "Data loaded - PCoA and alpha diversity computed!",
          type = "message"
        )
      }, error = function(e) {
        shiny::showNotification(paste("Error:", e$message), type = "error")
      })
    })

    # -- Dynamic sidebar controls (shown after load) -------------------
    output$post_load_ui <- shiny::renderUI({
      shiny::req(rv$obj)
      meta_cols <- rv$meta_cols
      cfg       <- rv$obj$config

      shiny::tagList(

        # -- Data scale detected for beta diversity (visible immediately) --
        shiny::uiOutput("otu_scale_info_sidebar_ui"),
        shiny::tags$hr(),

        # -- Step 2: Alpha diversity LMM -----------------------------
        shiny::tags$b("Step 2 - Alpha diversity LMM"),
        shiny::br(), shiny::br(),

        shiny::selectizeInput(
          "lmm_fixed_effects",
          "Fixed effects (covariates):",
          choices  = setdiff(
            colnames(rv$meta_enriched),
            c("SampleID", cfg$is_baseline_col)
          ),
          selected = intersect(
            c(cfg$timepoint_col,
              cfg$sequence_col,
              cfg$treatment_col,
              paste0("Baseline_", cfg$alpha_metrics)),
            colnames(rv$meta_enriched)
          ),
          multiple = TRUE,
          options  = list(
            placeholder  = "Click to select fixed effects...",
            plugins      = list("remove_button"),  # shows  on each chip
            maxOptions   = 200
          )
        ),
        shiny::selectInput(
          "lmm_random_effect",
          "Random effect (grouping / subject ID):",
          choices  = colnames(rv$meta_enriched),
          selected = if (cfg$subject_id_col %in% colnames(rv$meta_enriched))
            cfg$subject_id_col else colnames(rv$meta_enriched)[1]
        ),

        shiny::selectInput(
          "treat_ref", "Treatment reference level:",
          choices = sort(unique(as.character(
            rv$meta_enriched[
              rv$meta_enriched[[cfg$is_baseline_col]] == FALSE,
              cfg$treatment_col
            ]
          )))
        ),
        shiny::actionButton(
          "run_alpha_lmm", "Run alpha diversity LMM",
          class = "btn-warning btn-block"
        ),
        shiny::tags$div(
          class = "info-box",
          shiny::uiOutput("lmm_formula_info_ui")
        ),
        shiny::tags$hr(),

        # -- Step 3: PCoA display options ------------------------------
        shiny::tags$b("Step 3 - PCoA display options"),
        shiny::br(), shiny::br(),
        shiny::selectInput(
          "color_var", "Colour by:",
          choices  = meta_cols,
          selected = if (cfg$treatment_col %in% meta_cols)
            cfg$treatment_col else meta_cols[1]
        ),
        shiny::checkboxInput(
          "is_cont",
          "Treat as continuous (numeric gradient)?",
          value = FALSE
        ),
        shiny::selectInput(
          "palette", "Colour palette:",
          choices = CAT_PALETTES
        ),
        shiny::checkboxInput("rev_pal", "Reverse palette", value = FALSE),
        shiny::checkboxInput(
          "connect_pairs",
          "Connect paired samples (within-subject lines)",
          value = identical(cfg$study_design, "crossover")
        ),
        shiny::tags$div(
          class = "info-box",
          sprintf(
            "Draws a line between all samples sharing the same %s.
             Most useful for crossover / repeated-measures designs.",
            cfg$subject_id_col
          )
        ),
        shiny::tags$hr(),

        # -- 2D plot axes ---------------------------------------------
        shiny::tags$b("2D plot axes"),
        shiny::selectInput(
          "pc_x", "X axis:",
          choices  = setNames(1:3, c("PC1", "PC2", "PC3")),
          selected = 1
        ),
        shiny::selectInput(
          "pc_y", "Y axis:",
          choices  = setNames(1:3, c("PC1", "PC2", "PC3")),
          selected = 2
        ),
        shiny::tags$hr(),

        # -- Downloads ------------------------------------------------
        shiny::tags$b("Downloads"),
        shiny::br(), shiny::br(),
        shiny::downloadButton("dl_2d", "Download 2D (TIFF)"),
        shiny::br(), shiny::br(),
        shiny::downloadButton("dl_3d", "Download 3D (TIFF)"),
        shiny::tags$hr(),

        # -- Step 4: Statistical analysis ----------------------------
        shiny::tags$b("Step 4 - Beta diversity statistics"),
        shiny::br(), shiny::br(),
        shiny::selectInput(
          "stat_var", "Primary grouping variable:",
          choices  = meta_cols,
          selected = if (cfg$treatment_col %in% meta_cols)
            cfg$treatment_col else meta_cols[1]
        ),
        shiny::selectInput(
          "stat_method", "Test:",
          choices = c("PERMANOVA" = "permanova",
                      "ANOSIM"   = "anosim")
        ),
        shiny::radioButtons(
          "stat_scope", "Scope:",
          choices = c("Global"   = "global",
                      "Pairwise" = "pairwise"),
          inline  = TRUE
        ),

        # -- Covariates (PERMANOVA only) ------------------------------
        shiny::conditionalPanel(
          condition = "input.stat_method == 'permanova'",
          shiny::selectizeInput(
            "covariates",
            "Covariates (added before primary variable in model):",
            choices  = setdiff(
              meta_cols,
              c("SampleID", "PC1", "PC2", "PC3")
            ),
            selected = intersect(
              c(cfg$default_covariates,
                paste0("Baseline_", cfg$alpha_metrics)),
              meta_cols
            ),
            multiple = TRUE,
            options  = list(
              placeholder  = "Click to add covariates...",
              plugins      = list("remove_button"),
              maxOptions   = 200
            )
          ),
          shiny::tags$div(
            class = "info-box",
            sprintf(
              "For crossover design: include %s (period effect),
               %s (carryover), and a Baseline_* numeric covariate.",
              cfg$timepoint_col, cfg$sequence_col
            )
          )
        ),

        # -- ANOSIM covariate warning ---------------------------------
        shiny::conditionalPanel(
          condition = "input.stat_method == 'anosim'",
          shiny::tags$div(
            class = "warn-box",
            "\u26a0 ANOSIM does not support covariates. Use PERMANOVA
             for covariate-adjusted analysis."
          )
        ),

        # -- Strata (restricted permutations) ------------------------
        shiny::selectInput(
          "strata_var",
          "Restrict permutations within (strata):",
          choices  = c("None" = "", meta_cols),
          selected = if (!is.null(cfg$default_strata) &&
                         cfg$default_strata %in% meta_cols)
            cfg$default_strata else ""
        ),
        shiny::tags$div(
          class = "info-box",
          sprintf(
            "For crossover design: select %s to restrict permutations
             within each subject.", cfg$subject_id_col
          )
        ),
        shiny::numericInput(
          "stat_perm", "Permutations:",
          value = 999, min = 99, max = 9999, step = 100
        ),
        shiny::selectInput(
          "padj_method", "P-value adjustment:",
          choices  = PADJ_METHODS,
          selected = cfg$padj_method
        ),
        shiny::actionButton(
          "run_stat", "Run beta diversity test",
          class = "btn-success btn-block"
        ),
        shiny::tags$hr(),

        # -- Step 5: Taxonomy bar plot --------------------------------
        shiny::tags$b("Step 5 - Taxonomy bar plot"),
        shiny::br(), shiny::br(),
        shiny::selectInput(
          "taxon_level", "Taxonomic level:",
          choices  = TAXA_LEVELS,
          selected = "Family"
        ),
        shiny::selectInput(
          "taxon_group_var", "Facet / group by:",
          choices  = meta_cols,
          selected = if (cfg$treatment_col %in% meta_cols)
            cfg$treatment_col else meta_cols[1]
        ),
        shiny::numericInput(
          "taxon_top_n", "Top N taxa:",
          value = 14L, min = 3L, max = 30L, step = 1L
        ),
        shiny::checkboxInput(
          "taxon_exclude_baseline",
          "Exclude baseline samples",
          value = TRUE
        ),
        shiny::selectInput(
          "taxon_x_label_var", "X-axis label:",
          choices  = c("Sample ID" = "SampleID", meta_cols),
          selected = "SampleID"
        ),
        shiny::selectInput(
          "taxon_x_order_var", "Order bars by:",
          choices  = c("Sample ID (alphabetical)" = "SampleID", meta_cols),
          selected = "SampleID"
        ),
        shiny::tags$div(
          class = "info-box",
          "Label and order apply live \u2014 no need to re-click",
          shiny::tags$i("Build taxonomy plot"),
          "unless you change the taxonomic level or Top N."
        ),
        shiny::actionButton(
          "run_taxon", "Build taxonomy plot",
          class = "btn-info btn-block"
        ),
        shiny::tags$hr(),

        # -- Step 6: MaAsLin3 -----------------------------------------
        shiny::tags$b("Step 6 - MaAsLin3"),
        shiny::br(), shiny::br(),
        shiny::selectInput(
          "maaslin_level", "Taxonomic level:",
          choices  = TAXA_LEVELS,
          selected = "Genus"
        ),
        shiny::textInput(
          "maaslin_formula", "Fixed effects (formula, no leading '~'):",
          value = paste(
            intersect(
              c(cfg$timepoint_col, cfg$sequence_col, cfg$treatment_col),
              meta_cols
            ),
            collapse = " + "
          )
        ),
        shiny::tags$div(
          class = "info-box",
          "e.g. \"Treatment + Timepoint + Baseline_Shannon\". Do not
           include the grouping factor or random effect here \u2014
           set those separately below."
        ),
        shiny::selectInput(
          "maaslin_group_var",
          "Grouping factor (optional, ANOVA-style group() term):",
          choices  = c("None" = "", colnames(rv$meta_enriched)),
          selected = ""
        ),
        shiny::checkboxInput(
          "maaslin_use_random",
          "Include random effect",
          value = identical(cfg$study_design, "crossover")
        ),
        shiny::selectInput(
          "maaslin_random_effect", "Random effect (grouping / subject ID):",
          choices  = colnames(rv$meta_enriched),
          selected = if (cfg$subject_id_col %in% colnames(rv$meta_enriched))
            cfg$subject_id_col else colnames(rv$meta_enriched)[1]
        ),
        shiny::checkboxInput(
          "maaslin_small_random",
          "Random effect group sizes are small",
          value = FALSE
        ),
        shiny::tags$div(
          class = "info-box",
          "When checked (small_random_effects = TRUE), MaAsLin3
           automatically replaces the random effect with a fixed effect
           in the logistic prevalence model for groups with few
           observations."
        ),
        shiny::selectInput(
          "maaslin_normalization", "Normalization:",
          choices  = c("TSS", "CLR", "NONE"),
          selected = "TSS"
        ),
        shiny::selectInput(
          "maaslin_transform", "Transform:",
          choices  = c("LOG", "PLOG", "NONE"),
          selected = "LOG"
        ),
        shiny::numericInput(
          "maaslin_min_prevalence", "Min prevalence:",
          value = 0.1, min = 0, max = 1, step = 0.05
        ),
        shiny::numericInput(
          "maaslin_min_abundance", "Min abundance:",
          value = 0, min = 0, step = 0.001
        ),
        shiny::selectInput(
          "maaslin_padj_method", "P-value adjustment:",
          choices  = PADJ_METHODS,
          selected = cfg$padj_method
        ),
        shiny::numericInput(
          "maaslin_max_sig", "Max significance (q-value) reported:",
          value = 0.25, min = 0, max = 1, step = 0.05
        ),
        shiny::actionButton(
          "run_maaslin_btn", "Run MaAsLin3",
          class = "btn-danger btn-block"
        )
      )
    })

    # -- Keep palette dropdown in sync with continuous/categorical ------
    shiny::observeEvent(input$is_cont, {
      shiny::req(input$palette)
      new_choices <- if (isTRUE(input$is_cont)) CONT_PALETTES else CAT_PALETTES
      shiny::updateSelectInput(
        session, "palette",
        choices  = new_choices,
        selected = new_choices[1]
      )
    }, ignoreInit = TRUE)

    # -- Prevent same axis on X and Y ----------------------------------
    shiny::observeEvent(input$pc_x, {
      shiny::req(input$pc_y)
      if (input$pc_x == input$pc_y) {
        alt <- setdiff(c("1", "2", "3"), input$pc_x)[1]
        shiny::updateSelectInput(session, "pc_y", selected = alt)
      }
    })

    # -- Resolve "connect paired samples" column ------------------------
    connect_col <- shiny::reactive({
      shiny::req(rv$obj)
      if (isTRUE(input$connect_pairs)) rv$obj$config$subject_id_col else NULL
    })

    # -- 3D plot --------------------------------------------------------
    output$plot3d <- plotly::renderPlotly({
      shiny::req(rv$obj, input$color_var)
      build_3d(
        rv$obj,
        v          = input$color_var,
        palette    = input$palette,
        is_cont    = input$is_cont,
        rev_pal    = input$rev_pal,
        connect_by = connect_col()
      )
    })

    # -- 2D plot --------------------------------------------------------
    output$plot2d <- shiny::renderPlot({
      shiny::req(rv$obj, input$color_var, input$pc_x, input$pc_y)
      shiny::validate(
        shiny::need(
          input$pc_x != input$pc_y,
          "Please choose two different PC axes."
        )
      )
      build_2d(
        rv$obj,
        v          = input$color_var,
        ax_x       = as.integer(input$pc_x),
        ax_y       = as.integer(input$pc_y),
        palette    = input$palette,
        is_cont    = input$is_cont,
        rev_pal    = input$rev_pal,
        connect_by = connect_col()
      )
    })

    # -- Beta diversity statistics --------------------------------------
    # Results and parameters are frozen on button click so the table,
    # the info text, and PERMDISP always describe the same run - even if
    # the user changes a dropdown without re-clicking "Run".
    shiny::observeEvent(input$run_stat, {
      shiny::req(rv$obj, input$stat_var)

      obj_use                    <- rv$obj
      obj_use$config$padj_method <- input$padj_method

      covs <- if (input$stat_method == "permanova" &&
                  length(input$covariates) > 0)
        input$covariates else NULL

      strata <- if (nchar(input$strata_var) > 0) input$strata_var else NULL

      tryCatch({
        res <- shiny::withProgress(
          message = "Running permutation test...", value = 0.5, {
            if (input$stat_scope == "global") {
              run_global_test(
                obj          = obj_use,
                v            = input$stat_var,
                method       = input$stat_method,
                covariates   = covs,
                strata       = strata,
                permutations = input$stat_perm
              )
            } else {
              run_pairwise_test(
                obj          = obj_use,
                v            = input$stat_var,
                method       = input$stat_method,
                covariates   = covs,
                strata       = strata,
                permutations = input$stat_perm
              )
            }
          }
        )

        stat_store$result <- res
        stat_store$params <- list(
          method      = input$stat_method,
          scope       = input$stat_scope,
          var         = input$stat_var,
          covariates  = covs,
          strata      = strata,
          padj_method = input$padj_method,
          perms       = input$stat_perm
        )

      }, error = function(e) {
        shiny::showNotification(
          paste("Stat test failed:", e$message),
          type = "error"
        )
      })
    })

    # -- Stat info bar (reads FROZEN params, not live inputs) -----------
    output$stat_info <- shiny::renderText({
      shiny::req(stat_store$params)
      p <- stat_store$params

      covs_txt   <- if (length(p$covariates) > 0)
        paste(p$covariates, collapse = " + ") else "none"
      strata_txt <- if (!is.null(p$strata)) p$strata else "none"

      sprintf(
        paste0(
          "Test      : %s\n",
          "Scope     : %s\n",
          "Variable  : %s\n",
          "Covariates: %s\n",
          "Strata    : %s\n",
          "p-adjust  : %s\n",
          "Perms     : %s"
        ),
        toupper(p$method), p$scope, p$var,
        covs_txt, strata_txt, p$padj_method, p$perms
      )
    })

    # -- Stat results table ---------------------------------------------
    output$stat_table <- DT::renderDT({
      shiny::req(stat_store$result)
      tbl <- stat_store$result

      DT::datatable(
        tbl,
        rownames = FALSE,
        options  = list(
          scrollX    = TRUE,
          pageLength = 10,
          dom        = "tip"
        )
      ) %>%
        DT::formatRound(
          columns = intersect(
            c("Statistic_R2", "Statistic_R", "Statistic",
              "F_value", "p_value", "p_adjusted"),
            colnames(tbl)
          ),
          digits = 4
        ) %>%
        DT::formatStyle(
          "p_value",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.1),
            c("#d4edda", "#fff3cd", "white")
          )
        )
    })

    # -- PERMDISP output (reads FROZEN params) --------------------------
    output$permdisp_out <- shiny::renderPrint({
      shiny::req(stat_store$result, stat_store$params, rv$obj)
      p <- stat_store$params
      shiny::req(p$method == "permanova")

      tryCatch({
        obj <- rv$obj
        cfg <- obj$config
        md  <- obj$meta

        # Filter to non-baseline samples only
        if (cfg$is_baseline_col %in% colnames(md)) {
          md <- md[md[[cfg$is_baseline_col]] == FALSE, , drop = FALSE]
        }

        keep <- !is.na(md[[p$var]])
        md   <- md[keep, , drop = FALSE]

        d   <- stats::as.dist(
          as.matrix(obj$dist)[rownames(md), rownames(md)]
        )
        grp <- factor(md[[p$var]])

        if (nlevels(grp) < 2) {
          cat("Not enough levels for PERMDISP.\n")
          return(invisible(NULL))
        }

        bd <- vegan::betadisper(d, grp)
        pt <- vegan::permutest(bd, permutations = p$perms)

        cat("-- PERMDISP: Homogeneity of Multivariate Dispersion --\n")
        cat(sprintf("Grouping variable : %s\n", p$var))
        cat(sprintf("Permutations      : %s\n\n", p$perms))
        print(pt)
        cat("\nGroup dispersions (average distance to centroid):\n")
        print(round(tapply(bd$distances, bd$group, mean), 4))

      }, error = function(e) {
        cat("PERMDISP could not be computed:\n", e$message, "\n")
      })
    })

    # -- Alpha diversity controls (dynamic, shown after load) -----------
    output$alpha_controls_ui <- shiny::renderUI({
      shiny::req(rv$meta_enriched, rv$obj)
      cfg <- rv$obj$config

      available_metrics <- intersect(ALPHA_METRICS, colnames(rv$meta_enriched))

      shiny::tagList(
        shiny::fluidRow(
          shiny::column(
            4,
            shiny::selectInput(
              "alpha_metric",
              "Alpha diversity metric:",
              choices  = available_metrics,
              selected = if ("Shannon" %in% available_metrics)
                "Shannon" else available_metrics[1]
            )
          ),
          shiny::column(
            4,
            shiny::selectInput(
              "alpha_x_var",
              "X axis grouping:",
              choices  = rv$meta_cols,
              selected = if (cfg$treatment_col %in% rv$meta_cols)
                cfg$treatment_col else rv$meta_cols[1]
            )
          ),
          shiny::column(
            4,
            shiny::selectInput(
              "alpha_color_var",
              "Colour by:",
              choices  = rv$meta_cols,
              selected = if (cfg$timepoint_col %in% rv$meta_cols)
                cfg$timepoint_col else rv$meta_cols[1]
            )
          )
        ),
        shiny::checkboxInput(
          "alpha_exclude_baseline",
          "Exclude baseline samples from plot",
          value = TRUE
        ),
        shiny::downloadButton("dl_alpha", "Download alpha plot (TIFF)")
      )
    })

    # -- Alpha diversity plot -------------------------------------------
    output$alpha_plot <- shiny::renderPlot({
      shiny::req(
        rv$meta_enriched, rv$obj,
        input$alpha_metric, input$alpha_x_var,
        input$alpha_color_var, input$palette
      )
      cfg <- rv$obj$config

      shiny::validate(
        shiny::need(
          input$alpha_metric %in% colnames(rv$meta_enriched),
          paste("Metric", input$alpha_metric,
                "not found in enriched metadata.")
        )
      )

      build_alpha_scatter(
        meta_enriched    = rv$meta_enriched,
        metric           = input$alpha_metric,
        x_var            = input$alpha_x_var,
        color_var        = input$alpha_color_var,
        palette          = input$palette,
        rev_pal          = input$rev_pal,
        connect_by       = connect_col(),
        exclude_baseline = isTRUE(input$alpha_exclude_baseline),
        is_baseline_col  = cfg$is_baseline_col
      )
    })

    # -- Alpha diversity LMM --------------------------------------------
    alpha_lmm_result <- shiny::eventReactive(input$run_alpha_lmm, {
      shiny::req(rv$meta_enriched, rv$obj, input$treat_ref)
      cfg <- rv$obj$config

      shiny::validate(
        shiny::need(
          input$lmm_random_effect %in% colnames(rv$meta_enriched),
          sprintf("'%s' column required for LMM random effect.",
                  input$lmm_random_effect)
        ),
        shiny::need(
          cfg$treatment_col %in% colnames(rv$meta_enriched),
          sprintf("'%s' column required.", cfg$treatment_col)
        )
      )

      tryCatch({
        shiny::withProgress(
          message = "Fitting linear mixed models...", value = 0.3, {
            run_alpha_lmm_all_metrics(
              meta_enriched   = rv$meta_enriched,
              fixed_effects   = input$lmm_fixed_effects,
              random_effect   = input$lmm_random_effect,
              treatment_ref   = input$treat_ref,
              treatment_col   = cfg$treatment_col,
              is_baseline_col = cfg$is_baseline_col,
              padj_method     = input$padj_method
            )
          }
        )
      }, error = function(e) {
        shiny::showNotification(
          paste("Alpha LMM failed:", e$message),
          type = "error"
        )
        NULL
      })
    })

    # -- Alpha LMM formula info -----------------------------------------
    output$alpha_lmm_info_ui <- shiny::renderUI({
      shiny::req(rv$obj, input$lmm_fixed_effects, input$lmm_random_effect)

      shiny::tags$div(
        class = "info-box",
        shiny::HTML(sprintf(
          paste0(
            "Model: Alpha ~ Baseline_&lt;metric&gt; + %s + (1|%s).  ",
            "Reference: %s.<br>",
            "Baseline <b>samples</b> are excluded from the model rows, but ",
            "each subject's baseline <b>value</b> (e.g. Baseline_Shannon for ",
            "the Shannon model, Baseline_Observed for Observed_ASVs) is ",
            "automatically added as the first covariate when available, so ",
            "the Treatment effect is estimated adjusting for where each ",
            "subject started."
          ),
          paste(input$lmm_fixed_effects, collapse = " + "),
          input$lmm_random_effect,
          if (!is.null(input$treat_ref)) input$treat_ref else "(not set)"
        ))
      )
    })

    # -- Alpha LMM p-value column explanation (static per-run) -----------
    output$alpha_lmm_padj_info_ui <- shiny::renderUI({
      shiny::req(alpha_lmm_result())

      shiny::tags$div(
        class = "info-box",
        shiny::HTML(sprintf(
          paste0(
            "<b>p_value</b>: raw p-value for that term, from that metric's ",
            "own model.<br>",
            "<b>p_adjusted</b>: %s-adjusted <i>within</i> one metric's model ",
            "\u2014 across all its terms (Baseline, Timepoint, Sequence, ",
            "Treatment, etc.).<br>",
            "<b>p_adjusted_cross_metric</b>: %s-adjusted <i>across metrics</i> ",
            "\u2014 pools the raw Treatment p-value from every alpha metric's ",
            "model (Shannon, Simpson, Observed_ASVs, ...) and adjusts that set. ",
            "Shown for Treatment rows only; blank for other terms."
          ),
          input$padj_method, input$padj_method
        ))
      )
    })

    # -- Alpha LMM results table ----------------------------------------
    output$alpha_lmm_table <- DT::renderDT({
      shiny::req(alpha_lmm_result())
      tbl <- alpha_lmm_result()

      dt <- DT::datatable(
        tbl,
        rownames = FALSE,
        options  = list(
          scrollX    = TRUE,
          pageLength = 20,
          dom        = "tip"
        )
      ) %>%
        DT::formatRound(
          columns = intersect(
            c("Estimate", "SE", "df", "t_value",
              "p_value", "p_adjusted", "p_adjusted_cross_metric"),
            colnames(tbl)
          ),
          digits = 4
        )

      if ("p_value" %in% colnames(tbl)) {
        dt <- dt %>% DT::formatStyle(
          "p_value",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.1), c("#d4edda", "#fff3cd", "white")
          )
        )
      }

      if ("p_adjusted_cross_metric" %in% colnames(tbl)) {
        dt <- dt %>% DT::formatStyle(
          "p_adjusted_cross_metric",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.25), c("#d4edda", "#fff3cd", "white")
          )
        )
      } else if ("p_adjusted" %in% colnames(tbl)) {
        dt <- dt %>% DT::formatStyle(
          "p_adjusted",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.25), c("#d4edda", "#fff3cd", "white")
          )
        )
      }

      dt
    })

    # -- Enriched metadata table ----------------------------------------
    output$meta_table <- DT::renderDT({
      shiny::req(rv$meta_enriched)
      DT::datatable(
        rv$meta_enriched,
        rownames = FALSE,
        filter   = "top",
        options  = list(
          scrollX    = TRUE,
          pageLength = 15
        )
      )
    })

    # -- Download alpha diversity TIFF ----------------------------------
    output$dl_alpha <- shiny::downloadHandler(
      filename = function() {
        sprintf(
          "alpha_%s_%s_by_%s_%s.tiff",
          input$alpha_metric,
          gsub("[^A-Za-z0-9]", "_", input$alpha_x_var),
          gsub("[^A-Za-z0-9]", "_", input$alpha_color_var),
          format(Sys.Date(), "%Y%m%d")
        )
      },
      content = function(file) {
        shiny::req(rv$meta_enriched, rv$obj)
        cfg <- rv$obj$config

        export_alpha_tiff(
          meta_enriched    = rv$meta_enriched,
          outfile          = file,
          cfg              = cfg,
          metric           = input$alpha_metric,
          x_var            = input$alpha_x_var,
          color_var        = input$alpha_color_var,
          palette          = input$palette,
          rev_pal          = input$rev_pal,
          connect_by       = connect_col(),
          exclude_baseline = input$alpha_exclude_baseline,
          is_baseline_col  = cfg$is_baseline_col,
          facet_var        = NULL
        )
      }
    )

    # -- Download alpha diversity LMM results ----------------------------
    output$dl_alpha_lmm <- shiny::downloadHandler(
      filename = function() {
        paste0("alpha_lmm_results_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        shiny::req(alpha_lmm_result())
        write.csv(alpha_lmm_result(), file, row.names = FALSE)
      }
    )

    # -- Download beta diversity stat results ----------------------------
    output$dl_stat <- shiny::downloadHandler(
      filename = function() {
        p <- stat_store$params
        if (is.null(p)) return("beta_diversity_stats.csv")
        sprintf(
          "beta_%s_%s_%s_%s.csv",
          p$method,
          p$scope,
          gsub("[^A-Za-z0-9]", "_", p$var),
          format(Sys.Date(), "%Y%m%d")
        )
      },
      content = function(file) {
        shiny::req(stat_store$result)
        write.csv(stat_store$result, file, row.names = FALSE)
      }
    )

    # -- Download enriched metadata -------------------------------------
    output$dl_meta <- shiny::downloadHandler(
      filename = function() {
        paste0("enriched_metadata_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        shiny::req(rv$meta_enriched)
        write.csv(rv$meta_enriched, file, row.names = TRUE)
      }
    )

    # -- Taxonomy bar plot data (computed on button click) --------------
    taxon_data_r <- shiny::eventReactive(input$run_taxon, {
      shiny::req(rv$obj, input$taxon_level, input$taxon_group_var)
      cfg <- rv$obj$config

      shiny::withProgress(
        message = sprintf("Parsing %s abundances...", input$taxon_level),
        value   = 0.4, {

          tryCatch(
            compute_taxonomy_barplot_data(
              otu_path  = rv$obj$otu_path,    # see note below
              meta_path = rv$obj$meta_path,
              level     = input$taxon_level,
              group_var = input$taxon_group_var,
              top_n     = input$taxon_top_n,
              config    = cfg
            ),
            error = function(e) {
              shiny::showNotification(
                paste("Taxonomy plot failed:", e$message),
                type = "error"
              )
              NULL
            }
          )
        }
      )
    })

    # -- Taxonomy plot render -------------------------------------------
    output$taxon_plot <- shiny::renderPlot({
      shiny::req(taxon_data_r(), rv$obj, input$taxon_group_var)
      cfg <- rv$obj$config
      td  <- taxon_data_r()

      # taxon_data_r() only recomputes on "Build taxonomy plot" clicks.
      # If the level or Top N selectors have changed since the last
      # click, the cached data no longer matches what's shown in the
      # controls - warn instead of silently plotting stale data under
      # a title that would describe the *new* (unapplied) settings.
      shiny::validate(
        shiny::need(
          identical(attr(td, "level_used"), input$taxon_level) &&
            identical(attr(td, "top_n_used"), input$taxon_top_n),
          paste(
            "Taxonomic level or Top N has changed since the plot was",
            "last built. Click \"Build taxonomy plot\" to apply the",
            "new settings."
          )
        )
      )

      build_taxon_barplot(
        taxon_data        = td,
        level             = input$taxon_level,
        group_var         = input$taxon_group_var,
        palette           = input$palette,
        rev_pal           = isTRUE(input$rev_pal),
        top_n             = input$taxon_top_n,
        exclude_baseline  = isTRUE(input$taxon_exclude_baseline),
        is_baseline_col   = cfg$is_baseline_col,
        x_label_var       = input$taxon_x_label_var,
        x_order_var       = input$taxon_x_order_var
      )
    })

    # -- Download taxonomy TIFF -----------------------------------------
    output$dl_taxon <- shiny::downloadHandler(
      filename = function() {
        sprintf(
          "taxonomy_%s_by_%s_%s.tiff",
          input$taxon_level,
          gsub("[^A-Za-z0-9]", "_", input$taxon_group_var),
          format(Sys.Date(), "%Y%m%d")
        )
      },
      content = function(file) {
        shiny::req(taxon_data_r(), rv$obj)
        cfg <- rv$obj$config
        td  <- taxon_data_r()

        if (!identical(attr(td, "level_used"), input$taxon_level) ||
            !identical(attr(td, "top_n_used"), input$taxon_top_n)) {
          shiny::showNotification(
            paste(
              "Taxonomic level or Top N has changed since the plot was",
              "last built. Click \"Build taxonomy plot\" before downloading."
            ),
            type = "error"
          )
          shiny::req(FALSE)
        }

        export_taxon_barplot_tiff(
          taxon_data        = td,
          outfile           = file,
          level             = input$taxon_level,
          group_var         = input$taxon_group_var,
          palette           = input$palette,
          rev_pal           = isTRUE(input$rev_pal),
          top_n             = input$taxon_top_n,
          exclude_baseline  = isTRUE(input$taxon_exclude_baseline),
          is_baseline_col   = cfg$is_baseline_col,
          x_label_var       = input$taxon_x_label_var,
          x_order_var       = input$taxon_x_order_var,
          width             = 14,
          height            = 8,
          dpi               = 300
        )
      }
    )

    # -- MaAsLin3 -----------------------------------------------------
    shiny::observeEvent(input$run_maaslin_btn, {
      shiny::req(rv$obj, rv$meta_enriched, input$maaslin_formula)

      cfg <- rv$obj$config

      tryCatch({

        # Validate the fixed-effects portion up front for a clear error
        # before ever touching the OTU table / calling MaAsLin3.
        validate_maaslin_formula(input$maaslin_formula, rv$meta_enriched)

        group_var <- if (nchar(trimws(input$maaslin_group_var)) > 0)
          input$maaslin_group_var else NULL

        prepped <- NULL

        res <- shiny::withProgress(
          message = "Running MaAsLin3 (this can take a while)...",
          value = 0.2, {

            prepped <- prepare_maaslin_input(
              otu_path       = rv$obj$otu_path,
              meta_enriched  = rv$meta_enriched,
              level          = input$maaslin_level,
              config         = cfg,
              min_prevalence = input$maaslin_min_prevalence,
              min_abundance  = input$maaslin_min_abundance
            )

            shiny::incProgress(0.3, detail = "Fitting models...")

            run_maaslin(
              features              = prepped$features,
              metadata              = prepped$metadata,
              formula_str           = input$maaslin_formula,
              group_effects_var     = group_var,
              random_effect         = input$maaslin_random_effect,
              use_random            = isTRUE(input$maaslin_use_random),
              small_random_effects  = isTRUE(input$maaslin_small_random),
              normalization         = input$maaslin_normalization,
              transform             = input$maaslin_transform,
              min_prevalence        = input$maaslin_min_prevalence,
              min_abundance         = input$maaslin_min_abundance,
              padj_method           = input$maaslin_padj_method,
              max_significance      = input$maaslin_max_sig
            )
          }
        )

        maaslin_store$result <- res
        maaslin_store$params <- list(
          level             = input$maaslin_level,
          formula_used      = attr(res, "formula_used"),
          small_random      = isTRUE(input$maaslin_small_random),
          normalization     = input$maaslin_normalization,
          transform         = input$maaslin_transform,
          min_prevalence    = input$maaslin_min_prevalence,
          min_abundance     = input$maaslin_min_abundance,
          padj_method       = input$maaslin_padj_method,
          max_significance  = input$maaslin_max_sig,
          n_features        = ncol(prepped$features),
          n_samples         = nrow(prepped$features),
          otu_scale_info    = prepped$otu_scale_info
        )

        shiny::showNotification("MaAsLin3 finished!", type = "message")

      }, error = function(e) {
        shiny::showNotification(
          paste("MaAsLin3 failed:", e$message),
          type = "error", duration = NULL
        )
      })
    })

    # -- MaAsLin3 info bar (reads FROZEN params) -------------------------
    output$maaslin_info <- shiny::renderText({
      shiny::req(maaslin_store$params)
      p <- maaslin_store$params

      scale_msg <- .format_scale_info(p$otu_scale_info)
      if (is.null(scale_msg)) scale_msg <- "(not available)"

      sprintf(
        paste0(
          "Formula            : %s\n",
          "Taxonomic level    : %s\n",
          "Small random fx    : %s\n",
          "Normalization      : %s\n",
          "Transform          : %s\n",
          "Min prevalence     : %s\n",
          "Min abundance      : %s\n",
          "p-adjust           : %s\n",
          "Max significance   : %s\n",
          "Features x Samples : %d x %d\n",
          "Input scale        : %s"
        ),
        p$formula_used, p$level, p$small_random,
        p$normalization, p$transform,
        p$min_prevalence, p$min_abundance,
        p$padj_method, p$max_significance,
        p$n_features, p$n_samples,
        scale_msg
      )
    })

    # -- MaAsLin3 full results table --------------------------------------
    output$maaslin_table <- DT::renderDT({
      shiny::req(maaslin_store$result)
      tbl <- maaslin_store$result

      dt <- DT::datatable(
        tbl,
        rownames = FALSE,
        filter   = "top",
        options  = list(
          scrollX    = TRUE,
          pageLength = 15,
          dom        = "tip"
        )
      ) %>%
        DT::formatRound(
          columns = intersect(
            c("coef", "stderr", "pval", "qval"),
            colnames(tbl)
          ),
          digits = 4
        )

      if ("qval" %in% colnames(tbl)) {
        dt <- dt %>% DT::formatStyle(
          "qval",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.25), c("#d4edda", "#fff3cd", "white")
          )
        )
      }

      dt
    })

    # -- MaAsLin3 significant-hits table (uses FROZEN threshold) ---------
    output$maaslin_sig_table <- DT::renderDT({
      shiny::req(maaslin_store$result, maaslin_store$params)

      sig <- summarise_maaslin_results(
        maaslin_store$result,
        threshold = maaslin_store$params$max_significance
      )

      DT::datatable(
        sig,
        rownames = FALSE,
        options  = list(
          scrollX    = TRUE,
          pageLength = 15,
          dom        = "tip"
        )
      )
    })

    # -- Download MaAsLin3 results ----------------------------------------
    output$dl_maaslin <- shiny::downloadHandler(
      filename = function() {
        paste0("maaslin3_results_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        shiny::req(maaslin_store$result)
        write.csv(maaslin_store$result, file, row.names = FALSE)
      }
    )

    # -- Download 2D TIFF -----------------------------------------------
    output$dl_2d <- shiny::downloadHandler(
      filename = function() {
        sprintf(
          "pcoa_2D_%s_PC%s-PC%s.tiff",
          gsub("[^A-Za-z0-9]", "_", input$color_var),
          input$pc_x, input$pc_y
        )
      },
      content = function(file) {
        export_2d_tiff(
          rv$obj,
          v          = input$color_var,
          outfile    = file,
          ax_x       = as.integer(input$pc_x),
          ax_y       = as.integer(input$pc_y),
          palette    = input$palette,
          is_cont    = input$is_cont,
          rev_pal    = input$rev_pal,
          connect_by = connect_col()
        )
      }
    )

    # -- Download 3D TIFF -----------------------------------------------
    output$dl_3d <- shiny::downloadHandler(
      filename = function() {
        sprintf(
          "pcoa_3D_%s.tiff",
          gsub("[^A-Za-z0-9]", "_", input$color_var)
        )
      },
      content = function(file) {
        export_3d_tiff(
          rv$obj,
          v          = input$color_var,
          outfile    = file,
          palette    = input$palette,
          is_cont    = input$is_cont,
          rev_pal    = input$rev_pal,
          connect_by = connect_col()
        )
      }
    )

  }  # end server

  shiny::shinyApp(ui = ui, server = server)

}  # end run_app
