# export.R
# -----------------------------------------------------------------------------
# Export functions for publication-quality TIFF output:
#   1. export_2d_tiff()     - 2D PCoA ggplot -> TIFF
#   2. export_3d_tiff()     - 3D PCoA plotly -> TIFF (with fallback)
#   3. export_alpha_tiff()  - Alpha diversity plot -> TIFF
#
# All dimension and DPI settings are read from the pcoa_config() object
# (config.R), so changing config values propagates automatically.
#
# Compatible with:
#   config.R -> export_dpi, export_w_in, export_h_in,
#              export_3d_w_px, export_3d_h_px,
#              alpha_export_w_in, alpha_export_h_in
#   plots.R  -> build_2d(), build_3d(), build_alpha_scatter()
#   app.R    -> downloadHandler content functions
#   data.R   -> pcoa_obj$config, pcoa_obj$meta_enriched
#   mapping  -> CatID for connect_by
# -----------------------------------------------------------------------------

# -- 1. export_2d_tiff() -------------------------------------------------------

#' Export a 2D PCoA ggplot as a high-resolution TIFF
#'
#' Renders the plot via \code{build_2d()} from plots.R and saves it as
#' a LZW-compressed TIFF at the resolution specified in
#' \code{pcoa_config()} (config.R).
#'
#' For the crossover design, set \code{connect_by = "CatID"} to
#' draw within-subject paired lines between treatment observations.
#'
#' @param obj        A \code{pcoa_obj} from \code{prepare_pcoa()}.
#' @param v          Colour variable name (metadata column).
#' @param outfile    Output TIFF file path (provided by
#'   \code{downloadHandler}).
#' @param ax_x       Integer; PC axis for x (default 1).
#' @param ax_y       Integer; PC axis for y (default 2).
#' @param palette    RColorBrewer palette name (default \code{"Set2"}).
#' @param is_cont    Logical; use continuous colour scale?
#' @param rev_pal    Logical; reverse the palette?
#' @param connect_by Optional metadata column to draw within-subject
#'   connecting lines (e.g. \code{"CatID"} for crossover design).
#'   \code{NULL} = no lines (default).
#' @return Invisibly returns \code{outfile}.
#' @details Output is LZW-compressed TIFF (lossless). Dimensions and DPI
#'   are controlled by \code{export_w_in}, \code{export_h_in}, and
#'   \code{export_dpi} in \code{pcoa_config()} (config.R).
#' @importFrom ggplot2 ggsave
#' @export
export_2d_tiff <- function(obj, v, outfile,
                           ax_x       = 1,
                           ax_y       = 2,
                           palette    = "Set2",
                           is_cont    = FALSE,
                           rev_pal    = FALSE,
                           connect_by = NULL) {

  cfg <- obj$config

  # Build ggplot via plots.R - connect_by adds paired lines
  p <- build_2d(
    obj        = obj,
    v          = v,
    ax_x       = ax_x,
    ax_y       = ax_y,
    palette    = palette,
    is_cont    = is_cont,
    rev_pal    = rev_pal,
    connect_by = connect_by
  )

  ggplot2::ggsave(
    filename    = outfile,
    plot        = p,
    device      = "tiff",
    width       = cfg$export_w_in,
    height      = cfg$export_h_in,
    units       = "in",
    dpi         = cfg$export_dpi,
    compression = "lzw"          # lossless - required for publication
  )

  message(sprintf(
    "2D TIFF saved: %s (%.0f dpi, %g x %g in)",
    outfile, cfg$export_dpi, cfg$export_w_in, cfg$export_h_in
  ))

  invisible(outfile)
}

# -- 2. export_3d_tiff() -------------------------------------------------------

#' Export a 3D PCoA plotly figure as a high-resolution TIFF
#'
#' Attempts to use \code{plotly::save_image()} (requires the
#' \code{kaleido} Python package via \code{reticulate}) to render the
#' true 3D plotly figure. If \code{kaleido} is unavailable, the function
#' falls back gracefully to rendering the equivalent 2D (PC1  PC2) plot
#' via \code{build_2d()} and saves that as a TIFF with a caption noting
#' the fallback.
#'
#' @param obj        A \code{pcoa_obj} from \code{prepare_pcoa()}.
#' @param v          Colour variable name.
#' @param outfile    Output TIFF file path.
#' @param palette    RColorBrewer palette name.
#' @param is_cont    Logical; continuous colour scale?
#' @param rev_pal    Logical; reverse palette?
#' @param connect_by Optional column for within-subject lines
#'   (e.g. \code{"CatID"}). Passed to fallback \code{build_2d()}.
#' @return Invisibly returns \code{outfile}.
#' @details
#' \strong{Primary path:} Requires Python \code{kaleido} and
#'   \code{reticulate}. Install with:
#'   \preformatted{
#'     install.packages("reticulate")
#'     reticulate::install_miniconda()     # if needed
#'     reticulate::py_install("kaleido")
#'   }
#' \strong{Fallback path:} If \code{kaleido} is unavailable, the PC1  PC2
#'   ggplot is saved instead. A caption is added indicating this is a
#'   fallback. Users can also use the camera icon in the interactive
#'   Shiny panel for a PNG snapshot.
#'
#' Output dimensions are controlled by \code{export_3d_w_px} and
#' \code{export_3d_h_px} in \code{pcoa_config()} (config.R).
#' @importFrom ggplot2 ggsave labs theme element_text
#' @export
export_3d_tiff <- function(obj, v, outfile,
                           palette    = "Set2",
                           is_cont    = FALSE,
                           rev_pal    = FALSE,
                           connect_by = NULL) {

  cfg <- obj$config
  dpi <- cfg$export_dpi

  # Derived dimensions in inches from pixel config (config.R)
  w_in <- cfg$export_3d_w_px / dpi
  h_in <- cfg$export_3d_h_px / dpi

  # -- Try plotly/kaleido path ---------------------------------
  kaleido_ok <- requireNamespace("reticulate", quietly = TRUE) &&
    tryCatch({
      reticulate::import("kaleido")
      TRUE
    }, error = function(e) FALSE)

  if (kaleido_ok) {

    message("kaleido detected - rendering true 3D TIFF via plotly.")

    fig <- build_3d(
      obj        = obj,
      v          = v,
      palette    = palette,
      is_cont    = is_cont,
      rev_pal    = rev_pal,
      connect_by = connect_by
    )

    # Render to PNG first via kaleido
    tmp_png <- tempfile(fileext = ".png")

    tryCatch({
      plotly::save_image(
        fig,
        file   = tmp_png,
        width  = cfg$export_3d_w_px,
        height = cfg$export_3d_h_px,
        scale  = 1
      )

      # Convert PNG -> TIFF with correct DPI metadata via magick
      if (requireNamespace("magick", quietly = TRUE)) {
        img <- magick::image_read(tmp_png)
        img <- magick::image_set_density(img,
                                         paste0(dpi, "x", dpi))
        magick::image_write(
          img,
          path    = outfile,
          format  = "tiff"
        )
      } else {
        # magick unavailable - write PNG data directly as TIFF
        warning(paste0(
          "'magick' package not available. ",
          "TIFF DPI metadata may not be set correctly. ",
          "Install with: install.packages('magick')"
        ))
        file.copy(tmp_png, outfile, overwrite = TRUE)
      }

      if (file.exists(tmp_png)) file.remove(tmp_png)

      message(sprintf(
        "3D TIFF saved (kaleido): %s (%d x %d px, %.0f dpi)",
        outfile,
        cfg$export_3d_w_px,
        cfg$export_3d_h_px,
        dpi
      ))

    }, error = function(e) {
      warning(sprintf(
        "plotly::save_image() failed: %s\nFalling back to 2D export.",
        e$message
      ))
      if (file.exists(tmp_png)) file.remove(tmp_png)
      .export_3d_fallback(obj, v, outfile, palette, is_cont,
                          rev_pal, connect_by, cfg, w_in, h_in, dpi)
    })

  } else {
    # -- Fallback: render 2D equivalent --------------------------
    warning(paste0(
      "kaleido/reticulate not available - falling back to 2D TIFF export.\n",
      "To enable true 3D TIFF export:\n",
      "  install.packages('reticulate')\n",
      "  reticulate::py_install('kaleido')\n",
      "Alternatively, use the camera icon in the Shiny 3D panel for ",
      "an in-browser PNG snapshot."
    ))
    .export_3d_fallback(obj, v, outfile, palette, is_cont,
                        rev_pal, connect_by, cfg, w_in, h_in, dpi)
  }

  invisible(outfile)
}

# -- Internal fallback for export_3d_tiff() -----------------------------------
.export_3d_fallback <- function(obj, v, outfile,
                                palette, is_cont, rev_pal,
                                connect_by, cfg, w_in, h_in, dpi) {

  p <- build_2d(
    obj        = obj,
    v          = v,
    ax_x       = 1,
    ax_y       = 2,
    palette    = palette,
    is_cont    = is_cont,
    rev_pal    = rev_pal,
    connect_by = connect_by
  )

  # Add fallback caption so reader knows this is not the true 3D figure
  p <- p +
    ggplot2::labs(
      caption = paste0(
        "Note: True 3D export unavailable (kaleido not installed).\n",
        "Showing PC1 \u00d7 PC2 (2D). Use Shiny camera icon for 3D PNG."
      )
    ) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        colour = "grey50", size = 9, hjust = 0.5
      )
    )

  ggplot2::ggsave(
    filename    = outfile,
    plot        = p,
    device      = "tiff",
    width       = w_in,
    height      = h_in,
    units       = "in",
    dpi         = dpi,
    compression = "lzw"
  )

  message(sprintf(
    "3D TIFF (fallback 2D): %s (%.1f x %.1f in, %.0f dpi)",
    outfile, w_in, h_in, dpi
  ))
}

# -- 3. export_alpha_tiff() ----------------------------------------------------

#' Export an alpha diversity scatter/boxplot as a high-resolution TIFF
#'
#' Renders the plot via \code{build_alpha_scatter()} from plots.R.
#' Dimensions are controlled by \code{alpha_export_w_in} and
#' \code{alpha_export_h_in} in \code{pcoa_config()} (config.R).
#'
#' Designed for the crossover design:
#' \itemize{
#'   \item Set \code{connect_by = "CatID"} for within-cat paired lines.
#'   \item Set \code{exclude_baseline = TRUE} to remove
#'     \code{IsBaseline == TRUE} rows from the plot.
#'   \item Set \code{facet_var = "Timepoint"} to show
#'     Phase 1 and Phase 2 in separate panels.
#' }
#'
#' @param meta_enriched  Data frame from \code{pcoa_obj$meta_enriched}
#'   (produced by \code{prepare_all()} in data.R).
#' @param outfile        Output TIFF file path.
#' @param cfg            A \code{pcoa_config()} object (config.R).
#'   Used for DPI and export dimensions.
#' @param metric         Alpha metric column to plot
#'   (default \code{"Shannon"}).
#' @param x_var          X axis grouping column
#'   (default \code{"Treatment"}).
#' @param color_var      Point colour column
#'   (default \code{"Timepoint"}).
#' @param palette        RColorBrewer categorical palette.
#' @param rev_pal        Logical; reverse palette?
#' @param connect_by     Column for within-subject connecting lines.
#'   Default \code{"CatID"} for crossover design.
#' @param exclude_baseline Logical; exclude baseline rows (as flagged
#'   by \code{is_baseline_col})? Default \code{TRUE}.
#' @param is_baseline_col Metadata column flagging baseline rows
#'   (default \code{"IsBaseline"}). See \code{\link{pcoa_config}}.
#' @param facet_var      Optional column to facet the plot by
#'   (e.g. \code{"Timepoint"} for separate phase panels).
#'   \code{NULL} = no faceting (default).
#' @return Invisibly returns \code{outfile}.
#' @importFrom ggplot2 ggsave
#' @export
export_alpha_tiff <- function(meta_enriched,
                              outfile,
                              cfg              = pcoa_config(),
                              metric           = "Shannon",
                              x_var            = "Treatment",
                              color_var        = "Timepoint",
                              palette          = "Set2",
                              rev_pal          = FALSE,
                              connect_by       = "CatID",
                              exclude_baseline = TRUE,
                              is_baseline_col  = "IsBaseline",
                              facet_var        = NULL) {

  # -- Validate inputs -------------------------------------------
  if (!inherits(cfg, "pcoa_config"))
    stop("'cfg' must be a pcoa_config object from pcoa_config().")

  if (!metric %in% colnames(meta_enriched))
    stop(sprintf(
      "'%s' not found in meta_enriched. Available: %s",
      metric,
      paste(intersect(ALPHA_METRICS, colnames(meta_enriched)),
            collapse = ", ")
    ))

  # -- Build plot via plots.R ------------------------------------
  p <- build_alpha_scatter(
    meta_enriched    = meta_enriched,
    metric           = metric,
    x_var            = x_var,
    color_var        = color_var,
    palette          = palette,
    rev_pal          = rev_pal,
    connect_by       = connect_by,
    exclude_baseline = exclude_baseline,
    is_baseline_col  = is_baseline_col,
    facet_var        = facet_var
  )

  # -- Save as TIFF ----------------------------------------------
  # Dimensions from alpha_export_w_in / alpha_export_h_in (config.R)
  ggplot2::ggsave(
    filename    = outfile,
    plot        = p,
    device      = "tiff",
    width       = cfg$alpha_export_w_in,
    height      = cfg$alpha_export_h_in,
    units       = "in",
    dpi         = cfg$export_dpi,
    compression = "lzw"
  )

  message(sprintf(
    "Alpha TIFF saved: %s | metric: %s | x: %s | colour: %s\n  (%.0f dpi, %g x %g in)",
    outfile, metric, x_var, color_var,
    cfg$export_dpi, cfg$alpha_export_w_in, cfg$alpha_export_h_in
  ))

  invisible(outfile)
}

# -- export_taxon_barplot_tiff() -----------------------------------------------
#' Save a taxonomy bar plot to a TIFF file
#'
#' @param taxon_data   Data frame from \code{compute_taxonomy_barplot_data()}.
#' @param outfile      Output file path.
#' @param level        Taxonomic level label.
#' @param group_var    Metadata column used for faceting.
#' @param palette      Colour palette name.
#' @param rev_pal      Logical - reverse palette?
#' @param top_n        Number of top taxa.
#' @param exclude_baseline  Logical - remove baseline samples?
#' @param is_baseline_col   Baseline flag column name.
#' @param width,height Plot dimensions in inches (defaults: 14  8).
#' @param dpi          Resolution (default 300).
#' @param x_label_var  Metadata column shown as x-axis tick labels
#'   (default \code{"SampleID"}). See \code{\link{build_taxon_barplot}}.
#' @param x_order_var  Metadata column used to order bars along the
#'   x-axis (default \code{"SampleID"}). See
#'   \code{\link{build_taxon_barplot}}.
#' @export
export_taxon_barplot_tiff <- function(taxon_data,
                                      outfile,
                                      level             = "Family",
                                      group_var,
                                      palette           = "Set3",
                                      rev_pal           = FALSE,
                                      top_n             = TAXA_DEFAULT_TOP_N,
                                      exclude_baseline  = TRUE,
                                      is_baseline_col   = "IsBaseline",
                                      x_label_var       = "SampleID",
                                      x_order_var       = "SampleID",
                                      width             = 14,
                                      height            = 8,
                                      dpi               = 300) {

  p <- build_taxon_barplot(
    taxon_data        = taxon_data,
    level             = level,
    group_var         = group_var,
    palette           = palette,
    rev_pal           = rev_pal,
    top_n             = top_n,
    exclude_baseline  = exclude_baseline,
    is_baseline_col   = is_baseline_col,
    x_label_var       = x_label_var,
    x_order_var       = x_order_var
  )

  ggplot2::ggsave(
    filename = outfile,
    plot     = p,
    device   = "tiff",
    width    = width,
    height   = height,
    dpi      = dpi,
    compression = "lzw"
  )

  invisible(outfile)
}
