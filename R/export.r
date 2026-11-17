# export.R
# -----------------------------------------------------------------------------
# Publication TIFF export (LZW-compressed). Sizes and DPI come from
# pcoa_config(); every function builds its plot with plots.R and writes it
# with the shared .save_tiff() helper (utils.R).
# -----------------------------------------------------------------------------

#' Export a 2D PCoA plot as a TIFF
#'
#' @inheritParams build_2d
#' @param outfile Output file path.
#' @return Invisibly, \code{outfile}.
#' @export
export_2d_tiff <- function(obj, v, outfile, ax_x = 1, ax_y = 2, palette = "Set2",
                           is_cont = FALSE, rev_pal = FALSE, connect_by = NULL) {
  cfg <- obj$config
  .save_tiff(build_2d(obj, v, ax_x, ax_y, palette, is_cont, rev_pal, connect_by),
             outfile, cfg$export_w_in, cfg$export_h_in, cfg$export_dpi)
}

#' Export the 3D PCoA plot as a TIFF
#'
#' Renders the real 3D figure when Python \code{kaleido} (via
#' \code{reticulate}) and the \code{magick} package are available. Otherwise
#' it saves the PC1 x PC2 plot with a caption saying so. The camera icon in
#' the app always gives a 3D PNG snapshot.
#'
#' @inheritParams export_2d_tiff
#' @return Invisibly, \code{outfile}.
#' @export
export_3d_tiff <- function(obj, v, outfile, palette = "Set2", is_cont = FALSE,
                           rev_pal = FALSE, connect_by = NULL) {
  cfg <- obj$config
  dpi <- cfg$export_dpi
  can_3d <- requireNamespace("magick", quietly = TRUE) &&
    requireNamespace("reticulate", quietly = TRUE) &&
    tryCatch({ reticulate::import("kaleido"); TRUE }, error = function(e) FALSE)

  if (can_3d) {
    ok <- tryCatch({
      tmp <- tempfile(fileext = ".png"); on.exit(unlink(tmp), add = TRUE)
      plotly::save_image(build_3d(obj, v, palette, is_cont, rev_pal, connect_by),
                         file = tmp, width = cfg$export_3d_w_px,
                         height = cfg$export_3d_h_px, scale = 1)
      magick::image_write(magick::image_read(tmp), path = outfile, format = "tiff",
                          density = paste0(dpi, "x", dpi))
      TRUE
    }, error = function(e) {
      warning("3D render failed (", conditionMessage(e), "); saving 2D instead.",
              call. = FALSE)
      FALSE
    })
    if (ok) return(invisible(outfile))
  } else {
    warning("True 3D TIFF needs 'magick', 'reticulate' and Python 'kaleido'; ",
            "saving the PC1 x PC2 plot instead.", call. = FALSE)
  }

  p <- build_2d(obj, v, 1, 2, palette, is_cont, rev_pal, connect_by) +
    ggplot2::labs(caption = paste0("Note: 3D export unavailable - showing PC1 \u00d7 PC2. ",
                                   "Use the camera icon in the app for a 3D PNG.")) +
    ggplot2::theme(plot.caption = ggplot2::element_text(colour = "grey50", size = 9,
                                                        hjust = 0.5))
  .save_tiff(p, outfile, cfg$export_3d_w_px / dpi, cfg$export_3d_h_px / dpi, dpi)
}

#' Export an alpha diversity plot as a TIFF
#'
#' @inheritParams build_alpha_scatter
#' @param outfile Output file path.
#' @param cfg A \code{pcoa_config()} object (size and DPI).
#' @return Invisibly, \code{outfile}.
#' @export
export_alpha_tiff <- function(meta_enriched, outfile, cfg = pcoa_config(),
                              metric = "Shannon", x_var = "Treatment",
                              color_var = "Timepoint", palette = "Set2",
                              rev_pal = FALSE, connect_by = "CatID",
                              exclude_baseline = TRUE,
                              is_baseline_col = "IsBaseline", facet_var = NULL) {
  if (!inherits(cfg, "pcoa_config"))
    stop("'cfg' must be a pcoa_config object.", call. = FALSE)
  p <- build_alpha_scatter(meta_enriched, metric, x_var, color_var, palette, rev_pal,
                           connect_by, exclude_baseline, is_baseline_col, facet_var)
  .save_tiff(p, outfile, cfg$alpha_export_w_in, cfg$alpha_export_h_in, cfg$export_dpi)
}

#' Export a taxonomy bar plot as a TIFF
#'
#' @inheritParams build_taxon_barplot
#' @param outfile Output file path.
#' @param width,height Size in inches.
#' @param dpi Resolution.
#' @return Invisibly, \code{outfile}.
#' @export
export_taxon_barplot_tiff <- function(taxon_data, outfile, level = "Family", group_var,
                                      palette = "Set3", rev_pal = FALSE,
                                      top_n = TAXA_DEFAULT_TOP_N,
                                      exclude_baseline = TRUE,
                                      is_baseline_col = "IsBaseline",
                                      x_label_var = "SampleID",
                                      x_order_var = "SampleID",
                                      width = 14, height = 8, dpi = 300) {
  p <- build_taxon_barplot(taxon_data, level, group_var, palette, rev_pal, top_n,
                           exclude_baseline, is_baseline_col, x_label_var, x_order_var)
  .save_tiff(p, outfile, width, height, dpi)
}
