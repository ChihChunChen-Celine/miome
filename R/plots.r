# plots.R
# ─────────────────────────────────────────────────────────────────────────────
# Plotting functions for:
#   - PCoA 3D interactive (plotly)
#   - PCoA 2D static (ggplot2)
#   - Alpha diversity scatter/boxplot (ggplot2)
#
# TIFF export wrappers (export_2d_tiff, export_3d_tiff, export_alpha_tiff)
# live ONLY in export.R — they are not duplicated here.
#
# Compatible with:
#   config.R  → pcoa_config() fields: dist_method_label, marker_size,
#               marker_opacity, export_*, subject_id_col
#   data.R    → pcoa_obj$pcoa_df, $var_pct, $config, $meta_enriched
#   app.R     → build_3d(), build_2d(), build_alpha_scatter()
#   export.R  → export_2d_tiff(), export_3d_tiff(), export_alpha_tiff()
#               wrap build_2d()/build_3d()/build_alpha_scatter() from here
# ─────────────────────────────────────────────────────────────────────────────

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices colorRampPalette
get_cat_pal <- function(pal, n, reverse = FALSE) {
  max_cols <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
  cols <- if (n <= max_cols)
    RColorBrewer::brewer.pal(max(3, n), pal)[seq_len(n)]
  else
    grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(max_cols, pal))(n)
  if (reverse) rev(cols) else cols
}

make_colorscale <- function(pal_name, reverse = FALSE, n_steps = 11) {
  # Some CONT_PALETTES (e.g. "Blues", "YlOrRd") max out below 11 colours —
  # brewer.pal(11, pal_name) errors in that case. Request at most the
  # palette's real maximum, then ramp it up to n_steps.
  max_cols <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
  base_n   <- min(11, max_cols)
  hex <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(base_n, pal_name))(n_steps)
  if (reverse) hex <- rev(hex)
  pos <- seq(0, 1, length.out = n_steps)
  lapply(seq_along(hex), function(i) list(pos[i], hex[i]))
}

# Axis label: "PC1 (23.4%)"
.ax_lab <- function(var_pct, i) {
  sprintf("PC%d (%.1f%%)", i, var_pct[i])
}

# ── 1. build_3d() ─────────────────────────────────────────────────────────────

#' Build a 3D interactive PCoA plotly figure
#'
#' @param obj     A \code{pcoa_obj} from \code{prepare_pcoa()}.
#' @param v       Name of the colour variable (metadata column).
#' @param palette RColorBrewer palette name.
#' @param is_cont Logical; treat \code{v} as continuous?
#' @param rev_pal Logical; reverse the palette?
#' @param connect_by Optional metadata column to draw connecting lines
#'   between paired samples (e.g. \code{"CatID"} for crossover design).
#'   Only applies to categorical mode.
#' @return A \code{plotly} object.
#' @importFrom plotly plot_ly add_trace layout config
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @export
build_3d <- function(obj, v,
                     palette    = "Set2",
                     is_cont    = FALSE,
                     rev_pal    = FALSE,
                     connect_by = NULL) {

  cfg <- obj$config
  df  <- obj$pcoa_df[!is.na(obj$pcoa_df[[v]]), ]

  # Hover text
  make_hover <- function(sub_df) {
    paste0("<b>", sub_df$SampleID, "</b><br>", v, ": ", sub_df[[v]])
  }

  fig <- plotly::plot_ly()

  if (!is_cont) {
    # ── Categorical ───────────────────────────────────────────────
    lvls <- sort(unique(as.character(df[[v]])))
    cols <- stats::setNames(
      get_cat_pal(palette, length(lvls), rev_pal), lvls)

    for (lvl in lvls) {
      sub <- df[as.character(df[[v]]) == lvl, ]
      fig <- plotly::add_trace(
        fig,
        type       = "scatter3d",
        mode       = "markers",
        x = sub$PC1, y = sub$PC2, z = sub$PC3,
        name       = lvl,
        text       = make_hover(sub),
        hoverinfo  = "text",
        marker     = list(
          size    = cfg$marker_size,
          color   = cols[[lvl]],
          opacity = cfg$marker_opacity,
          line    = list(width = 0.5, color = "white")
        )
      )
    }

    # Optional: paired connecting lines (crossover design)
    if (!is.null(connect_by) && connect_by %in% colnames(df)) {
      subjects <- unique(df[[connect_by]])
      for (subj in subjects) {
        s <- df[df[[connect_by]] == subj, ]
        if (nrow(s) >= 2) {
          fig <- plotly::add_trace(
            fig,
            type      = "scatter3d",
            mode      = "lines",
            x = s$PC1, y = s$PC2, z = s$PC3,
            name      = subj,
            showlegend = FALSE,
            hoverinfo  = "skip",
            line       = list(color = "grey70", width = 1.5)
          )
        }
      }
    }

  } else {
    # ── Continuous ────────────────────────────────────────────────
    vals <- as.numeric(df[[v]])
    fig  <- plotly::add_trace(
      fig,
      type       = "scatter3d",
      mode       = "markers",
      x = df$PC1, y = df$PC2, z = df$PC3,
      name       = v,
      text       = make_hover(df),
      hoverinfo  = "text",
      showlegend = FALSE,
      marker     = list(
        size       = cfg$marker_size,
        opacity    = cfg$marker_opacity,
        color      = vals,
        colorscale = make_colorscale(palette, rev_pal),
        cmin       = min(vals, na.rm = TRUE),
        cmax       = max(vals, na.rm = TRUE),
        showscale  = TRUE,
        colorbar   = list(
          title     = list(text = v, side = "right"),
          thickness = 15, len = 0.5, x = 1.02
        ),
        line       = list(width = 0.5, color = "white")
      )
    )
  }

  # ── Layout ───────────────────────────────────────────────────────
  fig <- plotly::layout(
    fig,
    title  = list(
      # FIX [1]: was cfg$dist_method — now correctly cfg$dist_method_label
      text = paste0(cfg$plot_title,
                    " | ", cfg$dist_method_label,
                    " (n=", nrow(df), ")"),
      font = list(size = 14)
    ),
    scene  = list(
      xaxis  = list(title = .ax_lab(obj$var_pct, 1),
                    showgrid = TRUE, zeroline = FALSE),
      yaxis  = list(title = .ax_lab(obj$var_pct, 2),
                    showgrid = TRUE, zeroline = FALSE),
      zaxis  = list(title = .ax_lab(obj$var_pct, 3),
                    showgrid = TRUE, zeroline = FALSE),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.0))
    ),
    legend = list(
      title      = list(text = paste0("<b>", v, "</b>")),
      itemsizing = "constant"
    ),
    margin = list(t = 60, b = 10, l = 10, r = 10)
  )

  fig %>% plotly::config(
    toImageButtonOptions = list(
      format   = "png",
      filename = "pcoa_3D_snapshot",
      width    = cfg$export_3d_w_px,
      height   = cfg$export_3d_h_px,
      scale    = 1
    ),
    displaylogo = FALSE
  )
}

# ── 2. build_2d() ─────────────────────────────────────────────────────────────

#' Build a 2D PCoA ggplot
#'
#' @param obj        A \code{pcoa_obj} from \code{prepare_pcoa()}.
#' @param v          Name of the colour variable.
#' @param ax_x       Integer; PC axis for x (default 1).
#' @param ax_y       Integer; PC axis for y (default 2).
#' @param palette    RColorBrewer palette name.
#' @param is_cont    Logical; treat \code{v} as continuous?
#' @param rev_pal    Logical; reverse the palette?
#' @param connect_by Optional metadata column for paired connecting lines.
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom stats setNames
#' @export
build_2d <- function(obj, v, ax_x = 1, ax_y = 2,
                     palette    = "Set2",
                     is_cont    = FALSE,
                     rev_pal    = FALSE,
                     connect_by = NULL) {

  cfg  <- obj$config
  ax_x <- as.integer(ax_x)
  ax_y <- as.integer(ax_y)
  pc_x <- paste0("PC", ax_x)
  pc_y <- paste0("PC", ax_y)

  df <- obj$pcoa_df[!is.na(obj$pcoa_df[[v]]), ]

  # Subtitle grouping text
  if (is_cont) {
    grp_txt <- paste0(
      "range ",
      round(min(as.numeric(df[[v]]), na.rm = TRUE), 1), "\u2013",
      round(max(as.numeric(df[[v]]), na.rm = TRUE), 1)
    )
  } else {
    grp_txt <- paste(sort(unique(as.character(df[[v]]))), collapse = ", ")
  }

  title_txt    <- sprintf("PCoA \u2014 coloured by: %s", v)
  # FIX [1]: was cfg$dist_method — now correctly cfg$dist_method_label
  subtitle_txt <- sprintf("%s | %s | n = %d",
                          cfg$dist_method_label, grp_txt, nrow(df))

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[pc_x]], y = .data[[pc_y]])
  )

  # Optional paired lines — drawn BEFORE points so lines go under dots
  if (!is_cont && !is.null(connect_by) &&
      connect_by %in% colnames(df)) {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(group = .data[[connect_by]]),
        colour    = "grey75",
        linewidth = 0.4,
        alpha     = 0.7
      )
  }

  if (is_cont) {
    # Cap the request at the palette's real max (e.g. Blues/YlOrRd max
    # out at 9), then ramp up — mirrors make_colorscale()'s fix for
    # build_3d() so build_2d() no longer crashes on those palettes.
    max_cols  <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    ramp_cols <- RColorBrewer::brewer.pal(min(11, max_cols), palette)
    if (rev_pal) ramp_cols <- rev(ramp_cols)
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = as.numeric(.data[[v]])),
        size  = cfg$marker_size / 2,
        alpha = cfg$marker_opacity
      ) +
      ggplot2::scale_color_gradientn(colours = ramp_cols, name = v)

  } else {
    lvls     <- sort(unique(as.character(df[[v]])))
    cat_cols <- stats::setNames(
      get_cat_pal(palette, length(lvls), rev_pal), lvls)
    df[[v]]  <- factor(as.character(df[[v]]), levels = lvls)

    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data[[v]]),
        size   = cfg$marker_size / 2,
        alpha  = cfg$marker_opacity,
        shape  = 21,
        colour = "white",
        stroke = 0.4
      ) +
      ggplot2::scale_fill_manual(values = cat_cols, name = v)
  }

  p +
    ggplot2::labs(
      title    = title_txt,
      subtitle = subtitle_txt,
      x        = .ax_lab(obj$var_pct, ax_x),
      y        = .ax_lab(obj$var_pct, ax_y)
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face  = "bold",
                                               hjust = 0.5, size = 18),
      plot.subtitle    = ggplot2::element_text(hjust  = 0.5,
                                               colour = "grey40", size = 12),
      axis.title       = ggplot2::element_text(face = "bold"),
      legend.title     = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ── 3. build_alpha_scatter() ──────────────────────────────────────────────────

#' Build an alpha diversity boxplot + jitter scatter with optional
#' paired lines and grouping options
#'
#' Designed to work directly with the \code{meta_enriched} data frame
#' produced by \code{prepare_all()} in data.R.
#'
#' @param meta_enriched Data frame from \code{pcoa_obj$meta_enriched}.
#'   Must contain alpha diversity columns (Shannon, Simpson, etc.) (data.R).
#' @param metric        Alpha diversity column to plot (e.g. \code{"Shannon"}).
#' @param x_var         Metadata column for x axis grouping
#'   (e.g. \code{"Treatment"}).
#' @param color_var     Metadata column for point colour
#'   (e.g. \code{"Timepoint"}).
#' @param palette       RColorBrewer categorical palette.
#' @param rev_pal       Logical; reverse palette?
#' @param connect_by    Optional column to draw within-subject lines
#'   (e.g. \code{"CatID"} for crossover). Set \code{NULL} to disable.
#' @param exclude_baseline Logical; if \code{TRUE} (default) and the
#'   \code{is_baseline_col} column exists, baseline rows are removed
#'   before plotting.
#' @param is_baseline_col Metadata column flagging baseline rows
#'   (default \code{"IsBaseline"}). See \code{\link{pcoa_config}}.
#'   #' @param facet_var     Optional column to facet by
#'   (e.g. \code{"Timepoint"} for separate Phase panels).
#' @param marker_size   Point size (default 3).
#' @param marker_alpha  Point opacity (default 0.85).
#' @return A \code{ggplot} object. Can be passed to \code{export_alpha_tiff()}.
#' @import ggplot2
#' @importFrom stats setNames
#' @param facet_var Optional metadata column to facet by (or NULL).
#' @export
build_alpha_scatter <- function(meta_enriched,
                                metric             = "Shannon",
                                x_var              = "Treatment",
                                color_var          = "Timepoint",
                                palette            = "Set2",
                                rev_pal            = FALSE,
                                connect_by         = "CatID",
                                exclude_baseline   = TRUE,
                                is_baseline_col    = "IsBaseline",
                                facet_var          = NULL,
                                marker_size        = 3,
                                marker_alpha       = 0.85) {

  # ── Validate inputs ───────────────────────────────────────────
  if (!metric %in% colnames(meta_enriched))
    stop(sprintf(
      "'%s' not found in meta_enriched. Available metrics: %s",
      metric,
      paste(intersect(ALPHA_METRICS, colnames(meta_enriched)),
            collapse = ", ")
    ))

  if (!x_var %in% colnames(meta_enriched))
    stop(sprintf("x_var '%s' not found in meta_enriched.", x_var))

  if (!color_var %in% colnames(meta_enriched))
    stop(sprintf("color_var '%s' not found in meta_enriched.", color_var))

  # ── Filter baseline rows ──────────────────────────────────
  df <- meta_enriched
  if (isTRUE(exclude_baseline) && is_baseline_col %in% colnames(df)) {
    n_before <- nrow(df)
    df <- df[df[[is_baseline_col]] == FALSE, , drop = FALSE]
    message(sprintf(
      "Excluded %d baseline sample(s) from alpha diversity plot.",
      n_before - nrow(df)
    ))
  }

  if (nrow(df) == 0)
    stop("No samples remaining after baseline exclusion.")

  # Remove rows where metric or x_var is NA
  df <- df[!is.na(df[[metric]]) & !is.na(df[[x_var]]), , drop = FALSE]
  df[[metric]] <- as.numeric(df[[metric]])


  # ── Colour levels (continued) ─────────────────────────────────
  col_lvls <- sort(unique(as.character(df[[color_var]])))
  col_vals <- stats::setNames(
    get_cat_pal(palette, length(col_lvls), rev_pal),
    col_lvls
  )
  df[[color_var]] <- factor(as.character(df[[color_var]]),
                            levels = col_lvls)
  df[[x_var]]     <- factor(as.character(df[[x_var]]))

  # ── Build plot ────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x      = .data[[x_var]],
      y      = .data[[metric]],
      colour = .data[[color_var]]
    )
  )

  # ── Paired connecting lines (drawn UNDER points) ──────────────
  # Uses CatID to connect each cat's samples across x groups
  if (!is.null(connect_by) && connect_by %in% colnames(df)) {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(group = .data[[connect_by]]),
        colour    = "grey70",
        linewidth = 0.4,
        alpha     = 0.6
      )
  }

  # ── Boxplot (no outlier points — jitter handles them) ─────────
  p <- p +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = .data[[x_var]]),
      outlier.shape = NA,
      colour        = "grey45",
      fill          = NA,
      width         = 0.35,
      linewidth     = 0.5
    )

  # ── Jittered individual points ────────────────────────────────
  p <- p +
    ggplot2::geom_jitter(
      size   = marker_size,
      width  = 0.08,
      alpha  = marker_alpha
    ) +
    ggplot2::scale_colour_manual(values = col_vals,
                                 name   = color_var)

  # ── Optional faceting (e.g. by Timepoint / Phase) ─────────────
  if (!is.null(facet_var) && facet_var %in% colnames(df)) {
    p <- p +
      ggplot2::facet_wrap(
        stats::as.formula(paste("~", facet_var)),
        scales = "free_x"
      )
  }

  # ── Labels and theme ──────────────────────────────────────────
  # Build informative subtitle
  n_subjects <- if (!is.null(connect_by) && connect_by %in% colnames(df))
    length(unique(df[[connect_by]])) else nrow(df)

  subtitle_txt <- sprintf(
    "%d samples | %d subjects | grouped by: %s | coloured by: %s",
    nrow(df),
    n_subjects,
    x_var,
    color_var
  )

  p <- p +
    ggplot2::labs(
      title    = sprintf("Alpha Diversity \u2014 %s", metric),
      subtitle = subtitle_txt,
      x        = x_var,
      y        = metric
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face  = "bold",
                                               hjust = 0.5,
                                               size  = 16),
      plot.subtitle    = ggplot2::element_text(hjust  = 0.5,
                                               colour = "grey40",
                                               size   = 11),
      axis.title       = ggplot2::element_text(face = "bold"),
      legend.title     = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92"),
      strip.text       = ggplot2::element_text(face = "bold")
    )

  p
}

# ── build_taxon_barplot() ─────────────────────────────────────────────────────
#' Build a stacked taxonomy bar plot
#'
#' @param taxon_data  Data frame returned by
#'   \code{compute_taxonomy_barplot_data()}.
#' @param level       Taxonomic level label used in the plot title.
#' @param group_var   Metadata column used for faceting.
#' @param palette     RColorBrewer palette name (categorical).
#' @param rev_pal     Logical — reverse the palette?
#' @param top_n       Number of top taxa shown (for subtitle).
#' @param exclude_baseline  Logical — remove baseline samples?
#' @param is_baseline_col   Column name of the baseline flag.
#' @param x_label_var Metadata column whose values are shown as the
#'   x-axis tick labels (default \code{"SampleID"} — the raw sample
#'   ID, matching the previous fixed behaviour). Bars themselves are
#'   always keyed by the unique \code{SampleID}, so choosing a column
#'   that isn't unique per sample (e.g. \code{"Timepoint"}) will just
#'   show repeated tick text for different bars, not merge them.
#' @param x_order_var Metadata column used to order samples along the
#'   x-axis within each facet (ties broken by \code{SampleID} for a
#'   stable order). Default \code{"SampleID"} sorts alphabetically,
#'   matching the previous fixed behaviour.
#' @return A \code{ggplot} object.
#' @export
build_taxon_barplot <- function(taxon_data,
                                level             = "Family",
                                group_var,
                                palette           = "Set3",
                                rev_pal           = FALSE,
                                top_n             = TAXA_DEFAULT_TOP_N,
                                exclude_baseline  = TRUE,
                                is_baseline_col   = "IsBaseline",
                                x_label_var       = "SampleID",
                                x_order_var       = "SampleID") {

  df <- taxon_data

  # ── Optionally remove baseline samples ───────────────────────
  if (isTRUE(exclude_baseline) && is_baseline_col %in% colnames(df)) {
    df <- df[df[[is_baseline_col]] != TRUE | is.na(df[[is_baseline_col]]),
             , drop = FALSE]
  }

  if (nrow(df) == 0)
    stop("No data remaining after filtering. Check exclude_baseline setting.")

  if (!x_label_var %in% colnames(df))
    stop(sprintf("x_label_var '%s' not found in taxon_data.", x_label_var))
  if (!x_order_var %in% colnames(df))
    stop(sprintf("x_order_var '%s' not found in taxon_data.", x_order_var))

  # ── X-axis ordering ───────────────────────────────────────────
  # Bars are always keyed by the unique SampleID (so different
  # samples never merge), but the DISPLAY ORDER along the axis is
  # driven by x_order_var, with SampleID as a stable tie-breaker.
  sample_id_chr <- as.character(df$SampleID)
  sample_order  <- unique(
    sample_id_chr[order(df[[x_order_var]], sample_id_chr)]
  )
  df$SampleID <- factor(sample_id_chr, levels = sample_order)

  # ── X-axis labels ──────────────────────────────────────────────
  # Look up each sample's chosen label column and display that text
  # in place of the raw SampleID, without changing which bar is
  # which (bar identity/position still keyed by SampleID above).
  label_map <- stats::setNames(
    as.character(df[[x_label_var]])[match(sample_order, sample_id_chr)],
    sample_order
  )

  # ── Colour palette ─────────────────────────────────────────────
  # "Others" is always grey; named taxa get palette colours
  taxa_named <- sort(setdiff(unique(df$Taxon), "Others"))
  n_col      <- max(3L, length(taxa_named))   # RColorBrewer min = 3

  pal_cols <- tryCatch(
    RColorBrewer::brewer.pal(
      min(n_col, RColorBrewer::brewer.pal.info[palette, "maxcolors"]),
      palette
    ),
    error = function(e) grDevices::rainbow(n_col)
  )
  if (isTRUE(rev_pal)) pal_cols <- rev(pal_cols)

  # Recycle if more taxa than colours
  colour_map        <- stats::setNames(
    rep_len(pal_cols, length(taxa_named)),
    taxa_named
  )
  colour_map["Others"] <- "#CCCCCC"

  # Fix factor order: top taxa first (alphabetical), Others last
  taxon_order    <- c(taxa_named, "Others")
  df$Taxon       <- factor(df$Taxon, levels = taxon_order)

  # ── Build plot ─────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x    = SampleID,
      y    = RelAbundance,
      fill = Taxon
    )
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.85) +
    ggplot2::scale_fill_manual(values = colour_map) +
    ggplot2::scale_x_discrete(labels = label_map) +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", group_var)),
      scales = "free_x"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5, size = 8
      ),
      panel.grid       = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#f0f0f0"),
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      title    = sprintf("%s composition (Top %d taxa) grouped by %s",
                         level, top_n, group_var),
      x        = if (identical(x_label_var, "SampleID")) "Sample"
      else x_label_var,
      y        = "Relative abundance (%)",
      fill     = "Taxon"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    # coord_cartesian() only zooms the displayed range — unlike
    # scale_y_continuous(limits = ...), it never censors/drops data
    # whose stacked cumulative position creeps slightly outside
    # [0, 100] due to floating-point summation of many small
    # percentages. Using scale limits here was silently deleting the
    # topmost stacked segment for whichever sample's running total
    # happened to round a hair above 100 at a given taxonomic level —
    # which taxon totals are being summed (and hence which samples
    # are affected) differs by level, matching the pattern of
    # different "short" samples at each rank.
    ggplot2::coord_cartesian(ylim = c(0, 100))

  p
}
