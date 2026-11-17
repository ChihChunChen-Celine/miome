# plots.R
# -----------------------------------------------------------------------------
# Plot builders: 3D PCoA (plotly), 2D PCoA, alpha diversity, taxonomy bars.
# TIFF export wrappers live in export.R.
# -----------------------------------------------------------------------------

# -- Internal helpers ----------------------------------------------------------

# -- Custom colours ---------------------------------------------------------------
# Wherever a plot takes `palette`, it may be
#   * a palette name ("Miome", "Set2", ...; for gradients "RdBu", ...),
#   * a NAMED vector of colours, one per group you want to fix, e.g.
#     c(Control = "#1B9E77", MUO = "#D95F02"). Groups not named keep the colour of
#     the base palette (attribute "base", default "Miome"),
#   * for gradients, an UNNAMED vector of two or more colours, low to high.
.is_colour <- function(x) is.character(x) & grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", x)

# One colour per level: base palette first, then the user's fixed colours on top.
.level_colours <- function(palette, lvls, reverse = FALSE) {
  custom <- if (!is.null(names(palette))) palette[nzchar(names(palette))] else character(0)
  base   <- attr(palette, "base")
  if (is.null(base)) base <- if (length(custom) || length(palette) != 1) "Miome" else palette
  cols <- stats::setNames(get_cat_pal(base, length(lvls), reverse), lvls)
  hit  <- intersect(names(custom)[.is_colour(custom)], lvls)
  cols[hit] <- custom[hit]
  cols
}

# n distinct colours from a palette (interpolates beyond its size).
get_cat_pal <- function(pal, n, reverse = FALSE) {
  if (identical(pal, "Miome")) {
    cols <- if (n <= length(MIOME_PALETTE)) MIOME_PALETTE[seq_len(n)]
            else grDevices::colorRampPalette(MIOME_PALETTE)(n)
    return(if (reverse) rev(cols) else cols)
  }
  max_cols <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
  base     <- RColorBrewer::brewer.pal(max(3, min(n, max_cols)), pal)
  cols     <- if (n <= max_cols) base[seq_len(n)]
              else grDevices::colorRampPalette(base)(n)
  if (reverse) rev(cols) else cols
}

# Continuous ramp (used by both 2D and 3D plots).
.cont_ramp <- function(pal, reverse = FALSE, n = 11) {
  if (length(pal) >= 2 && all(.is_colour(pal))) {             # custom gradient, low -> high
    hex <- grDevices::colorRampPalette(unname(pal))(n)
    return(if (reverse) rev(hex) else hex)
  }
  max_cols <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
  hex <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(11, max_cols), pal))(n)
  if (reverse) rev(hex) else hex
}

make_colorscale <- function(pal_name, reverse = FALSE, n_steps = 11) {
  hex <- .cont_ramp(pal_name, reverse, n_steps)
  pos <- seq(0, 1, length.out = n_steps)
  lapply(seq_along(hex), function(i) list(pos[i], hex[i]))
}

.ax_lab <- function(var_pct, i) sprintf("PC%d (%.1f%%)", i, var_pct[i])

# Rows ordered so paired lines run subject by subject in time order.
.order_for_paths <- function(df, connect_by, time_col) {
  keys <- list(df[[connect_by]])
  if (!is.null(time_col) && time_col %in% colnames(df)) keys <- c(keys, list(df[[time_col]]))
  keys <- c(keys, list(df$SampleID))
  df[do.call(order, keys), , drop = FALSE]
}

.has_col <- function(df, col) !is.null(col) && length(col) == 1 && col %in% colnames(df)

# -- 1. build_3d() -------------------------------------------------------------
#' Build a 3D interactive PCoA plotly figure
#'
#' @param obj A \code{pcoa_obj} from \code{prepare_all()}.
#' @param v Colour variable (metadata column).
#' @param palette Palette name, or your own colours: a named vector such as
#'   \code{c(Control = "#1B9E77", MUO = "#D95F02")} fixes the colour of those
#'   groups (the others keep the Miome palette); for a numeric variable an
#'   unnamed vector such as \code{c("#2F6FB0", "#F2F2F2", "#B8456B")} is a
#'   gradient from low to high.
#' @param is_cont Logical; treat \code{v} as continuous?
#' @param rev_pal Logical; reverse the palette?
#' @param connect_by Optional column for within-subject lines (categorical
#'   mode only). Lines follow the configured timepoint column.
#' @return A \code{plotly} object.
#' @export
build_3d <- function(obj, v, palette = "Set2", is_cont = FALSE,
                     rev_pal = FALSE, connect_by = NULL) {
  cfg <- obj$config
  df  <- obj$pcoa_df[!is.na(obj$pcoa_df[[v]]), ]
  hover <- function(d) paste0("<b>", d$SampleID, "</b><br>", v, ": ", d[[v]])
  marker <- function(...) list(size = cfg$marker_size, opacity = cfg$marker_opacity,
                               line = list(width = 0.5, color = "white"), ...)
  fig <- plotly::plot_ly()

  if (!is_cont) {
    lvls <- sort(unique(as.character(df[[v]])))
    cols <- .level_colours(palette, lvls, rev_pal)
    for (lvl in lvls) {
      s <- df[as.character(df[[v]]) == lvl, ]
      fig <- plotly::add_trace(fig, type = "scatter3d", mode = "markers",
                               x = s$PC1, y = s$PC2, z = s$PC3, name = lvl,
                               text = hover(s), hoverinfo = "text",
                               marker = marker(color = cols[[lvl]]))
    }
    if (.has_col(df, connect_by)) {
      dp <- .order_for_paths(df, connect_by, cfg$timepoint_col)
      for (subj in unique(dp[[connect_by]])) {
        s <- dp[dp[[connect_by]] %in% subj, ]
        if (nrow(s) >= 2)
          fig <- plotly::add_trace(fig, type = "scatter3d", mode = "lines",
                                   x = s$PC1, y = s$PC2, z = s$PC3,
                                   name = as.character(subj), showlegend = FALSE,
                                   hoverinfo = "skip",
                                   line = list(color = "grey70", width = 1.5))
      }
    }
  } else {
    vals <- as.numeric(df[[v]])
    fig <- plotly::add_trace(fig, type = "scatter3d", mode = "markers",
      x = df$PC1, y = df$PC2, z = df$PC3, name = v, text = hover(df),
      hoverinfo = "text", showlegend = FALSE,
      marker = marker(color = vals, colorscale = make_colorscale(palette, rev_pal),
                      cmin = min(vals, na.rm = TRUE), cmax = max(vals, na.rm = TRUE),
                      showscale = TRUE,
                      colorbar = list(title = list(text = v, side = "right"),
                                      thickness = 15, len = 0.5, x = 1.02)))
  }

  ax <- function(i) list(title = .ax_lab(obj$var_pct, i), showgrid = TRUE, zeroline = FALSE)
  fig <- plotly::layout(fig,
    title  = list(text = sprintf("%s | %s (n=%d)", cfg$plot_title,
                                 cfg$dist_method_label, nrow(df)),
                  font = list(size = 14)),
    scene  = list(xaxis = ax(1), yaxis = ax(2), zaxis = ax(3),
                  camera = list(eye = list(x = 1.5, y = 1.5, z = 1.0))),
    legend = list(title = list(text = paste0("<b>", v, "</b>")), itemsizing = "constant"),
    margin = list(t = 60, b = 10, l = 10, r = 10))

  plotly::config(fig,
    toImageButtonOptions = list(format = "png", filename = "pcoa_3D_snapshot",
                                width = cfg$export_3d_w_px,
                                height = cfg$export_3d_h_px, scale = 1),
    displaylogo = FALSE)
}

# -- 2. build_2d() -------------------------------------------------------------
#' Build a 2D PCoA ggplot
#'
#' @inheritParams build_3d
#' @param ax_x,ax_y PC axes for x and y (1-3).
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @export
build_2d <- function(obj, v, ax_x = 1, ax_y = 2, palette = "Set2",
                     is_cont = FALSE, rev_pal = FALSE, connect_by = NULL) {
  cfg  <- obj$config
  pc_x <- paste0("PC", as.integer(ax_x))
  pc_y <- paste0("PC", as.integer(ax_y))
  df   <- obj$pcoa_df[!is.na(obj$pcoa_df[[v]]), ]

  grp_txt <- if (is_cont) {
    r <- range(as.numeric(df[[v]]), na.rm = TRUE)
    sprintf("range %s\u2013%s", round(r[1], 1), round(r[2], 1))
  } else paste(sort(unique(as.character(df[[v]]))), collapse = ", ")

  p <- ggplot(df, aes(x = .data[[pc_x]], y = .data[[pc_y]]))

  if (!is_cont && .has_col(df, connect_by))
    p <- p + geom_path(data = .order_for_paths(df, connect_by, cfg$timepoint_col),
                       aes(group = .data[[connect_by]]),
                       colour = "grey75", linewidth = 0.4, alpha = 0.7)

  if (is_cont) {
    p <- p +
      geom_point(aes(color = as.numeric(.data[[v]])),
                 size = cfg$marker_size / 2, alpha = cfg$marker_opacity) +
      scale_color_gradientn(colours = .cont_ramp(palette, rev_pal), name = v)
  } else {
    lvls <- sort(unique(as.character(df[[v]])))
    df[[v]] <- factor(as.character(df[[v]]), levels = lvls)
    p <- p +
      geom_point(data = df, aes(fill = .data[[v]]), size = cfg$marker_size / 2,
                 alpha = cfg$marker_opacity, shape = 21, colour = "white",
                 stroke = 0.4) +
      scale_fill_manual(values = .level_colours(palette, lvls, rev_pal), name = v)
  }

  p +
    labs(title = sprintf("PCoA \u2014 coloured by: %s", v),
         subtitle = sprintf("%s | %s | n = %d", cfg$dist_method_label, grp_txt, nrow(df)),
         x = .ax_lab(obj$var_pct, ax_x), y = .ax_lab(obj$var_pct, ax_y)) +
    .miome_theme(14)
}

.miome_theme <- function(base_size) {
  theme_bw(base_size = base_size) +
    theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = base_size + 3),
          plot.subtitle    = element_text(hjust = 0.5, colour = "grey40", size = base_size - 2),
          axis.title       = element_text(face = "bold"),
          legend.title     = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text       = element_text(face = "bold"))
}

# -- 3. build_alpha_scatter() --------------------------------------------------
#' Alpha diversity boxplot with jittered points and optional paired lines
#'
#' @param meta_enriched Data frame from \code{pcoa_obj$meta_enriched}.
#' @param metric Alpha diversity column to plot.
#' @param x_var Column for x-axis groups.
#' @param color_var Column for point colour.
#' @param palette Palette name, or your own colours: a named vector such as
#'   \code{c(Control = "#1B9E77", MUO = "#D95F02")} fixes the colour of those
#'   groups (the others keep the Miome palette); for a numeric variable an
#'   unnamed vector such as \code{c("#2F6FB0", "#F2F2F2", "#B8456B")} is a
#'   gradient from low to high.
#' @param rev_pal Logical; reverse palette?
#' @param connect_by Optional column for within-subject lines, or \code{NULL}.
#' @param exclude_baseline Logical; drop baseline rows first.
#' @param is_baseline_col Baseline flag column.
#' @param facet_var Optional column to facet by, or \code{NULL}.
#' @param marker_size Point size.
#' @param marker_alpha Point opacity.
#' @param title,y_label Optional plot title and y-axis label (defaults describe
#'   an alpha-diversity metric; set them when plotting another measurement).
#' @return A \code{ggplot} object.
#' @export
build_alpha_scatter <- function(meta_enriched,
                                metric           = "Shannon",
                                x_var            = "Treatment",
                                color_var        = "Timepoint",
                                palette          = "Set2",
                                rev_pal          = FALSE,
                                connect_by       = "CatID",
                                exclude_baseline = TRUE,
                                is_baseline_col  = "IsBaseline",
                                facet_var        = NULL,
                                marker_size      = 3,
                                marker_alpha     = 0.85,
                                title            = NULL,
                                y_label          = NULL) {
  for (col in c(metric, x_var, color_var))
    if (!col %in% colnames(meta_enriched))
      stop(sprintf("Column '%s' not found in the metadata.", col), call. = FALSE)

  df <- meta_enriched
  if (isTRUE(exclude_baseline)) df <- df[!.is_baseline(df, is_baseline_col), , drop = FALSE]
  df <- df[!is.na(df[[metric]]) & !is.na(df[[x_var]]), , drop = FALSE]
  if (!nrow(df)) stop("No samples left to plot.", call. = FALSE)

  df[[metric]]    <- as.numeric(df[[metric]])
  col_lvls        <- sort(unique(as.character(df[[color_var]])))
  df[[color_var]] <- factor(as.character(df[[color_var]]), levels = col_lvls)
  df[[x_var]]     <- factor(as.character(df[[x_var]]))
  df$SampleID     <- rownames(df)

  # One jitter shared by points and lines, so lines end on the points;
  # no vertical jitter, so integer metrics keep their true values.
  jit   <- position_jitter(width = 0.08, height = 0, seed = 1)
  pairs <- .has_col(df, connect_by)
  if (pairs) df <- .order_for_paths(df, connect_by, x_var)

  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[metric]]))
  if (pairs)
    p <- p + geom_path(aes(group = .data[[connect_by]]), position = jit,
                       colour = "grey70", linewidth = 0.4, alpha = 0.6)
  p <- p +
    geom_boxplot(aes(group = .data[[x_var]]), outlier.shape = NA,
                 colour = "grey45", fill = NA, width = 0.35, linewidth = 0.5) +
    geom_point(aes(colour = .data[[color_var]]), position = jit,
               size = marker_size, alpha = marker_alpha) +
    scale_colour_manual(values = .level_colours(palette, col_lvls, rev_pal), name = color_var)

  if (.has_col(df, facet_var))
    p <- p + facet_wrap(stats::as.formula(paste0("~`", facet_var, "`")), scales = "free_x")

  n_subj <- if (pairs) length(unique(df[[connect_by]])) else nrow(df)
  p + labs(title = if (is.null(title)) sprintf("Alpha Diversity \u2014 %s", metric) else title,
           subtitle = sprintf("%d samples | %d subjects | grouped by: %s | coloured by: %s",
                              nrow(df), n_subj, x_var, color_var),
           x = x_var, y = if (is.null(y_label)) metric else y_label) +
    .miome_theme(13)
}

# -- 4. build_taxon_barplot() --------------------------------------------------
#' Build a stacked taxonomy bar plot
#'
#' @param taxon_data Data frame from \code{compute_taxonomy_barplot_data()}.
#' @param level Taxonomic level label for the title.
#' @param group_var Metadata column used for faceting.
#' @param palette Palette name, or your own colours: a named vector such as
#'   \code{c(Control = "#1B9E77", MUO = "#D95F02")} fixes the colour of those
#'   groups (the others keep the Miome palette); for a numeric variable an
#'   unnamed vector such as \code{c("#2F6FB0", "#F2F2F2", "#B8456B")} is a
#'   gradient from low to high.
#' @param rev_pal Logical; reverse the palette?
#' @param top_n Number of top taxa (for the title).
#' @param exclude_baseline Logical; remove baseline samples?
#' @param is_baseline_col Baseline flag column.
#' @param x_label_var Column shown as x-axis tick labels (bars stay keyed by
#'   \code{SampleID}).
#' @param x_order_var Column used to order bars (ties broken by SampleID).
#' @return A \code{ggplot} object.
#' @export
build_taxon_barplot <- function(taxon_data,
                                level            = "Family",
                                group_var,
                                palette          = "Set3",
                                rev_pal          = FALSE,
                                top_n            = TAXA_DEFAULT_TOP_N,
                                exclude_baseline = TRUE,
                                is_baseline_col  = "IsBaseline",
                                x_label_var      = "SampleID",
                                x_order_var      = "SampleID") {
  df <- taxon_data
  if (isTRUE(exclude_baseline)) df <- df[!.is_baseline(df, is_baseline_col), , drop = FALSE]
  if (!nrow(df)) stop("No data left after filtering. Check 'Exclude baseline'.", call. = FALSE)
  for (col in c(x_label_var, x_order_var, group_var))
    if (!col %in% colnames(df)) stop(sprintf("Column '%s' not found.", col), call. = FALSE)

  ids   <- as.character(df$SampleID)
  order <- unique(ids[order(df[[x_order_var]], ids)])
  df$SampleID <- factor(ids, levels = order)
  labels <- stats::setNames(as.character(df[[x_label_var]])[match(order, ids)], order)

  taxa   <- sort(setdiff(unique(df$Taxon), "Others"))
  colmap <- c(.level_colours(palette, taxa, isTRUE(rev_pal)), Others = "#CCCCCC")
  if ("Others" %in% names(palette) && .is_colour(palette[["Others"]])) colmap["Others"] <- palette[["Others"]]
  df$Taxon <- factor(df$Taxon, levels = c(taxa, "Others"))

  ggplot(df, aes(x = SampleID, y = RelAbundance, fill = Taxon)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = colmap) +
    scale_x_discrete(labels = labels) +
    facet_wrap(stats::as.formula(paste0("~`", group_var, "`")), scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    # coord_cartesian zooms without dropping bars whose stacked total rounds
    # a hair above 100 (scale limits would delete the top segment).
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = sprintf("%s composition (Top %d taxa) grouped by %s", level, top_n, group_var),
         x = if (identical(x_label_var, "SampleID")) "Sample" else x_label_var,
         y = "Relative abundance (%)", fill = "Taxon") +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "#f0f0f0"),
          legend.text = element_text(size = 9))
}
