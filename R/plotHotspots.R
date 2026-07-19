# Internal helper: resolve a row index from optional ligand/receptor args.
# Returns `index` unchanged when neither ligand nor receptor is supplied.
# Multi-subunit complexes match on the full subunit set, regardless of the
# order/separator the user supplies (e.g. "TGFBR2,TGFBR1", "TGFBR1_TGFBR2",
# and "TGFBR1|TGFBR2" all match a stored "TGFBR2, TGFBR1").
.resolve_lr_index <- function(LR_results, index, ligand, receptor) {
  if (!is.null(ligand) && !is.null(receptor)) {
    want_l <- .norm_units(ligand)
    want_r <- .norm_units(receptor)
    same <- function(col, want)
      vapply(col, function(s) identical(.norm_units(s), want), logical(1))
    matches <- which(same(LR_results$ligand.symbol, want_l) &
                       same(LR_results$receptor.symbol, want_r))
    if (length(matches) == 0)
      stop("No LR pair found with ligand = '", ligand,
           "' and receptor = '", receptor, "'.")
    if (length(matches) > 1)
      warning("Multiple rows match ligand = '", ligand, "' and receptor = '",
              receptor, "'. Using the first match (index ", matches[1], ").")
    return(matches[1])
  }
  if (xor(!is.null(ligand), !is.null(receptor)))
    stop("'ligand' and 'receptor' must be provided together.")
  index
}


#' Spatial hotspot map for one ligand-receptor pair
#'
#' Generic function. Plots each bin coloured by significance status: empty,
#' non-significant, or significant hotspot (continuous gradient of
#' -log10 p-value or 1 - p-value).
#'
#' @param x A \code{blisa} object.
#' @param ... Additional arguments passed to the method.
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' # Continuing from the blisa() example:
#' # result <- blisa(spe, bin_size = 50, group = "cell_type")
#' plotHotspots(result, index = 1)
#' plotHotspots(result, index = 1, log_pval = FALSE)
#' }
#' @export
plotHotspots <- function(x, ...) UseMethod("plotHotspots")


#' @describeIn plotHotspots Method for a \code{blisa} object.
#'
#' @param index Integer. Row index into \code{LR_results} selecting the
#'   ligand-receptor pair to visualise. Ignored when both \code{ligand} and
#'   \code{receptor} are supplied. Default \code{1} (top-ranked pair).
#' @param ligand Character. Ligand gene symbol. When both \code{ligand} and
#'   \code{receptor} are provided the matching LR pair is located automatically
#'   and \code{index} is ignored. Must be supplied together with
#'   \code{receptor}.
#' @param receptor Character. Receptor gene symbol. Must be supplied together
#'   with \code{ligand}.
#' @param as_points Logical. If \code{TRUE}, draw each bin as a dot at its
#'   centroid instead of its polygon. Useful for Visium, where each bin is a
#'   single spot. Default \code{FALSE} (draw polygons).
#' @param size Numeric. Point size when \code{as_points = TRUE} (or when
#'   \code{background} is supplied). Default \code{1.5}.
#' @param background A \code{ggplot} object to draw the hotspots on top of, e.g.
#'   \code{scider::plotImage(spe)} to place the H&E image behind the spots.
#'   When supplied, the hotspots are always rendered as dots (the raster
#'   background is incompatible with \code{geom_sf} coordinates), and the bins
#'   must share the background plot's coordinate frame (true for objects read by
#'   \code{scider::readVisium()}). Default \code{NULL}.
#' @param log_pval Logical. If \code{TRUE} (default), colour significant bins
#'   by -log10(p-value). If \code{FALSE}, use 1 - p-value.
#' @param p_cutoff Numeric or \code{NULL}. When \code{NULL} (default), the
#'   pre-computed hotspot bins stored in the \code{blisa} object are used,
#'   reflecting the \code{p_cutoff} and High-High quadrant classification
#'   applied during \code{\link{blisa}}. When a numeric value is supplied,
#'   bins are re-defined on the fly as those with \code{all_pval <= p_cutoff}
#'   and quadrant label \code{"High-High"} (from the stored
#'   \code{all_quadrant}), giving an exact re-threshold consistent with the
#'   original classification.
#'
#' @export
plotHotspots.blisa <- function(x, index = 1, ligand = NULL, receptor = NULL,
                               as_points = FALSE, size = 1.5, background = NULL,
                               log_pval = TRUE, p_cutoff = NULL, ...) {
  bins       <- x$bins
  LR_results <- x$LR_results

  # Resolve which LR pair to plot
  index <- .resolve_lr_index(LR_results, index, ligand, receptor)

  interaction <- rownames(LR_results)[index]

  # Resolve significant bin indices and their p-values
  if (!is.null(p_cutoff)) {
    all_pval      <- LR_results$all_pval[[index]]
    all_quadrant  <- LR_results$all_quadrant[[index]]
    if (is.null(all_pval) || is.null(all_quadrant))
      stop("'all_pval' and 'all_quadrant' must be present in LR_results to ",
           "use 'p_cutoff'. Ensure blisa() was run to completion.")
    sig_indices <- which(all_pval <= p_cutoff & all_quadrant == "High-High")
    p_values    <- all_pval[sig_indices]
  } else {
    sig_indices <- LR_results$sig_index[[index]]
    p_values    <- LR_results$sig_pval[[index]]
  }

  sig_num <- length(sig_indices)

  if (log_pval) {
    plot_vals <- -log10(p_values)
    lgd_title <- "-log10(pval)"
  } else {
    plot_vals <- 1 - p_values
    lgd_title <- "1-pval"
  }

  # Pre-compute a fill colour for every bin. Bins not included in this LR
  # pair's LISA test (empty, isolated, low-cell, or low total counts) are
  # treated the same as empty bins; only tested bins can be "non-significant".
  tested    <- !is.na(LR_results$all_quadrant[[index]])
  fill_cols <- rep("#FFFFFF", nrow(bins))    # empty / excluded from testing (white)
  fill_cols[tested] <- "#D3D3D3"             # tested but non-significant (light grey)

  if (length(sig_indices) > 0) {
    pval_range <- range(plot_vals)
    if (pval_range[1] == pval_range[2]) pval_range[2] <- pval_range[1] + 1
    pval_norm  <- (plot_vals - pval_range[1]) / diff(pval_range)
    pal        <- colorRampPalette(LRI_spatial_colors)(100)
    fill_cols[sig_indices] <- pal[as.integer(pval_norm * 99) + 1]
  }

  bins$fill_col <- fill_cols

  # Invisible dummy points to drive the continuous colorbar legend
  val_range <- if (length(plot_vals) > 0) range(plot_vals) else c(0, 1)
  bbox <- sf::st_bbox(bins)
  cx   <- unname((bbox["xmin"] + bbox["xmax"]) / 2)
  cy   <- unname((bbox["ymin"] + bbox["ymax"]) / 2)
  dummy <- data.frame(
    x   = cx,
    y   = cy,
    val = seq(val_range[1], val_range[2], length.out = 100)
  )

  title_suffix <- if (!is.null(p_cutoff)) paste0(" (p \u2264 ", p_cutoff, ")") else ""

  # Overlay on an existing background plot (e.g. scider::plotImage(spe)), which
  # already places the H&E image in the (micron) coordinate frame of the spots.
  # We add the hotspot dots on top -- the same pattern scider::plotSpatial uses.
  # Requires the bins to share that coordinate frame (true for objects read by
  # scider::readVisium()).
  if (!is.null(background)) {
    # Prefer the raw image-registered coordinates (stored by visiumSpotBins) so
    # the dots align with the H&E; fall back to bin centroids otherwise (which
    # may be in a de-tilted/rescaled analysis frame and not match the image).
    if (!is.null(bins$img_x) && !is.null(bins$img_y)) {
      pts_df <- data.frame(x = bins$img_x, y = bins$img_y,
                           fill_col = bins$fill_col)
    } else {
      cent   <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(bins)))
      pts_df <- data.frame(x = cent[, 1], y = cent[, 2],
                           fill_col = bins$fill_col)
    }
    p <- background +
      geom_point(data = pts_df, aes(x = x, y = y, fill = fill_col),
                 shape = 21, color = "transparent", size = size) +
      scale_fill_identity() +
      geom_point(data = dummy, aes(x = x, y = y, color = val), alpha = 0) +
      scale_color_gradientn(colours = LRI_spatial_colors, name = lgd_title) +
      guides(color = guide_colorbar(title.position = "top")) +
      labs(
        title    = paste0(interaction, " - ", sig_num, " hotspots", title_suffix),
        subtitle = "White: Empty/untested bins | Light grey: Non-significant bins"
      )
    return(p)
  }

  # Draw the bins either as polygons (geom_sf) or as dots at their centroids.
  # Dots use shape 21 so they carry the same `fill` scale as the polygons,
  # leaving the `color` scale free for the p-value colorbar legend.
  if (as_points) {
    cent   <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(bins)))
    pts_df <- data.frame(x = cent[, 1], y = cent[, 2], fill_col = bins$fill_col)
    # shape 21 = filled circle carrying the `fill` scale; use a transparent
    # (not NA) border, since colour = NA drops every point via remove_missing.
    bin_layer   <- geom_point(data = pts_df, aes(x = x, y = y, fill = fill_col),
                              shape = 21, color = "transparent", size = size)
    coord_layer <- ggplot2::coord_equal()
  } else {
    bin_layer   <- geom_sf(aes(fill = fill_col), color = NA)
    coord_layer <- NULL
  }

  p <- ggplot(bins) +
    bin_layer +
    coord_layer +
    scale_fill_identity() +
    geom_point(data = dummy, aes(x = x, y = y, color = val), alpha = 0) +
    scale_color_gradientn(colours = LRI_spatial_colors, name = lgd_title) +
    guides(color = guide_colorbar(title.position = "top")) +
    labs(
      title    = paste0(interaction, " - ", sig_num, " hotspots", title_suffix),
      subtitle = "White: Empty/untested bins | Light grey: Non-significant bins"
    ) +
    theme_void() +
    theme(legend.position = "right")

  return(p)
}
