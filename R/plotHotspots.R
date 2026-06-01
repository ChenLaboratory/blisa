# Internal helper: resolve a row index from optional ligand/receptor args.
# Returns `index` unchanged when neither ligand nor receptor is supplied.
.resolve_lr_index <- function(LR_results, index, ligand, receptor) {
  if (!is.null(ligand) && !is.null(receptor)) {
    matches <- which(LR_results$ligand.symbol == ligand &
                       LR_results$receptor.symbol == receptor)
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

  # Pre-compute a fill colour for every bin
  fill_cols <- rep("#F0F0F0", nrow(bins))    # empty
  fill_cols[bins$n_cells > 0] <- "#D3D3D3"  # non-significant

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

  p <- ggplot(bins) +
    geom_sf(aes(fill = fill_col), color = NA) +
    scale_fill_identity() +
    geom_point(data = dummy, aes(x = x, y = y, color = val), alpha = 0) +
    scale_color_gradientn(colours = LRI_spatial_colors, name = lgd_title) +
    guides(color = guide_colorbar(title.position = "top")) +
    labs(
      title    = paste0(interaction, " - ", sig_num, " hotspots", title_suffix),
      subtitle = "Grey scale: Light (Empty bins) | Medium (Non-sig bins)"
    ) +
    theme_void() +
    theme(legend.position = "right")

  return(p)
}
